// M6-11 research accel mode ('ode'): cross-path consistency, boundary
// continuity, analytic accuracy, lookahead-range validity, and sleep
// parity under the new acceleration definition. Default 'hermite'
// behavior is covered by the existing suites (bitwise unchanged).
import { describe, it, expect } from 'vitest'
import { singleToneScore, type PerformanceScore } from '../sequence'
import { TransientPlate, hermiteQdRange } from '../plateTransient'
import { createTransientSim, calibratedForceN, type TransientSandSim } from '../transientSand'
import { G } from '../surface'

function excitedPlate(mode: 'hermite' | 'ode', freq = 152, forceN = 0): TransientPlate {
  const plate = new TransientPlate(singleToneScore(freq, forceN, 30), [0.05, 0.03])
  plate.accelMode = mode
  let seed = 222
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  return plate
}

describe('M6-11 ode accel mode (plate)', () => {
  it('surfAt / heightAccelAt / accelAt agree bitwise in ode mode', () => {
    const plate = excitedPlate('ode', 649, calibratedForceN())
    const dt = 1 / (44 * 649)
    plate.beginStep(0, dt)
    const out = new Float64Array(5); const out2 = new Float64Array(2)
    let seed = 99
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    for (let k = 0; k < 500; k++) {
      const x = (rnd() * 2 - 1) * plate.a * 0.95
      const y = (rnd() * 2 - 1) * plate.a * 0.95
      const tau = rnd() * dt
      plate.surfAt(x, y, tau, out)
      plate.heightAccelAt(x, y, tau, out2)
      expect(out2[0]).toBe(out[0])
      expect(out2[1]).toBe(out[4])
      expect(plate.accelAt(x, y, tau)).toBe(out[4])
    }
  })

  it('ode wtt is continuous across segment boundaries where hermite jumps', () => {
    const dt = 1 / (44 * 152)
    const jump = (mode: 'hermite' | 'ode'): number => {
      const plate = excitedPlate(mode)
      const out = new Float64Array(5)
      plate.beginStep(0, dt)
      plate.surfAt(0.04, -0.03, dt * (1 - 1e-9), out)
      const wL = out[4]
      plate.beginStep(dt, 2 * dt)
      plate.surfAt(0.04, -0.03, dt * (1 + 1e-9), out)
      return Math.abs(out[4] - wL)
    }
    const jH = jump('hermite'); const jO = jump('ode')
    expect(jH).toBeGreaterThan(10) // excited case: hermite wtt really jumps
    expect(jO).toBeLessThan(jH * 0.01)
  })

  it('ode wtt is far closer to the analytic free-decay solution', () => {
    const err = (mode: 'hermite' | 'ode'): number => {
      const plate = excitedPlate(mode)
      const zeta = plate.zeta
      const q0s = Float64Array.from(plate.q); const qd0s = Float64Array.from(plate.qd)
      const phi = new Float64Array(plate.ndof)
      const qAna = (w: number, q0: number, qd0: number, t: number): number => {
        const wd = w * Math.sqrt(1 - zeta * zeta)
        return Math.exp(-zeta * w * t) * (q0 * Math.cos(wd * t) + ((qd0 + zeta * w * q0) / wd) * Math.sin(wd * t))
      }
      const qdAna = (w: number, q0: number, qd0: number, t: number): number => {
        const wd = w * Math.sqrt(1 - zeta * zeta)
        const c = (qd0 + zeta * w * q0) / wd
        const e = Math.exp(-zeta * w * t)
        return e * (-zeta * w * (q0 * Math.cos(wd * t) + c * Math.sin(wd * t)) + wd * (-q0 * Math.sin(wd * t) + c * Math.cos(wd * t)))
      }
      const dt = 1 / (44 * 152)
      const out = new Float64Array(5)
      let segT1 = 0; let maxE = 0
      for (let t = dt / 4; t <= 0.005; t += dt / 4) {
        while (segT1 < t) { plate.beginStep(segT1, segT1 + dt); segT1 += dt }
        for (const [x, y] of [[0.03, -0.02], [-0.05, 0.04]] as const) {
          plate.surfAt(x, y, t, out)
          plate.phiAll(x, y, phi)
          let wA = 0
          for (let j = 0; j < plate.ndof; j++) {
            const w = plate.omega[j]
            wA += (-2 * zeta * w * qdAna(w, q0s[j], qd0s[j], t) - w * w * qAna(w, q0s[j], qd0s[j], t)) * phi[j]
          }
          maxE = Math.max(maxE, Math.abs(out[4] - wA))
        }
      }
      return maxE
    }
    const eH = err('hermite'); const eO = err('ode')
    expect(eH).toBeGreaterThan(100)
    expect(eO).toBeLessThan(eH / 5)
  })

  it('hermiteQdRange contains dense velocity samples', () => {
    let seed = 4242
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let k = 0; k < 200; k++) {
      const q0 = rnd() * 1e-3; const qd0 = rnd(); const q1 = rnd() * 1e-3; const qd1 = rnd()
      const h = 1e-4 + Math.abs(rnd()) * 1e-3
      const [lo, hi] = hermiteQdRange(q0, qd0, q1, qd1, h)
      for (let sN = 0; sN <= 20; sN++) {
        const s = sN / 20; const s2 = s * s
        const qd = ((6 * s2 - 6 * s) * q0 + (3 * s2 - 4 * s + 1) * h * qd0
          + (-6 * s2 + 6 * s) * q1 + (3 * s2 - 2 * s) * h * qd1) / h
        expect(qd).toBeGreaterThanOrEqual(lo - 1e-12)
        expect(qd).toBeLessThanOrEqual(hi + 1e-12)
      }
    }
  })

  it('ode lookahead ranges bound the implemented acceleration', () => {
    const plate = excitedPlate('ode', 152, calibratedForceN())
    const dt = 1 / (44 * 152)
    const lc = plate.lookaheadRanges(0, 0.02, () => dt)
    const phi = new Float64Array(plate.ndof)
    let segT1 = 0
    let seed = 1234
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const times = Array.from({ length: 300 }, () => rnd() * 0.019).sort((a, b) => a - b)
    for (const t of times) {
      while (segT1 < t) { plate.beginStep(segT1, segT1 + dt); segT1 += dt }
      const x = (rnd() * 2 - 1) * plate.a * 0.95
      const y = (rnd() * 2 - 1) * plate.a * 0.95
      plate.phiAll(x, y, phi)
      let lb = 0; let ub = 0
      for (let j = 0; j < plate.ndof; j++) {
        lb += Math.min(phi[j] * lc.aMin[j], phi[j] * lc.aMax[j])
        ub += Math.max(phi[j] * lc.aMin[j], phi[j] * lc.aMax[j])
      }
      const w = plate.accelAt(x, y, t)
      expect(w).toBeGreaterThanOrEqual(lb - 1e-6)
      expect(w).toBeLessThanOrEqual(ub + 1e-6)
    }
  })
})

describe('M6-11 ode accel mode (Codex review regressions)', () => {
  it('directly forced DOFs are included even with zero Hermite endpoints', () => {
    // Short force event inside one segment, plate at rest: RK4 substeps
    // need not sample the event, so Hermite endpoints can stay exactly
    // zero while the ODE acceleration at mid-event is nonzero. The
    // live-DOF mask must not drop those DOFs from accelAt.
    const score: PerformanceScore = {
      schema: 'physicymatics.score.v1', duration_s: 30, initialPhase_rad: Math.PI / 2, audioGain: 0.2,
      events: [{ start_s: 1e-6, end_s: 2e-6, frequency_hz: 152, forcePeak_N: 10, attack_s: 0, release_s: 0 }],
    }
    const plate = new TransientPlate(score, [0, 0])
    plate.accelMode = 'ode'
    const dt = 1 / (44 * 152)
    plate.beginStep(0, dt)
    const out = new Float64Array(5)
    plate.surfAt(0.02, 0.01, 1.5e-6, out)
    const wA = plate.accelAt(0.02, 0.01, 1.5e-6)
    expect(Math.abs(out[4])).toBeGreaterThan(1)
    expect(wA).toBe(out[4])
    // Same when the mode is switched only after beginStep: the setter
    // must rebuild the current segment's live mask, not just caches.
    const plate2 = new TransientPlate(score, [0, 0])
    plate2.beginStep(0, dt)
    plate2.accelMode = 'ode'
    plate2.surfAt(0.02, 0.01, 1.5e-6, out)
    expect(Math.abs(out[4])).toBeGreaterThan(1)
    expect(plate2.accelAt(0.02, 0.01, 1.5e-6)).toBe(out[4])
  })

  it('switching accelMode mid-segment invalidates cached coefficients', () => {
    const plate = excitedPlate('hermite')
    const dt = 1 / (44 * 152)
    plate.beginStep(0, dt)
    const out = new Float64Array(5)
    plate.surfAt(0.02, 0.01, 1e-5, out)
    const wH = out[4]
    plate.accelMode = 'ode'
    plate.surfAt(0.02, 0.01, 1e-5, out)
    expect(out[4]).not.toBe(wH)
    expect(plate.accelAt(0.02, 0.01, 1e-5)).toBe(out[4])
    // Switching back restores the hermite definition consistently.
    plate.accelMode = 'hermite'
    plate.surfAt(0.02, 0.01, 1e-5, out)
    expect(out[4]).toBe(wH)
    expect(plate.accelAt(0.02, 0.01, 1e-5)).toBe(wH)
  })
})

describe('M6-11 ode accel mode (sand)', () => {
  const runMode = (enabled: boolean, tEnd = 0.3): TransientSandSim => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 60, 7, [0, 0], { enabled })
    sim.plate.accelMode = 'ode'
    sim.advanceTo(tEnd)
    return sim
  }

  it('sleep on/off reproduce the same physics in ode mode', () => {
    const off = runMode(false); const on = runMode(true)
    expect(on.nBounce).toBe(off.nBounce)
    expect(on.nDetach).toBe(off.nDetach)
    expect(on.nEventCap).toBe(off.nEventCap)
    let m = 0
    for (let i = 0; i < on.state.id.length; i++) {
      m = Math.max(m, Math.abs(on.state.x[i] - off.state.x[i]), Math.abs(on.state.y[i] - off.state.y[i]))
    }
    // Same tolerance as the established sleep-parity tests (counts must
    // match exactly; positions to <5e-6). Observed diff here: ~1e-13.
    expect(m).toBeLessThan(5e-6)
    expect([...on.state.sitting]).toEqual([...off.state.sitting])
  })

  it('ode mode detaches seated grains and stays finite', () => {
    const sim = runMode(true, 0.3)
    expect(sim.nDetach).toBeGreaterThan(0)
    expect(Number.isFinite(sim.state.x[0])).toBe(true)
    expect(sim.t).toBeCloseTo(0.3, 12)
    void G
  })
})
