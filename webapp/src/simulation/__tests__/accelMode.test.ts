// M6-11 research accel mode ('ode'): cross-path consistency, boundary
// continuity, analytic accuracy, lookahead-range validity, and sleep
// parity under the new acceleration definition. Default 'hermite'
// behavior is covered by the existing suites (bitwise unchanged).
import { describe, it, expect } from 'vitest'
import { singleToneScore, type PerformanceScore } from '../sequence'
import { TransientPlate, hermiteQdRange } from '../plateTransient'
import { TransientSandSim, createTransientSim, calibratedForceN } from '../transientSand'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'
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

describe('M6-11R degenerate tiny final step (production bugfix)', () => {
  it('advanceTo on the step grid matches pure step() driving (500Hz)', () => {
    // Regression: advanceTo(0.1) at 500Hz used to create a ~1e-15 s
    // final segment whose Hermite acceleration is cancellation
    // garbage, spuriously detaching ~51 of 60 seated grains.
    const mkSim = () => {
      const plate = new TransientPlate(singleToneScore(500, calibratedForceN(), 30), [0, 0])
      return new TransientSandSim(plate, initialState(plate.a, 60, 12345), { ...DEFAULT_PARAMS, spp: 44 }, { enabled: true })
    }
    const stepped = mkSim()
    const nSteps = Math.round(0.1 / stepped.dtNow())
    for (let k = 0; k < nSteps; k++) stepped.step()
    const adv = mkSim()
    adv.advanceTo(0.1)
    expect(adv.t).toBe(0.1)
    expect(adv.nDetach).toBe(stepped.nDetach)
    expect(adv.nBounce).toBe(stepped.nBounce)
  })

  it('M6-13R: consecutive tiny advanceTo targets keep plate/sim clocks synchronized', () => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 5, 7)
    sim.advanceTo(0.75e-12)
    expect(sim.t).toBe(0.75e-12)
    expect(sim.plate.t).toBe(sim.t)
    sim.advanceTo(1.5e-12)
    expect(sim.plate.t).toBe(sim.t)
    sim.advanceTo(1e-12) // reversed target: no-op
    expect(sim.t).toBe(1.5e-12)
    expect(sim.plate.t).toBe(sim.t)
    sim.step()
    expect(sim.plate.t).toBe(sim.t)
    for (const z of sim.state.z) expect(Number.isFinite(z)).toBe(true)
  })

  it('tiny segments evaluate sane acceleration (no 1/h^2 garbage)', () => {
    const plate = excitedPlate('hermite')
    plate.beginStep(0, 1e-12)
    const out = new Float64Array(5)
    plate.surfAt(0.02, 0.01, 5e-13, out)
    expect(Number.isFinite(out[4])).toBe(true)
    expect(Math.abs(out[4])).toBeLessThan(1e5)
    expect(plate.accelAt(0.02, 0.01, 5e-13)).toBe(out[4])
  })
})

describe('M6-11R interval-assisted detach search', () => {
  it('wttLowerBound bounds dense accelAt samples (ode mode)', () => {
    const plate = excitedPlate('ode', 152, calibratedForceN())
    const dt = 1 / (44 * 152)
    plate.beginStep(0, dt)
    let seed = 77
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    for (let k = 0; k < 60; k++) {
      const x = (rnd() * 2 - 1) * plate.a * 0.9
      const y = (rnd() * 2 - 1) * plate.a * 0.9
      const ta = rnd() * dt * 0.8
      const tb = ta + (0.2 + 0.8 * rnd()) * (dt - ta)
      const lb = plate.wttLowerBound(x, y, ta, tb)
      for (let sN = 0; sN <= 40; sN++) {
        const w = plate.accelAt(x, y, ta + ((tb - ta) * sN) / 40)
        expect(w).toBeGreaterThanOrEqual(lb - 1e-6)
      }
    }
  })

  it('firstDetachBelow matches a dense-scan first root (ode mode)', () => {
    const plate = excitedPlate('ode', 152, calibratedForceN())
    const dt = 1 / (44 * 152)
    const denseFirst = (x: number, y: number, t0: number): number | null => {
      let tp = t0; let fp = plate.accelAt(x, y, t0) + 9.81
      if (fp < 0) return t0
      for (let k = 1; k <= 2000; k++) {
        const t = t0 + (dt * k) / 2000
        const ft = plate.accelAt(x, y, t) + 9.81
        if (ft < 0) {
          let lo = tp; let hi = t
          for (let it = 0; it < 50; it++) {
            const m = 0.5 * (lo + hi)
            if (plate.accelAt(x, y, m) + 9.81 >= 0) lo = m; else hi = m
          }
          return 0.5 * (lo + hi)
        }
        tp = t; fp = ft
      }
      return null
    }
    let segT = 0; let rootsFound = 0
    for (let s = 0; s < 30 && rootsFound < 5; s++) {
      plate.beginStep(segT, segT + dt)
      for (const [x, y] of [[0.02, 0.01], [-0.04, 0.03], [0.0, 0.0]] as const) {
        const got = plate.firstDetachBelow(x, y, segT, segT + dt)
        const want = denseFirst(x, y, segT)
        expect(got.kind).not.toBe('unknown')
        if (want === null) expect(got.kind).toBe('none')
        else {
          rootsFound++
          expect(got.kind).toBe('found')
          if (got.kind === 'found') expect(Math.abs(got.t - want)).toBeLessThan(2e-6)
        }
      }
      segT += dt
    }
    expect(rootsFound).toBeGreaterThan(0)
  })

  // M6-12R review counterexamples (PR #14): envelope peak, tiny-segment
  // mask, first-root ordering. Each is a from-rest plate so the only
  // acceleration is the direct forcing term.
  const evScore = (events: PerformanceScore['events'], phase: number): PerformanceScore => ({
    schema: 'physicymatics.score.v1', duration_s: 30, initialPhase_rad: phase, audioGain: 0.2, events,
  })

  it('M6-12R: wttLowerBound includes the overlapping attack/release envelope peak', () => {
    // Event 1-2 us, F=-0.445 N, attack 0.8 us, release 0.7 us (ramps
    // overlap; envelope peaks at their crossing, below 1).
    const score = evScore([{
      start_s: 1e-6, end_s: 2e-6, frequency_hz: 152, forcePeak_N: -0.445,
      attack_s: 0.8e-6, release_s: 0.7e-6,
    }], Math.PI / 2)
    const plate = new TransientPlate(score, [0, 0])
    plate.accelMode = 'ode'
    const dt = 1 / (44 * 152)
    plate.beginStep(0, dt)
    const ta = 1.5e-6; const tb = 1.5625e-6
    const lb = plate.wttLowerBound(0.02, 0.01, ta, tb)
    let trueMin = Infinity
    for (let k = 0; k <= 400; k++) {
      trueMin = Math.min(trueMin, plate.accelAt(0.02, 0.01, ta + ((tb - ta) * k) / 400))
    }
    expect(trueMin).toBeLessThan(-G) // the counterexample really crosses
    expect(lb).toBeLessThanOrEqual(trueMin + 1e-6)
    const r = plate.firstDetachBelow(0.02, 0.01, 0, dt)
    expect(r.kind).toBe('found')
    if (r.kind === 'found') {
      expect(r.t).toBeGreaterThan(1.4e-6)
      expect(r.t).toBeLessThan(1.7e-6)
    }
  })

  it('M6-12R: tiny segments keep directly forced DOFs live in hermite mode', () => {
    // 10 N event spanning 10-20 ns inside a 90 ns segment: Hermite
    // endpoints stay zero, but the tiny-segment ODE-form acceleration
    // must see the forcing in accelAt as well as surfAt.
    const score = evScore([{
      start_s: 10e-9, end_s: 20e-9, frequency_hz: 152, forcePeak_N: 10,
      attack_s: 0, release_s: 0,
    }], Math.PI / 2)
    const plate = new TransientPlate(score, [0, 0])
    expect(plate.accelMode).toBe('hermite')
    plate.beginStep(0, 9e-8)
    expect(plate.odeFormActive()).toBe(true)
    const out = new Float64Array(5)
    plate.surfAt(0.02, 0.01, 15e-9, out)
    const a = plate.accelAt(0.02, 0.01, 15e-9)
    expect(Math.abs(out[4])).toBeGreaterThan(100)
    expect(a).toBe(out[4])
    const out2 = new Float64Array(2)
    plate.heightAccelAt(0.02, 0.01, 15e-9, out2)
    expect(out2[1]).toBe(out[4])
  })

  it('M6-12R: firstDetachBelow returns the FIRST root across two short events', () => {
    // Two 10 N events (1-2 us and 18-19 us) inside one nominal step.
    // The first episode sits between scan points; the old fast path
    // returned the second event's crossing (18 us).
    const score = evScore([
      { start_s: 1e-6, end_s: 2e-6, frequency_hz: 152, forcePeak_N: 10, attack_s: 0, release_s: 0 },
      { start_s: 18e-6, end_s: 19e-6, frequency_hz: 152, forcePeak_N: 10, attack_s: 0, release_s: 0 },
    ], -Math.PI / 2)
    const plate = new TransientPlate(score, [0, 0])
    plate.accelMode = 'ode'
    const dt = 1 / (44 * 152)
    plate.beginStep(0, dt)
    expect(plate.accelAt(0.02, 0.01, 1.5e-6)).toBeLessThan(-G)
    const r = plate.firstDetachBelow(0.02, 0.01, 0, dt)
    expect(r.kind).toBe('found')
    if (r.kind === 'found') {
      expect(r.t).toBeGreaterThanOrEqual(1e-6)
      expect(r.t).toBeLessThan(2e-6)
    }
    // 'unknown' is distinguishable and only appears when the budget
    // is deliberately starved; it is never reported as 'none'.
    const starved = plate.firstDetachBelow(0.02, 0.01, 0, dt, 0)
    expect(starved.kind).toBe('unknown')
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

describe('M6-14R setDrive cache invalidation', () => {
  it('newly forced DOFs enter ode acceleration after drive move', () => {
    // Zero state, center drive: only m=0 DOFs forced. Move drive
    // eccentrically mid-segment: m>0 DOFs become directly forced and
    // must appear in accelAt/surfAt identically (stale live mask
    // dropped them before the fix).
    const plate = new TransientPlate(singleToneScore(152, 1.0, 30), [0, 0])
    plate.accelMode = 'ode'
    plate.beginStep(0, 0.002)
    const before = plate.accelAt(0.05, 0.02, 0.0011)
    plate.setDrive(0.06, 0.0)
    const after = plate.accelAt(0.05, 0.02, 0.0011)
    const out = new Float64Array(5)
    plate.surfAt(0.05, 0.02, 0.0011, out)
    expect(after).toBe(out[4])
    expect(after).not.toBe(before)
    expect(Math.abs(after)).toBeGreaterThan(0)
    // All evaluation paths on the same plate agree after the move.
    const out2 = new Float64Array(2)
    plate.heightAccelAt(0.05, 0.02, 0.0011, out2)
    expect(out2[1]).toBe(after)
  })
})
