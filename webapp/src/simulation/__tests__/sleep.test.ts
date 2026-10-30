// M6-2 gates: sleep on/off must reproduce the same physics (positions,
// counts), sleep must actually engage, input changes must invalidate it.
import { describe, it, expect } from 'vitest'
import { singleToneScore, type PerformanceScore } from '../sequence'
import { TransientSandSim, createTransientSim, calibratedForceN } from '../transientSand'
import { TransientPlate, hermiteQRange } from '../plateTransient'
import { initialState } from '../sand'

function run(enabled: boolean, freq = 152, tEnd = 0.6, n = 60, forceScale = 1): TransientSandSim {
  const sim = createTransientSim(singleToneScore(freq, calibratedForceN() * forceScale, 30), n, 7, [0, 0], { enabled })
  sim.advanceTo(tEnd)
  return sim
}

function maxPosDiff(a: TransientSandSim, b: TransientSandSim): number {
  let m = 0
  for (let i = 0; i < a.state.id.length; i++) {
    m = Math.max(m, Math.abs(a.state.x[i] - b.state.x[i]), Math.abs(a.state.y[i] - b.state.y[i]))
  }
  return m
}

describe('M6-2 sleep parity vs baseline', () => {
  it('strong drive 152Hz: same trajectories and counts, sleep engages', () => {
    const off = run(false)
    const on = run(true)
    expect(on.skippedUpdates).toBeGreaterThan(0)
    expect(Math.abs(on.t - off.t)).toBeLessThan(1e-12)
    expect(on.nBounce).toBe(off.nBounce)
    expect(on.nDetach).toBe(off.nDetach)
    expect(on.nAlive).toBe(off.nAlive)
    // Position tolerance: ballistic reassociation (chunked materialization)
    // differs at ~1e-12 per wake and is chaotically amplified through
    // collisions; event counts above are the primary gate.
    expect(maxPosDiff(on, off)).toBeLessThan(5e-6)
    expect([...on.state.sitting]).toEqual([...off.state.sitting])
  })

  it('weak drive: sitting sleep dominates and still matches', () => {
    const off = run(false, 152, 0.6, 60, 0.02)
    const on = run(true, 152, 0.6, 60, 0.02)
    expect(on.skippedUpdates).toBeGreaterThan(0)
    expect(maxPosDiff(on, off)).toBeLessThan(1e-6)
    expect(on.nAlive).toBe(off.nAlive)
  })

  it('frequency sweep samples match baseline', () => {
    for (const f of [70, 173.31, 368.6, 649]) {
      const off = run(false, f, 0.3, 40)
      const on = run(true, f, 0.3, 40)
      expect(on.nAlive, `f=${f}`).toBe(off.nAlive)
      // Higher frequencies amplify reassociation differences faster
      // (more collisions per second); counts remain the primary gate.
      expect(maxPosDiff(on, off), `f=${f}`).toBeLessThan(f >= 600 ? 1e-4 : 5e-6)
      expect(on.nBounce, `f=${f}`).toBe(off.nBounce)
    }
  })

  it('frequency change invalidates sleep and keeps state', () => {
    const score: PerformanceScore = singleToneScore(152, calibratedForceN() * 0.02, 30)
    const sim = createTransientSim(score, 40, 3, [0, 0], { enabled: true })
    sim.advanceTo(0.4)
    const sleepingBefore = sim.nSleeping
    const xBefore = Float64Array.from(sim.state.x)
    sim.wakeAll()
    sim.plate.setFrequency(649, sim.t)
    sim.refreshForceBounds()
    expect(sim.nSleeping).toBe(0)
    if (sleepingBefore > 0) {
      // Wake materialization must not teleport sitting particles in xy.
      let jump = 0
      for (let i = 0; i < 40; i++) {
        if (sim.state.sitting[i]) jump = Math.max(jump, Math.abs(sim.state.x[i] - xBefore[i]))
      }
      expect(jump).toBe(0)
    }
    sim.advanceTo(0.7)
    expect(Number.isFinite(sim.state.x[0])).toBe(true)
  })

  it('undriven plate: dropped particle settles identically', () => {
    const off = run(false, 152, 0.4, 10, 0)
    const on = run(true, 152, 0.4, 10, 0)
    expect(maxPosDiff(on, off)).toBeLessThan(1e-9)
    expect(on.nBounce).toBe(off.nBounce)
  })
})

describe('M6-2R regressions (Codex review)', () => {
  it('hermiteQRange reaches 4/27 (old 0.105 bound was unsafe)', () => {
    const [, hi] = hermiteQRange(0, 1, 0, 0, 1)
    expect(Math.abs(hi - 4 / 27)).toBeLessThan(1e-12)
    expect(hi).toBeGreaterThan(0.105)
  })

  it('residual 649Hz mode under 70Hz score: no missed detach (Codex counterexample)', () => {
    const mk = (enabled: boolean) => {
      const plate = new TransientPlate(singleToneScore(70, 0, 30))
      const st = initialState(plate.a, 0, 0, { x: [0], y: [0] })
      const sim = new TransientSandSim(plate, st, undefined, { enabled })
      const i649 = plate.modes.findIndex((m) => m.m === 0 && m.root_idx === 1)
      const j = plate.dofs.findIndex((d) => d.mode === i649)
      const w = plate.omega[j]
      plate.q[j] = 9.4 / (w * w)
      plate.qd[j] = 0
      // Seat at the center on the surface at t=0 (h = q * Phi(0), Phi(0)=R(0)).
      const S5 = new Float64Array(5)
      plate.beginStep(0, 1e-9)
      plate.surfAt(0, 0, 0, S5)
      st.z[0] = S5[0]; st.sitting[0] = 1
      sim.t = plate.t
      return sim
    }
    const off = mk(false); const on = mk(true)
    off.advanceTo(0.05); on.advanceTo(0.05)
    expect(off.nDetach).toBeGreaterThan(0) // baseline really detaches
    expect(on.nDetach).toBe(off.nDetach) // sleep must not miss it
    expect(Math.abs(on.state.z[0] - off.state.z[0])).toBeLessThan(1e-6)
  })

  it('event cap does not trigger in normal scenarios', () => {
    const on = run(true)
    expect(on.nEventCap).toBe(0)
  })

  it('fallback substeps never accept a below-surface endpoint (Codex re-review)', () => {
    // Unforced plate, one long step: many bounces + collapse inside the
    // step. The event-cap fallback must keep checking contacts over each
    // substep remainder instead of ballistically accepting z < 0.
    const sim = createTransientSim(singleToneScore(152, 0, 30), 1, 7, [0, 0], { enabled: false })
    const st = sim.state
    st.x[0] = 0; st.y[0] = 0; st.z[0] = 0.001
    st.vx[0] = 0; st.vy[0] = 0; st.vz[0] = -0.1; st.sitting[0] = 0
    sim.step(1)
    expect(sim.nEventCap).toBeGreaterThan(0)
    expect(st.z[0]).toBeGreaterThanOrEqual(-1e-12)
  })

  it('lost particles are excluded from the sleeping fraction (Codex re-review)', () => {
    const sim = createTransientSim(singleToneScore(152, 0, 30), 2, 7, [0, 0], { enabled: true })
    const st = sim.state; const a = sim.plate.a
    st.x[0] = a - 0.001; st.y[0] = 0; st.z[0] = 0.01
    st.vx[0] = 20; st.vy[0] = 0; st.vz[0] = 0; st.sitting[0] = 0
    st.x[1] = 0; st.y[1] = 0; st.z[1] = 0
    st.vx[1] = 0; st.vy[1] = 0; st.vz[1] = 0; st.sitting[1] = 1
    sim.step()
    expect(st.alive[0]).toBe(0)
    expect(sim.sleeping[0]).toBe(0)
    expect(sim.stats().sleepingFrac).toBeLessThanOrEqual(1)
  })
})
