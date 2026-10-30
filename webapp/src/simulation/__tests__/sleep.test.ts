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

describe('M6-3C differential frame transfer', () => {
  it('reconstructed display positions match snapshot at every frame', () => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN() * 0.05, 30), 60, 7, [0, 0], { enabled: true })
    const n = sim.state.id.length
    const cacheXy = new Float32Array(n * 2)
    const cacheAlive = new Uint8Array(n)
    let lastChanged = -1
    for (let f = 0; f < 12; f++) {
      sim.advanceTo(f / 30)
      const frame = sim.frameDelta()
      if (f === 0) expect(frame.full).toBe(true)
      for (let k = 0; k < frame.idx.length; k++) {
        cacheXy[2 * frame.idx[k]] = frame.xy[2 * k]
        cacheXy[2 * frame.idx[k] + 1] = frame.xy[2 * k + 1]
        cacheAlive[frame.idx[k]] = 1
      }
      for (const dIdx of frame.dead) cacheAlive[dIdx] = 0
      lastChanged = frame.idx.length + frame.dead.length
      // Reference: full snapshot at the same physics time.
      const snap = sim.snapshot()
      expect(frame.t).toBe(snap.t)
      expect(frame.aliveCount).toBe(snap.aliveCount)
      let k2 = 0
      for (let i = 0; i < n; i++) {
        if (!sim.state.alive[i]) { expect(cacheAlive[i]).toBe(0); continue }
        expect(cacheAlive[i]).toBe(1)
        expect(cacheXy[2 * i]).toBe(snap.xy[k2])
        expect(cacheXy[2 * i + 1]).toBe(snap.xy[k2 + 1])
        k2 += 2
      }
    }
    expect(lastChanged).toBeGreaterThanOrEqual(0)
  })
})

describe('M6-3B region-shared sitting bounds', () => {
  it('cellPhiRanges contains phiAll at every sampled point', () => {
    const plate = new TransientPlate(singleToneScore(152, 0.05, 30))
    const nd = plate.ndof
    const lo = new Float64Array(nd); const hi = new Float64Array(nd); const v = new Float64Array(nd)
    let seed = 12345
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const a = plate.a
    const cells: [number, number, number, number][] = [
      [-a, -a / 2, -a, -a / 2], [0, a / 3, 0, a / 3], [-0.01, 0.01, -0.01, 0.01],
      [a * 0.7, a, -a * 0.2, 0], [-a, 0, a * 0.5, a], [0.001, 0.01, -0.05, -0.04],
    ]
    for (const [x0, x1, y0, y1] of cells) {
      plate.cellPhiRanges(x0, x1, y0, y1, lo, hi)
      for (let k = 0; k < 2000; k++) {
        const x = x0 + rnd() * (x1 - x0); const y = y0 + rnd() * (y1 - y0)
        plate.phiAll(x, y, v)
        for (let j = 0; j < nd; j++) {
          expect(v[j]).toBeGreaterThanOrEqual(lo[j] - 1e-12)
          expect(v[j]).toBeLessThanOrEqual(hi[j] + 1e-12)
        }
      }
    }
  })

  it('regions on/off reproduce the same physics and regions engage', () => {
    const mk = (regions: boolean) => createTransientSim(
      singleToneScore(152, calibratedForceN() * 0.05, 30), 60, 7, [0, 0],
      { enabled: true, scheduler: true, regions })
    const off = mk(false); const on = mk(true)
    off.advanceTo(0.6); on.advanceTo(0.6)
    expect(on.regionChecks).toBeGreaterThan(0)
    expect(on.nBounce).toBe(off.nBounce)
    expect(on.nDetach).toBe(off.nDetach)
    expect(on.nAlive).toBe(off.nAlive)
    expect(maxPosDiff(on, off)).toBeLessThan(5e-6)
  })
})

describe('M6-3A scheduler equivalence vs legacy full scan', () => {
  it('same counts, same trajectories, same skipped count', () => {
    const legacy = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 60, 7, [0, 0], { enabled: true, scheduler: false })
    const sched = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 60, 7, [0, 0], { enabled: true, scheduler: true })
    legacy.advanceTo(0.6); sched.advanceTo(0.6)
    expect(sched.skippedUpdates).toBeGreaterThan(0)
    expect(sched.skippedUpdates).toBe(legacy.skippedUpdates)
    expect(sched.nBounce).toBe(legacy.nBounce)
    expect(sched.nDetach).toBe(legacy.nDetach)
    expect(sched.nEventCap).toBe(legacy.nEventCap)
    expect(sched.wakeCount).toBe(legacy.wakeCount)
    expect(sched.nAlive).toBe(legacy.nAlive)
    // Per-particle math is independent of processing order, so results
    // agree to floating-point identity level, not just tolerance level.
    expect(maxPosDiff(sched, legacy)).toBeLessThan(1e-9)
    expect([...sched.state.sitting]).toEqual([...legacy.state.sitting])
  })

  it('scheduler survives input invalidation and edge loss', () => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN() * 0.05, 30), 60, 7, [0, 0], { enabled: true, scheduler: true })
    sim.advanceTo(0.3)
    sim.wakeAll()
    sim.plate.setFrequency(368.6, sim.t)
    sim.advanceTo(0.5)
    expect(Math.abs(sim.t - 0.5)).toBeLessThan(1e-12)
    expect(sim.stats().sleepingFrac).toBeLessThanOrEqual(1)
    for (let i = 0; i < sim.state.id.length; i++) {
      expect(Number.isFinite(sim.state.x[i])).toBe(true)
      expect(Number.isFinite(sim.state.z[i])).toBe(true)
    }
  })
})

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
