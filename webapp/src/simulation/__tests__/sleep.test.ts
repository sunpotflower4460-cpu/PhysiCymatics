// M6-2 gates: sleep on/off must reproduce the same physics (positions,
// counts), sleep must actually engage, input changes must invalidate it.
import { describe, it, expect } from 'vitest'
import { singleToneScore, type PerformanceScore } from '../sequence'
import { TransientSandSim, createTransientSim, createTransientSimShard, shardRanges, calibratedForceN } from '../transientSand'
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

describe('M6-4B shard equivalence (in-process replicas)', () => {
  it('K shards stepped in lockstep reproduce the single sim bitwise', () => {
    const score = singleToneScore(152, calibratedForceN(), 30)
    const N = 90
    for (const K of [2, 3, 4]) {
      const whole = createTransientSim(score, N, 7, [0, 0], { enabled: true })
      const shards = Array.from({ length: K }, (_, k) =>
        createTransientSimShard(score, N, 7, [0, 0], { enabled: true }, k, K))
      // Same ensemble: merged initial positions equal the single sim's.
      const ranges = shardRanges(N, K)
      shards.forEach((s, k) => {
        const [a] = ranges[k]
        for (let i = 0; i < s.state.id.length; i++) {
          expect(s.state.id[i]).toBe(a + i)
          expect(s.state.x[i]).toBe(whole.state.x[a + i])
          expect(s.state.y[i]).toBe(whole.state.y[a + i])
        }
      })
      const STEPS = 1500
      for (let s = 0; s < STEPS; s++) {
        whole.step()
        for (const sh of shards) sh.step()
        if (s === 749) { // mid-run input change at the same step boundary
          whole.wakeAll(); whole.plate.setFrequency(649, whole.t)
          for (const sh of shards) { sh.wakeAll(); sh.plate.setFrequency(649, sh.t) }
        }
      }
      whole.materializeAll()
      let bSum = 0; let dSum = 0
      shards.forEach((s, k) => {
        s.materializeAll()
        bSum += s.nBounce; dSum += s.nDetach
        const [a] = ranges[k]
        for (let i = 0; i < s.state.id.length; i++) {
          const g = a + i
          expect(s.state.x[i]).toBe(whole.state.x[g])
          expect(s.state.y[i]).toBe(whole.state.y[g])
          expect(s.state.z[i]).toBe(whole.state.z[g])
          expect(s.state.vz[i]).toBe(whole.state.vz[g])
          expect(s.state.sitting[i]).toBe(whole.state.sitting[g])
          expect(s.state.alive[i]).toBe(whole.state.alive[g])
        }
      })
      expect(bSum).toBe(whole.nBounce)
      expect(dSum).toBe(whole.nDetach)
    }
  })
})

describe('M6-4A surfAt optimization equivalence', () => {
  it('optimized surfAt is bitwise identical to the reference', () => {
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30))
    const o1 = new Float64Array(5); const o2 = new Float64Array(5)
    let seed = 4242
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const check = (x: number, y: number, tau: number) => {
      plate.surfAt(x, y, tau, o1)
      plate.surfAtReference(x, y, tau, o2)
      for (let k = 0; k < 5; k++) expect(o1[k]).toBe(o2[k])
    }
    let t = 0
    const dt = 1 / (44 * 152)
    for (let step = 0; step < 60; step++) {
      plate.beginStep(t, t + dt)
      // Deterministic boundary points + random points, shared & unique taus.
      check(0, 0, t); check(plate.a, 0, t + dt); check(-plate.a * 0.99, 0.001, t + dt / 2)
      check(0.03, -0.04, t + dt * 0.37)
      for (let k = 0; k < 40; k++) {
        check((rnd() * 2 - 1) * plate.a, (rnd() * 2 - 1) * plate.a, t + rnd() * dt)
      }
      // Repeated same-tau calls exercise the coefficient cache path.
      for (let k = 0; k < 10; k++) check(rnd() * 0.1 - 0.05, rnd() * 0.1 - 0.05, t + dt)
      t += dt
      if (step === 29) plate.setFrequency(649, t)
    }
  })

  it('full-sim results unchanged with optimized surfAt (counts + positions)', () => {
    // The optimized path is the only surfAt now; equivalence to the
    // reference above plus the existing parity gates cover this. Here we
    // assert the standard strong-drive counts stay at their pinned values
    // from the M6-3R benchmark (500 grains, 1s, sleep OFF).
    const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 500, 7, [0, 0], { enabled: false })
    sim.advanceTo(1)
    expect(sim.nBounce).toBe(19357)
    expect(sim.nDetach).toBe(10947)
  })
})

describe('M6-3R audit: scheduler invariants', () => {
  const internals = (sim: TransientSandSim) => sim as unknown as {
    activeList: Int32Array; activePos: Int32Array; activeCount: number; heapSize: number
  }

  function checkInvariant(sim: TransientSandSim): void {
    const s = internals(sim); const st = sim.state
    const seen = new Set<number>()
    for (let k = 0; k < s.activeCount; k++) {
      const i = s.activeList[k]
      expect(seen.has(i)).toBe(false) // no double registration
      seen.add(i)
      expect(s.activePos[i]).toBe(k)
      expect(st.alive[i]).toBe(1) // no lost particle in the active set
      expect(sim.sleeping[i]).toBe(0) // no sleeper in the active set
    }
    let awakeAlive = 0
    for (let i = 0; i < st.id.length; i++) {
      if (st.alive[i] && !sim.sleeping[i]) { awakeAlive++; expect(seen.has(i)).toBe(true) }
      if (sim.sleeping[i]) expect(sim.wakeT[i]).toBeGreaterThan(sim.t - 1e-12) // never overdue
    }
    expect(s.activeCount).toBe(awakeAlive) // nobody skipped by swap-remove
  }

  it('active-set invariant holds at every step incl. invalidation', () => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN() * 0.3, 30), 80, 7, [0, 0], { enabled: true, scheduler: true })
    for (let k = 0; k < 400; k++) { sim.step(); checkInvariant(sim) }
    sim.wakeAll(); sim.plate.setFrequency(649, sim.t); checkInvariant(sim)
    for (let k = 0; k < 400; k++) { sim.step(); checkInvariant(sim) }
  })

  it('stale heap entries after wakeAll never wake anything', () => {
    const mk = (scheduler: boolean) => createTransientSim(
      singleToneScore(152, calibratedForceN() * 0.05, 30), 60, 7, [0, 0], { enabled: true, scheduler })
    const legacy = mk(false); const sched = mk(true)
    for (const s of [legacy, sched]) {
      s.advanceTo(0.25)
      s.wakeAll(); s.plate.setFrequency(368.6, s.t)
      s.advanceTo(0.45)
      s.wakeAll(); s.plate.setFrequency(152, s.t)
      s.advanceTo(0.6)
    }
    expect(sched.nBounce).toBe(legacy.nBounce)
    expect(sched.nDetach).toBe(legacy.nDetach)
    expect(sched.wakeCount).toBe(legacy.wakeCount)
    expect(maxPosDiff(sched, legacy)).toBeLessThan(1e-9)
  })

  it('lost particles never reactivate and leave the active set', () => {
    const sim = createTransientSim(singleToneScore(152, 0, 30), 3, 7, [0, 0], { enabled: true, scheduler: true })
    const st = sim.state; const a = sim.plate.a
    st.x[0] = a - 0.0005; st.y[0] = 0; st.z[0] = 0.01; st.vx[0] = 20; st.vy[0] = 0; st.vz[0] = 5; st.sitting[0] = 0
    for (let k = 0; k < 50; k++) sim.step()
    expect(st.alive[0]).toBe(0)
    expect(sim.sleeping[0]).toBe(0)
    expect(internals(sim).activePos[0]).toBe(-1)
    for (let k = 0; k < 50; k++) sim.step()
    expect(st.alive[0]).toBe(0)
    checkInvariant(sim)
  })
})

describe('M6-3R audit: region bounds at boundaries', () => {
  it('deterministic grid containment incl. origin/angular-wrap/bin-edge cells', () => {
    const plate = new TransientPlate(singleToneScore(152, 0.05, 30))
    const nd = plate.ndof
    const lo = new Float64Array(nd); const hi = new Float64Array(nd); const v = new Float64Array(nd)
    const a = plate.a
    const ngrid = plate.modes[0].R.length
    const binW = a / (ngrid - 1)
    const rEdge = 40 * binW // a radial-bin boundary
    const cells: [number, number, number, number][] = [
      [-0.002, 0.002, -0.002, 0.002], // contains origin
      [-0.11, -0.09, -0.003, 0.003], // angular wrap at theta = +/-pi
      [rEdge - 0.0015, rEdge + 0.0015, -0.001, 0.001], // radial bin edge
      [a - 0.004, a, -0.002, 0.002], // outer rim
      [-a, -a + 0.004, a - 0.004, a], // corner outside disk partially
      [0.02, 0.05, 0.02, 0.05], // generic interior
      [-0.05, -0.02, 0.001, 0.004], // near negative x-axis
    ]
    for (const [x0, x1, y0, y1] of cells) {
      plate.cellPhiRanges(x0, x1, y0, y1, lo, hi)
      for (let kx = 0; kx <= 20; kx++) {
        for (let ky = 0; ky <= 20; ky++) {
          const x = x0 + (kx / 20) * (x1 - x0); const y = y0 + (ky / 20) * (y1 - y0)
          plate.phiAll(x, y, v)
          for (let j = 0; j < nd; j++) {
            expect(v[j]).toBeGreaterThanOrEqual(lo[j] - 1e-12)
            expect(v[j]).toBeLessThanOrEqual(hi[j] + 1e-12)
          }
        }
      }
    }
  })

  it('every grid cell of the actual 64x64 layout contains its sample points', () => {
    const plate = new TransientPlate(singleToneScore(152, 0.05, 30))
    const nd = plate.ndof
    const lo = new Float64Array(nd); const hi = new Float64Array(nd); const v = new Float64Array(nd)
    const a = plate.a; const G = 64; const w = (2 * a) / G
    for (let gy = 0; gy < G; gy++) {
      for (let gx = 0; gx < G; gx++) {
        const x0 = -a + gx * w; const y0 = -a + gy * w
        plate.cellPhiRanges(x0, x0 + w, y0, y0 + w, lo, hi)
        for (const fx of [0, 0.25, 0.5, 0.75, 1]) {
          for (const fy of [0, 0.25, 0.5, 0.75, 1]) {
            plate.phiAll(x0 + fx * w, y0 + fy * w, v)
            for (let j = 0; j < nd; j++) {
              if (v[j] < lo[j] - 1e-12 || v[j] > hi[j] + 1e-12) {
                throw new Error(`cell ${gx},${gy} dof ${j}: ${v[j]} outside [${lo[j]}, ${hi[j]}]`)
              }
            }
          }
        }
      }
    }
  })

  it('cell lower bound <= exact lower bound <= cell upper bound', () => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 10, 7, [0, 0], { enabled: true, regions: true })
    sim.step()
    const s = sim as unknown as {
      lookCache: { aMin: Float64Array; aMax: Float64Array }
      cellLb: Float64Array; cellLbHi: Float64Array
      ensureLookahead: (t: number) => void; ensureCellLb: () => void; cellIndex: (x: number, y: number) => number
    }
    s.ensureLookahead(sim.t); s.ensureCellLb()
    const lc = s.lookCache
    const phi = new Float64Array(sim.plate.ndof)
    let seed = 999
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    for (let k = 0; k < 3000; k++) {
      const x = (rnd() * 2 - 1) * sim.plate.a; const y = (rnd() * 2 - 1) * sim.plate.a
      sim.plate.phiAll(x, y, phi)
      let exact = 0
      for (let j = 0; j < phi.length; j++) exact += Math.min(phi[j] * lc.aMin[j], phi[j] * lc.aMax[j])
      const c = s.cellIndex(x, y)
      expect(s.cellLb[c]).toBeLessThanOrEqual(exact + 1e-9)
      expect(exact).toBeLessThanOrEqual(s.cellLbHi[c] + 1e-9)
    }
  })
})

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

  it('reset resends a full frame; deaths are reported exactly once', () => {
    const sim = createTransientSim(singleToneScore(152, 0, 30), 3, 7, [0, 0], { enabled: true })
    const st = sim.state; const a = sim.plate.a
    const f0 = sim.frameDelta()
    expect(f0.full).toBe(true)
    expect(f0.idx.length).toBe(3)
    st.x[0] = a - 0.0005; st.y[0] = 0; st.z[0] = 0.01; st.vx[0] = 20; st.vz[0] = 5; st.sitting[0] = 0
    let deadSeen = 0
    for (let k = 0; k < 30; k++) {
      sim.step()
      const f = sim.frameDelta()
      deadSeen += [...f.dead].filter((d) => d === 0).length
      // A dead slot is never sent as a position update afterwards.
      if (st.alive[0] === 0) expect([...f.idx]).not.toContain(0)
    }
    expect(deadSeen).toBe(1)
    sim.resetFrameDelta()
    const fr = sim.frameDelta()
    expect(fr.full).toBe(true)
    expect(fr.idx.length).toBe(sim.nAlive)
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
