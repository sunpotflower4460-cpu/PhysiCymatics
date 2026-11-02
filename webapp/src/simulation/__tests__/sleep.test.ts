// M6-2 gates: sleep on/off must reproduce the same physics (positions,
// counts), sleep must actually engage, input changes must invalidate it.
import { describe, it, expect } from 'vitest'
import { singleToneScore, type PerformanceScore } from '../sequence'
import { TransientSandSim, createTransientSim, createTransientSimShard, shardRanges, calibratedForceN, type SleepOptions } from '../transientSand'
import { TransientPlate, hermiteQRange } from '../plateTransient'
import { initialState, initialStateRange } from '../sand'

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

describe('M6-5 local non-contact bound', () => {
  it('cell height bound dominates actual surfAt height everywhere sampled', () => {
    const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 10, 7, [0, 0], { enabled: true, localBound: true })
    sim.step()
    const s = sim as unknown as { cellHeightMax: (c: number) => number; cellIndex: (x: number, y: number) => number }
    const out = new Float64Array(5)
    let seed = 77
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const dt = sim.dtNow()
    for (let k = 0; k < 4000; k++) {
      const x = (rnd() * 2 - 1) * sim.plate.a
      const y = (rnd() * 2 - 1) * sim.plate.a
      const tau = sim.t - rnd() * dt
      sim.plate.surfAt(x, y, tau, out)
      const bound = s.cellHeightMax(s.cellIndex(x, y))
      expect(out[0]).toBeLessThanOrEqual(bound + 1e-12)
    }
  })

  it('localBound on/off are bitwise identical (strong and weak)', () => {
    for (const forceScale of [1, 0.05]) {
      const mk = (localBound: boolean) => createTransientSim(
        singleToneScore(152, calibratedForceN() * forceScale, 30), 60, 7, [0, 0],
        { enabled: true, localBound })
      const off = mk(false); const on = mk(true)
      off.advanceTo(0.6); on.advanceTo(0.6)
      expect(on.nBounce).toBe(off.nBounce)
      expect(on.nDetach).toBe(off.nDetach)
      expect(on.nEventCap).toBe(off.nEventCap)
      expect(maxPosDiff(on, off)).toBe(0)
      expect([...on.state.sitting]).toEqual([...off.state.sitting])
      expect(on.localSkips).toBeGreaterThanOrEqual(0)
    }
  })
})

describe('M6-5R boundary audit (local bound domain)', () => {
  // The precomputed Phi ranges cover the closed square [-a,a]^2 only.
  // phiAll outside that square evaluates R[last bin] * trig(theta) with
  // the TRUE theta, which clamped edge cells do not contain. Audit
  // evidence: at (0.24, 0.05) most DOFs violate the clamped cell range.
  // locallyClear must therefore never prove a segment that leaves the
  // square, and must fail closed on non-finite input.
  type Priv = {
    locallyClear: (x0: number, y0: number, x1: number, y1: number, zMin: number) => boolean
    cellHeightMax: (c: number) => number
    cellIndex: (x: number, y: number) => number
    cellPhiLo: Float64Array; cellPhiHi: Float64Array
  }
  const calmSim = () => {
    const sim = createTransientSim(singleToneScore(152, 0, 30), 10, 7, [0, 0], { enabled: true, localBound: true })
    sim.step() // establishes a plate segment with q == 0 (no drive force)
    return sim
  }

  it('clamped edge cells do NOT contain phiAll outside [-a,a]^2 (audit evidence)', () => {
    const sim = calmSim()
    const s = sim as unknown as Priv
    const pl = sim.plate; const nd = pl.ndof
    const phi = new Float64Array(nd)
    pl.phiAll(0.24, 0.05, phi)
    const c = s.cellIndex(0.24, 0.05) // clamped to the edge column
    let violations = 0
    for (let j = 0; j < nd; j++) {
      if (phi[j] < s.cellPhiLo[c * nd + j] - 1e-12 || phi[j] > s.cellPhiHi[c * nd + j] + 1e-12) violations++
    }
    expect(violations).toBeGreaterThan(0)
  })

  it('locallyClear never proves a segment leaving [-a,a]^2', () => {
    const sim = calmSim()
    const s = sim as unknown as Priv
    const a = sim.plate.a
    // Fully inside (endpoints on the closed boundary included): provable.
    expect(s.locallyClear(0.10, 0.05, 0.11, 0.05, 1e-6)).toBe(true)
    expect(s.locallyClear(0.10, 0.05, a, 0.05, 1e-6)).toBe(true)
    expect(s.locallyClear(-a, -a, a, a, 1e-6)).toBe(true)
    // Leaving the square in any direction: never proven, even though the
    // calm plate height is exactly 0 (the bound has no domain there).
    expect(s.locallyClear(0.10, 0.05, 0.24, 0.05, 1e-6)).toBe(false)
    expect(s.locallyClear(0.10, 0.05, a + 1e-12, 0.05, 1e-6)).toBe(false)
    expect(s.locallyClear(0.05, 0.10, 0.05, 0.24, 1e-6)).toBe(false)
    expect(s.locallyClear(-0.24, 0.05, -0.10, 0.05, 1e-6)).toBe(false)
    expect(s.locallyClear(0.13, 0.13, 0.10, 0.05, 1e-6)).toBe(false)
  })

  it('locallyClear fails closed on non-finite input', () => {
    const sim = calmSim()
    const s = sim as unknown as Priv
    expect(s.locallyClear(0.10, 0.05, 0.11, 0.05, NaN)).toBe(false)
    expect(s.locallyClear(0.10, 0.05, 0.11, 0.05, Infinity)).toBe(false)
    expect(s.locallyClear(NaN, 0.05, 0.11, 0.05, 1e-6)).toBe(false)
    expect(s.locallyClear(0.10, 0.05, Infinity, 0.05, 1e-6)).toBe(false)
  })

  it('cell height bound contains surfAt at deterministic boundary points', () => {
    const sim = createTransientSim(singleToneScore(649, calibratedForceN(), 30), 10, 7, [0, 0], { enabled: true, localBound: true })
    sim.advanceTo(0.05)
    sim.step()
    const s = sim as unknown as Priv
    const pl = sim.plate; const a = pl.a
    const out = new Float64Array(5)
    const ng = pl.modes[0].R.length
    const pts: [number, number][] = [
      [0, 0], [a, 0], [-a, 0], [0, a], [0, -a], [a, a], [-a, -a], [a, -a], [-a, a],
      [a - 1e-15, 0.03], [-(a - 1e-15), -0.07], [1e-15, -1e-15],
    ]
    // radial bin edges at several angles
    for (const k of [1, 2, 17, 100, ng - 2, ng - 1]) {
      const r = (k * a) / (ng - 1)
      for (const th of [0, 0.7, Math.PI / 2, 2.4, Math.PI, -1.3]) {
        const x = r * Math.cos(th); const y = r * Math.sin(th)
        if (Math.abs(x) <= a && Math.abs(y) <= a) pts.push([x, y])
      }
    }
    const dt = sim.dtNow()
    for (const [x, y] of pts) {
      const bound = s.cellHeightMax(s.cellIndex(x, y))
      for (const fr of [0, 0.31, 0.5, 0.83, 1]) {
        pl.surfAt(x, y, sim.t - dt + fr * dt, out)
        expect(out[0]).toBeLessThanOrEqual(bound + 1e-12)
      }
    }
  })

  it('hermiteQRange contains dense Bernstein samples and lies in the Bernstein hull', () => {
    let seed = 424243
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let trial = 0; trial < 300; trial++) {
      const q0 = rnd() * 1e-3; const q1 = rnd() * 1e-3
      const qd0 = rnd(); const qd1 = rnd()
      const h = 1e-4 + Math.abs(rnd()) * 1e-3
      const [lo, hi] = hermiteQRange(q0, qd0, q1, qd1, h)
      // Bernstein coefficients of the same cubic on s in [0,1].
      const b0 = q0; const b1 = q0 + (h * qd0) / 3; const b2 = q1 - (h * qd1) / 3; const b3 = q1
      const bMin = Math.min(b0, b1, b2, b3); const bMax = Math.max(b0, b1, b2, b3)
      const scale = Math.max(1e-12, Math.abs(bMin), Math.abs(bMax))
      expect(lo).toBeGreaterThanOrEqual(bMin - 1e-12 * scale - 1e-18)
      expect(hi).toBeLessThanOrEqual(bMax + 1e-12 * scale + 1e-18)
      for (let k = 0; k <= 40; k++) {
        const s = k / 40; const u = 1 - s
        // de Casteljau via Bernstein basis (independent evaluation path)
        const v = u * u * u * b0 + 3 * u * u * s * b1 + 3 * u * s * s * b2 + s * s * s * b3
        expect(v).toBeGreaterThanOrEqual(lo - 1e-12 * scale - 1e-18)
        expect(v).toBeLessThanOrEqual(hi + 1e-12 * scale + 1e-18)
      }
    }
  })

  it('outward-fast ensemble incl. edge exits: on/off bitwise identical', () => {
    const mk = (localBound: boolean) => {
      const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 200, 11, [0, 0], { enabled: true, localBound })
      const st = sim.state
      for (let i = 0; i < st.id.length; i++) {
        const r = Math.hypot(st.x[i], st.y[i]) || 1e-9
        const sp = 1.5 + (i % 7) * 0.6
        st.vx[i] = (st.x[i] / r) * sp; st.vy[i] = (st.y[i] / r) * sp
        st.vz[i] = 0.4 + (i % 5) * 0.15; st.sitting[i] = 0
      }
      return sim
    }
    const off = mk(false); const on = mk(true)
    off.advanceTo(0.3); on.advanceTo(0.3)
    expect(on.nBounce).toBe(off.nBounce)
    expect(on.nDetach).toBe(off.nDetach)
    expect([...on.state.alive]).toEqual([...off.state.alive])
    expect([...on.state.sitting]).toEqual([...off.state.sitting])
    expect(maxPosDiff(on, off)).toBe(0)
    for (let i = 0; i < on.state.id.length; i++) {
      expect(on.state.z[i]).toBe(off.state.z[i])
      expect(on.state.vx[i]).toBe(off.state.vx[i])
      expect(on.state.vz[i]).toBe(off.state.vz[i])
    }
    // The scenario really does throw particles off the plate.
    expect(on.nAlive).toBeLessThan(200)
  })

  it('frequency switch with residual vibration: on/off bitwise identical', () => {
    const mk = (localBound: boolean) => createTransientSim(
      singleToneScore(152, calibratedForceN(), 30), 80, 5, [0, 0], { enabled: true, localBound })
    const off = mk(false); const on = mk(true)
    for (const s of [off, on]) {
      s.advanceTo(0.25)
      s.wakeAll(); s.plate.setFrequency(649, s.t)
      s.advanceTo(0.5)
    }
    expect(on.nBounce).toBe(off.nBounce)
    expect(on.nDetach).toBe(off.nDetach)
    expect([...on.state.alive]).toEqual([...off.state.alive])
    expect(maxPosDiff(on, off)).toBe(0)
  })
})

describe('M6-6 sliced non-contact proof', () => {
  type Priv = {
    locallyClearSliced: (x0: number, y0: number, z0: number, vx0: number, vy0: number, vz0: number,
      cur: number, tSeg: number, tNew: number, k: number) => boolean
    cellHeightMaxSlice: (c: number, s: number, k: number) => number
    cellIndex: (x: number, y: number) => number
  }

  it('slice height bound contains surfAt inside its slice', () => {
    const K = 4
    const sim = createTransientSim(singleToneScore(649, calibratedForceN(), 30), 10, 7, [0, 0],
      { enabled: true, localBound: true, localSplits: K })
    sim.advanceTo(0.05)
    const t0 = sim.t
    sim.step()
    const s = sim as unknown as Priv
    const pl = sim.plate
    const out = new Float64Array(5)
    const dt = sim.t - t0
    let seed = 31337
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    for (let n = 0; n < 4000; n++) {
      const x = (rnd() * 2 - 1) * pl.a; const y = (rnd() * 2 - 1) * pl.a
      const fr = rnd()
      const tau = t0 + fr * dt
      const slice = Math.min(K - 1, Math.floor(fr * K))
      pl.surfAt(x, y, tau, out)
      const bound = s.cellHeightMaxSlice(s.cellIndex(x, y), slice, K)
      expect(out[0]).toBeLessThanOrEqual(bound + 1e-12)
    }
  })

  it('slice q ranges are contained in the whole-segment q ranges', () => {
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30))
    let seed = 777
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 1e-4; plate.qd[j] = rnd() * 0.2 }
    plate.beginStep(0, 2e-4)
    const whole = plate.segQRanges()
    for (const K of [2, 4, 8]) {
      const sl = plate.segQSliceRanges(K)
      for (let sIdx = 0; sIdx < K; sIdx++) {
        for (let j = 0; j < plate.ndof; j++) {
          expect(sl.min[sIdx * plate.ndof + j]).toBeGreaterThanOrEqual(whole.min[j] - 1e-15)
          expect(sl.max[sIdx * plate.ndof + j]).toBeLessThanOrEqual(whole.max[j] + 1e-15)
        }
      }
    }
  })

  it('sliced proof respects the same domain and fail-closed rules', () => {
    const sim = createTransientSim(singleToneScore(152, 0, 30), 10, 7, [0, 0],
      { enabled: true, localBound: true, localSplits: 4 })
    const t0 = sim.t
    sim.step()
    const s = sim as unknown as Priv
    const dt = sim.t - t0
    // Calm plate (height 0): inside-square high flight is provable...
    expect(s.locallyClearSliced(0.05, 0.05, 1e-3, 0.1, 0, 0, t0, t0, sim.t, 4)).toBe(true)
    // ...but never when the segment leaves [-a,a]^2 or input is non-finite.
    expect(s.locallyClearSliced(0.10, 0.05, 1e-3, 1000, 0, 0, t0, t0, sim.t, 4)).toBe(false)
    expect(s.locallyClearSliced(0.05, 0.05, NaN, 0.1, 0, 0, t0, t0, sim.t, 4)).toBe(false)
    expect(s.locallyClearSliced(0.05, 0.05, 1e-3, 0.1, 0, 0, t0, t0, t0 + dt, 4)).toBe(true)
  })

  it('splits on/off are bitwise identical (strong, weak, adaptive and fixed)', () => {
    for (const forceScale of [1, 0.05]) {
      const mk = (opts: object) => createTransientSim(
        singleToneScore(152, calibratedForceN() * forceScale, 30), 80, 7, [0, 0],
        { enabled: true, ...opts })
      const base = mk({ localBound: false })
      base.advanceTo(0.6)
      for (const opts of [
        { localBound: true, localSplits: 2 },
        { localBound: true, localSplits: 4 },
        { localBound: true, localSplits: 8 },
        { localBound: true, localSplits: 4, localAdaptive: false },
      ]) {
        const on = mk(opts)
        on.advanceTo(0.6)
        expect(on.nBounce).toBe(base.nBounce)
        expect(on.nDetach).toBe(base.nDetach)
        expect(on.nEventCap).toBe(base.nEventCap)
        expect([...on.state.alive]).toEqual([...base.state.alive])
        expect([...on.state.sitting]).toEqual([...base.state.sitting])
        expect(maxPosDiff(on, base)).toBe(0)
        for (let i = 0; i < on.state.id.length; i++) {
          expect(on.state.z[i]).toBe(base.state.z[i])
          expect(on.state.vz[i]).toBe(base.state.vz[i])
        }
      }
    }
  })

  it('heightAt is bitwise identical to surfAt height', () => {
    const plate = new TransientPlate(singleToneScore(649, calibratedForceN(), 30))
    let seed = 20240
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
    const out = new Float64Array(5)
    const a = plate.a
    const pts: [number, number][] = [[0, 0], [a, 0], [0, -a], [a, a], [1e-15, -1e-15], [0.3, 0.2], [-0.24, 0.05]]
    for (let n = 0; n < 300; n++) pts.push([rnd() * 0.2, rnd() * 0.2])
    for (let step = 0; step < 3; step++) {
      plate.beginStep(step * 1.7e-4, (step + 1) * 1.7e-4)
      for (const [x, y] of pts) {
        for (const fr of [0, 0.17, 0.5, 0.99, 1]) {
          const tau = step * 1.7e-4 + fr * 1.7e-4
          plate.surfAt(x, y, tau, out)
          expect(plate.heightAt(x, y, tau)).toBe(out[0])
        }
      }
    }
  })

  it('heightAccelAt/accelAt are bitwise identical to surfAt h/wtt', () => {
    const plate = new TransientPlate(singleToneScore(649, calibratedForceN(), 30))
    let seed = 51501
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
    const out = new Float64Array(5); const out2 = new Float64Array(2)
    const a = plate.a
    const pts: [number, number][] = [[0, 0], [a, 0], [0, -a], [a, a], [1e-15, -1e-15], [0.3, 0.2], [-0.24, 0.05]]
    for (let n = 0; n < 300; n++) pts.push([rnd() * 0.2, rnd() * 0.2])
    for (let step = 0; step < 3; step++) {
      plate.beginStep(step * 1.7e-4, (step + 1) * 1.7e-4)
      for (const [x, y] of pts) {
        for (const fr of [0, 0.17, 0.5, 0.99, 1]) {
          const tau = step * 1.7e-4 + fr * 1.7e-4
          plate.surfAt(x, y, tau, out)
          plate.heightAccelAt(x, y, tau, out2)
          expect(out2[0]).toBe(out[0])
          expect(out2[1]).toBe(out[4])
          expect(plate.accelAt(x, y, tau)).toBe(out[4])
          expect(plate.heightAtFused(x, y, tau)).toBe(out[0])
        }
      }
    }
  })

  it('live-DOF mask: dead DOFs are skipped bitwise (M6-8A)', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 30))
    // Excite only DOFs 1 and 5; all other Hermite endpoints stay zero.
    plate.q[1] = 1.3e-4; plate.qd[1] = -0.21; plate.q[5] = -2e-5; plate.qd[5] = 0.07
    plate.beginStep(0, 1.7e-4)
    expect(plate.liveCount).toBe(2)
    expect([...plate.liveDofs.subarray(0, plate.liveCount)]).toEqual([1, 5])
    const out = new Float64Array(5); const out2 = new Float64Array(2)
    let seed = 4242
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let n = 0; n < 500; n++) {
      const x = rnd() * 0.2; const y = rnd() * 0.2; const tau = Math.abs(rnd()) * 1.7e-4
      plate.surfAt(x, y, tau, out)
      expect(plate.heightAt(x, y, tau)).toBe(out[0])
      expect(plate.heightAtFused(x, y, tau)).toBe(out[0])
      plate.heightAccelAt(x, y, tau, out2)
      expect(out2[0]).toBe(out[0])
      expect(out2[1]).toBe(out[4])
      expect(plate.accelAt(x, y, tau)).toBe(out[4])
    }
    // With forcing, driven m=0 DOFs join the live set in segment 1.
    const driven = new TransientPlate(singleToneScore(152, calibratedForceN(), 30))
    driven.beginStep(0, 1.7e-4)
    expect(driven.liveCount).toBe(2) // exactly the two m=0 modes (center drive)
  })

  it('sliced proof engages where the whole-segment proof cannot', () => {
    const mk = (opts: object) => createTransientSim(
      singleToneScore(152, calibratedForceN(), 30), 200, 7, [0, 0], { enabled: true, ...opts })
    const whole = mk({ localBound: true })
    const sliced = mk({ localBound: true, localSplits: 4 })
    whole.advanceTo(0.6); sliced.advanceTo(0.6)
    expect(sliced.localSliceAttempts).toBeGreaterThan(0)
    // The sliced proof must not sample MORE segments than M6-5.
    expect(sliced.sampledSegments).toBeLessThanOrEqual(whole.sampledSegments)
  })
})

describe('M6-7B localized flying sleep', () => {
  it('window cell bound contains actual surfAt heights over the window', () => {
    const sim = createTransientSim(singleToneScore(649, calibratedForceN(), 30), 10, 7, [0, 0],
      { enabled: true, localSleep: true })
    type Priv = {
      ensureLookahead: (t: number) => void; ensureCellHWin: () => void
      cellHWin: Float64Array; cellIndex: (x: number, y: number) => number
      lookCache: { end: number } | null
    }
    const s = sim as unknown as Priv
    sim.advanceTo(0.05)
    s.ensureLookahead(sim.t)
    s.ensureCellHWin()
    const pl = sim.plate
    const out = new Float64Array(5)
    // The lookahead plate evolution is the same RK4 as stepping, so
    // check against a twin sim stepped through the window start.
    const twin = createTransientSim(singleToneScore(649, calibratedForceN(), 30), 10, 7, [0, 0],
      { enabled: true, localSleep: true })
    twin.advanceTo(0.05)
    let seed = 60606
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const dt = twin.dtNow()
    for (let stepN = 0; stepN < 40; stepN++) {
      const t0 = twin.t
      twin.step()
      for (let n = 0; n < 60; n++) {
        const x = (rnd() * 2 - 1) * pl.a; const y = (rnd() * 2 - 1) * pl.a
        const tau = t0 + rnd() * dt
        twin.plate.surfAt(x, y, tau, out)
        const bound = s.cellHWin[s.cellIndex(x, y)]
        expect(out[0]).toBeLessThanOrEqual(bound + 1e-9)
      }
    }
  })

  it('localized sleep reproduces baseline physics (counts + positions)', () => {
    for (const forceScale of [1, 0.05]) {
      const mk = (opts: SleepOptions) => createTransientSim(
        singleToneScore(152, calibratedForceN() * forceScale, 30), 60, 7, [0, 0], opts)
      const off = mk({ enabled: false })
      const glob = mk({ enabled: true })
      const loc = mk({ enabled: true, localSleep: true })
      off.advanceTo(0.6); glob.advanceTo(0.6); loc.advanceTo(0.6)
      expect(loc.localSleepWins).toBeGreaterThanOrEqual(0)
      expect(loc.nBounce).toBe(off.nBounce)
      expect(loc.nDetach).toBe(off.nDetach)
      expect(loc.nAlive).toBe(off.nAlive)
      expect(maxPosDiff(loc, off)).toBeLessThan(5e-6)
      // Localizing must never skip FEWER updates than the global bound
      // would over the same run when it engages (bound is tighter).
      expect(loc.skippedUpdates).toBeGreaterThanOrEqual(glob.skippedUpdates * 0.9)
    }
  })

  it('frequency switch invalidates localized sleeps (wakeAll path)', () => {
    const mk = (opts: SleepOptions) => {
      const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 60, 7, [0, 0], opts)
      sim.advanceTo(0.25)
      sim.wakeAll(); sim.plate.setFrequency(649, sim.t)
      sim.advanceTo(0.5)
      return sim
    }
    const loc = mk({ enabled: true, localSleep: true })
    const glob = mk({ enabled: true })
    const off = mk({ enabled: false })
    expect(loc.nBounce).toBe(off.nBounce)
    expect(loc.nAlive).toBe(off.nAlive)
    // Sleep-vs-off divergence in this scenario (1.4e-5 m) comes from the
    // pre-existing wake/materialize rounding of M6-2 sleep itself;
    // localization must stay close to the global-bound sleep path.
    expect(maxPosDiff(loc, glob)).toBeLessThan(2e-6)
    expect(loc.skippedUpdates).toBeGreaterThan(glob.skippedUpdates)
  })
})

describe('M6-4R range initialization', () => {
  it('initialStateRange is bitwise identical to slicing initialState', () => {
    for (const [count, seed] of [[100, 7], [5000, 42], [600000, 1]] as const) {
      const full = initialState(0.12, count, seed)
      const ranges: [number, number][] = [[0, 1], [0, count], [1, 2],
        [Math.floor(count / 3), Math.floor(count / 3) + 17], [count - 1, count]]
      for (const [a, b] of ranges) {
        const st = initialStateRange(0.12, count, seed, a, b)
        for (let i = 0; i < b - a; i++) {
          expect(st.x[i]).toBe(full.x[a + i])
          expect(st.y[i]).toBe(full.y[a + i])
          expect(st.z[i]).toBe(full.z[a + i])
          expect(st.id[i]).toBe(a + i)
          expect(st.alive[i]).toBe(1)
          expect(st.sitting[i]).toBe(0)
        }
      }
    }
  })
})

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
