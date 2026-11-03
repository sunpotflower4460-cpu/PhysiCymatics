// M6-9R corrected re-validation. The M6-9 A2 experiment compared
// surfAt(tEnd - 0.37*dt_spp), i.e. DIFFERENT absolute times per spp
// (up to ~48us apart at 152Hz). This file redoes the comparison on a
// single shared absolute-time grid and adds an independent analytic
// free-decay reference that separates RK4 endpoint error from cubic
// Hermite interpolation error. Bench bodies are skipped unless
// PHYSI_BENCH=1; the grid/offset regression test always runs.
import { describe, expect, it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientPlate } from '../plateTransient'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function excited(freq: number, forceN: number, seedIn: number): TransientPlate {
  const plate = new TransientPlate(singleToneScore(freq, forceN, 30))
  let seed = seedIn
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  return plate
}

const PTS_SEED = 333
function samplePoints(plate: TransientPlate, n = 40): [number, number][] {
  let seed = PTS_SEED
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  const pts: [number, number][] = []
  for (let i = 0; i < n; i++) pts.push([(rnd() * 2 - 1) * plate.a * 0.95, (rnd() * 2 - 1) * plate.a * 0.95])
  return pts
}

/** Shared absolute sample times: k * (dt44/4), k = 1..K. dt44/4 is the
 *  spp=176 segment length, so every spp in {44,88,176,352} evaluates
 *  the SAME Float64 values (interior Hermite points or boundaries). */
function commonTimes(freq: number, durationS: number): Float64Array {
  const dT = 1 / (176 * freq)
  const K = Math.round(durationS / dT)
  const ts = new Float64Array(K)
  for (let k = 0; k < K; k++) ts[k] = (k + 1) * dT
  return ts
}

class Runner {
  segT1 = 0
  lastTau = NaN
  constructor(readonly plate: TransientPlate, readonly dt: number) {}
  evalAt(x: number, y: number, t: number, out: Float64Array): void {
    while (this.segT1 < t) { this.plate.beginStep(this.segT1, this.segT1 + this.dt); this.segT1 += this.dt }
    this.lastTau = t
    this.plate.surfAt(x, y, t, out)
  }
}

describe('m69c comparison-time regression', () => {
  it('common grid is bit-identical across spp and offsets are detectable', () => {
    const freq = 152
    const times = commonTimes(freq, 0.005)
    // Every spp runner must receive exactly the same Float64 times.
    for (const spp of [44, 88, 176, 352]) {
      const r = new Runner(excited(freq, 0, 222), 1 / (spp * freq))
      const out = new Float64Array(5)
      for (const t of times) { r.evalAt(0.03, -0.02, t, out); expect(r.lastTau).toBe(t) }
    }
    // An intentional offset must change evaluated heights by a clearly
    // measurable amount somewhere (otherwise alignment is untestable).
    const rA = new Runner(excited(freq, 0, 222), 1 / (44 * freq))
    const rB = new Runner(excited(freq, 0, 222), 1 / (44 * freq))
    const oA = new Float64Array(5); const oB = new Float64Array(5)
    let maxShift = 0
    for (const t of times) {
      rA.evalAt(0.03, -0.02, t, oA)
      rB.evalAt(0.03, -0.02, t + 48e-6, oB)
      maxShift = Math.max(maxShift, Math.abs(oA[0] - oB[0]))
    }
    expect(maxShift).toBeGreaterThan(1e-6)
  })
})

benchIt('M6-9R corrected A3 + analytic separation', () => {
  // ---- A3: forced+excited plate, shared absolute times ----
  for (const freq of [152, 649, 1241.6]) {
    const times = commonTimes(freq, 0.01)
    const sample = (spp: number): Float64Array => {
      const plate = excited(freq, 3.6, 222)
      const pts = samplePoints(plate)
      const r = new Runner(plate, 1 / (spp * freq))
      const out = new Float64Array(5)
      const res = new Float64Array(times.length * pts.length * 2)
      let w = 0
      for (const t of times) for (const [x, y] of pts) { r.evalAt(x, y, t, out); res[w++] = out[0]; res[w++] = out[4] }
      return res
    }
    const ref = sample(352)
    for (const spp of [44, 88, 176]) {
      const got = sample(spp)
      let maxH = 0; let maxW = 0
      for (let i = 0; i < ref.length; i += 2) {
        maxH = Math.max(maxH, Math.abs(ref[i] - got[i]))
        maxW = Math.max(maxW, Math.abs(ref[i + 1] - got[i + 1]))
      }
      console.log(`A3-common freq=${freq} spp=${spp}: max|dh|=${maxH.toExponential(3)} max|dwtt|=${maxW.toExponential(3)} (vs spp=352, same absolute times)`)
    }
    // Flawed M6-9 A2 replica (tEnd - 0.37*dt_spp) for side-by-side record.
    const flawed = (spp: number): Float64Array => {
      const plate = excited(freq, 3.6, 222)
      const pts = samplePoints(plate)
      const dt = 1 / (spp * freq); const dt44 = 1 / (44 * freq)
      const ratio = Math.round(dt44 / dt)
      const total = Math.round(0.01 / dt44)
      const out = new Float64Array(5); const res: number[] = []
      let seg = 0
      for (let s = 0; s < total; s++) {
        for (let u = 0; u < ratio; u++) { plate.beginStep(seg * dt, (seg + 1) * dt); seg++ }
        const tEnd = (s + 1) * dt44
        for (const [x, y] of pts) { plate.surfAt(x, y, tEnd - 0.37 * dt, out); res.push(out[0]) }
      }
      return Float64Array.from(res)
    }
    const fRef = flawed(352)
    for (const spp of [44, 88, 176]) {
      const got = flawed(spp)
      let maxH = 0
      for (let i = 0; i < fRef.length; i++) maxH = Math.max(maxH, Math.abs(fRef[i] - got[i]))
      console.log(`A2-flawed-replica freq=${freq} spp=${spp}: max|dh|=${maxH.toExponential(3)} (times differ by up to ${(0.37 * (1 / (spp * freq) - 1 / (352 * freq)) * 1e6).toFixed(1)}us)`)
    }
  }

  // ---- Analytic free decay: separate RK4 endpoint vs Hermite error ----
  // q'' + 2*zeta*w*q' + w^2 q = 0, zeta<1:
  // q(t)=e^{-aw t}[q0 cos(wd t) + ((qd0+aw q0)/wd) sin(wd t)], a=zeta.
  const base = excited(152, 0, 222)
  const zeta = base.zeta
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
  const q0s = Float64Array.from(base.q); const qd0s = Float64Array.from(base.qd)
  const pts = samplePoints(base)
  const phi = new Float64Array(base.ndof)
  // (1) RK4 endpoint error + total surface error vs analytic, per spp
  // (segment length set by "drive freq" 152Hz, as in the sim).
  for (const spp of [44, 88, 176, 352]) {
    const plate = excited(152, 0, 222)
    const dt = 1 / (spp * 152)
    const times = commonTimes(152, 0.01)
    const r = new Runner(plate, dt)
    const out = new Float64Array(5)
    let maxSurfErr = 0; let maxQEndErr = 0
    for (const t of times) {
      for (const [x, y] of pts) {
        r.evalAt(x, y, t, out)
        plate.phiAll(x, y, phi)
        let hAna = 0
        for (let j = 0; j < plate.ndof; j++) hAna += qAna(plate.omega[j], q0s[j], qd0s[j], t) * phi[j]
        maxSurfErr = Math.max(maxSurfErr, Math.abs(out[0] - hAna))
      }
      // endpoint error at the runner's current segment end
      const tEnd = r.segT1
      for (let j = 0; j < plate.ndof; j++) {
        maxQEndErr = Math.max(maxQEndErr, Math.abs(plate.q[j] - qAna(plate.omega[j], q0s[j], qd0s[j], tEnd)))
      }
    }
    console.log(`ANA spp=${spp} (seg=${(dt * 1e6).toFixed(1)}us): RK4 endpoint max|dq|=${maxQEndErr.toExponential(3)} ; total surface max|dh vs analytic|=${maxSurfErr.toExponential(3)}`)
  }
  // (2) Pure Hermite interpolation error: analytic endpoints in,
  // cubic Hermite (replicated formula) vs analytic interior values.
  const hermite = (q0: number, qd0: number, q1: number, qd1: number, s: number, h: number): number => {
    const h00 = 2 * s ** 3 - 3 * s ** 2 + 1; const h10 = s ** 3 - 2 * s ** 2 + s
    const h01 = -2 * s ** 3 + 3 * s ** 2; const h11 = s ** 3 - s ** 2
    return h00 * q0 + h10 * h * qd0 + h01 * q1 + h11 * h * qd1
  }
  for (const spp of [44, 88, 176, 352]) {
    const h = 1 / (spp * 152)
    let maxErr = 0
    for (let seg = 0; seg < 20; seg++) {
      const t0 = seg * h
      for (let j = 0; j < base.ndof; j++) {
        const w = base.omega[j]
        const qA = qAna(w, q0s[j], qd0s[j], t0); const qdA = qdAna(w, q0s[j], qd0s[j], t0)
        const qB = qAna(w, q0s[j], qd0s[j], t0 + h); const qdB = qdAna(w, q0s[j], qd0s[j], t0 + h)
        for (const s of [0.13, 0.37, 0.5, 0.71, 0.9]) {
          const got = hermite(qA, qdA, qB, qdB, s, h)
          const want = qAna(w, q0s[j], qd0s[j], t0 + s * h)
          maxErr = Math.max(maxErr, Math.abs(got - want))
        }
      }
    }
    console.log(`ANA pure-Hermite seg=${(h * 1e6).toFixed(1)}us (spp=${spp}@152Hz): max|dq|=${maxErr.toExponential(3)} (analytic endpoints, all DOFs incl. 1241Hz mode)`)
  }
  // (3) Forced convergence by halving (same RK4 family — only a
  // self-consistency check, anchored by the analytic case above).
  for (const freq of [649, 1241.6]) {
    const times = commonTimes(freq, 0.01)
    const samp = (spp: number): Float64Array => {
      const plate = excited(freq, 3.6, 222)
      const pp = samplePoints(plate)
      const r = new Runner(plate, 1 / (spp * freq))
      const out = new Float64Array(5); const res: number[] = []
      for (const t of times) for (const [x, y] of pp) { r.evalAt(x, y, t, out); res.push(out[0]) }
      return Float64Array.from(res)
    }
    const a = samp(44); const b = samp(88); const c = samp(176)
    const d = (u: Float64Array, v: Float64Array) => { let m = 0; for (let i = 0; i < u.length; i++) m = Math.max(m, Math.abs(u[i] - v[i])); return m }
    console.log(`FORCED-CONV freq=${freq}: |44-88|=${d(a, b).toExponential(3)} |88-176|=${d(b, c).toExponential(3)} ratio=${(d(a, b) / Math.max(d(b, c), 1e-300)).toFixed(2)}`)
  }
}, 600000)
