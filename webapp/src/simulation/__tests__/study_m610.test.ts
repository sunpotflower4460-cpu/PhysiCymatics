// M6-10 plate-only study: decouple the plate time representation from
// the drive frequency. Candidates are evaluated WITHOUT touching
// production code:
//  CUR : current representation (segment = 1/(44*f_drive), Hermite wtt)
//  A   : sub-segmented Hermite — segment capped so w_max*h <= C,
//        endpoints from the plate's own RK4 (advanceState)
//  B   : ODE acceleration — qdd from Hermite q/qd via the modal ODE
//        (including the drive term when forced)
//  C   : finer global segment (spp176-like) as an upper reference
// Truth for the free-decay part is the analytic damped-oscillator
// solution (independent of RK4+Hermite). The forced part uses the
// production surfAt at spp=352 as a self-reference (labelled as such;
// it shares the RK4 family, so it only checks representation
// differences, anchored by the analytic part).
// Metrics: max |dh|, max |dwtt| on common absolute times;
// detach-threshold (wtt + g = 0) root-time errors; cost indications.
// Skipped unless PHYSI_BENCH=1.
import { it } from 'vitest'
import { ScoreSignal, singleToneScore } from '../sequence'
import { TransientPlate } from '../plateTransient'
import { calibratedForceN } from '../transientSand'
import { G } from '../surface'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function excited(freq: number, forceN: number, seedIn: number): TransientPlate {
  const plate = new TransientPlate(singleToneScore(freq, forceN, 30))
  let seed = seedIn
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  return plate
}

function hermite(q0: number, qd0: number, q1: number, qd1: number, h: number, s: number): [number, number, number] {
  const s2 = s * s; const s3 = s2 * s
  const q = (2 * s3 - 3 * s2 + 1) * q0 + (s3 - 2 * s2 + s) * h * qd0 + (-2 * s3 + 3 * s2) * q1 + (s3 - s2) * h * qd1
  const dqds = (6 * s2 - 6 * s) * q0 + (3 * s2 - 4 * s + 1) * h * qd0 + (-6 * s2 + 6 * s) * q1 + (3 * s2 - 2 * s) * h * qd1
  const d2qds2 = (12 * s - 6) * q0 + (6 * s - 4) * h * qd0 + (-12 * s + 6) * q1 + (6 * s - 2) * h * qd1
  return [q, dqds / h, d2qds2 / (h * h)]
}

interface Seg { t0: number; q0: Float64Array; qd0: Float64Array; q1: Float64Array; qd1: Float64Array }

function endpoints(plate: TransientPlate, hSeg: number, T: number): Seg[] {
  let q: Float64Array = Float64Array.from(plate.q); let qd: Float64Array = Float64Array.from(plate.qd)
  const segs: Seg[] = []
  for (let t0 = 0; t0 < T; t0 += hSeg) {
    const [q1, qd1] = plate.advanceState(q, qd, t0, t0 + hSeg)
    segs.push({ t0, q0: q, qd0: qd, q1, qd1 })
    q = q1; qd = qd1
  }
  return segs
}

benchIt('M6-10 plate representation candidates', () => {
  const base0 = excited(152, 0, 222)
  const zeta = base0.zeta
  const q0s = Float64Array.from(base0.q); const qd0s = Float64Array.from(base0.qd)
  const phi = new Float64Array(base0.ndof)
  let seed = 333
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  const pts: [number, number][] = []
  for (let i = 0; i < 20; i++) pts.push([(rnd() * 2 - 1) * base0.a * 0.95, (rnd() * 2 - 1) * base0.a * 0.95])
  const T = 0.01

  const evalSurf = (plate: TransientPlate, segs: Seg[], hSeg: number, x: number, y: number, t: number,
    odeAccel: boolean, forceAt: ((t: number) => number) | null): [number, number] => {
    // Production convention: tau belongs to the segment [t0, t1] with
    // t1 >= tau (surfAt evaluates inside the closing segment). Hermite
    // qdd is discontinuous across segment boundaries, so this matters.
    const segIdx = Math.max(0, Math.min(segs.length - 1, Math.ceil(t / hSeg) - 1))
    const seg = segs[segIdx]
    const s = (t - seg.t0) / hSeg
    plate.phiAll(x, y, phi)
    const F = odeAccel && forceAt ? forceAt(t) : 0
    let h = 0; let w = 0
    for (let j = 0; j < plate.ndof; j++) {
      const [q, qd, qddH] = hermite(seg.q0[j], seg.qd0[j], seg.q1[j], seg.qd1[j], hSeg, s)
      const om = plate.omega[j]
      const qdd = odeAccel ? (plate.phiDrive[j] * F) / plate.mass[j] - 2 * zeta * om * qd - om * om * q : qddH
      h += q * phi[j]; w += qdd * phi[j]
    }
    return [h, w]
  }

  // ===== Part 1: free decay vs analytic truth (segment driver 152Hz) =====
  const FREQ = 152
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
  const qddAna = (w: number, q0: number, qd0: number, t: number): number =>
    -2 * zeta * w * qdAna(w, q0, qd0, t) - w * w * qAna(w, q0, qd0, t)
  const surfAna = (x: number, y: number, t: number): [number, number] => {
    base0.phiAll(x, y, phi)
    let h = 0; let w = 0
    for (let j = 0; j < base0.ndof; j++) { h += qAna(base0.omega[j], q0s[j], qd0s[j], t) * phi[j]; w += qddAna(base0.omega[j], q0s[j], qd0s[j], t) * phi[j] }
    return [h, w]
  }
  const times: number[] = []
  for (let t = 1 / (176 * FREQ); t <= T; t += 1 / (176 * FREQ)) times.push(t)
  const wMax = Math.max(...base0.omega)
  const dt44 = 1 / (44 * FREQ)
  const variants: [string, number, boolean][] = [
    ['CUR seg=dt44 Hermite-wtt', dt44, false],
    ['B ODE-wtt seg=dt44', dt44, true],
    ['A subseg C=0.6', Math.min(dt44, 0.6 / wMax), false],
    ['A subseg C=0.3', Math.min(dt44, 0.3 / wMax), false],
    ['A subseg C=0.15', Math.min(dt44, 0.15 / wMax), false],
    ['A+B C=0.3 ODE-wtt', Math.min(dt44, 0.3 / wMax), true],
    ['C seg=dt44/4', dt44 / 4, false],
  ]
  for (const [name, hSeg, ode] of variants) {
    const plate = excited(FREQ, 0, 222)
    const segs = endpoints(plate, hSeg, T)
    let maxH = 0; let maxW = 0
    for (const t of times) for (const [x, y] of pts) {
      const [ha, wa] = surfAna(x, y, t)
      const [h, w] = evalSurf(plate, segs, hSeg, x, y, t, ode, null)
      maxH = Math.max(maxH, Math.abs(h - ha)); maxW = Math.max(maxW, Math.abs(w - wa))
    }
    console.log(`M610 free ${name}: hSeg=${(hSeg * 1e6).toFixed(1)}us subsegsPerStep44=${(dt44 / hSeg).toFixed(1)} max|dh|=${maxH.toExponential(3)} max|dwtt|=${maxW.toExponential(3)}`)
  }
  // Cross-check: test-side CUR must match production surfAt.
  {
    const plate = excited(FREQ, 0, 222)
    const segs = endpoints(plate, dt44, T)
    const out = new Float64Array(5)
    let maxH = 0; let maxW = 0; let segT1 = 0; let nBadW = 0; let nTot = 0
    for (const tRaw of times) {
      const t = tRaw + 1e-7 // stay off segment boundaries (wtt is discontinuous there)
      while (segT1 < t) { plate.beginStep(segT1, segT1 + dt44); segT1 += dt44 }
      for (const [x, y] of pts) {
        plate.surfAt(x, y, t, out)
        const [h, w] = evalSurf(plate, segs, dt44, x, y, t, false, null)
        maxH = Math.max(maxH, Math.abs(h - out[0])); maxW = Math.max(maxW, Math.abs(w - out[4]))
        nTot++; if (Math.abs(w - out[4]) > 1e-6) nBadW++
      }
    }
    console.log(`M610 crosscheck test-side vs production surfAt (interior times): max|dh|=${maxH.toExponential(3)} max|dwtt|=${maxW.toExponential(3)} badW=${nBadW}/${nTot}`)
  }
  // Detach-threshold root times vs analytic.
  const rootsOf = (f: (t: number) => number, tMax: number): number[] => {
    const roots: number[] = []
    const n = 20000; let tp = 0; let fp = f(0)
    for (let i = 1; i <= n; i++) {
      const t = (i / n) * tMax; const ft = f(t)
      if (fp !== 0 && ft !== 0 && Math.sign(ft) !== Math.sign(fp)) {
        let lo = tp; let hi = t
        for (let k = 0; k < 60; k++) { const m = 0.5 * (lo + hi); if (Math.sign(f(m)) === Math.sign(fp)) lo = m; else hi = m }
        roots.push(0.5 * (lo + hi))
      }
      tp = t; fp = ft
    }
    return roots
  }
  for (const [name, hSeg, ode] of variants) {
    const plate = excited(FREQ, 0, 222)
    const segs = endpoints(plate, hSeg, T)
    let maxRootErr = 0; let matched = 0; let anaCount = 0; let varCount = 0
    for (const [x, y] of pts.slice(0, 6)) {
      const ra = rootsOf((t) => surfAna(x, y, t)[1] + G, T)
      const rv = rootsOf((t) => evalSurf(plate, segs, hSeg, x, y, t, ode, null)[1] + G, T)
      anaCount += ra.length; varCount += rv.length
      for (const r of ra) {
        let best = Infinity
        for (const v of rv) best = Math.min(best, Math.abs(v - r))
        if (best < 2e-4) { matched++; maxRootErr = Math.max(maxRootErr, best) }
      }
    }
    console.log(`M610 detach-roots ${name}: analytic=${anaCount} variant=${varCount} matched=${matched} maxRootTimeErr=${maxRootErr.toExponential(3)}s`)
  }

  // ===== Part 2: forced drive, production spp=352 as self-reference =====
  for (const fDrive of [152, 649]) {
    const FN = calibratedForceN()
    const score = singleToneScore(fDrive, FN, 30)
    const signal = new ScoreSignal(score)
    const forceAt = (t: number) => signal.forcePeakAt(t) * signal.u(t)
    const dtF = 1 / (44 * fDrive)
    const tGrid: number[] = []
    for (let t = 1 / (176 * fDrive); t <= T; t += 1 / (176 * fDrive)) tGrid.push(t)
    // Reference samples from production surfAt with dt44/8 segments.
    const refPlate = excited(fDrive, FN, 222)
    const refDt = dtF / 8
    const ref = new Float64Array(tGrid.length * pts.length * 2)
    {
      const out = new Float64Array(5)
      let segT1 = 0; let wI = 0
      for (const t of tGrid) {
        while (segT1 < t) { refPlate.beginStep(segT1, segT1 + refDt); segT1 += refDt }
        for (const [x, y] of pts) { refPlate.surfAt(x, y, t, out); ref[wI++] = out[0]; ref[wI++] = out[4] }
      }
    }
    const fVariants: [string, number, boolean][] = [
      ['CUR seg=dt44 Hermite-wtt', dtF, false],
      ['B ODE-wtt seg=dt44', dtF, true],
      ['A subseg C=0.3', Math.min(dtF, 0.3 / wMax), false],
      ['A+B C=0.3 ODE-wtt', Math.min(dtF, 0.3 / wMax), true],
      ['C seg=dt44/4', dtF / 4, false],
    ]
    for (const [name, hSeg, ode] of fVariants) {
      const plate = excited(fDrive, FN, 222)
      const segs = endpoints(plate, hSeg, T)
      let maxH = 0; let maxW = 0; let wI = 0
      for (const t of tGrid) for (const [x, y] of pts) {
        const [h, w] = evalSurf(plate, segs, hSeg, x, y, t, ode, forceAt)
        maxH = Math.max(maxH, Math.abs(h - ref[wI])); maxW = Math.max(maxW, Math.abs(w - ref[wI + 1])); wI += 2
      }
      console.log(`M610 forced f=${fDrive} ${name}: hSeg=${(hSeg * 1e6).toFixed(1)}us subsegsPerStep44=${(dtF / hSeg).toFixed(1)} max|dh|=${maxH.toExponential(3)} max|dwtt|=${maxW.toExponential(3)} (vs production spp352 self-ref)`)
    }
  }

  // ===== Cost indications (same session ratios only) =====
  for (const [name, hSeg] of [['dt44@152', dt44], ['C=0.3', Math.min(dt44, 0.3 / wMax)], ['C=0.15', Math.min(dt44, 0.15 / wMax)]] as [string, number][]) {
    const plate = excited(FREQ, 0, 222)
    const t0 = performance.now()
    for (let rep = 0; rep < 20; rep++) endpoints(plate, hSeg, T)
    console.log(`M610 cost endpoints ${name}: ${((performance.now() - t0) / 20).toFixed(2)}ms per 0.01s span`)
  }
  {
    const plate = excited(FREQ, 0, 222)
    plate.beginStep(0, dt44)
    const out = new Float64Array(5)
    const N = 200000
    let acc = 0
    let t0 = performance.now()
    for (let i = 0; i < N; i++) { const p = pts[i % pts.length]; plate.surfAt(p[0], p[1], ((i % 97) / 97) * dt44, out); acc += out[0] }
    const prodMs = performance.now() - t0
    console.log(`M610 cost eval production surfAt: ${((N / prodMs) * 1000).toExponential(3)} evals/s (acc=${acc.toFixed(3)})`)
    for (const [name, hSeg, ode] of variants.slice(0, 4)) {
      const p2 = excited(FREQ, 0, 222)
      const segs = endpoints(p2, hSeg, T)
      acc = 0; t0 = performance.now()
      for (let i = 0; i < N; i++) { const p = pts[i % pts.length]; const [h] = evalSurf(p2, segs, hSeg, p[0], p[1], ((i % 997) / 997) * T, ode, null); acc += h }
      const ms = performance.now() - t0
      console.log(`M610 cost eval test-side ${name}: ${((N / ms) * 1000).toExponential(3)} evals/s ratio-vs-prod=${(ms / prodMs).toFixed(2)} (acc=${acc.toFixed(3)})`)
    }
  }
}, 900000)
