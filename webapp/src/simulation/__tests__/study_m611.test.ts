// M6-11 plate-only study: kinematic/dynamic consistency of time
// representations, against INDEPENDENT analytic solutions for both
// free decay and sinusoidal forcing (attack=0 score, so the drive is
// exactly Fpk*sin(2*pi*f*t) and each modal ODE has a closed form).
// Candidates (all test-side, production untouched by this file):
//  A: cubic Hermite q/qd/qdd (current production representation)
//  B: Hermite q/qd + ODE qdd (M6-10 candidate B)
//  C: quintic Hermite using endpoint q/qd/qdd(ODE) — q, qd, qdd are
//     derivatives of one polynomial, and qdd is continuous across
//     segment boundaries by construction.
// Metrics: errors vs analytic (h, uz, wtt); kinematic defect
// (implemented wtt minus exact 2nd derivative of implemented h);
// dynamic (ODE) residual; boundary jumps of h/uz/wtt; detach-root
// finder miss rates for B/C. Skipped unless PHYSI_BENCH=1.
import { it } from 'vitest'
import { ScoreSignal, type PerformanceScore } from '../sequence'
import { TransientPlate } from '../plateTransient'
import { calibratedForceN } from '../transientSand'
import { G } from '../surface'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function pureToneScore(freq: number, forceN: number): PerformanceScore {
  return {
    schema: 'physicymatics.score.v1', duration_s: 30, initialPhase_rad: 0, audioGain: 0.2,
    events: [{ start_s: 0, end_s: 30, frequency_hz: freq, forcePeak_N: forceN, attack_s: 0, release_s: 0 }],
  }
}
function excited(score: PerformanceScore, driveXY: [number, number], seedIn: number, presetScale = 1): TransientPlate {
  const plate = new TransientPlate(score, driveXY)
  let seed = seedIn
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4 * presetScale; plate.qd[j] = rnd() * 0.3 * presetScale }
  return plate
}

function hermite3(q0: number, qd0: number, q1: number, qd1: number, h: number, s: number): [number, number, number] {
  const s2 = s * s; const s3 = s2 * s
  const q = (2 * s3 - 3 * s2 + 1) * q0 + (s3 - 2 * s2 + s) * h * qd0 + (-2 * s3 + 3 * s2) * q1 + (s3 - s2) * h * qd1
  const dq = (6 * s2 - 6 * s) * q0 + (3 * s2 - 4 * s + 1) * h * qd0 + (-6 * s2 + 6 * s) * q1 + (3 * s2 - 2 * s) * h * qd1
  const ddq = (12 * s - 6) * q0 + (6 * s - 4) * h * qd0 + (-12 * s + 6) * q1 + (6 * s - 2) * h * qd1
  return [q, dq / h, ddq / (h * h)]
}
function hermite5(q0: number, qd0: number, a0: number, q1: number, qd1: number, a1: number, h: number, s: number): [number, number, number] {
  const s2 = s * s; const s3 = s2 * s; const s4 = s3 * s; const s5 = s4 * s
  const H0 = 1 - 10 * s3 + 15 * s4 - 6 * s5; const H1 = s - 6 * s3 + 8 * s4 - 3 * s5
  const H2 = (s2 - 3 * s3 + 3 * s4 - s5) / 2; const H3 = 10 * s3 - 15 * s4 + 6 * s5
  const H4 = -4 * s3 + 7 * s4 - 3 * s5; const H5 = (s3 - 2 * s4 + s5) / 2
  const dH0 = -30 * s2 + 60 * s3 - 30 * s4; const dH1 = 1 - 18 * s2 + 32 * s3 - 15 * s4
  const dH2 = (2 * s - 9 * s2 + 12 * s3 - 5 * s4) / 2; const dH3 = 30 * s2 - 60 * s3 + 30 * s4
  const dH4 = -12 * s2 + 28 * s3 - 15 * s4; const dH5 = (3 * s2 - 8 * s3 + 5 * s4) / 2
  const ddH0 = -60 * s + 180 * s2 - 120 * s3; const ddH1 = -36 * s + 96 * s2 - 60 * s3
  const ddH2 = (2 - 18 * s + 36 * s2 - 20 * s3) / 2; const ddH3 = 60 * s - 180 * s2 + 120 * s3
  const ddH4 = -24 * s + 84 * s2 - 60 * s3; const ddH5 = (6 * s - 24 * s2 + 20 * s3) / 2
  const hh = h * h
  const q = H0 * q0 + H1 * h * qd0 + H2 * hh * a0 + H3 * q1 + H4 * h * qd1 + H5 * hh * a1
  const dq = dH0 * q0 + dH1 * h * qd0 + dH2 * hh * a0 + dH3 * q1 + dH4 * h * qd1 + dH5 * hh * a1
  const ddq = ddH0 * q0 + ddH1 * h * qd0 + ddH2 * hh * a0 + ddH3 * q1 + ddH4 * h * qd1 + ddH5 * hh * a1
  return [q, dq / h, ddq / hh]
}

interface SegE { t0: number; q0: Float64Array; qd0: Float64Array; a0: Float64Array; q1: Float64Array; qd1: Float64Array; a1: Float64Array }

benchIt('M6-11 consistency A/B/C vs analytic (free + forced)', () => {
  const FN = calibratedForceN()
  const T = 0.01
  const makePts = (plate: TransientPlate) => {
    let seed = 333
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const pts: [number, number][] = []
    for (let i = 0; i < 20; i++) pts.push([(rnd() * 2 - 1) * plate.a * 0.95, (rnd() * 2 - 1) * plate.a * 0.95])
    return pts
  }

  const runCase = (label: string, fDrive: number, forceN: number, driveXY: [number, number], presetScale = 1) => {
    const score = pureToneScore(fDrive, forceN)
    const signal = new ScoreSignal(score)
    const forceAt = (t: number) => signal.forcePeakAt(t) * signal.u(t)
    const plate = excited(score, driveXY, 222, presetScale)
    const zeta = plate.zeta
    const qInit = Float64Array.from(plate.q); const qdInit = Float64Array.from(plate.qd)
    const Om = 2 * Math.PI * fDrive
    // Analytic per-DOF solution for q''+2*zeta*w*q'+w^2 q = Aj*sin(Om t).
    const ana = (j: number, t: number): [number, number, number] => {
      const w = plate.omega[j]
      const A = (plate.phiDrive[j] * forceN) / plate.mass[j]
      const den = (w * w - Om * Om) ** 2 + (2 * zeta * w * Om) ** 2
      const al = den > 0 ? (A * (w * w - Om * Om)) / den : 0
      const be = den > 0 ? (-A * 2 * zeta * w * Om) / den : 0
      const qp = al * Math.sin(Om * t) + be * Math.cos(Om * t)
      const qdp = al * Om * Math.cos(Om * t) - be * Om * Math.sin(Om * t)
      const wd = w * Math.sqrt(1 - zeta * zeta)
      const h0 = qInit[j] - be; const hd0 = qdInit[j] - al * Om
      const c = (hd0 + zeta * w * h0) / wd
      const e = Math.exp(-zeta * w * t)
      const qh = e * (h0 * Math.cos(wd * t) + c * Math.sin(wd * t))
      const qdh = e * (-zeta * w * (h0 * Math.cos(wd * t) + c * Math.sin(wd * t)) + wd * (-h0 * Math.sin(wd * t) + c * Math.cos(wd * t)))
      const q = qh + qp; const qd = qdh + qdp
      const qdd = A * Math.sin(Om * t) - 2 * zeta * w * qd - w * w * q
      return [q, qd, qdd]
    }
    const phi = new Float64Array(plate.ndof)
    const surfAna = (x: number, y: number, t: number): [number, number, number] => {
      plate.phiAll(x, y, phi)
      let h = 0; let u = 0; let w = 0
      for (let j = 0; j < plate.ndof; j++) { const [q, qd, qdd] = ana(j, t); h += q * phi[j]; u += qd * phi[j]; w += qdd * phi[j] }
      return [h, u, w]
    }
    // Verify analytic derivation against plate RK4 endpoints first.
    {
      const [qE] = [plate.advanceState(qInit, qdInit, 0, T)]
      let maxE = 0
      for (let j = 0; j < plate.ndof; j++) maxE = Math.max(maxE, Math.abs(qE[0][j] - ana(j, T)[0]))
      console.log(`M611 ${label}: analytic-vs-RK4 endpoint max|dq| at T=${T} = ${maxE.toExponential(3)} (derivation check; RK4 err scale)`)
    }
    const odeAccel = (j: number, q: number, qd: number, t: number): number =>
      (plate.phiDrive[j] * forceAt(t)) / plate.mass[j] - 2 * zeta * plate.omega[j] * qd - plate.omega[j] ** 2 * q
    const hSeg = 1 / (44 * fDrive)
    const segs: SegE[] = []
    {
      let q: Float64Array = Float64Array.from(qInit); let qd: Float64Array = Float64Array.from(qdInit)
      for (let t0 = 0; t0 < T; t0 += hSeg) {
        const [q1, qd1] = plate.advanceState(q, qd, t0, t0 + hSeg)
        const a0 = new Float64Array(plate.ndof); const a1 = new Float64Array(plate.ndof)
        for (let j = 0; j < plate.ndof; j++) { a0[j] = odeAccel(j, q[j], qd[j], t0); a1[j] = odeAccel(j, q1[j], qd1[j], t0 + hSeg) }
        segs.push({ t0, q0: q, qd0: qd, a0, q1, qd1, a1 })
        q = q1; qd = qd1
      }
    }
    type Mode = 'A' | 'B' | 'C'
    const evalCand = (mode: Mode, x: number, y: number, t: number): [number, number, number, number, number] => {
      const idx = Math.max(0, Math.min(segs.length - 1, Math.ceil(t / hSeg) - 1))
      const seg = segs[idx]
      const s = (t - seg.t0) / hSeg
      plate.phiAll(x, y, phi)
      let h = 0; let u = 0; let w = 0; let kinDef = 0; let dynRes = 0
      for (let j = 0; j < plate.ndof; j++) {
        let q: number; let qd: number; let qdd: number; let qddOwn: number
        if (mode === 'C') {
          [q, qd, qdd] = hermite5(seg.q0[j], seg.qd0[j], seg.a0[j], seg.q1[j], seg.qd1[j], seg.a1[j], hSeg, s)
          qddOwn = qdd
        } else {
          [q, qd, qddOwn] = hermite3(seg.q0[j], seg.qd0[j], seg.q1[j], seg.qd1[j], hSeg, s)
          qdd = mode === 'B' ? odeAccel(j, q, qd, t) : qddOwn
        }
        h += q * phi[j]; u += qd * phi[j]; w += qdd * phi[j]
        kinDef += (qdd - qddOwn) * phi[j]
        dynRes += (qdd - odeAccel(j, q, qd, t)) * phi[j]
      }
      return [h, u, w, kinDef, dynRes]
    }
    const pts = makePts(plate)
    const times: number[] = []
    for (let t = 1 / (176 * fDrive); t <= T; t += 1 / (176 * fDrive)) times.push(t + 1e-7)
    for (const mode of ['A', 'B', 'C'] as Mode[]) {
      let eH = 0; let eU = 0; let eW = 0; let kD = 0; let dR = 0
      for (const t of times) for (const [x, y] of pts) {
        const [ha, ua, wa] = surfAna(x, y, t)
        const [h, u, w, kd, dr] = evalCand(mode, x, y, t)
        eH = Math.max(eH, Math.abs(h - ha)); eU = Math.max(eU, Math.abs(u - ua)); eW = Math.max(eW, Math.abs(w - wa))
        kD = Math.max(kD, Math.abs(kd)); dR = Math.max(dR, Math.abs(dr))
      }
      console.log(`M611 ${label} cand=${mode}: max|dh|=${eH.toExponential(3)} max|duz|=${eU.toExponential(3)} max|dwtt|=${eW.toExponential(3)} kinDefect=${kD.toExponential(3)} dynResidual=${dR.toExponential(3)}`)
    }
    for (const mode of ['A', 'B', 'C'] as Mode[]) {
      let jH = 0; let jU = 0; let jW = 0
      const evalAt = (seg: SegE, s: number, x: number, y: number): [number, number, number] => {
        plate.phiAll(x, y, phi)
        let h = 0; let u = 0; let w = 0
        const t = seg.t0 + s * hSeg
        for (let j = 0; j < plate.ndof; j++) {
          let q: number; let qd: number; let qdd: number
          if (mode === 'C') [q, qd, qdd] = hermite5(seg.q0[j], seg.qd0[j], seg.a0[j], seg.q1[j], seg.qd1[j], seg.a1[j], hSeg, s)
          else { [q, qd, qdd] = hermite3(seg.q0[j], seg.qd0[j], seg.q1[j], seg.qd1[j], hSeg, s); if (mode === 'B') qdd = odeAccel(j, q, qd, t) }
          h += q * phi[j]; u += qd * phi[j]; w += qdd * phi[j]
        }
        return [h, u, w]
      }
      for (let k = 0; k + 1 < segs.length; k++) for (const [x, y] of pts) {
        const L = evalAt(segs[k], 1, x, y); const R = evalAt(segs[k + 1], 0, x, y)
        jH = Math.max(jH, Math.abs(L[0] - R[0])); jU = Math.max(jU, Math.abs(L[1] - R[1])); jW = Math.max(jW, Math.abs(L[2] - R[2]))
      }
      console.log(`M611 ${label} cand=${mode} boundaryJumps: h=${jH.toExponential(3)} uz=${jU.toExponential(3)} wtt=${jW.toExponential(3)}`)
    }
    return { plate, segs, hSeg, evalCand, surfAna, pts }
  }

  const free = runCase('free152', 152, 0, [0, 0])
  runCase('forced152-ecc', 152, FN, [0.05, 0.03])
  runCase('forced649-ecc', 649, FN, [0.05, 0.03])
  runCase('forced1241-ecc', 1241.6, FN, [0.05, 0.03])
  runCase('forcedPure152-ecc', 152, FN, [0.05, 0.03], 0)
  runCase('forcedPure649-ecc', 649, FN, [0.05, 0.03], 0)

  // Detach-root finder for B/C (free decay): reference = dense scan of
  // the candidate's own wtt; truth comparison = analytic roots.
  {
    const { segs, hSeg, evalCand, surfAna, pts } = free
    const findRoots = (f: (t: number) => number, perSeg: number): number[] => {
      const roots: number[] = []
      for (const seg of segs) {
        let tp = seg.t0; let fp = f(tp)
        for (let kS = 1; kS <= perSeg; kS++) {
          const t = seg.t0 + (hSeg * kS) / perSeg; const ft = f(t)
          if (fp === 0) roots.push(tp)
          else if (ft !== 0 && Math.sign(ft) !== Math.sign(fp)) {
            let lo = tp; let hi = t
            for (let it2 = 0; it2 < 50; it2++) { const m = 0.5 * (lo + hi); if (Math.sign(f(m)) === Math.sign(fp)) lo = m; else hi = m }
            roots.push(0.5 * (lo + hi))
          }
          tp = t; fp = ft
        }
      }
      return roots
    }
    for (const mode of ['B', 'C'] as const) {
      for (const K of [2, 4, 8, 16]) {
        let miss = 0; let extra = 0; let refN = 0; let maxAnaErr = 0; let anaMiss = 0; let anaN = 0
        for (const [x, y] of pts.slice(0, 6)) {
          const f = (t: number) => evalCand(mode, x, y, Math.min(t, T - 1e-12))[2] + G
          const ref = findRoots(f, 2000)
          const got = findRoots(f, K)
          refN += ref.length
          for (const r of ref) if (!got.some((g2) => Math.abs(g2 - r) < 1e-6)) miss++
          for (const g2 of got) if (!ref.some((r) => Math.abs(g2 - r) < 1e-6)) extra++
          const ra = findRoots((t) => surfAna(x, y, Math.min(t, T - 1e-12))[2] + G, 2000)
          anaN += ra.length
          for (const r of ra) {
            let best = Infinity
            for (const g2 of got) best = Math.min(best, Math.abs(g2 - r))
            if (best < 2e-4) maxAnaErr = Math.max(maxAnaErr, best); else anaMiss++
          }
        }
        console.log(`M611 detach-finder mode=${mode} K=${K}: denseRefRoots=${refN} missedVsDense=${miss} extraVsDense=${extra} | vs analytic roots=${anaN} unmatched=${anaMiss} maxTimeErr=${maxAnaErr.toExponential(3)}s`)
      }
    }
  }
}, 900000)
