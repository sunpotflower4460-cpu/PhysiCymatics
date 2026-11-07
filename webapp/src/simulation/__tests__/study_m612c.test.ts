// M6-12 candidate C (quintic Hermite) range-bound derivation study
// (skipped unless PHYSI_BENCH=1). Research only: nothing here is wired
// into the simulator. Goal: conservative per-DOF ranges of the quintic
// q(s) and qdd(s) over a segment — the missing ingredient every
// existing safety proof (non-contact, seated/flying sleep, detach
// bounds) would need before C could ever be connected to sand.
//
// Method: Bernstein form. A degree-n polynomial on [0,1] lies inside
// [min, max] of its Bernstein coefficients (convex-hull property,
// exact arithmetic). de Casteljau subdivision at s=1/2 tightens the
// inclusion quadratically. Floating-point widening hair is applied
// the same way as the existing cubic bounds (documented, measured).
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

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

/** Power coefficients (in s) of the quintic, from the H-basis rows. */
function quinticPower(q0: number, qd0: number, a0: number, q1: number, qd1: number, a1: number, h: number): number[] {
  // basis power rows for [H0,H1,H2,H3,H4,H5]
  const rows = [
    [1, 0, 0, -10, 15, -6],
    [0, 1, 0, -6, 8, -3],
    [0, 0, 0.5, -1.5, 1.5, -0.5],
    [0, 0, 0, 10, -15, 6],
    [0, 0, 0, -4, 7, -3],
    [0, 0, 0, 0.5, -1, 0.5],
  ]
  const w = [q0, h * qd0, h * h * a0, q1, h * qd1, h * h * a1]
  const c = [0, 0, 0, 0, 0, 0]
  for (let b = 0; b < 6; b++) for (let i = 0; i < 6; i++) c[i] += w[b] * rows[b][i]
  return c
}

function binom(n: number, k: number): number {
  let r = 1
  for (let i = 0; i < k; i++) r = (r * (n - i)) / (i + 1)
  return r
}

/** Bernstein coefficients of a power-basis polynomial on [0,1]. */
function toBernstein(c: number[]): number[] {
  const n = c.length - 1
  const b = new Array(n + 1).fill(0)
  for (let k = 0; k <= n; k++) {
    let s = 0
    for (let i = 0; i <= k; i++) s += (binom(k, i) / binom(n, i)) * c[i]
    b[k] = s
  }
  return b
}

function splitBernstein(b: number[]): [number[], number[]] {
  const n = b.length - 1
  const left: number[] = [b[0]]; const rightRev: number[] = [b[n]]
  let cur = b.slice()
  for (let r = 1; r <= n; r++) {
    const nxt: number[] = []
    for (let i = 0; i < cur.length - 1; i++) nxt.push((cur[i] + cur[i + 1]) / 2)
    left.push(nxt[0]); rightRev.push(nxt[nxt.length - 1])
    cur = nxt
  }
  return [left, rightRev.reverse()]
}

/** Conservative [lo,hi] of a power-basis polynomial on [0,1] via
 *  Bernstein inclusion with `depth` bisection levels. */
function polyRange(c: number[], depth: number): [number, number] {
  let polys = [toBernstein(c)]
  for (let d = 0; d < depth; d++) {
    const next: number[][] = []
    for (const b of polys) { const [l, r] = splitBernstein(b); next.push(l, r) }
    polys = next
  }
  let lo = Infinity; let hi = -Infinity
  for (const b of polys) for (const v of b) { lo = Math.min(lo, v); hi = Math.max(hi, v) }
  return [lo, hi]
}

benchIt('M6-12C quintic range bounds: inclusion + tightness on real driven segments', () => {
  let seed = 777
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  for (const freq of [152, 649, 1241.6]) {
    const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), [0.03, -0.02])
    const signal = plate.signal
    const forceAt = (t: number) => signal.forcePeakAt(t - plate.epoch) * signal.u(t - plate.epoch)
    const odeA = (j: number, q: number, qd: number, t: number) =>
      (plate.phiDrive[j] * forceAt(t)) / plate.mass[j] - 2 * plate.zeta * plate.omega[j] * qd - plate.omega[j] ** 2 * q
    const dt = 1 / (44 * freq)
    // Random nonzero modal excitation so every DOF is live.
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = (rnd() - 0.5) * 2e-4; plate.qd[j] = (rnd() - 0.5) * 2 }
    type Seg = { t0: number; t1: number; q0: Float64Array; qd0: Float64Array; q1: Float64Array; qd1: Float64Array }
    let checked = 0; let violations = 0
    const slackQ: Record<number, number> = { 0: 0, 2: 0, 4: 0 }
    const slackA: Record<number, number> = { 0: 0, 2: 0, 4: 0 }
    let minSlackQ4 = Infinity; let minSlackA4 = Infinity
    for (let k = 0; k < 200; k++) {
      plate.beginStep(k * dt, (k + 1) * dt)
      const sg = (plate as unknown as { seg: Seg }).seg
      const h = sg.t1 - sg.t0
      for (let j = 0; j < plate.ndof; j++) {
        const a0 = odeA(j, sg.q0[j], sg.qd0[j], sg.t0)
        const a1 = odeA(j, sg.q1[j], sg.qd1[j], sg.t1)
        const c = quinticPower(sg.q0[j], sg.qd0[j], a0, sg.q1[j], sg.qd1[j], a1, h)
        // qdd(s)/h^2 power coeffs: [2c2, 6c3, 12c4, 20c5] / h^2
        const ca = [2 * c[2] / (h * h), 6 * c[3] / (h * h), 12 * c[4] / (h * h), 20 * c[5] / (h * h)]
        // Dense truth.
        let qLo = Infinity; let qHi = -Infinity; let aLo = Infinity; let aHi = -Infinity
        const NS = 1000
        for (let si = 0; si <= NS; si++) {
          const [q, , qdd] = hermite5(sg.q0[j], sg.qd0[j], a0, sg.q1[j], sg.qd1[j], a1, h, si / NS)
          qLo = Math.min(qLo, q); qHi = Math.max(qHi, q)
          aLo = Math.min(aLo, qdd); aHi = Math.max(aHi, qdd)
        }
        for (const depth of [0, 2, 4]) {
          const [lo, hi] = polyRange(c, depth)
          const [alo, ahi] = polyRange(ca, depth)
          const hairQ = 1e-9 * (1 + Math.abs(lo) + Math.abs(hi))
          const hairA = 1e-9 * (1 + Math.abs(alo) + Math.abs(ahi))
          if (qLo < lo - hairQ || qHi > hi + hairQ) violations++
          if (aLo < alo - hairA || aHi > ahi + hairA) violations++
          const trueW = qHi - qLo; const trueWA = aHi - aLo
          if (trueW > 0) slackQ[depth] = Math.max(slackQ[depth], (hi - lo) / trueW)
          if (trueWA > 0) slackA[depth] = Math.max(slackA[depth], (ahi - alo) / trueWA)
          if (depth === 4) {
            if (trueW > 0) minSlackQ4 = Math.min(minSlackQ4, (hi - lo) / trueW)
            if (trueWA > 0) minSlackA4 = Math.min(minSlackA4, (ahi - alo) / trueWA)
          }
        }
        checked++
      }
    }
    console.log(`M612C f=${freq}: segments*DOF=${checked} inclusion violations=${violations} (hair 1e-9 rel) maxWidthRatio q[d0=${slackQ[0].toFixed(2)} d2=${slackQ[2].toFixed(3)} d4=${slackQ[4].toFixed(4)}] qdd[d0=${slackA[0].toFixed(2)} d2=${slackA[2].toFixed(3)} d4=${slackA[4].toFixed(4)}] minWidthRatio@d4 q=${minSlackQ4.toFixed(4)} qdd=${minSlackA4.toFixed(4)}`)
  }
}, 1_800_000)
