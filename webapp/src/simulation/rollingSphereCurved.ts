// M6-16 RESEARCH ONLY — not wired into TransientSandSim.
// Rolling of a rigid sphere on a STATIC curved surface z = h(x,y)
// with Coulomb friction (docs/m6_16_design.md). The center rides
// the offset surface C = S + R*n; its second partials need third
// derivatives of h, supplied by the SurfaceFn3 input.
//
// PUT-IN: sphere (I = beta m R^2), static graph surface, Coulomb
// friction, no spinning friction.
// EMERGED: N = n.Hc + g*n_z (normal equation decouples because
// the offset-surface Jacobian columns are tangential); rolling
// (qdd, f) from the coupled 4x4 system of tangential equations
// and slip-acceleration zero; slip onset direction = slip
// acceleration at f = 0.
// claim-tier: established (classical rigid-body mechanics).
import { G } from './surface'
import type { RollingParams, RollingState, RollingMode } from './rollingSphere'

export interface SurfaceField3 {
  h: number
  hx: number; hy: number
  hxx: number; hxy: number; hyy: number
  hxxx: number; hxxy: number; hxyy: number; hyyy: number
}
export type SurfaceFn3 = (x: number, y: number, out: SurfaceField3) => void

export interface CurvedRollingAccel {
  ax: number; ay: number
  dwx: number; dwy: number; dwz: number
  nPerMass: number
  slipSpeed: number
  mode: RollingMode
}

type V3 = [number, number, number]
const dot = (a: V3, b: V3): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
const cross = (a: V3, b: V3): V3 => [
  a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0],
]
const add = (a: V3, b: V3): V3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
const scale = (a: V3, k: number): V3 => [a[0] * k, a[1] * k, a[2] * k]
const sub = (a: V3, b: V3): V3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

const sf3: SurfaceField3 = {
  h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0,
}

export interface CenterPartials {
  C: V3
  /** Jacobian columns [C_x, C_y]. */
  J: [V3, V3]
  /** Second partials [C_xx, C_xy, C_yy]. */
  H: [V3, V3, V3]
  n: V3
}

/** Offset-center position and partials from the surface field.
 *  Exported for finite-difference validation in tests. */
export function centerPartials(x: number, y: number, surf: SurfaceFn3, radius: number): CenterPartials {
  surf(x, y, sf3)
  return partialsFromField(sf3, x, y, radius)
}

function partialsFromField(s: SurfaceField3, x: number, y: number, R: number): CenterPartials {
  const m: V3 = [-s.hx, -s.hy, 1]
  const sn = Math.hypot(m[0], m[1], m[2])
  const n: V3 = [m[0] / sn, m[1] / sn, m[2] / sn]
  // Gradient component derivatives: G[d] = d(hx,hy)/dd (Hessian rows),
  // G2[d1][d2] = d^2(hx,hy)/dd1 dd2.
  const G: [number, number][] = [[s.hxx, s.hxy], [s.hxy, s.hyy]]
  const G2: [number, number][][] = [
    [[s.hxxx, s.hxxy], [s.hxxy, s.hxyy]],
    [[s.hxxy, s.hxyy], [s.hxyy, s.hyyy]],
  ]
  const md: V3[] = [0, 1].map((d): V3 => [-G[d][0], -G[d][1], 0])
  const snd = md.map((v) => dot(m, v) / sn)
  const nd: V3[] = [0, 1].map((d): V3 => sub(scale(md[d], 1 / sn), scale(n, snd[d] / sn)))
  const ndd: V3[] = [] // [xx, xy, yy]
  const pairs: [number, number][] = [[0, 0], [0, 1], [1, 1]]
  for (const [d1, d2] of pairs) {
    const mdd: V3 = [-G2[d1][d2][0], -G2[d1][d2][1], 0]
    const sndd = (dot(md[d1], md[d2]) + dot(m, mdd) - snd[d1] * snd[d2]) / sn
    // n_d1d2 = mdd/sn - md1*snd2/sn^2 - nd2*snd1/sn - n*sndd/sn
    //           + n*snd1*snd2/sn^2
    ndd.push(add(
      add(scale(mdd, 1 / sn), scale(md[d1], -snd[d2] / (sn * sn))),
      add(scale(nd[d2], -snd[d1] / sn),
        scale(n, -sndd / sn + (snd[d1] * snd[d2]) / (sn * sn))),
    ))
  }
  const S: V3 = [x, y, s.h]
  const Sd: V3[] = [[1, 0, s.hx], [0, 1, s.hy]]
  const Sdd: V3[] = [[0, 0, s.hxx], [0, 0, s.hxy], [0, 0, s.hyy]]
  return {
    C: add(S, scale(n, R)),
    J: [add(Sd[0], scale(nd[0], R)), add(Sd[1], scale(nd[1], R))],
    H: [
      add(Sdd[0], scale(ndd[0], R)),
      add(Sdd[1], scale(ndd[1], R)),
      add(Sdd[2], scale(ndd[2], R)),
    ],
    n,
  }
}

interface Frame {
  cp: CenterPartials
  t1: V3; t2: V3
  /** Hc(v,v) vector. */
  Hc: V3
  ndot: V3
  MQ: [number, number, number] // J^T J symmetric entries
  nPerMass: number
}

function buildFrame(st: RollingState, surf: SurfaceFn3, p: RollingParams): Frame {
  surf(st.x, st.y, sf3)
  const cp = partialsFromField(sf3, st.x, st.y, p.radius)
  const Sx: V3 = [1, 0, sf3.hx]
  const sxl = Math.hypot(Sx[0], Sx[1], Sx[2])
  const t1: V3 = [Sx[0] / sxl, Sx[1] / sxl, Sx[2] / sxl]
  const t2 = cross(cp.n, t1)
  const Hc = add(
    add(scale(cp.H[0], st.vx * st.vx), scale(cp.H[1], 2 * st.vx * st.vy)),
    scale(cp.H[2], st.vy * st.vy),
  )
  // ndot recovered from Cp partials: nd = (J - Sd)/R per column.
  const ndx: V3 = [
    (cp.J[0][0] - 1) / p.radius, cp.J[0][1] / p.radius,
    (cp.J[0][2] - sf3.hx) / p.radius,
  ]
  const ndy: V3 = [
    cp.J[1][0] / p.radius, (cp.J[1][1] - 1) / p.radius,
    (cp.J[1][2] - sf3.hy) / p.radius,
  ]
  const ndot = add(scale(ndx, st.vx), scale(ndy, st.vy))
  const MQ: [number, number, number] = [
    dot(cp.J[0], cp.J[0]), dot(cp.J[0], cp.J[1]), dot(cp.J[1], cp.J[1]),
  ]
  return { cp, t1, t2, Hc, ndot, MQ, nPerMass: dot(cp.n, Hc) + G * cp.n[2] }
}

function slipVec(st: RollingState, fr: Frame, p: RollingParams): V3 {
  const Cdot = add(scale(fr.cp.J[0], st.vx), scale(fr.cp.J[1], st.vy))
  const wxn = cross([st.wx, st.wy, st.wz], fr.cp.n)
  return sub(Cdot, scale(wxn, p.radius))
}

/** Tangential slip acceleration given friction vector f (uses the
 *  tangential equations to eliminate qdd). */
function slipAccelTangential(
  st: RollingState, fr: Frame, p: RollingParams, f: V3,
): { qdd: [number, number]; udot1: number; udot2: number } {
  const { cp, Hc, MQ } = fr
  const rhs: [number, number] = [
    dot(cp.J[0], add([0, 0, -G], f)) - dot(cp.J[0], Hc),
    dot(cp.J[1], add([0, 0, -G], f)) - dot(cp.J[1], Hc),
  ]
  const det = MQ[0] * MQ[2] - MQ[1] * MQ[1]
  const qdd: [number, number] = [
    (MQ[2] * rhs[0] - MQ[1] * rhs[1]) / det,
    (-MQ[1] * rhs[0] + MQ[0] * rhs[1]) / det,
  ]
  // udot_t = J qdd + Hc + f/beta - R (w x ndot), tangential part:
  // J qdd and f are tangential; project Hc and the last term.
  const wxnd = cross([st.wx, st.wy, st.wz], fr.ndot)
  const rest = sub(Hc, scale(wxnd, p.radius))
  const Jq = add(scale(cp.J[0], qdd[0]), scale(cp.J[1], qdd[1]))
  return {
    qdd,
    udot1: dot(fr.t1, Jq) + dot(fr.t1, rest) + dot(fr.t1, f) / p.beta,
    udot2: dot(fr.t2, Jq) + dot(fr.t2, rest) + dot(fr.t2, f) / p.beta,
  }
}

function spinDeriv(n: V3, f: V3, p: RollingParams): V3 {
  const c = cross(n, f)
  const k = -1 / (p.beta * p.radius)
  return [k * c[0], k * c[1], k * c[2]]
}

function solve4(A: number[][], b: number[]): number[] {
  const M = A.map((row, i) => [...row, b[i]])
  for (let c = 0; c < 4; c++) {
    let piv = c
    for (let r = c + 1; r < 4; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r
    const tmp = M[c]; M[c] = M[piv]; M[piv] = tmp
    const d = M[c][c]
    for (let r = 0; r < 4; r++) {
      if (r === c) continue
      const f = M[r][c] / d
      for (let k = c; k < 5; k++) M[r][k] -= f * M[c][k]
    }
  }
  return [M[0][4] / M[0][0], M[1][4] / M[1][1], M[2][4] / M[2][2], M[3][4] / M[3][3]]
}

/** Rolling solution (qdd, f) at the current state. */
function rollSolution(
  st: RollingState, fr: Frame, p: RollingParams,
): { qdd: [number, number]; f: V3 } {
  const { cp, Hc, MQ, t1, t2 } = fr
  const wxnd = cross([st.wx, st.wy, st.wz], fr.ndot)
  const rest = sub(Hc, scale(wxnd, p.radius))
  const B = (v: V3): [number, number] => [dot(cp.J[0], v), dot(cp.J[1], v)]
  const [b1x, b1y] = B(t1); const [b2x, b2y] = B(t2)
  const grav: [number, number] = B([0, 0, -G])
  const hc: [number, number] = B(Hc)
  // Rows: [Mq | -J^T T] (qdd,fc) = grav - J^T Hc ; [T^T J | I/beta] = -T^T rest
  const sol = solve4([
    [MQ[0], MQ[1], -b1x, -b2x],
    [MQ[1], MQ[2], -b1y, -b2y],
    [dot(t1, cp.J[0]), dot(t1, cp.J[1]), 1 / p.beta, 0],
    [dot(t2, cp.J[0]), dot(t2, cp.J[1]), 0, 1 / p.beta],
  ], [grav[0] - hc[0], grav[1] - hc[1], -dot(t1, rest), -dot(t2, rest)])
  return {
    qdd: [sol[0], sol[1]],
    f: add(scale(t1, sol[2]), scale(t2, sol[3])),
  }
}

export function curvedRollingAccel(
  st: RollingState, surf: SurfaceFn3, p: RollingParams,
): CurvedRollingAccel {
  const fr = buildFrame(st, surf, p)
  const nPerMass = fr.nPerMass
  if (nPerMass <= 0) {
    return { ax: 0, ay: 0, dwx: 0, dwy: 0, dwz: 0, nPerMass, slipSpeed: 0, mode: 'detach' }
  }
  const u = slipVec(st, fr, p)
  const us = Math.hypot(u[0], u[1], u[2])
  if (us < 1e-12) {
    const sol = rollSolution(st, fr, p)
    const fMag = Math.hypot(sol.f[0], sol.f[1], sol.f[2])
    if (fMag <= p.muS * nPerMass) {
      const dw = spinDeriv(fr.cp.n, sol.f, p)
      return {
        ax: sol.qdd[0], ay: sol.qdd[1], dwx: dw[0], dwy: dw[1], dwz: dw[2],
        nPerMass, slipSpeed: 0, mode: 'roll',
      }
    }
    // Impending slip direction: slip acceleration at f = 0.
    const sa = slipAccelTangential(st, fr, p, [0, 0, 0])
    const dir3 = add(scale(fr.t1, sa.udot1), scale(fr.t2, sa.udot2))
    const dm = Math.hypot(dir3[0], dir3[1], dir3[2]) || 1
    const dirU: V3 = [dir3[0] / dm, dir3[1] / dm, dir3[2] / dm]
    return kineticCurved(st, fr, p, dirU)
  }
  const dirU: V3 = [u[0] / us, u[1] / us, u[2] / us]
  return kineticCurved(st, fr, p, dirU)
}

function kineticCurved(
  st: RollingState, fr: Frame, p: RollingParams, dirU: V3,
): CurvedRollingAccel {
  const fm = p.muK * fr.nPerMass
  const f: V3 = [-fm * dirU[0], -fm * dirU[1], -fm * dirU[2]]
  const sa = slipAccelTangential(st, fr, p, f)
  const dw = spinDeriv(fr.cp.n, f, p)
  const u = slipVec(st, fr, p)
  return {
    ax: sa.qdd[0], ay: sa.qdd[1], dwx: dw[0], dwy: dw[1], dwz: dw[2],
    nPerMass: fr.nPerMass, slipSpeed: Math.hypot(u[0], u[1], u[2]), mode: 'slip',
  }
}

/** Kinetic energy per unit mass on a static surface. */
export function curvedKineticEnergy(st: RollingState, surf: SurfaceFn3, p: RollingParams): number {
  const fr = buildFrame(st, surf, p)
  const Cdot = add(scale(fr.cp.J[0], st.vx), scale(fr.cp.J[1], st.vy))
  const w2 = st.wx * st.wx + st.wy * st.wy + st.wz * st.wz
  return 0.5 * dot(Cdot, Cdot) + 0.5 * p.beta * p.radius * p.radius * w2
}

export function curvedCenterHeight(st: RollingState, surf: SurfaceFn3, p: RollingParams): number {
  return centerPartials(st.x, st.y, surf, p.radius).C[2]
}

export interface CurvedStepResult { mode: RollingMode; advanced: number }

type Deriv = (st: RollingState) => { d: Float64Array; detach: boolean }

function rk4Curved(s0: RollingState, h: number, deriv: Deriv): { state: RollingState; sawDetach: boolean } {
  let sawDetach = false
  const D = (st: RollingState): Float64Array => {
    const r = deriv(st)
    if (r.detach) sawDetach = true
    return r.d
  }
  const mk = (base: RollingState, k: Float64Array, f: number): RollingState => ({
    x: base.x + k[0] * f, y: base.y + k[1] * f,
    vx: base.vx + k[2] * f, vy: base.vy + k[3] * f,
    wx: base.wx + k[4] * f, wy: base.wy + k[5] * f, wz: base.wz + k[6] * f,
  })
  const k1 = D(s0)
  const k2 = D(mk(s0, k1, h / 2))
  const k3 = D(mk(s0, k2, h / 2))
  const k4 = D(mk(s0, k3, h))
  const comb = new Float64Array(7)
  for (let i = 0; i < 7; i++) comb[i] = (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
  return { state: mk(s0, comb, 1), sawDetach }
}

function rollDeriv(surf: SurfaceFn3, p: RollingParams): Deriv {
  return (st) => {
    const fr = buildFrame(st, surf, p)
    if (fr.nPerMass <= 0) {
      return { d: Float64Array.from([st.vx, st.vy, 0, 0, 0, 0, 0]), detach: true }
    }
    const sol = rollSolution(st, fr, p)
    const dw = spinDeriv(fr.cp.n, sol.f, p)
    return {
      d: Float64Array.from([st.vx, st.vy, sol.qdd[0], sol.qdd[1], dw[0], dw[1], dw[2]]),
      detach: false,
    }
  }
}

function slipDeriv(surf: SurfaceFn3, p: RollingParams, dirHint: V3): Deriv {
  return (st) => {
    const fr = buildFrame(st, surf, p)
    if (fr.nPerMass <= 0) {
      return { d: Float64Array.from([st.vx, st.vy, 0, 0, 0, 0, 0]), detach: true }
    }
    const u = slipVec(st, fr, p)
    const us = Math.hypot(u[0], u[1], u[2])
    let dir: V3 = dirHint
    if (us >= 1e-12) {
      const cand: V3 = [u[0] / us, u[1] / us, u[2] / us]
      if (dot(cand, dirHint) > 0) dir = cand
    }
    const a = kineticCurved(st, fr, p, dir)
    return {
      d: Float64Array.from([st.vx, st.vy, a.ax, a.ay, a.dwx, a.dwy, a.dwz]),
      detach: false,
    }
  }
}

function rollMarginAt(st: RollingState, surf: SurfaceFn3, p: RollingParams): { margin: number; nPerMass: number } {
  const fr = buildFrame(st, surf, p)
  if (fr.nPerMass <= 0) return { margin: -Infinity, nPerMass: fr.nPerMass }
  const sol = rollSolution(st, fr, p)
  const fm = Math.hypot(sol.f[0], sol.f[1], sol.f[2])
  return { margin: p.muS * fr.nPerMass - fm, nPerMass: fr.nPerMass }
}

/** Advance one step on a static curved surface with event
 *  handling (roll-break, slip-zero, detach), earliest event first.
 *  Roll-break detection evaluates the stick margin at integrated
 *  trial states (bisection on the step fraction) — a sampling
 *  check along the trajectory, not a proof. */
export function stepCurvedRolling(
  st: RollingState, surf: SurfaceFn3, p: RollingParams, dt: number,
): CurvedStepResult {
  let advanced = 0
  for (let iter = 0; iter < 16; iter++) {
    const rem = dt - advanced
    if (rem <= 0) break
    const a0 = curvedRollingAccel(st, surf, p)
    if (a0.mode === 'detach') return { mode: 'detach', advanced }
    const start: RollingState = { ...st }

    if (a0.mode === 'roll') {
      const deriv = rollDeriv(surf, p)
      const trial = rk4Curved(start, rem, deriv)
      const mEnd = rollMarginAt(trial.state, surf, p)
      const bad = trial.sawDetach || mEnd.nPerMass <= 0 || mEnd.margin < 0
      if (!bad) {
        Object.assign(st, trial.state)
        return { mode: 'roll', advanced: dt }
      }
      let lo = 0; let hi = 1
      let detachSide = mEnd.nPerMass <= 0
      for (let it = 0; it < 50; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Curved(start, rem * mid, deriv)
        const mm = rollMarginAt(r.state, surf, p)
        if (r.sawDetach || mm.nPerMass <= 0 || mm.margin < 0) {
          hi = mid; detachSide = r.sawDetach || mm.nPerMass <= 0
        } else lo = mid
      }
      const r = rk4Curved(start, rem * hi, deriv)
      Object.assign(st, r.state)
      advanced += rem * hi
      if (detachSide) return { mode: 'detach', advanced }
      continue
    }

    // Slip.
    const fr0 = buildFrame(st, surf, p)
    const u0 = slipVec(st, fr0, p)
    const us0 = Math.hypot(u0[0], u0[1], u0[2])
    let dirHint: V3
    if (us0 >= 1e-12) dirHint = [u0[0] / us0, u0[1] / us0, u0[2] / us0]
    else {
      const sa = slipAccelTangential(st, fr0, p, [0, 0, 0])
      const dir3 = add(scale(fr0.t1, sa.udot1), scale(fr0.t2, sa.udot2))
      const dm = Math.hypot(dir3[0], dir3[1], dir3[2]) || 1
      dirHint = [dir3[0] / dm, dir3[1] / dm, dir3[2] / dm]
    }
    const deriv = slipDeriv(surf, p, dirHint)
    const trial = rk4Curved(start, rem, deriv)
    const endA = curvedRollingAccel(trial.state, surf, p)
    const detachInside = trial.sawDetach || endA.mode === 'detach'
    const signedSlip = (s2: RollingState): number => {
      const fr = buildFrame(s2, surf, p)
      return dot(slipVec(s2, fr, p), dirHint)
    }
    const zeroInside = us0 >= 1e-12 && signedSlip(trial.state) <= 0
    if (!detachInside && !zeroInside) {
      Object.assign(st, trial.state)
      return { mode: 'slip', advanced: dt }
    }
    let fracDetach = Infinity
    if (detachInside) {
      let lo = 0; let hi = 1
      for (let it = 0; it < 50; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Curved(start, rem * mid, deriv)
        const aa = curvedRollingAccel(r.state, surf, p)
        if (r.sawDetach || aa.mode === 'detach') hi = mid; else lo = mid
      }
      fracDetach = hi
    }
    let fracZero = Infinity
    if (zeroInside) {
      let lo = 0; let hi = 1
      for (let it = 0; it < 60; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Curved(start, rem * mid, deriv)
        if (signedSlip(r.state) > 0) lo = mid; else hi = mid
      }
      fracZero = hi
    }
    if (fracDetach <= fracZero) {
      const r = rk4Curved(start, rem * fracDetach, deriv)
      Object.assign(st, r.state)
      return { mode: 'detach', advanced: advanced + rem * fracDetach }
    }
    const r = rk4Curved(start, rem * fracZero, deriv)
    Object.assign(st, r.state)
    const frE = buildFrame(st, surf, p)
    const uE = slipVec(st, frE, p)
    if (Math.hypot(uE[0], uE[1], uE[2]) > 1e-9) {
      Object.assign(st, trial.state)
      return { mode: 'slip', advanced: dt }
    }
    curvedProjectToRolling(st, surf, p)
    advanced += rem * fracZero
  }
  const aEnd = curvedRollingAccel(st, surf, p)
  return { mode: aEnd.mode, advanced }
}

/** Project omega so slip is exactly zero (keeps normal spin). */
export function curvedProjectToRolling(st: RollingState, surf: SurfaceFn3, p: RollingParams): void {
  const fr = buildFrame(st, surf, p)
  const u = slipVec(st, fr, p)
  const d = cross(fr.cp.n, u)
  st.wx += d[0] / p.radius; st.wy += d[1] / p.radius; st.wz += d[2] / p.radius
}
