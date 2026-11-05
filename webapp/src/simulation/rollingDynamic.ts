// M6-17 RESEARCH ONLY — not wired into TransientSandSim.
// Rolling of a rigid sphere on a TIME-DEPENDENT curved surface
// z = h(x,y,t) with Coulomb friction (docs/m6_17_design.md).
// Extends the static M6-16 formulation: the offset center
// C = S + R*n needs partials in (x,y,t), i.e. spatial derivatives
// of h to 3rd order, of ht to 2nd, of htt to 1st.
//
// PUT-IN: sphere (I = beta m R^2), graph surface whose material
// points move vertically (0,0,ht), Coulomb friction, no spinning
// friction.
// EMERGED: N = n.A0 + g*n_z with A0 the known part of the center
// acceleration (normal equation still decouples); rolling from
// the coupled 4x4 system; slip u = J qdot + R n_t - R (w x n).
// claim-tier: established (classical rigid-body mechanics).
import { G } from './surface'
import type { RollingParams, RollingState, RollingMode } from './rollingSphere'

export interface SurfaceFieldDyn {
  h: number
  hx: number; hy: number
  hxx: number; hxy: number; hyy: number
  hxxx: number; hxxy: number; hxyy: number; hyyy: number
  ht: number; htt: number
  htx: number; hty: number
  htxx: number; htxy: number; htyy: number
  httx: number; htty: number
}
export type SurfaceFnDyn = (x: number, y: number, t: number, out: SurfaceFieldDyn) => void

export interface DynRollingParams extends RollingParams {
  /** M6-25: Kirchhoff top-surface offset H/2 (m). Material
   *  surface point S = (x - ht2 hx, y - ht2 hy, h + ht2) with
   *  the mid-surface graph normal (linearized plate theory;
   *  docs/m6_25_design.md). 0/undefined = pre-M6-25 model. */
  plateHalfThickness?: number
}

export interface DynRollingAccel {
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

const sf: SurfaceFieldDyn = {
  h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0,
  hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0,
  ht: 0, htt: 0, htx: 0, hty: 0, htxx: 0, htxy: 0, htyy: 0, httx: 0, htty: 0,
}

// Derivative directions 0=x, 1=y, 2=t. D1[d] = d(hx,hy)/dd,
// D2[d1][d2] = d^2(hx,hy)/dd1 dd2 (symmetric in d1,d2 here).
function derivTables(s: SurfaceFieldDyn): { D1: [number, number][]; D2: [number, number][][] } {
  const D1: [number, number][] = [[s.hxx, s.hxy], [s.hxy, s.hyy], [s.htx, s.hty]]
  const D2: [number, number][][] = [
    [[s.hxxx, s.hxxy], [s.hxxy, s.hxyy], [s.htxx, s.htxy]],
    [[s.hxxy, s.hxyy], [s.hxyy, s.hyyy], [s.htxy, s.htyy]],
    [[s.htxx, s.htxy], [s.htxy, s.htyy], [s.httx, s.htty]],
  ]
  return { D1, D2 }
}

export interface DynCenterPartials {
  C: V3
  /** First partials [C_x, C_y, C_t]. */
  D: [V3, V3, V3]
  /** Second partials [xx, xy, yy, xt, yt, tt]. */
  H: [V3, V3, V3, V3, V3, V3]
  n: V3
  /** Normal partials [n_x, n_y, n_t]. */
  nD: [V3, V3, V3]
}

export function dynCenterPartials(
  x: number, y: number, t: number, surf: SurfaceFnDyn, radius: number,
  halfThick = 0,
): DynCenterPartials {
  surf(x, y, t, sf)
  return partialsFromField(sf, x, y, radius, halfThick)
}

function partialsFromField(
  s: SurfaceFieldDyn, x: number, y: number, R: number, ht2 = 0,
): DynCenterPartials {
  const { D1, D2 } = derivTables(s)
  const m: V3 = [-s.hx, -s.hy, 1]
  const sn = Math.hypot(m[0], m[1], m[2])
  const n: V3 = [m[0] / sn, m[1] / sn, m[2] / sn]
  const md: V3[] = [0, 1, 2].map((d): V3 => [-D1[d][0], -D1[d][1], 0])
  const snd = md.map((v) => dot(m, v) / sn)
  const nd: V3[] = [0, 1, 2].map((d): V3 => sub(scale(md[d], 1 / sn), scale(n, snd[d] / sn)))
  const pairIdx: [number, number][] = [[0, 0], [0, 1], [1, 1], [0, 2], [1, 2], [2, 2]]
  const ndd: V3[] = []
  for (const [d1, d2] of pairIdx) {
    const mdd: V3 = [-D2[d1][d2][0], -D2[d1][d2][1], 0]
    const sndd = (dot(md[d1], md[d2]) + dot(m, mdd) - snd[d1] * snd[d2]) / sn
    ndd.push(add(
      add(scale(mdd, 1 / sn), scale(md[d1], -snd[d2] / (sn * sn))),
      add(scale(nd[d2], -snd[d1] / sn),
        scale(n, -sndd / sn + (snd[d1] * snd[d2]) / (sn * sn))),
    ))
  }
  // Material surface point: mid-surface graph (ht2 = 0) or
  // Kirchhoff top surface S = (x - ht2 hx, y - ht2 hy, h + ht2).
  const S: V3 = [x - ht2 * s.hx, y - ht2 * s.hy, s.h + ht2]
  const Sd: V3[] = [
    [1 - ht2 * s.hxx, -ht2 * s.hxy, s.hx],
    [-ht2 * s.hxy, 1 - ht2 * s.hyy, s.hy],
    [-ht2 * s.htx, -ht2 * s.hty, s.ht],
  ]
  const second: [number, number, number][] = [
    // [d2(hx), d2(hy), d2(h)] per pair [xx, xy, yy, xt, yt, tt]
    [s.hxxx, s.hxxy, s.hxx], [s.hxxy, s.hxyy, s.hxy],
    [s.hxyy, s.hyyy, s.hyy], [s.htxx, s.htxy, s.htx],
    [s.htxy, s.htyy, s.hty], [s.httx, s.htty, s.htt],
  ]
  const Sdd: V3[] = second.map(([a, b, c]): V3 => [-ht2 * a, -ht2 * b, c])
  return {
    C: add(S, scale(n, R)),
    D: [
      add(Sd[0], scale(nd[0], R)),
      add(Sd[1], scale(nd[1], R)),
      add(Sd[2], scale(nd[2], R)),
    ],
    H: [
      add(Sdd[0], scale(ndd[0], R)), add(Sdd[1], scale(ndd[1], R)),
      add(Sdd[2], scale(ndd[2], R)), add(Sdd[3], scale(ndd[3], R)),
      add(Sdd[4], scale(ndd[4], R)), add(Sdd[5], scale(ndd[5], R)),
    ] as [V3, V3, V3, V3, V3, V3],
    n,
    nD: [nd[0], nd[1], nd[2]],
  }
}

interface Frame {
  cp: DynCenterPartials
  t1: V3; t2: V3
  MQ: [number, number, number]
  A0: V3
  K: V3
  nPerMass: number
}

function buildFrame(st: RollingState, surf: SurfaceFnDyn, p: DynRollingParams, t: number): Frame {
  surf(st.x, st.y, t, sf)
  const cp = partialsFromField(sf, st.x, st.y, p.radius, p.plateHalfThickness ?? 0)
  const Sx: V3 = [1, 0, sf.hx]
  const sxl = Math.hypot(Sx[0], Sx[1], Sx[2])
  const t1: V3 = [Sx[0] / sxl, Sx[1] / sxl, Sx[2] / sxl]
  const t2 = cross(cp.n, t1)
  const Hc = add(
    add(scale(cp.H[0], st.vx * st.vx), scale(cp.H[1], 2 * st.vx * st.vy)),
    scale(cp.H[2], st.vy * st.vy),
  )
  const X = add(scale(cp.H[3], st.vx), scale(cp.H[4], st.vy)) // C_xt vx + C_yt vy
  const A0 = add(add(Hc, scale(X, 2)), cp.H[5])
  const ndot = add(add(scale(cp.nD[0], st.vx), scale(cp.nD[1], st.vy)), cp.nD[2])
  const wxnd = cross([st.wx, st.wy, st.wz], ndot)
  // K0 parts: Hc + (0,0,grad ht . v) + 2R(n_xt vx + n_yt vy) + R n_tt
  // Recover n_xt etc from C partials: C_xt = (0,0,htx) + R n_xt.
  // Recover n_xt, n_yt, n_tt from C partials minus the material
  // surface partials (top-surface aware, M6-25).
  const ht2 = p.plateHalfThickness ?? 0
  const nXt: V3 = [
    (cp.H[3][0] + ht2 * sf.htxx) / p.radius,
    (cp.H[3][1] + ht2 * sf.htxy) / p.radius,
    (cp.H[3][2] - sf.htx) / p.radius,
  ]
  const nYt: V3 = [
    (cp.H[4][0] + ht2 * sf.htxy) / p.radius,
    (cp.H[4][1] + ht2 * sf.htyy) / p.radius,
    (cp.H[4][2] - sf.hty) / p.radius,
  ]
  const nTt: V3 = [
    (cp.H[5][0] + ht2 * sf.httx) / p.radius,
    (cp.H[5][1] + ht2 * sf.htty) / p.radius,
    (cp.H[5][2] - sf.htt) / p.radius,
  ]
  const K = sub(add(
    add(Hc, [0, 0, sf.htx * st.vx + sf.hty * st.vy]),
    add(scale(add(scale(nXt, st.vx), scale(nYt, st.vy)), 2 * p.radius),
      scale(nTt, p.radius)),
  ), scale(wxnd, p.radius))
  const MQ: [number, number, number] = [
    dot(cp.D[0], cp.D[0]), dot(cp.D[0], cp.D[1]), dot(cp.D[1], cp.D[1]),
  ]
  return { cp, t1, t2, MQ, A0, K, nPerMass: dot(cp.n, A0) + G * cp.n[2] }
}

function slipVec(st: RollingState, fr: Frame, p: DynRollingParams): V3 {
  const Jq = add(scale(fr.cp.D[0], st.vx), scale(fr.cp.D[1], st.vy))
  const wxn = cross([st.wx, st.wy, st.wz], fr.cp.n)
  return sub(add(Jq, scale(fr.cp.nD[2], p.radius)), scale(wxn, p.radius))
}

function qddForFriction(fr: Frame, f: V3): [number, number] {
  const { cp, MQ, A0 } = fr
  const rhs: [number, number] = [
    dot(cp.D[0], sub(add([0, 0, -G], f), A0)),
    dot(cp.D[1], sub(add([0, 0, -G], f), A0)),
  ]
  const det = MQ[0] * MQ[2] - MQ[1] * MQ[1]
  return [
    (MQ[2] * rhs[0] - MQ[1] * rhs[1]) / det,
    (-MQ[1] * rhs[0] + MQ[0] * rhs[1]) / det,
  ]
}

function slipAccelTangential(fr: Frame, p: DynRollingParams, f: V3): { qdd: [number, number]; udot1: number; udot2: number } {
  const qdd = qddForFriction(fr, f)
  const Jq = add(scale(fr.cp.D[0], qdd[0]), scale(fr.cp.D[1], qdd[1]))
  return {
    qdd,
    udot1: dot(fr.t1, Jq) + dot(fr.t1, fr.K) + dot(fr.t1, f) / p.beta,
    udot2: dot(fr.t2, Jq) + dot(fr.t2, fr.K) + dot(fr.t2, f) / p.beta,
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

function rollSolution(fr: Frame, p: DynRollingParams): { qdd: [number, number]; f: V3 } {
  const { cp, MQ, A0, K, t1, t2 } = fr
  const B = (v: V3): [number, number] => [dot(cp.D[0], v), dot(cp.D[1], v)]
  const [b1x, b1y] = B(t1); const [b2x, b2y] = B(t2)
  const grav: [number, number] = B([0, 0, -G])
  const a0c: [number, number] = B(A0)
  const sol = solve4([
    [MQ[0], MQ[1], -b1x, -b2x],
    [MQ[1], MQ[2], -b1y, -b2y],
    [dot(t1, cp.D[0]), dot(t1, cp.D[1]), 1 / p.beta, 0],
    [dot(t2, cp.D[0]), dot(t2, cp.D[1]), 0, 1 / p.beta],
  ], [grav[0] - a0c[0], grav[1] - a0c[1], -dot(t1, K), -dot(t2, K)])
  return { qdd: [sol[0], sol[1]], f: add(scale(t1, sol[2]), scale(t2, sol[3])) }
}

function kineticDyn(fr: Frame, p: DynRollingParams, dirU: V3, st: RollingState): DynRollingAccel {
  const fm = p.muK * fr.nPerMass
  const f: V3 = [-fm * dirU[0], -fm * dirU[1], -fm * dirU[2]]
  const qdd = qddForFriction(fr, f)
  const dw = spinDeriv(fr.cp.n, f, p)
  const u = slipVec(st, fr, p)
  return {
    ax: qdd[0], ay: qdd[1], dwx: dw[0], dwy: dw[1], dwz: dw[2],
    nPerMass: fr.nPerMass, slipSpeed: Math.hypot(u[0], u[1], u[2]), mode: 'slip',
  }
}

export function dynRollingAccel(
  st: RollingState, surf: SurfaceFnDyn, p: DynRollingParams, t: number,
): DynRollingAccel {
  const fr = buildFrame(st, surf, p, t)
  if (fr.nPerMass <= 0) {
    return { ax: 0, ay: 0, dwx: 0, dwy: 0, dwz: 0, nPerMass: fr.nPerMass, slipSpeed: 0, mode: 'detach' }
  }
  const u = slipVec(st, fr, p)
  const us = Math.hypot(u[0], u[1], u[2])
  if (us < 1e-12) {
    const sol = rollSolution(fr, p)
    const fMag = Math.hypot(sol.f[0], sol.f[1], sol.f[2])
    if (fMag <= p.muS * fr.nPerMass) {
      const dw = spinDeriv(fr.cp.n, sol.f, p)
      return {
        ax: sol.qdd[0], ay: sol.qdd[1], dwx: dw[0], dwy: dw[1], dwz: dw[2],
        nPerMass: fr.nPerMass, slipSpeed: 0, mode: 'roll',
      }
    }
    const sa = slipAccelTangential(fr, p, [0, 0, 0])
    const dir3 = add(scale(fr.t1, sa.udot1), scale(fr.t2, sa.udot2))
    const dm = Math.hypot(dir3[0], dir3[1], dir3[2]) || 1
    return kineticDyn(fr, p, [dir3[0] / dm, dir3[1] / dm, dir3[2] / dm], st)
  }
  return kineticDyn(fr, p, [u[0] / us, u[1] / us, u[2] / us], st)
}

/** Center velocity (for flight transitions / energy). */
export function dynCenterVelocity(st: RollingState, surf: SurfaceFnDyn, p: DynRollingParams, t: number, out: V3): void {
  const fr = buildFrame(st, surf, p, t)
  const v = add(add(scale(fr.cp.D[0], st.vx), scale(fr.cp.D[1], st.vy)), fr.cp.D[2])
  out[0] = v[0]; out[1] = v[1]; out[2] = v[2]
}

export interface DynStepResult { mode: RollingMode; advanced: number }

type Deriv = (st: RollingState, t: number) => { d: Float64Array; detach: boolean }

function rk4Dyn(s0: RollingState, t: number, h: number, deriv: Deriv): { state: RollingState; sawDetach: boolean } {
  let sawDetach = false
  const D = (st: RollingState, tt: number): Float64Array => {
    const r = deriv(st, tt)
    if (r.detach) sawDetach = true
    return r.d
  }
  const mk = (base: RollingState, k: Float64Array, f: number): RollingState => ({
    x: base.x + k[0] * f, y: base.y + k[1] * f,
    vx: base.vx + k[2] * f, vy: base.vy + k[3] * f,
    wx: base.wx + k[4] * f, wy: base.wy + k[5] * f, wz: base.wz + k[6] * f,
  })
  const k1 = D(s0, t)
  const k2 = D(mk(s0, k1, h / 2), t + h / 2)
  const k3 = D(mk(s0, k2, h / 2), t + h / 2)
  const k4 = D(mk(s0, k3, h), t + h)
  const comb = new Float64Array(7)
  for (let i = 0; i < 7; i++) comb[i] = (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
  return { state: mk(s0, comb, 1), sawDetach }
}

function rollDeriv(surf: SurfaceFnDyn, p: DynRollingParams): Deriv {
  return (st, tt) => {
    const fr = buildFrame(st, surf, p, tt)
    if (fr.nPerMass <= 0) return { d: Float64Array.from([st.vx, st.vy, 0, 0, 0, 0, 0]), detach: true }
    const sol = rollSolution(fr, p)
    const dw = spinDeriv(fr.cp.n, sol.f, p)
    return {
      d: Float64Array.from([st.vx, st.vy, sol.qdd[0], sol.qdd[1], dw[0], dw[1], dw[2]]),
      detach: false,
    }
  }
}

function slipDeriv(surf: SurfaceFnDyn, p: DynRollingParams, dirHint: V3): Deriv {
  return (st, tt) => {
    const fr = buildFrame(st, surf, p, tt)
    if (fr.nPerMass <= 0) return { d: Float64Array.from([st.vx, st.vy, 0, 0, 0, 0, 0]), detach: true }
    const u = slipVec(st, fr, p)
    const us = Math.hypot(u[0], u[1], u[2])
    let dir: V3 = dirHint
    if (us >= 1e-12) {
      const cand: V3 = [u[0] / us, u[1] / us, u[2] / us]
      if (dot(cand, dirHint) > 0) dir = cand
    }
    const a = kineticDyn(fr, p, dir, st)
    return {
      d: Float64Array.from([st.vx, st.vy, a.ax, a.ay, a.dwx, a.dwy, a.dwz]),
      detach: false,
    }
  }
}

function rollMarginAt(st: RollingState, surf: SurfaceFnDyn, p: DynRollingParams, t: number): { margin: number; nPerMass: number } {
  const fr = buildFrame(st, surf, p, t)
  if (fr.nPerMass <= 0) return { margin: -Infinity, nPerMass: fr.nPerMass }
  const sol = rollSolution(fr, p)
  const fm = Math.hypot(sol.f[0], sol.f[1], sol.f[2])
  return { margin: p.muS * fr.nPerMass - fm, nPerMass: fr.nPerMass }
}

/** Advance one step on a time-dependent curved surface with
 *  event handling (roll-break, slip-zero, detach), earliest
 *  event first. Same sampling caveat as M6-16. */
export function stepDynRolling(
  st: RollingState, surf: SurfaceFnDyn, p: DynRollingParams, t: number, dt: number,
): DynStepResult {
  let advanced = 0
  for (let iter = 0; iter < 16; iter++) {
    const rem = dt - advanced
    if (rem <= 0) break
    const tc = t + advanced
    const a0 = dynRollingAccel(st, surf, p, tc)
    if (a0.mode === 'detach') return { mode: 'detach', advanced }
    const start: RollingState = { ...st }

    if (a0.mode === 'roll') {
      const deriv = rollDeriv(surf, p)
      const trial = rk4Dyn(start, tc, rem, deriv)
      const mEnd = rollMarginAt(trial.state, surf, p, tc + rem)
      const bad = trial.sawDetach || mEnd.nPerMass <= 0 || mEnd.margin < 0
      if (!bad) {
        Object.assign(st, trial.state)
        return { mode: 'roll', advanced: dt }
      }
      let lo = 0; let hi = 1
      let detachSide = mEnd.nPerMass <= 0
      for (let it = 0; it < 50; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Dyn(start, tc, rem * mid, deriv)
        const mm = rollMarginAt(r.state, surf, p, tc + rem * mid)
        if (r.sawDetach || mm.nPerMass <= 0 || mm.margin < 0) {
          hi = mid; detachSide = r.sawDetach || mm.nPerMass <= 0
        } else lo = mid
      }
      const r = rk4Dyn(start, tc, rem * hi, deriv)
      Object.assign(st, r.state)
      advanced += rem * hi
      if (detachSide) return { mode: 'detach', advanced }
      continue
    }

    const fr0 = buildFrame(st, surf, p, tc)
    const u0 = slipVec(st, fr0, p)
    const us0 = Math.hypot(u0[0], u0[1], u0[2])
    let dirHint: V3
    if (us0 >= 1e-12) dirHint = [u0[0] / us0, u0[1] / us0, u0[2] / us0]
    else {
      const sa = slipAccelTangential(fr0, p, [0, 0, 0])
      const dir3 = add(scale(fr0.t1, sa.udot1), scale(fr0.t2, sa.udot2))
      const dm = Math.hypot(dir3[0], dir3[1], dir3[2]) || 1
      dirHint = [dir3[0] / dm, dir3[1] / dm, dir3[2] / dm]
    }
    const deriv = slipDeriv(surf, p, dirHint)
    const trial = rk4Dyn(start, tc, rem, deriv)
    const endA = dynRollingAccel(trial.state, surf, p, tc + rem)
    const detachInside = trial.sawDetach || endA.mode === 'detach'
    const signedSlip = (s2: RollingState, tt: number): number => {
      const fr = buildFrame(s2, surf, p, tt)
      return dot(slipVec(s2, fr, p), dirHint)
    }
    const zeroInside = us0 >= 1e-12 && signedSlip(trial.state, tc + rem) <= 0
    if (!detachInside && !zeroInside) {
      Object.assign(st, trial.state)
      return { mode: 'slip', advanced: dt }
    }
    let fracDetach = Infinity
    if (detachInside) {
      let lo = 0; let hi = 1
      for (let it = 0; it < 50; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Dyn(start, tc, rem * mid, deriv)
        const aa = dynRollingAccel(r.state, surf, p, tc + rem * mid)
        if (r.sawDetach || aa.mode === 'detach') hi = mid; else lo = mid
      }
      fracDetach = hi
    }
    let fracZero = Infinity
    if (zeroInside) {
      let lo = 0; let hi = 1
      for (let it = 0; it < 60; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Dyn(start, tc, rem * mid, deriv)
        if (signedSlip(r.state, tc + rem * mid) > 0) lo = mid; else hi = mid
      }
      fracZero = hi
    }
    if (fracDetach <= fracZero) {
      const r = rk4Dyn(start, tc, rem * fracDetach, deriv)
      Object.assign(st, r.state)
      return { mode: 'detach', advanced: advanced + rem * fracDetach }
    }
    const r = rk4Dyn(start, tc, rem * fracZero, deriv)
    Object.assign(st, r.state)
    const frE = buildFrame(st, surf, p, tc + rem * fracZero)
    const uE = slipVec(st, frE, p)
    if (Math.hypot(uE[0], uE[1], uE[2]) > 1e-9) {
      Object.assign(st, trial.state)
      return { mode: 'slip', advanced: dt }
    }
    dynProjectToRolling(st, surf, p, tc + rem * fracZero)
    advanced += rem * fracZero
  }
  const aEnd = dynRollingAccel(st, surf, p, t + advanced)
  return { mode: aEnd.mode, advanced }
}

export function dynProjectToRolling(st: RollingState, surf: SurfaceFnDyn, p: DynRollingParams, t: number): void {
  const fr = buildFrame(st, surf, p, t)
  const u = slipVec(st, fr, p)
  const d = cross(fr.cp.n, u)
  st.wx += d[0] / p.radius; st.wy += d[1] / p.radius; st.wz += d[2] / p.radius
}
