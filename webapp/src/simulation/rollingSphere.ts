// M6-15 RESEARCH ONLY — not wired into TransientSandSim.
// Rolling of a rigid sphere (radius R, I = beta*m*R^2 per unit mass)
// on a PLANE surface h(x,y,t) = hx*x + hy*y + s(t) with Coulomb
// friction. Validated standalone against analytic results
// (docs/m6_15_design.md): incline acceleration g*sin/(1+beta),
// stick boundary tan = muS*(1+beta)/beta, flat slip-to-roll
// transition v0/(1+beta), energy conservation.
//
// PUT-IN: sphere + plane + Coulomb friction; plate material points
// move vertically (0,0,ht). Curved or time-varying-slope surfaces
// need third derivatives of h for the offset-surface constraint
// and are rejected (throw) rather than approximated.
// EMERGED: no-slip acceleration a = -(g+htt)/(1+beta) M^-1 grad h
// from M a = Qf - grad h (g+htt) with f = -beta * dv_t/dt;
// N = (g+htt)/sqrt(1+|grad h|^2); spin about the normal is
// conserved (no spinning friction in this model).
// claim-tier: established (classical rigid-body mechanics).
import { G } from './surface'
import type { SurfaceField, SurfaceFn, FrictionParams } from './seatedFriction'

export interface RollingParams extends FrictionParams {
  radius: number
  /** I/(m R^2): 2/5 solid sphere, 2/3 hollow shell. */
  beta: number
}

export interface RollingState {
  x: number; y: number; vx: number; vy: number
  /** Angular velocity (rad/s), world frame. */
  wx: number; wy: number; wz: number
}

export type RollingMode = 'roll' | 'slip' | 'detach'

export interface RollingAccel {
  ax: number; ay: number
  dwx: number; dwy: number; dwz: number
  nPerMass: number
  /** Tangential slip speed |u| (0 in roll mode). */
  slipSpeed: number
  mode: RollingMode
}

const sf: SurfaceField = {
  h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0,
}

function checkPlane(s: SurfaceField): void {
  if (s.hxx !== 0 || s.hxy !== 0 || s.hyy !== 0 || s.htx !== 0 || s.hty !== 0) {
    throw new Error('rollingSphere: only plane surfaces (constant gradient) are supported')
  }
}

type V3 = [number, number, number]
const cross = (a: V3, b: V3): V3 => [
  a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0],
]

function frame(s: SurfaceField): { n: V3; sn: number; e1: V3; e2: V3 } {
  const sn = Math.sqrt(1 + s.hx * s.hx + s.hy * s.hy)
  return {
    n: [-s.hx / sn, -s.hy / sn, 1 / sn], sn,
    e1: [1, 0, s.hx], e2: [0, 1, s.hy],
  }
}

/** Solve M a = rhs with M = I + grad grad^T (det = sn^2). */
function solveM(s: SurfaceField, sn: number, rx: number, ry: number): [number, number] {
  const d = sn * sn
  return [
    ((1 + s.hy * s.hy) * rx - s.hx * s.hy * ry) / d,
    (-s.hx * s.hy * rx + (1 + s.hx * s.hx) * ry) / d,
  ]
}

/** Slip velocity u = v_t - R (w x n) as a 3-vector (tangential). */
export function slipVelocity(st: RollingState, surf: SurfaceFn, p: RollingParams, t: number, out: V3): void {
  surf(st.x, st.y, t, sf); checkPlane(sf)
  const { n, e1, e2 } = frame(sf)
  const wxn = cross([st.wx, st.wy, st.wz], n)
  for (let i = 0; i < 3; i++) out[i] = st.vx * e1[i] + st.vy * e2[i] - p.radius * wxn[i]
}

/** Required static-friction magnitude and N for rolling at time t
 *  (plane: depends on t only through htt). */
function rollRequirement(s: SurfaceField, p: RollingParams): { fReq: V3; nPerMass: number; ax: number; ay: number } {
  const { sn, e1, e2 } = frame(s)
  const nPerMass = (G + s.htt) / sn
  const [gx, gy] = solveM(s, sn, s.hx, s.hy)
  const c = -(G + s.htt) / (1 + p.beta)
  const ax = c * gx; const ay = c * gy
  const fReq: V3 = [0, 0, 0]
  for (let i = 0; i < 3; i++) fReq[i] = -p.beta * (ax * e1[i] + ay * e2[i])
  return { fReq, nPerMass, ax, ay }
}

function accelWithFriction(s: SurfaceField, f: V3): [number, number] {
  const { sn, e1, e2 } = frame(s)
  const qx = e1[0] * f[0] + e1[1] * f[1] + e1[2] * f[2]
  const qy = e2[0] * f[0] + e2[1] * f[1] + e2[2] * f[2]
  return solveM(s, sn, qx - s.hx * (G + s.htt), qy - s.hy * (G + s.htt))
}

function spinDeriv(n: V3, f: V3, p: RollingParams): V3 {
  const c = cross(n, f)
  const k = -1 / (p.beta * p.radius)
  return [k * c[0], k * c[1], k * c[2]]
}

/** Acceleration under a prescribed kinetic-friction direction
 *  (unit tangential vector dirU = slip direction). */
function kineticRollingAccel(
  st: RollingState, surf: SurfaceFn, p: RollingParams, t: number, dirU: V3,
): RollingAccel {
  surf(st.x, st.y, t, sf); checkPlane(sf)
  const s = sf
  const { n, sn } = frame(s)
  const nPerMass = (G + s.htt) / sn
  if (nPerMass <= 0) return { ax: 0, ay: 0, dwx: 0, dwy: 0, dwz: 0, nPerMass, slipSpeed: 0, mode: 'detach' }
  const u: V3 = [0, 0, 0]
  slipVelocity(st, surf, p, t, u)
  const fm = p.muK * nPerMass
  const f: V3 = [-fm * dirU[0], -fm * dirU[1], -fm * dirU[2]]
  const [ax, ay] = accelWithFriction(s, f)
  const dw = spinDeriv(n, f, p)
  return {
    ax, ay, dwx: dw[0], dwy: dw[1], dwz: dw[2],
    nPerMass, slipSpeed: Math.hypot(u[0], u[1], u[2]), mode: 'slip',
  }
}

/** Pointwise mode + acceleration at the current state. */
export function rollingAccel(
  st: RollingState, surf: SurfaceFn, p: RollingParams, t: number,
): RollingAccel {
  surf(st.x, st.y, t, sf); checkPlane(sf)
  const s = sf
  const { n } = frame(s)
  const req = rollRequirement(s, p)
  if (req.nPerMass <= 0) {
    return { ax: 0, ay: 0, dwx: 0, dwy: 0, dwz: 0, nPerMass: req.nPerMass, slipSpeed: 0, mode: 'detach' }
  }
  const u: V3 = [0, 0, 0]
  slipVelocity(st, surf, p, t, u)
  const us = Math.hypot(u[0], u[1], u[2])
  if (us < 1e-12) {
    const fMag = Math.hypot(req.fReq[0], req.fReq[1], req.fReq[2])
    if (fMag <= p.muS * req.nPerMass) {
      const dw = spinDeriv(n, req.fReq, p)
      return {
        ax: req.ax, ay: req.ay, dwx: dw[0], dwy: dw[1], dwz: dw[2],
        nPerMass: req.nPerMass, slipSpeed: 0, mode: 'roll',
      }
    }
    // Impending slip: slip grows opposite the required static
    // friction, i.e. along -fReq (the would-be acceleration).
    const dirU: V3 = [-req.fReq[0] / fMag, -req.fReq[1] / fMag, -req.fReq[2] / fMag]
    return kineticRollingAccel(st, surf, p, t, dirU)
  }
  const dirU: V3 = [u[0] / us, u[1] / us, u[2] / us]
  return kineticRollingAccel(st, surf, p, t, dirU)
}

/** Project omega so the slip velocity is exactly zero (preserves
 *  the normal spin component). */
export function projectToRolling(st: RollingState, surf: SurfaceFn, p: RollingParams, t: number): void {
  surf(st.x, st.y, t, sf); checkPlane(sf)
  const { n } = frame(sf)
  const u: V3 = [0, 0, 0]
  slipVelocity(st, surf, p, t, u)
  const d = cross(n, u)
  st.wx += d[0] / p.radius; st.wy += d[1] / p.radius; st.wz += d[2] / p.radius
}

/** Kinetic energy per unit mass (translation of center + spin). */
export function kineticEnergyPerMass(st: RollingState, surf: SurfaceFn, p: RollingParams, t: number): number {
  surf(st.x, st.y, t, sf); checkPlane(sf)
  const vzC = sf.hx * st.vx + sf.hy * st.vy + sf.ht
  const v2 = st.vx * st.vx + st.vy * st.vy + vzC * vzC
  const w2 = st.wx * st.wx + st.wy * st.wy + st.wz * st.wz
  return 0.5 * v2 + 0.5 * p.beta * p.radius * p.radius * w2
}

export interface RollingStepResult { mode: RollingMode; advanced: number }

type Deriv = (st: RollingState, t: number) => { d: Float64Array; detach: boolean }

function rk4Rolling(
  s0: RollingState, t: number, h: number, deriv: Deriv,
): { state: RollingState; sawDetach: boolean } {
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

function rollDeriv(surf: SurfaceFn, p: RollingParams): Deriv {
  return (st, tt) => {
    surf(st.x, st.y, tt, sf); checkPlane(sf)
    const req = rollRequirement(sf, p)
    const { n } = frame(sf)
    const dw = spinDeriv(n, req.fReq, p)
    return {
      d: Float64Array.from([st.vx, st.vy, req.ax, req.ay, dw[0], dw[1], dw[2]]),
      detach: req.nPerMass <= 0,
    }
  }
}

function slipDeriv(surf: SurfaceFn, p: RollingParams, dirHint: V3): Deriv {
  return (st, tt) => {
    // Direction continuation: a stage past the slip-zero point must
    // not flip friction inside the same RK4 step (event handling in
    // stepRolling localizes the crossing instead).
    const u: V3 = [0, 0, 0]
    slipVelocity(st, surf, p, tt, u)
    const us = Math.hypot(u[0], u[1], u[2])
    let dir: V3 = dirHint
    if (us >= 1e-12) {
      const cand: V3 = [u[0] / us, u[1] / us, u[2] / us]
      if (cand[0] * dirHint[0] + cand[1] * dirHint[1] + cand[2] * dirHint[2] > 0) dir = cand
    }
    const a = kineticRollingAccel(st, surf, p, tt, dir)
    return {
      d: Float64Array.from([st.vx, st.vy, a.ax, a.ay, a.dwx, a.dwy, a.dwz]),
      detach: a.mode === 'detach',
    }
  }
}

/** Advance one step with event handling (roll-break, slip-zero,
 *  detach), mirroring the seated-friction module's discipline:
 *  the earliest event in the step is processed first, and the
 *  remainder continues in the new mode. Roll-break scanning uses
 *  8 samples of the (plane: time-only) stick margin — a sampling
 *  check, not a proof against narrower episodes. */
export function stepRolling(
  st: RollingState, surf: SurfaceFn, p: RollingParams, t: number, dt: number,
): RollingStepResult {
  let advanced = 0
  for (let iter = 0; iter < 16; iter++) {
    const rem = dt - advanced
    if (rem <= 0) break
    const tc = t + advanced
    const a0 = rollingAccel(st, surf, p, tc)
    if (a0.mode === 'detach') return { mode: 'detach', advanced }

    if (a0.mode === 'roll') {
      const marginAt = (tt: number): { margin: number; nPerMass: number } => {
        surf(st.x, st.y, tt, sf); checkPlane(sf)
        const req = rollRequirement(sf, p)
        const fm = Math.hypot(req.fReq[0], req.fReq[1], req.fReq[2])
        return { margin: p.muS * req.nPerMass - fm, nPerMass: req.nPerMass }
      }
      let evT = -1; let evDetach = false
      let prevT = tc
      for (let k = 1; k <= 8; k++) {
        const tt = tc + (rem * k) / 8
        const m = marginAt(tt)
        if (m.nPerMass <= 0 || m.margin < 0) {
          let lo = prevT; let hi = tt
          let detachSide = m.nPerMass <= 0
          for (let it = 0; it < 50; it++) {
            const mid = 0.5 * (lo + hi)
            const mm = marginAt(mid)
            if (mm.nPerMass <= 0 || mm.margin < 0) { hi = mid; detachSide = mm.nPerMass <= 0 } else lo = mid
          }
          evT = hi; evDetach = detachSide
          break
        }
        prevT = tt
      }
      const target = evT < 0 ? rem : evT - tc
      const r = rk4Rolling(st, tc, target, rollDeriv(surf, p))
      Object.assign(st, r.state)
      advanced += target
      if (evT < 0) return { mode: 'roll', advanced: dt }
      if (evDetach) return { mode: 'detach', advanced }
      continue
    }

    // Slip.
    const u0: V3 = [0, 0, 0]
    slipVelocity(st, surf, p, tc, u0)
    const us0 = Math.hypot(u0[0], u0[1], u0[2])
    let dirHint: V3
    if (us0 >= 1e-12) dirHint = [u0[0] / us0, u0[1] / us0, u0[2] / us0]
    else {
      // Impending slip from rest-of-slip: direction of -fReq.
      surf(st.x, st.y, tc, sf)
      const req = rollRequirement(sf, p)
      const fm = Math.hypot(req.fReq[0], req.fReq[1], req.fReq[2]) || 1
      dirHint = [-req.fReq[0] / fm, -req.fReq[1] / fm, -req.fReq[2] / fm]
    }
    const deriv = slipDeriv(surf, p, dirHint)
    const start: RollingState = { ...st }
    const signedSlip = (s2: RollingState, tt: number): number => {
      const u: V3 = [0, 0, 0]
      slipVelocity(s2, surf, p, tt, u)
      return u[0] * dirHint[0] + u[1] * dirHint[1] + u[2] * dirHint[2]
    }
    // M6-26: scan integrated sub-states (8 sub-fractions) for the
    // first detach / slip-direction crossing. Endpoint + stage
    // checks miss a negative-normal interval between stage times,
    // and a rejected projected zero must not discard a later
    // detach in the same step.
    let prevFrac = 0; let badFrac = Infinity
    let detachInside = false; let zeroInside = false
    let lastGood: { state: RollingState; sawDetach: boolean } | null = null
    for (let k = 1; k <= 8; k++) {
      const f = k / 8
      const r = rk4Rolling(start, tc, rem * f, deriv)
      const aa = rollingAccel(r.state, surf, p, tc + rem * f)
      const dBad = r.sawDetach || aa.mode === 'detach'
      const zBad = us0 >= 1e-12 && signedSlip(r.state, tc + rem * f) <= 0
      if (dBad || zBad) {
        badFrac = f; detachInside = dBad; zeroInside = zBad
        break
      }
      prevFrac = f; lastGood = r
    }
    if (!Number.isFinite(badFrac)) {
      Object.assign(st, (lastGood as { state: RollingState }).state)
      return { mode: 'slip', advanced: dt }
    }
    let fracDetach = Infinity
    if (detachInside) {
      let lo = prevFrac; let hi = badFrac
      for (let it = 0; it < 50; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Rolling(start, tc, rem * mid, deriv)
        const aa = rollingAccel(r.state, surf, p, tc + rem * mid)
        if (r.sawDetach || aa.mode === 'detach') hi = mid; else lo = mid
      }
      fracDetach = hi
    }
    let fracZero = Infinity
    if (zeroInside) {
      let lo = prevFrac; let hi = badFrac
      for (let it = 0; it < 60; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Rolling(start, tc, rem * mid, deriv)
        if (signedSlip(r.state, tc + rem * mid) > 0) lo = mid; else hi = mid
      }
      fracZero = hi
    }
    if (fracDetach <= fracZero) {
      const r = rk4Rolling(start, tc, rem * fracDetach, deriv)
      Object.assign(st, r.state)
      return { mode: 'detach', advanced: advanced + rem * fracDetach }
    }
    const r = rk4Rolling(start, tc, rem * fracZero, deriv)
    Object.assign(st, r.state)
    const uEnd: V3 = [0, 0, 0]
    slipVelocity(st, surf, p, tc + rem * fracZero, uEnd)
    if (Math.hypot(uEnd[0], uEnd[1], uEnd[2]) > 1e-9) {
      // Direction rotation, not a rolling transition: subdivide
      // here and re-integrate the remainder with the new slip
      // direction (a pending detach is re-detected there).
      advanced += rem * fracZero
      continue
    }
    projectToRolling(st, surf, p, tc + rem * fracZero)
    advanced += rem * fracZero
    // Loop: roll check at the transition point decides the rest.
  }
  const aEnd = rollingAccel(st, surf, p, t + advanced)
  return { mode: aEnd.mode, advanced }
}

/** Center height z_C = h + R/sn (offset along the normal). */
export function centerHeight(st: RollingState, surf: SurfaceFn, p: RollingParams, t: number): number {
  surf(st.x, st.y, t, sf); checkPlane(sf)
  const sn = Math.sqrt(1 + sf.hx * sf.hx + sf.hy * sf.hy)
  return sf.h + p.radius / sn
}
