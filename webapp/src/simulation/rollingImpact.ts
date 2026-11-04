// M6-18 RESEARCH ONLY — not wired into TransientSandSim.
// Frictional rigid-sphere impact and flight<->rolling state
// conversion on the M6-17 dynamic surface (docs/m6_18_design.md).
//
// PUT-IN: instantaneous impact at a contact the caller has
// localized; normal restitution e; Coulomb friction impulse
// with fixed slip direction during the impulse, guarded so the
// slip never reverses inside one impact.
// EMERGED: delta_u_t = Jt (1 + 1/beta); stick impulse
// Jt = -beta u_t/(1+beta) when inside the friction cone;
// angular momentum about the contact point is conserved by
// construction (impulses act at the contact).
// claim-tier: established (rigid impact mechanics).
import {
  dynCenterPartials, dynCenterVelocity,
  type SurfaceFnDyn, type SurfaceFieldDyn,
} from './rollingDynamic'
import type { RollingParams, RollingState } from './rollingSphere'

export interface ImpactParams extends RollingParams {
  restitution: number
}

/** Free-flight state of the sphere center. */
export interface FlightState {
  cx: number; cy: number; cz: number
  vx: number; vy: number; vz: number
  wx: number; wy: number; wz: number
}

export type ImpactMode = 'stick' | 'slide' | 'separating'

export interface ImpactResult {
  mode: ImpactMode
  normalImpulse: number
  frictionImpulse: [number, number, number]
  /** Contact parameters (x,y) used. */
  contact: [number, number]
}

type V3 = [number, number, number]
const dot = (a: V3, b: V3): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
const cross = (a: V3, b: V3): V3 => [
  a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0],
]

const sf: SurfaceFieldDyn = {
  h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0,
  hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0,
  ht: 0, htt: 0, htx: 0, hty: 0, htxx: 0, htxy: 0, htyy: 0, httx: 0, htty: 0,
}

/** Solve contact parameters (x,y) whose offset center equals C.
 *  Newton iteration with the horizontal Jacobian of C(q). */
export function solveContactParams(
  cx: number, cy: number, surf: SurfaceFnDyn, radius: number, t: number,
  guess?: [number, number],
): [number, number] {
  let x = guess ? guess[0] : cx
  let y = guess ? guess[1] : cy
  for (let it = 0; it < 50; it++) {
    const cp = dynCenterPartials(x, y, t, surf, radius)
    const ex = cp.C[0] - cx; const ey = cp.C[1] - cy
    if (Math.hypot(ex, ey) < 1e-13) break
    const a = cp.D[0][0]; const b = cp.D[1][0]
    const c = cp.D[0][1]; const d = cp.D[1][1]
    const det = a * d - b * c
    x -= (d * ex - b * ey) / det
    y -= (-c * ex + a * ey) / det
  }
  return [x, y]
}

/** Apply the frictional impact to a flight state at contact
 *  (mutates velocities and spin in place; position unchanged). */
export function sphereImpact(
  st: FlightState, surf: SurfaceFnDyn, p: ImpactParams, t: number,
): ImpactResult {
  const [qx, qy] = solveContactParams(st.cx, st.cy, surf, p.radius, t)
  surf(qx, qy, t, sf)
  const sn = Math.hypot(sf.hx, sf.hy, 1)
  const n: V3 = [-sf.hx / sn, -sf.hy / sn, 1 / sn]
  const w: V3 = [st.wx, st.wy, st.wz]
  const wxn = cross(w, n)
  const u: V3 = [
    st.vx - p.radius * wxn[0],
    st.vy - p.radius * wxn[1],
    st.vz - p.radius * wxn[2] - sf.ht,
  ]
  const un = dot(u, n)
  if (un >= 0) {
    return { mode: 'separating', normalImpulse: 0, frictionImpulse: [0, 0, 0], contact: [qx, qy] }
  }
  const Jn = -(1 + p.restitution) * un
  const ut: V3 = [u[0] - un * n[0], u[1] - un * n[1], u[2] - un * n[2]]
  const uts = Math.hypot(ut[0], ut[1], ut[2])
  const stickMag = (p.beta / (1 + p.beta)) * uts
  let Jt: V3
  let mode: ImpactMode
  if (uts === 0 || stickMag <= p.muS * Jn) {
    Jt = [
      -ut[0] * p.beta / (1 + p.beta),
      -ut[1] * p.beta / (1 + p.beta),
      -ut[2] * p.beta / (1 + p.beta),
    ]
    mode = 'stick'
  } else if (p.muK * Jn >= stickMag) {
    // Guard: kinetic impulse would reverse the slip inside the
    // impact (only possible when muK > muS in effect); cap at stick.
    Jt = [
      -ut[0] * p.beta / (1 + p.beta),
      -ut[1] * p.beta / (1 + p.beta),
      -ut[2] * p.beta / (1 + p.beta),
    ]
    mode = 'stick'
  } else {
    const k = -p.muK * Jn / uts
    Jt = [k * ut[0], k * ut[1], k * ut[2]]
    mode = 'slide'
  }
  st.vx += Jn * n[0] + Jt[0]
  st.vy += Jn * n[1] + Jt[1]
  st.vz += Jn * n[2] + Jt[2]
  const nxJ = cross(n, Jt)
  const kw = -1 / (p.beta * p.radius)
  st.wx += kw * nxJ[0]; st.wy += kw * nxJ[1]; st.wz += kw * nxJ[2]
  return { mode, normalImpulse: Jn, frictionImpulse: Jt, contact: [qx, qy] }
}

/** Convert a flight state at contact into a rolling state. The
 *  parameter velocity solves J qdot = V - C_t in the tangential
 *  (least-squares / J^T) sense; any residual normal velocity is
 *  the caller's responsibility (impact first). */
export function flightToRolling(
  st: FlightState, surf: SurfaceFnDyn, p: RollingParams, t: number,
): RollingState {
  const [qx, qy] = solveContactParams(st.cx, st.cy, surf, p.radius, t)
  const cp = dynCenterPartials(qx, qy, t, surf, p.radius)
  const rhsV: V3 = [
    st.vx - cp.D[2][0], st.vy - cp.D[2][1], st.vz - cp.D[2][2],
  ]
  const m00 = dot(cp.D[0], cp.D[0]); const m01 = dot(cp.D[0], cp.D[1]); const m11 = dot(cp.D[1], cp.D[1])
  const r0 = dot(cp.D[0], rhsV); const r1 = dot(cp.D[1], rhsV)
  const det = m00 * m11 - m01 * m01
  return {
    x: qx, y: qy,
    vx: (m11 * r0 - m01 * r1) / det,
    vy: (-m01 * r0 + m00 * r1) / det,
    wx: st.wx, wy: st.wy, wz: st.wz,
  }
}

/** Convert a rolling state into a flight state (detach). */
export function rollingToFlight(
  st: RollingState, surf: SurfaceFnDyn, p: RollingParams, t: number,
): FlightState {
  const cp = dynCenterPartials(st.x, st.y, t, surf, p.radius)
  const v: [number, number, number] = [0, 0, 0]
  dynCenterVelocity(st, surf, p, t, v)
  return {
    cx: cp.C[0], cy: cp.C[1], cz: cp.C[2],
    vx: v[0], vy: v[1], vz: v[2],
    wx: st.wx, wy: st.wy, wz: st.wz,
  }
}

/** Angular momentum about the contact point, per unit mass:
 *  L = R n x V + beta R^2 w. Conserved across sphereImpact. */
export function contactAngularMomentum(
  st: FlightState, surf: SurfaceFnDyn, p: RollingParams, t: number,
  out: [number, number, number],
): void {
  const [qx, qy] = solveContactParams(st.cx, st.cy, surf, p.radius, t)
  surf(qx, qy, t, sf)
  const sn = Math.hypot(sf.hx, sf.hy, 1)
  const n: V3 = [-sf.hx / sn, -sf.hy / sn, 1 / sn]
  const nxV = cross(n, [st.vx, st.vy, st.vz])
  out[0] = p.radius * nxV[0] + p.beta * p.radius * p.radius * st.wx
  out[1] = p.radius * nxV[1] + p.beta * p.radius * p.radius * st.wy
  out[2] = p.radius * nxV[2] + p.beta * p.radius * p.radius * st.wz
}

/** Flight kinetic energy per unit mass. */
export function flightEnergyPerMass(st: FlightState, p: RollingParams): number {
  const v2 = st.vx * st.vx + st.vy * st.vy + st.vz * st.vz
  const w2 = st.wx * st.wx + st.wy * st.wy + st.wz * st.wz
  return 0.5 * v2 + 0.5 * p.beta * p.radius * p.radius * w2
}
