// M6-27 RESEARCH ONLY — not wired into any sim.
// Rigid dumbbell grain (two equal spheres rigidly joined):
// torque-free flight and single-contact frictional impact.
// docs/m6_27_design.md. Stage A only: no sustained contact.
//
// PUT-IN: lobes (radius r, mass 1/2 each per unit total mass)
// at C ± d*e; axisymmetric inertia I_par = beta*r^2,
// I_perp = beta*r^2 + d^2; instantaneous Coulomb impact at a
// caller-localized contact point.
// EMERGED: contact-point velocity change Du = K*J with
// K = I3 - [r]x I^-1 [r]x; stick impulse from the tangential
// 2x2 block including normal-tangential coupling; angular
// momentum about the contact point conserved by construction.
// claim-tier: established (rigid-body mechanics).
import { G } from './surface'

export interface DumbbellParams {
  lobeRadius: number
  /** Half the lobe-center separation (m). 0 = sphere control. */
  halfSeparation: number
  /** Lobe-sphere inertia factor I_lobe/(m_lobe r^2), 2/5 solid. */
  beta: number
  muS: number
  muK: number
  restitution: number
}

export interface DumbbellState {
  cx: number; cy: number; cz: number
  vx: number; vy: number; vz: number
  /** Unit axis vector (lobe offsets are ±halfSeparation * e). */
  ex: number; ey: number; ez: number
  wx: number; wy: number; wz: number
}

type V3 = [number, number, number]
const dot = (a: V3, b: V3): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
const cross = (a: V3, b: V3): V3 => [
  a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0],
]
const add = (a: V3, b: V3): V3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
const sub = (a: V3, b: V3): V3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
const scale = (a: V3, k: number): V3 => [a[0] * k, a[1] * k, a[2] * k]
const norm = (a: V3): V3 => {
  const m = Math.hypot(a[0], a[1], a[2]) || 1
  return [a[0] / m, a[1] / m, a[2] / m]
}

export function inertiaPerMass(p: DumbbellParams): { iPar: number; iPerp: number } {
  const base = p.beta * p.lobeRadius * p.lobeRadius
  return { iPar: base, iPerp: base + p.halfSeparation * p.halfSeparation }
}

function axisOf(st: DumbbellState): V3 {
  return norm([st.ex, st.ey, st.ez])
}

/** Apply the axisymmetric inertia tensor (per unit mass). */
function applyI(e: V3, v: V3, p: DumbbellParams): V3 {
  const { iPar, iPerp } = inertiaPerMass(p)
  const par = dot(v, e)
  return [
    iPerp * v[0] + (iPar - iPerp) * par * e[0],
    iPerp * v[1] + (iPar - iPerp) * par * e[1],
    iPerp * v[2] + (iPar - iPerp) * par * e[2],
  ]
}

function applyIInv(e: V3, v: V3, p: DumbbellParams): V3 {
  const { iPar, iPerp } = inertiaPerMass(p)
  const par = dot(v, e)
  return [
    v[0] / iPerp + (1 / iPar - 1 / iPerp) * par * e[0],
    v[1] / iPerp + (1 / iPar - 1 / iPerp) * par * e[1],
    v[2] / iPerp + (1 / iPar - 1 / iPerp) * par * e[2],
  ]
}

/** Angular momentum about the center of mass, per unit mass. */
export function angularMomentum(st: DumbbellState, p: DumbbellParams, out: V3): void {
  const L = applyI(axisOf(st), [st.wx, st.wy, st.wz], p)
  out[0] = L[0]; out[1] = L[1]; out[2] = L[2]
}

/** Kinetic energy per unit mass (translation + rotation). */
export function kineticEnergyPerMass(st: DumbbellState, p: DumbbellParams): number {
  const L: V3 = [0, 0, 0]
  angularMomentum(st, p, L)
  return 0.5 * (st.vx ** 2 + st.vy ** 2 + st.vz ** 2)
    + 0.5 * (st.wx * L[0] + st.wy * L[1] + st.wz * L[2])
}

/** Lobe centers [C + d*e, C - d*e]. */
export function lobeCenters(st: DumbbellState, p: DumbbellParams): [V3, V3] {
  const e = axisOf(st)
  const d = p.halfSeparation
  return [
    [st.cx + d * e[0], st.cy + d * e[1], st.cz + d * e[2]],
    [st.cx - d * e[0], st.cy - d * e[1], st.cz - d * e[2]],
  ]
}

/** Advance torque-free flight by dt (gravity on the COM;
 *  rotation with L fixed over the step, RK4 on e with stage
 *  renormalization — a unit-constraint projection, documented
 *  in docs/m6_27_design.md). Mutates st. */
export function dumbbellFlightStep(st: DumbbellState, p: DumbbellParams, dt: number): void {
  st.cx += st.vx * dt
  st.cy += st.vy * dt
  st.cz += st.vz * dt - 0.5 * G * dt * dt
  st.vz -= G * dt
  const L: V3 = [0, 0, 0]
  angularMomentum(st, p, L)
  const deriv = (eIn: V3): V3 => {
    const e = norm(eIn)
    const w = applyIInv(e, L, p)
    return cross(w, e)
  }
  const e0 = axisOf(st)
  const k1 = deriv(e0)
  const k2 = deriv(add(e0, scale(k1, dt / 2)))
  const k3 = deriv(add(e0, scale(k2, dt / 2)))
  const k4 = deriv(add(e0, scale(k3, dt)))
  const e1 = norm([
    e0[0] + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
    e0[1] + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
    e0[2] + (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]),
  ])
  st.ex = e1[0]; st.ey = e1[1]; st.ez = e1[2]
  const w = applyIInv(e1, L, p)
  st.wx = w[0]; st.wy = w[1]; st.wz = w[2]
}

export type DumbbellImpactMode = 'stick' | 'slide' | 'separating'

export interface DumbbellImpactResult {
  mode: DumbbellImpactMode
  normalImpulse: number
  impulse: V3
  /** Pre-impact contact-point relative velocity. */
  uPre: V3
}

/** Contact-point relative velocity u = V + w x r - Vs,
 *  r = point - C. */
export function contactVelocity(
  st: DumbbellState, point: V3, surfaceVel: V3,
): V3 {
  const r = sub(point, [st.cx, st.cy, st.cz])
  const wxr = cross([st.wx, st.wy, st.wz], r)
  return [
    st.vx + wxr[0] - surfaceVel[0],
    st.vy + wxr[1] - surfaceVel[1],
    st.vz + wxr[2] - surfaceVel[2],
  ]
}

/** Apply a frictional rigid-body impact at `point` with unit
 *  normal `normal` and plate material velocity `surfaceVel`.
 *  Mutates velocities/spin; position and axis unchanged. */
export function dumbbellImpact(
  st: DumbbellState, point: V3, normal: V3, surfaceVel: V3, p: DumbbellParams,
): DumbbellImpactResult {
  const e = axisOf(st)
  const n = norm(normal)
  const r = sub(point, [st.cx, st.cy, st.cz])
  const u = contactVelocity(st, point, surfaceVel)
  const un = dot(u, n)
  if (un >= 0) return { mode: 'separating', normalImpulse: 0, impulse: [0, 0, 0], uPre: u }
  // Effective-mass operator K (per unit mass): K(J) = J - [r]x I^-1 [r]x J.
  const applyK = (J: V3): V3 => {
    const t = applyIInv(e, cross(r, J), p)
    const c = cross(r, t)
    return [J[0] - c[0], J[1] - c[1], J[2] - c[2]]
  }
  const Kn = applyK(n)
  const Jn = -(1 + p.restitution) * un / dot(n, Kn)
  // Tangential basis.
  const helper: V3 = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]
  const t1 = norm(cross(helper, n))
  const t2 = cross(n, t1)
  const Kt1 = applyK(t1); const Kt2 = applyK(t2)
  const ut: V3 = sub(u, scale(n, un))
  const ut1 = dot(ut, t1); const ut2 = dot(ut, t2)
  // Stick solve: tangential Du = -ut given the normal impulse.
  const a11 = dot(t1, Kt1); const a12 = dot(t1, Kt2)
  const a21 = dot(t2, Kt1); const a22 = dot(t2, Kt2)
  const rhs1 = -ut1 - Jn * dot(t1, Kn)
  const rhs2 = -ut2 - Jn * dot(t2, Kn)
  const det = a11 * a22 - a12 * a21
  const stickA = (a22 * rhs1 - a12 * rhs2) / det
  const stickB = (-a21 * rhs1 + a11 * rhs2) / det
  const stickMag = Math.hypot(stickA, stickB)
  let Jt: V3
  let mode: DumbbellImpactMode
  if (stickMag <= p.muS * Jn) {
    Jt = add(scale(t1, stickA), scale(t2, stickB))
    mode = 'stick'
  } else {
    // Impending tangential motion after the normal impulse
    // (equals ut whenever the coupling Kt_n vanishes, e.g. a sphere).
    const impending = add(ut, sub(scale(Kn, Jn), scale(n, Jn * dot(n, Kn))))
    const im = Math.hypot(impending[0], impending[1], impending[2])
    if (im === 0 || p.muK * Jn >= stickMag) {
      // Reversal guard: kinetic impulse would overshoot stick.
      Jt = add(scale(t1, stickA), scale(t2, stickB))
      mode = 'stick'
    } else {
      const k = -p.muK * Jn / im
      Jt = scale(impending, k)
      mode = 'slide'
      // Guard against tangential reversal inside this impact.
      const du = applyK(add(scale(n, Jn), Jt))
      const uPost = add(u, du)
      const utPost = sub(uPost, scale(n, dot(uPost, n)))
      if (dot(utPost, impending) < 0) {
        Jt = add(scale(t1, stickA), scale(t2, stickB))
        mode = 'stick'
      }
    }
  }
  const J = add(scale(n, Jn), Jt)
  st.vx += J[0]; st.vy += J[1]; st.vz += J[2]
  const dw = applyIInv(e, cross(r, J), p)
  st.wx += dw[0]; st.wy += dw[1]; st.wz += dw[2]
  return { mode, normalImpulse: Jn, impulse: J, uPre: u }
}

/** Angular momentum about a fixed point P, per unit mass:
 *  H_P = (C - P) x V + I*w. Conserved across dumbbellImpact
 *  when P is the contact point (impulses act at P). */
export function angularMomentumAbout(
  st: DumbbellState, p: DumbbellParams, point: V3, out: V3,
): void {
  const L: V3 = [0, 0, 0]
  angularMomentum(st, p, L)
  const c = cross(sub([st.cx, st.cy, st.cz], point), [st.vx, st.vy, st.vz])
  out[0] = c[0] + L[0]; out[1] = c[1] + L[1]; out[2] = c[2] + L[2]
}
