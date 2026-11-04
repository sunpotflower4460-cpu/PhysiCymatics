// M6-13 RESEARCH ONLY — not wired into TransientSandSim.
// Seated-grain horizontal dynamics on a graph surface z = h(x, y, t)
// with Coulomb friction. Validated standalone in
// __tests__/seatedFriction.test.ts against analytic incline/vibration
// results (docs/m6_13_design.md). Connecting this to the sand sim
// invalidates seated-sleep proofs (fixed-position assumption) and
// must stay behind a default-off research switch until those are
// rebuilt.
//
// PUT-IN: point mass constrained to z=h, plate material points move
// vertically only (0,0,ht) — as in the modal plate model.
// EMERGED: normal force per unit mass (B+g)/s with
// B = v^T H v + 2 (grad ht).v + htt; detach iff B+g <= 0, which
// reduces to the existing wtt < -g at rest on a flat surface.
// claim-tier: established (classical constrained mechanics).
import { G } from './surface'

export interface SurfaceField {
  h: number
  hx: number; hy: number
  hxx: number; hxy: number; hyy: number
  ht: number; htt: number
  htx: number; hty: number
}

export type SurfaceFn = (x: number, y: number, t: number, out: SurfaceField) => void

export interface FrictionParams { muS: number; muK: number }

export interface SeatedState { x: number; y: number; vx: number; vy: number }

export type SeatedMode = 'stick' | 'slip' | 'detach'

export interface SeatedAccel {
  ax: number; ay: number
  /** Normal force per unit mass (m/s^2). <= 0 means detach. */
  nPerMass: number
  mode: SeatedMode
}

const F = (): SurfaceField => ({
  h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0,
})
const sf = F()

/** Tangential relative speed of the grain vs the surface point. */
export function tangentialSpeed(st: SeatedState, s: SurfaceField): number {
  const vzRel = s.hx * st.vx + s.hy * st.vy
  return Math.hypot(st.vx, st.vy, vzRel)
}

/** Vertical velocity of a seated grain (for flight transition). */
export function seatedVz(st: SeatedState, s: SurfaceField): number {
  return s.hx * st.vx + s.hy * st.vy + s.ht
}

function accelWithFrictionQ(
  st: SeatedState, s: SurfaceField, qfx: number, qfy: number,
): [number, number] {
  const { hx, hy, hxx, hxy, hyy, ht, htt, htx, hty } = s
  const vx = st.vx; const vy = st.vy
  // d(grad h)/dt along the grain trajectory
  const dhx = hxx * vx + hxy * vy + htx
  const dhy = hxy * vx + hyy * vy + hty
  // Mdot v, M = I + grad h grad h^T
  const mdvx = dhx * (hx * vx + hy * vy) + hx * (dhx * vx + dhy * vy)
  const mdvy = dhy * (hx * vx + hy * vy) + hy * (dhx * vx + dhy * vy)
  const hvx = hxx * vx + hxy * vy; const hvy = hxy * vx + hyy * vy
  const ghv = hx * vx + hy * vy
  const rx = -G * hx + qfx - mdvx - htt * hx - ht * (hvx + htx) + (ghv + ht) * (hvx + htx)
  const ry = -G * hy + qfy - mdvy - htt * hy - ht * (hvy + hty) + (ghv + ht) * (hvy + hty)
  // Solve M a = r (2x2 SPD, det = 1 + |grad h|^2)
  const det = 1 + hx * hx + hy * hy
  const ax = ((1 + hy * hy) * rx - hx * hy * ry) / det
  const ay = (-hx * hy * rx + (1 + hx * hx) * ry) / det
  return [ax, ay]
}

/** Generalized friction force (per unit mass) for a 3D tangential
 *  friction acceleration vector f (f . n = 0). */
function frictionQ(s: SurfaceField, fx: number, fy: number, fz: number): [number, number] {
  return [fx + s.hx * fz, fy + s.hy * fz]
}

/** Seated acceleration under Coulomb friction. Pure function of the
 *  state and the surface field at (x, y, t). */
export function seatedAccel(
  st: SeatedState, surf: SurfaceFn, p: FrictionParams, t: number,
): SeatedAccel {
  surf(st.x, st.y, t, sf)
  const s = sf
  const B = st.vx * (s.hxx * st.vx + s.hxy * st.vy)
    + st.vy * (s.hxy * st.vx + s.hyy * st.vy)
    + 2 * (s.htx * st.vx + s.hty * st.vy) + s.htt
  const sn = Math.sqrt(1 + s.hx * s.hx + s.hy * s.hy)
  const nPerMass = (B + G) / sn
  if (nPerMass <= 0) return { ax: 0, ay: 0, nPerMass, mode: 'detach' }

  const vt = tangentialSpeed(st, s)
  if (vt < 1e-12) {
    // Stick candidate: required generalized friction for a = 0.
    // At v = 0 the equation reduces to Qf = (g + htt) grad h.
    const qfx = (G + s.htt) * s.hx; const qfy = (G + s.htt) * s.hy
    // Invert Q = M (fx, fy) (with fz = hx fx + hy fy) for |f|.
    const det = 1 + s.hx * s.hx + s.hy * s.hy
    const fx = ((1 + s.hy * s.hy) * qfx - s.hx * s.hy * qfy) / det
    const fy = (-s.hx * s.hy * qfx + (1 + s.hx * s.hx) * qfy) / det
    const fz = s.hx * fx + s.hy * fy
    const fReq = Math.hypot(fx, fy, fz)
    if (fReq <= p.muS * nPerMass) return { ax: 0, ay: 0, nPerMass, mode: 'stick' }
    // Impending slip: the grain would move opposite the required
    // static friction, so kinetic friction acts along +f_req.
    const fm = (p.muK * nPerMass) / fReq
    const [qfx2, qfy2] = frictionQ(s, fm * fx, fm * fy, fm * fz)
    const [ax, ay] = accelWithFrictionQ(st, s, qfx2, qfy2)
    return { ax, ay, nPerMass, mode: 'slip' }
  }

  // Kinetic friction opposing the tangential relative velocity.
  // v_rel = (vx, vy, grad h . v) is already tangential.
  const vzRel = s.hx * st.vx + s.hy * st.vy
  const fm = p.muK * nPerMass / vt
  const fx = -fm * st.vx; const fy = -fm * st.vy; const fz = -fm * vzRel
  const [qfx, qfy] = frictionQ(s, fx, fy, fz)
  const [ax, ay] = accelWithFrictionQ(st, s, qfx, qfy)
  return { ax, ay, nPerMass, mode: 'slip' }
}

/** One RK4 step of seated motion. Returns 'detach' without advancing
 *  if the normal force vanishes at the step start; the caller then
 *  transitions to flight with vz = seatedVz. Sticking grains return
 *  'stick' and do not move. */
export function stepSeated(
  st: SeatedState, surf: SurfaceFn, p: FrictionParams, t: number, dt: number,
): SeatedMode {
  const a0 = seatedAccel(st, surf, p, t)
  if (a0.mode !== 'slip') return a0.mode
  const deriv = (s2: SeatedState, tt: number): [number, number, number, number] => {
    const a = seatedAccel(s2, surf, p, tt)
    return [s2.vx, s2.vy, a.ax, a.ay]
  }
  const s0 = { ...st }
  const k1 = deriv(s0, t)
  const s1 = { x: s0.x + k1[0] * dt / 2, y: s0.y + k1[1] * dt / 2, vx: s0.vx + k1[2] * dt / 2, vy: s0.vy + k1[3] * dt / 2 }
  const k2 = deriv(s1, t + dt / 2)
  const s2 = { x: s0.x + k2[0] * dt / 2, y: s0.y + k2[1] * dt / 2, vx: s0.vx + k2[2] * dt / 2, vy: s0.vy + k2[3] * dt / 2 }
  const k3 = deriv(s2, t + dt / 2)
  const s3 = { x: s0.x + k3[0] * dt, y: s0.y + k3[1] * dt, vx: s0.vx + k3[2] * dt, vy: s0.vy + k3[3] * dt }
  const k4 = deriv(s3, t + dt)
  st.x = s0.x + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
  st.y = s0.y + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
  st.vx = s0.vx + (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
  st.vy = s0.vy + (dt / 6) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
  return 'slip'
}
