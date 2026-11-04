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

/** Constraint acceleration term B = v^T H v + 2 (grad ht).v + htt.
 *  The normal force and the tangential equation share this B, so
 *  acceleration and normal force derive from one constraint. */
function constraintB(st: SeatedState, s: SurfaceField): number {
  return st.vx * (s.hxx * st.vx + s.hxy * st.vy)
    + st.vy * (s.hxy * st.vx + s.hyy * st.vy)
    + 2 * (s.htx * st.vx + s.hty * st.vy) + s.htt
}

function accelWithFrictionQ(
  st: SeatedState, s: SurfaceField, qfx: number, qfy: number,
): [number, number] {
  // Lagrange equations for T = (|v|^2 + u^2)/2 with u = grad h . v + ht
  // simplify (the u(Hv+grad ht) terms cancel against d/dt of the
  // momentum) to exactly:
  //   M a = Qf - grad h (g + B),   M = I + grad h grad h^T,
  // with B as in constraintB — including BOTH mixed space-time
  // contributions to u-dot (M6-13R review; the previous expansion
  // double-counted M-dot terms and dropped one grad ht . v).
  const { hx, hy } = s
  const gb = G + constraintB(st, s)
  const rx = qfx - hx * gb
  const ry = qfy - hy * gb
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
  const B = constraintB(st, s)
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

export interface SeatedStepResult {
  mode: SeatedMode
  /** Seconds of the requested step actually advanced (<= dt). On
   *  'detach' the state is at the detach crossing and `advanced`
   *  is the time to reach it. */
  advanced: number
}

/** Kinetic-slip acceleration with an explicit motion direction
 *  (unit tangential vector). Used inside RK4 stages so a stage state
 *  whose velocity is near or past zero does not flip the friction
 *  direction mid-step: stages with negligible speed keep the
 *  step-initial direction (the mathematical continuation whose
 *  signed speed crosses zero at the stop event). */
function kineticAccel(
  st: SeatedState, surf: SurfaceFn, p: FrictionParams, t: number,
  dirHint: [number, number, number],
): SeatedAccel {
  surf(st.x, st.y, t, sf)
  const s = sf
  const B = constraintB(st, s)
  const sn = Math.sqrt(1 + s.hx * s.hx + s.hy * s.hy)
  const nPerMass = (B + G) / sn
  if (nPerMass <= 0) return { ax: 0, ay: 0, nPerMass, mode: 'detach' }
  const vt = tangentialSpeed(st, s)
  let dx: number; let dy: number; let dz: number
  if (vt >= 1e-12) {
    const vzRel = s.hx * st.vx + s.hy * st.vy
    dx = st.vx / vt; dy = st.vy / vt; dz = vzRel / vt
    // A stage past the stop point has reversed velocity; using its
    // direction would flip friction inside the same RK4 step and
    // corrupt event localization. Keep the step-initial direction
    // (continuation) whenever the stage direction has reversed.
    if (dx * dirHint[0] + dy * dirHint[1] + dz * dirHint[2] <= 0) {
      [dx, dy, dz] = dirHint
    }
  } else {
    [dx, dy, dz] = dirHint
  }
  const fm = p.muK * nPerMass
  const [qfx, qfy] = frictionQ(s, -fm * dx, -fm * dy, -fm * dz)
  const [ax, ay] = accelWithFrictionQ(st, s, qfx, qfy)
  return { ax, ay, nPerMass, mode: 'slip' }
}

function rk4Advance(
  s0: SeatedState, surf: SurfaceFn, p: FrictionParams, t: number, h: number,
): { state: SeatedState; sawDetach: boolean } {
  surf(s0.x, s0.y, t, sf)
  const vz0 = sf.hx * s0.vx + sf.hy * s0.vy
  const sp0 = Math.hypot(s0.vx, s0.vy, vz0) || 1
  const dirHint: [number, number, number] = [s0.vx / sp0, s0.vy / sp0, vz0 / sp0]
  let sawDetach = false
  const deriv = (s2: SeatedState, tt: number): [number, number, number, number] => {
    const a = kineticAccel(s2, surf, p, tt, dirHint)
    if (a.mode === 'detach') sawDetach = true
    return [s2.vx, s2.vy, a.ax, a.ay]
  }
  const k1 = deriv(s0, t)
  const s1 = { x: s0.x + k1[0] * h / 2, y: s0.y + k1[1] * h / 2, vx: s0.vx + k1[2] * h / 2, vy: s0.vy + k1[3] * h / 2 }
  const k2 = deriv(s1, t + h / 2)
  const s2 = { x: s0.x + k2[0] * h / 2, y: s0.y + k2[1] * h / 2, vx: s0.vx + k2[2] * h / 2, vy: s0.vy + k2[3] * h / 2 }
  const k3 = deriv(s2, t + h / 2)
  const s3 = { x: s0.x + k3[0] * h, y: s0.y + k3[1] * h, vx: s0.vx + k3[2] * h, vy: s0.vy + k3[3] * h }
  const k4 = deriv(s3, t + h)
  return {
    state: {
      x: s0.x + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
      y: s0.y + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
      vx: s0.vx + (h / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]),
      vy: s0.vy + (h / 6) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]),
    },
    sawDetach,
  }
}

/** Static-friction margin at rest at fixed (x, y): muS*N - |f_req|
 *  per unit mass (>0 can stick). Also returns nPerMass (<=0 detach). */
function stickMargin(
  st: SeatedState, surf: SurfaceFn, p: FrictionParams, t: number,
): { margin: number; nPerMass: number } {
  const rest: SeatedState = { x: st.x, y: st.y, vx: 0, vy: 0 }
  surf(rest.x, rest.y, t, sf)
  const s = sf
  const sn = Math.sqrt(1 + s.hx * s.hx + s.hy * s.hy)
  const nPerMass = (s.htt + G) / sn
  if (nPerMass <= 0) return { margin: -Infinity, nPerMass }
  const qfx = (G + s.htt) * s.hx; const qfy = (G + s.htt) * s.hy
  const det = 1 + s.hx * s.hx + s.hy * s.hy
  const fx = ((1 + s.hy * s.hy) * qfx - s.hx * s.hy * qfy) / det
  const fy = (-s.hx * s.hy * qfx + (1 + s.hx * s.hx) * qfy) / det
  const fz = s.hx * fx + s.hy * fy
  return { margin: p.muS * nPerMass - Math.hypot(fx, fy, fz), nPerMass }
}

/** Advance seated motion by up to dt with event handling (M6-13R):
 *  - detach (normal force crossing zero) inside the step is located
 *    by bisection; the state is left at the crossing.
 *  - kinetic-friction stopping (tangential velocity reaching zero)
 *    inside the step is located by bisection; the grain is set
 *    exactly to rest there, then stick is re-evaluated and any
 *    remaining time continues in the resulting mode.
 *  - stick breaking mid-step (time-varying surface) is bracketed by
 *    sampling the stick margin at 8 sub-times, then bisected.
 *  Without this, RK4 stages straddling a stop reverse the friction
 *  direction and produce false oscillation/drift. Event brackets
 *  found only by the 8-point stick-margin sample are documented as
 *  sampling-based (not a proof of no narrower stick-break episode). */
export function stepSeated(
  st: SeatedState, surf: SurfaceFn, p: FrictionParams, t: number, dt: number,
): SeatedStepResult {
  let advanced = 0
  for (let iter = 0; iter < 16; iter++) {
    const rem = dt - advanced
    if (rem <= 0) break
    const tc = t + advanced
    const a0 = seatedAccel(st, surf, p, tc)
    if (a0.mode === 'detach') return { mode: 'detach', advanced }

    if (a0.mode === 'stick') {
      // Scan for a stick-break or detach inside the remaining step
      // at the fixed position (state does not move while stuck).
      let evT = -1; let evDetach = false
      let prevM = stickMargin(st, surf, p, tc)
      let prevT = tc
      for (let k = 1; k <= 8; k++) {
        const tt = tc + (rem * k) / 8
        const m = stickMargin(st, surf, p, tt)
        if (m.nPerMass <= 0 || m.margin < 0) {
          // Bisect the first sign change of (detach or margin).
          let lo = prevT; let hi = tt
          const badAt = (mm: { margin: number; nPerMass: number }) => mm.nPerMass <= 0 || mm.margin < 0
          let detachSide = m.nPerMass <= 0
          for (let it = 0; it < 50; it++) {
            const mid = 0.5 * (lo + hi)
            const mm = stickMargin(st, surf, p, mid)
            if (badAt(mm)) { hi = mid; detachSide = mm.nPerMass <= 0 } else lo = mid
          }
          evT = hi; evDetach = detachSide
          break
        }
        prevM = m; prevT = tt
      }
      void prevM
      if (evT < 0) return { mode: 'stick', advanced: dt }
      advanced = evT - t
      if (evDetach) return { mode: 'detach', advanced }
      continue // re-evaluate at the stick-break instant (now slip)
    }

    // Slip: trial RK4 over the remainder, watching for events.
    const start: SeatedState = { ...st }
    surf(start.x, start.y, tc, sf)
    const vt0 = tangentialSpeed(start, sf)
    const ezx = start.vx; const ezy = start.vy
    const trial = rk4Advance(start, surf, p, tc, rem)
    surf(trial.state.x, trial.state.y, tc + rem, sf)
    const endAccel = seatedAccel(trial.state, surf, p, tc + rem)
    const detachInside = trial.sawDetach || endAccel.mode === 'detach'
    // Signed tangential speed along the initial direction.
    const signedSpeed = (s2: SeatedState): number => {
      if (vt0 <= 0) return 0
      return (s2.vx * ezx + s2.vy * ezy) / Math.hypot(ezx, ezy)
    }
    const stopInside = vt0 > 0 && signedSpeed(trial.state) <= 0

    if (!detachInside && !stopInside) {
      Object.assign(st, trial.state)
      return { mode: 'slip', advanced: dt }
    }

    if (detachInside) {
      // Bisect the fraction whose integrated state first detaches.
      let lo = 0; let hi = 1
      for (let it = 0; it < 50; it++) {
        const mid = 0.5 * (lo + hi)
        const r = rk4Advance(start, surf, p, tc, rem * mid)
        const aa = seatedAccel(r.state, surf, p, tc + rem * mid)
        if (r.sawDetach || aa.mode === 'detach') hi = mid; else lo = mid
      }
      const r = rk4Advance(start, surf, p, tc, rem * hi)
      Object.assign(st, r.state)
      return { mode: 'detach', advanced: advanced + rem * hi }
    }

    // Stopping event: bisect the fraction where signed speed hits 0.
    let lo = 0; let hi = 1
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi)
      const r = rk4Advance(start, surf, p, tc, rem * mid)
      if (signedSpeed(r.state) > 0) lo = mid; else hi = mid
    }
    const r = rk4Advance(start, surf, p, tc, rem * hi)
    surf(r.state.x, r.state.y, tc + rem * hi, sf)
    if (tangentialSpeed(r.state, sf) > 1e-9) {
      // Direction rotated through perpendicular without stopping;
      // not a stop event. Accept the trial step.
      Object.assign(st, trial.state)
      return { mode: 'slip', advanced: dt }
    }
    Object.assign(st, r.state)
    st.vx = 0; st.vy = 0
    advanced += rem * hi
    // Loop: stick check at the stop point decides the remainder.
  }
  const aEnd = seatedAccel(st, surf, p, t + advanced)
  return { mode: aEnd.mode, advanced }
}
