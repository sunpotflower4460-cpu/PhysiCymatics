// M6-22 RESEARCH ONLY — not wired into TransientSandSim.
// Frictional impact between two equal rigid spheres
// (docs/m6_22_design.md). Instantaneous impulse per unit mass:
// Jn = -(1+e) u_n / 2; delta_u_t = 2 (1 + 1/beta) Jt; stick when
// |Jt_stick| <= muS |Jn|, else slide with muK (reversal guard).
// Total momentum and total angular momentum are conserved by
// construction; energy never gains for e <= 1.
// claim-tier: established (rigid impact mechanics).
import type { FlightState, ImpactParams } from './rollingImpact'

export type PairImpactMode = 'stick' | 'slide' | 'separating'

type V3 = [number, number, number]
const cross = (a: V3, b: V3): V3 => [
  a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0],
]

/** Apply the pair impact to two flight states (mutates both).
 *  Caller guarantees the spheres are in contact (distance ~2R). */
export function pairImpact(a: FlightState, b: FlightState, p: ImpactParams): PairImpactMode {
  const dx = b.cx - a.cx; const dy = b.cy - a.cy; const dz = b.cz - a.cz
  const dist = Math.hypot(dx, dy, dz)
  if (dist === 0) return 'separating'
  const n: V3 = [dx / dist, dy / dist, dz / dist]
  // Contact-point relative velocity u = (Va - Vb) + R (wa + wb) x n.
  const ws: V3 = [a.wx + b.wx, a.wy + b.wy, a.wz + b.wz]
  const wxn = cross(ws, n)
  const u: V3 = [
    (a.vx - b.vx) + p.radius * wxn[0],
    (a.vy - b.vy) + p.radius * wxn[1],
    (a.vz - b.vz) + p.radius * wxn[2],
  ]
  const un = u[0] * n[0] + u[1] * n[1] + u[2] * n[2]
  // Closing when un > 0 (centers approach along n from a to b).
  if (un <= 0) return 'separating'
  const Jn = -(1 + p.restitution) * un / 2
  const ut: V3 = [u[0] - un * n[0], u[1] - un * n[1], u[2] - un * n[2]]
  const uts = Math.hypot(ut[0], ut[1], ut[2])
  const stickScale = -p.beta / (2 * (1 + p.beta))
  const stickMag = Math.abs(stickScale) * uts
  let Jt: V3
  let mode: PairImpactMode
  if (uts === 0 || stickMag <= p.muS * Math.abs(Jn) || p.muK * Math.abs(Jn) >= stickMag) {
    Jt = [stickScale * ut[0], stickScale * ut[1], stickScale * ut[2]]
    mode = 'stick'
  } else {
    const k = -p.muK * Math.abs(Jn) / uts
    Jt = [k * ut[0], k * ut[1], k * ut[2]]
    mode = 'slide'
  }
  const J: V3 = [Jn * n[0] + Jt[0], Jn * n[1] + Jt[1], Jn * n[2] + Jt[2]]
  a.vx += J[0]; a.vy += J[1]; a.vz += J[2]
  b.vx -= J[0]; b.vy -= J[1]; b.vz -= J[2]
  const nxJ = cross(n, J)
  const kw = 1 / (p.beta * p.radius)
  a.wx += kw * nxJ[0]; a.wy += kw * nxJ[1]; a.wz += kw * nxJ[2]
  b.wx += kw * nxJ[0]; b.wy += kw * nxJ[1]; b.wz += kw * nxJ[2]
  return mode
}
