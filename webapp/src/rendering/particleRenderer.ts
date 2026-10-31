// particleRenderer.ts — reads particle snapshots, draws them. Nothing else.
// Top-down orthographic projection of the actual (x, y) state (M2 scope).
// No highlighting by position, no trails, no added points. Point size is a
// display choice, recorded here; it is not the physical grain diameter.
import type { SandSnapshot } from '../simulation/types'
import type { TransientFrame } from '../simulation/transientSand'

export const PARTICLE_POINT_PX = 2
const C_BG = '#000000'
const C_PLATE_EDGE = '#3a3a3a'
const C_GRAIN = '#f5e9c8'

/** Display-side cache for differential frames (M6-3C). Holds the last
 *  received xy per particle slot; physics state itself stays in the
 *  Worker — this cache only avoids resending unchanged coordinates. */
export class FrameCache {
  xy: Float32Array = new Float32Array(0)
  alive: Uint8Array = new Uint8Array(0)
  t = 0

  apply(frame: TransientFrame): void {
    const n = frame.initialCount
    if (frame.full || this.xy.length !== n * 2) {
      this.xy = new Float32Array(n * 2)
      this.alive = new Uint8Array(n)
    }
    for (let k = 0; k < frame.idx.length; k++) {
      const i = frame.idx[k]
      this.xy[2 * i] = frame.xy[2 * k]
      this.xy[2 * i + 1] = frame.xy[2 * k + 1]
      this.alive[i] = 1
    }
    for (let k = 0; k < frame.dead.length; k++) this.alive[frame.dead[k]] = 0
    this.t = frame.t
  }
}

export function drawParticles(
  canvas: HTMLCanvasElement, snap: SandSnapshot, plateRadiusM: number,
): void {
  const ctx = canvas.getContext('2d')
  if (!ctx) return
  const { width: w, height: h } = canvas
  ctx.fillStyle = C_BG
  ctx.fillRect(0, 0, w, h)
  const scale = (Math.min(w, h) / 2) * 0.92 / plateRadiusM
  const cx = w / 2; const cy = h / 2
  ctx.strokeStyle = C_PLATE_EDGE
  ctx.lineWidth = 1
  ctx.beginPath()
  ctx.arc(cx, cy, plateRadiusM * scale, 0, 2 * Math.PI)
  ctx.stroke()
  ctx.fillStyle = C_GRAIN
  const xy = snap.xy
  for (let i = 0; i < xy.length; i += 2) {
    ctx.fillRect(cx + xy[i] * scale, cy - xy[i + 1] * scale, PARTICLE_POINT_PX, PARTICLE_POINT_PX)
  }
}

/** Draw from the differential FrameCache (same projection/points). */
export function drawFrameCache(
  canvas: HTMLCanvasElement, cache: FrameCache, plateRadiusM: number,
): void {
  const ctx = canvas.getContext('2d')
  if (!ctx) return
  const { width: w, height: h } = canvas
  ctx.fillStyle = C_BG
  ctx.fillRect(0, 0, w, h)
  const scale = (Math.min(w, h) / 2) * 0.92 / plateRadiusM
  const cx = w / 2; const cy = h / 2
  ctx.strokeStyle = C_PLATE_EDGE
  ctx.lineWidth = 1
  ctx.beginPath()
  ctx.arc(cx, cy, plateRadiusM * scale, 0, 2 * Math.PI)
  ctx.stroke()
  ctx.fillStyle = C_GRAIN
  const { xy, alive } = cache
  for (let i = 0; i < alive.length; i++) {
    if (!alive[i]) continue
    ctx.fillRect(cx + xy[2 * i] * scale, cy - xy[2 * i + 1] * scale, PARTICLE_POINT_PX, PARTICLE_POINT_PX)
  }
}
