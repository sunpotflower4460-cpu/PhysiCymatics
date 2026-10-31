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
  // Changed slot span of the last applied frame(s), for partial GPU
  // buffer updates; -1 when nothing changed. Reset by consumeSpan().
  private spanLo = -1
  private spanHi = -1

  reset(n: number): void {
    this.xy = new Float32Array(n * 2)
    this.alive = new Uint8Array(n)
    this.t = 0
    this.spanLo = n > 0 ? 0 : -1
    this.spanHi = n > 0 ? n - 1 : -1
  }

  consumeSpan(): [number, number] {
    const r: [number, number] = [this.spanLo, this.spanHi]
    this.spanLo = -1; this.spanHi = -1
    return r
  }

  /** Apply a (possibly sharded) delta frame. offset is the shard's global
   *  start index; total is the global particle count. A full frame from a
   *  shard only rewrites its own range. */
  apply(frame: TransientFrame, offset = 0, total?: number): void {
    const n = total ?? frame.initialCount
    if (this.xy.length !== n * 2) this.reset(n)
    if (frame.full) {
      // Clear only this frame's range before applying its living slots.
      for (let i = offset; i < offset + frame.initialCount; i++) this.alive[i] = 0
      if (frame.initialCount > 0) this.touch(offset, offset + frame.initialCount - 1)
    }
    for (let k = 0; k < frame.idx.length; k++) {
      const i = offset + frame.idx[k]
      this.xy[2 * i] = frame.xy[2 * k]
      this.xy[2 * i + 1] = frame.xy[2 * k + 1]
      this.alive[i] = 1
      this.touch(i, i)
    }
    for (let k = 0; k < frame.dead.length; k++) {
      const i = offset + frame.dead[k]
      this.alive[i] = 0
      this.touch(i, i)
    }
    this.t = frame.t
  }

  private touch(lo: number, hi: number): void {
    if (this.spanLo < 0 || lo < this.spanLo) this.spanLo = lo
    if (hi > this.spanHi) this.spanHi = hi
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
