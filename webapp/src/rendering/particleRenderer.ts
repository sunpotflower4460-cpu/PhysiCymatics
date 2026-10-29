// particleRenderer.ts — reads particle snapshots, draws them. Nothing else.
// Top-down orthographic projection of the actual (x, y) state (M2 scope).
// No highlighting by position, no trails, no added points. Point size is a
// display choice, recorded here; it is not the physical grain diameter.
import type { SandSnapshot } from '../simulation/types'

export const PARTICLE_POINT_PX = 2
const C_BG = '#000000'
const C_PLATE_EDGE = '#3a3a3a'
const C_GRAIN = '#f5e9c8'

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
