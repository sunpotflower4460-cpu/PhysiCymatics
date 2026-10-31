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

const VS = `
attribute vec2 aPos;
attribute float aAlive;
uniform float uSx;
uniform float uSy;
uniform float uPointPx;
void main() {
  if (aAlive < 0.5) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
  } else {
    gl_Position = vec4(aPos.x * uSx, aPos.y * uSy, 0.0, 1.0);
  }
  gl_PointSize = uPointPx;
}`
const FS = `
precision mediump float;
uniform vec4 uColor;
void main() { gl_FragColor = uColor; }`

function hexRgb(hex: string): [number, number, number] {
  return [parseInt(hex.slice(1, 3), 16) / 255, parseInt(hex.slice(3, 5), 16) / 255, parseInt(hex.slice(5, 7), 16) / 255]
}

/** M6-4C WebGL point renderer over the same FrameCache data and the same
 *  projection as the Canvas 2D path. Positions are the physical xy slots;
 *  dead slots are clipped away in the vertex shader (no added points).
 *  Buffer updates use the cache's changed span: a partial update when the
 *  span is narrow, a full re-upload when it covers most slots. Returns
 *  false from draw() when the context is unavailable/lost so the caller
 *  can fall back to Canvas 2D. */
export class GlParticleRenderer {
  private gl: WebGLRenderingContext | null
  private prog: WebGLProgram | null = null
  private posBuf: WebGLBuffer | null = null
  private aliveBuf: WebGLBuffer | null = null
  private circleBuf: WebGLBuffer | null = null
  private bufN = 0
  private aliveF = new Float32Array(0)
  private scratchAlive = new Float32Array(0)
  lost = false

  constructor(readonly canvas: HTMLCanvasElement) {
    this.gl = (canvas.getContext('webgl2') || canvas.getContext('webgl')) as WebGLRenderingContext | null
    canvas.addEventListener('webglcontextlost', (e) => { e.preventDefault(); this.lost = true })
    if (this.gl) this.initGl()
  }

  get available(): boolean { return !!this.gl && !this.lost }

  private initGl(): void {
    const gl = this.gl!
    const compile = (type: number, src: string) => {
      const s = gl.createShader(type)!
      gl.shaderSource(s, src); gl.compileShader(s)
      return s
    }
    const prog = gl.createProgram()!
    gl.attachShader(prog, compile(gl.VERTEX_SHADER, VS))
    gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, FS))
    gl.linkProgram(prog)
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) { this.gl = null; return }
    this.prog = prog
    this.posBuf = gl.createBuffer()
    this.aliveBuf = gl.createBuffer()
    this.circleBuf = gl.createBuffer()
    const pts: number[] = []
    for (let k = 0; k < 128; k++) {
      const th = (k / 128) * 2 * Math.PI
      pts.push(Math.cos(th), Math.sin(th))
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, this.circleBuf)
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(pts), gl.STATIC_DRAW)
  }

  draw(cache: FrameCache, plateRadiusM: number): boolean {
    const gl = this.gl
    if (!gl || this.lost || !this.prog) return false
    const { width: w, height: h } = this.canvas
    const n = cache.alive.length
    gl.viewport(0, 0, w, h)
    gl.clearColor(0, 0, 0, 1)
    gl.clear(gl.COLOR_BUFFER_BIT)
    gl.useProgram(this.prog)
    const scale = (Math.min(w, h) / 2) * 0.92 / plateRadiusM
    const uSx = (scale * 2) / w
    const uSy = (scale * 2) / h
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uSx'), uSx)
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uSy'), uSy)
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uPointPx'), PARTICLE_POINT_PX)
    const aPos = gl.getAttribLocation(this.prog, 'aPos')
    const aAlive = gl.getAttribLocation(this.prog, 'aAlive')
    const uColor = gl.getUniformLocation(this.prog, 'uColor')
    // Plate rim (unit circle scaled by plateRadiusM via aPos transform:
    // feed meter coordinates on the rim).
    gl.uniform4f(uColor, ...hexRgb(C_PLATE_EDGE), 1)
    gl.bindBuffer(gl.ARRAY_BUFFER, this.circleBuf)
    gl.enableVertexAttribArray(aPos)
    gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0)
    gl.disableVertexAttribArray(aAlive)
    gl.vertexAttrib1f(aAlive, 1)
    // Scale rim points (unit circle) by drawing with a temporary uniform
    // multiplier folded into uSx/uSy would distort; instead rim buffer is
    // unitless and we scale uniforms by plateRadiusM for this draw.
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uSx'), uSx * plateRadiusM)
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uSy'), uSy * plateRadiusM)
    gl.drawArrays(gl.LINE_LOOP, 0, 128)
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uSx'), uSx)
    gl.uniform1f(gl.getUniformLocation(this.prog, 'uSy'), uSy)
    if (n === 0) return true
    // Particle buffers.
    if (this.bufN !== n) {
      this.bufN = n
      this.aliveF = new Float32Array(n)
      this.scratchAlive = new Float32Array(n)
      for (let i = 0; i < n; i++) this.aliveF[i] = cache.alive[i]
      gl.bindBuffer(gl.ARRAY_BUFFER, this.posBuf)
      gl.bufferData(gl.ARRAY_BUFFER, cache.xy, gl.DYNAMIC_DRAW)
      gl.bindBuffer(gl.ARRAY_BUFFER, this.aliveBuf)
      gl.bufferData(gl.ARRAY_BUFFER, this.aliveF, gl.DYNAMIC_DRAW)
      cache.consumeSpan()
    } else {
      const [lo, hi] = cache.consumeSpan()
      if (lo >= 0) {
        const wide = hi - lo + 1 > 0.6 * n
        gl.bindBuffer(gl.ARRAY_BUFFER, this.posBuf)
        if (wide) {
          gl.bufferData(gl.ARRAY_BUFFER, cache.xy, gl.DYNAMIC_DRAW)
          for (let i = 0; i < n; i++) this.aliveF[i] = cache.alive[i]
          gl.bindBuffer(gl.ARRAY_BUFFER, this.aliveBuf)
          gl.bufferData(gl.ARRAY_BUFFER, this.aliveF, gl.DYNAMIC_DRAW)
        } else {
          gl.bufferSubData(gl.ARRAY_BUFFER, lo * 2 * 4, cache.xy.subarray(lo * 2, (hi + 1) * 2))
          const sub = this.scratchAlive.subarray(0, hi - lo + 1)
          for (let i = lo; i <= hi; i++) sub[i - lo] = cache.alive[i]
          gl.bindBuffer(gl.ARRAY_BUFFER, this.aliveBuf)
          gl.bufferSubData(gl.ARRAY_BUFFER, lo * 4, sub)
        }
      }
    }
    gl.uniform4f(uColor, ...hexRgb(C_GRAIN), 1)
    gl.bindBuffer(gl.ARRAY_BUFFER, this.posBuf)
    gl.enableVertexAttribArray(aPos)
    gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0)
    gl.bindBuffer(gl.ARRAY_BUFFER, this.aliveBuf)
    gl.enableVertexAttribArray(aAlive)
    gl.vertexAttribPointer(aAlive, 1, gl.FLOAT, false, 0, 0)
    gl.drawArrays(gl.POINTS, 0, n)
    return true
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
