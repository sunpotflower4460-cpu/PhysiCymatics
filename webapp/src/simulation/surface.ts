// surface.ts — signed, SI, time-dependent plate surface for mode (0,1).
//
// PUT-IN     : 4096-point signed R / dR table generated from
//              reference/sand_bounce.make_mode(0,0) (scripts/gen_sand_web_data.py);
//              w = A*W*sin(w t) with A = Gamma*g/omega^2. Fast lookup identical
//              to c002b (index = min(int(r*(ngrid-1)/a), ngrid-1), no interpolation).
// EMERGED    : nothing here — this only evaluates the prescribed surface.
// CLAIM-TIER : measured (within the c002 numeric model).
// FLOORS     : single prescribed sinusoid; pack's normalized |field| preview
//              (physics/plate.ts) is NOT used for contact physics.
import modeData from '../data/sand_mode01.json'

export const G = 9.81

interface ModeData {
  m: number; f_hz: number; a_m: number; ngrid: number
  r_nodes_over_a: number[]; R: number[]; dR_drho: number[]
}
const MD = modeData as unknown as ModeData

export class SingleModeSurface {
  readonly a = MD.a_m
  readonly m = MD.m
  readonly fHz = MD.f_hz
  readonly omega = 2 * Math.PI * MD.f_hz
  readonly gamma: number
  readonly A: number
  readonly rNode = MD.r_nodes_over_a.length ? MD.r_nodes_over_a[0] * MD.a_m : NaN
  private readonly ng1 = MD.ngrid - 1
  private readonly R = Float64Array.from(MD.R)
  private readonly dR = Float64Array.from(MD.dR_drho)

  constructor(gamma: number) {
    this.gamma = gamma
    this.A = (gamma * G) / this.omega ** 2
  }

  /** (W, dW/dx, dW/dy) with gradients wrt metres (includes 1/a). */
  wGrad(xp: number, yp: number, out: [number, number, number]): void {
    const r = Math.hypot(xp, yp)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const R = this.R[idx]
    const dR = this.dR[idx] / this.a
    const rs = r < 1e-12 ? 1e-12 : r
    const c = xp / rs
    const s = yp / rs
    if (this.m === 0) {
      out[0] = R; out[1] = dR * c; out[2] = dR * s
      return
    }
    const th = Math.atan2(yp, xp)
    const cm = Math.cos(this.m * th)
    const sm = Math.sin(this.m * th)
    const wr = dR * cm
    const wt = (-this.m * R * sm) / rs
    out[0] = R * cm; out[1] = wr * c - wt * s; out[2] = wr * s + wt * c
  }

  /** (h, hx, hy, uz, W) at trig values s=sin(w t), c=cos(w t) — c002b order. */
  surf(xp: number, yp: number, s: number, c: number, out: Float64Array): void {
    const g: [number, number, number] = [0, 0, 0]
    this.wGrad(xp, yp, g)
    out[0] = this.A * g[0] * s
    out[1] = this.A * g[1] * s
    out[2] = this.A * g[2] * s
    out[3] = this.A * this.omega * g[0] * c
    out[4] = g[0]
  }

  sample(xM: number, yM: number, tS: number) {
    const g: [number, number, number] = [0, 0, 0]
    this.wGrad(xM, yM, g)
    const s = Math.sin(this.omega * tS)
    const c = Math.cos(this.omega * tS)
    return {
      inside: Math.hypot(xM, yM) <= this.a,
      w_m: this.A * g[0] * s,
      wt_ms: this.A * this.omega * g[0] * c,
      wtt_ms2: -this.A * this.omega ** 2 * g[0] * s,
      wx: this.A * g[1] * s,
      wy: this.A * g[2] * s,
    }
  }
}
