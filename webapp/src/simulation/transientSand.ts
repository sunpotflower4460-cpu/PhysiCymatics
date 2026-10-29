// transientSand.ts — sand on the transient plate (TS port of
// reference/plate_transient.TransientSandSim). Contact order as c002;
// detach when wtt < -g; stick when |vn2| < 0.3*vScale (reduces to c002).
import { G } from './surface'
import { TransientPlate } from './plateTransient'
import { initialState } from './sand'
import { DEFAULT_PARAMS, type ParticleState, type SandParams } from './types'
import type { PerformanceScore } from './sequence'

export interface TransientStats {
  t: number; alive: number; sitting: number; nBounce: number
  stepsPerSec: number; freqHz: number | null; ampUm: number; eventIndex: number
}

const S5 = new Float64Array(5)

export class TransientSandSim {
  t = 0
  nBounce = 0
  constructor(
    readonly plate: TransientPlate,
    readonly state: ParticleState,
    readonly p: SandParams = DEFAULT_PARAMS,
  ) {}

  get nAlive(): number { let s = 0; for (const v of this.state.alive) s += v; return s }

  dtNow(): number {
    const ts = this.plate.scoreTime(this.t)
    const f = this.plate.signal.freqAt(ts) ?? this.plate.signal.freqAt(ts + 1e-9) ?? 152
    return 1 / (this.p.spp * Math.max(f, 20))
  }

  step(): void {
    const { plate: pl, state: st, p } = this
    const dt = this.dtNow(); const t = this.t; const tNew = t + dt
    pl.beginStep(t, tNew)
    const vStick = p.vStickFrac * pl.vScale()
    const n = st.id.length
    for (let i = 0; i < n; i++) {
      if (!st.alive[i] || st.sitting[i]) continue
      const xf = st.x[i] + st.vx[i] * dt
      const yf = st.y[i] + st.vy[i] * dt
      const zf = st.z[i] + st.vz[i] * dt - 0.5 * G * dt * dt
      const vzf = st.vz[i] - G * dt
      pl.surfAt(xf, yf, tNew, S5)
      if (zf >= S5[0]) { st.x[i] = xf; st.y[i] = yf; st.z[i] = zf; st.vz[i] = vzf; continue }
      const x0 = st.x[i]; const y0 = st.y[i]; const z0 = st.z[i]
      const vx0 = st.vx[i]; const vy0 = st.vy[i]; const vz0 = st.vz[i]
      let t0 = t; let t1 = tNew
      for (let k = 0; k < p.bisectIters; k++) {
        const tm = 0.5 * (t0 + t1); const ddt = tm - t
        pl.surfAt(x0 + vx0 * ddt, y0 + vy0 * ddt, tm, S5)
        if (z0 + vz0 * ddt - 0.5 * G * ddt * ddt < S5[0]) t1 = tm; else t0 = tm
      }
      const tc = 0.5 * (t0 + t1); const ddt = tc - t
      const xc = x0 + vx0 * ddt; const yc = y0 + vy0 * ddt
      const vzc = vz0 - G * ddt
      pl.surfAt(xc, yc, tc, S5)
      const hc = S5[0]; const hx = S5[1]; const hy = S5[2]; const uz = S5[3]
      const nn = Math.sqrt(hx * hx + hy * hy + 1)
      const nx = -hx / nn; const ny = -hy / nn; const nz = 1 / nn
      const rvz = vzc - uz
      const vn = vx0 * nx + vy0 * ny + rvz * nz
      if (vn >= 0) continue
      const k1 = (1 + p.e) * vn
      const r2x = vx0 - k1 * nx; const r2y = vy0 - k1 * ny; const r2z = rvz - k1 * nz
      const vn2 = r2x * nx + r2y * ny + r2z * nz
      st.x[i] = xc; st.y[i] = yc; st.z[i] = hc + 1e-9
      st.vx[i] = vn2 * nx + p.beta * (r2x - vn2 * nx)
      st.vy[i] = vn2 * ny + p.beta * (r2y - vn2 * ny)
      st.vz[i] = vn2 * nz + p.beta * (r2z - vn2 * nz) + uz
      this.nBounce++
      if (Math.abs(vn2) < vStick) { st.sitting[i] = 1; st.vx[i] = 0; st.vy[i] = 0 }
    }
    for (let i = 0; i < n; i++) {
      if (!st.alive[i] || !st.sitting[i]) continue
      pl.surfAt(st.x[i], st.y[i], tNew, S5)
      st.z[i] = S5[0]
      if (S5[4] < -G) { st.sitting[i] = 0; st.vz[i] = S5[3]; st.z[i] = S5[0] + 1e-9 }
    }
    for (let i = 0; i < n; i++) {
      if (st.alive[i] && Math.hypot(st.x[i], st.y[i]) > pl.a) st.alive[i] = 0
    }
    this.t = tNew
  }

  stats(stepsPerSec = 0): TransientStats {
    let seated = 0
    for (let i = 0; i < this.state.id.length; i++) {
      if (this.state.alive[i] && this.state.sitting[i]) seated++
    }
    const alive = this.nAlive
    const ev = this.plate.signal.eventAt(this.plate.scoreTime(this.t))
    return {
      t: Math.round(this.t * 100) / 100, alive,
      sitting: alive ? seated / alive : NaN, nBounce: this.nBounce, stepsPerSec,
      freqHz: this.plate.signal.freqAt(this.plate.scoreTime(this.t)), ampUm: this.plate.ampScale() * 1e6,
      eventIndex: ev ? ev.index : -1,
    }
  }

  snapshot(stepsPerSec = 0) {
    const st = this.state
    const xy = new Float32Array(this.nAlive * 2)
    let k = 0
    for (let i = 0; i < st.id.length; i++) {
      if (!st.alive[i]) continue
      xy[k++] = st.x[i]; xy[k++] = st.y[i]
    }
    return { t: this.t, xy, aliveCount: k / 2, initialCount: st.id.length, stats: this.stats(stepsPerSec) }
  }
}

export function createTransientSim(
  score: PerformanceScore, count: number, seed: number, driveXY: [number, number] = [0, 0],
): TransientSandSim {
  const plate = new TransientPlate(score, driveXY)
  return new TransientSandSim(plate, initialState(plate.a, count, seed))
}

/** Force giving ~860 um steady amplitude at the 152 Hz mode (center drive). */
export function calibratedForceN(driveXY: [number, number] = [0, 0]): number {
  const probe = new TransientPlate(
    { schema: 'physicymatics.score.v1', events: [{ start_s: 0, end_s: 1, frequency_hz: 152, forcePeak_N: 1, attack_s: 0, release_s: 0 }], duration_s: 1, initialPhase_rad: 0, audioGain: 0.2 },
    driveXY)
  return probe.forceForAmplitude(probe.modeIndex152(), 860e-6)
}
