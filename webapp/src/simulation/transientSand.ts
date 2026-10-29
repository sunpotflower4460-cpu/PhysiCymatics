// transientSand.ts — sand on the transient plate (TS port of
// reference/plate_transient.TransientSandSim, M6-1 event-accurate version).
// Same algorithm/constants as Python: per-particle time cursor, first-contact
// search (8 gap samples + 30 bisections), residual time after bounce/detach,
// exact detach root (wtt is linear in a Hermite segment), event cap fallback.
// Detection floor: a dip entirely between adjacent samples can be missed.
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
export const CONTACT_SAMPLES = 8
export const EVENT_BISECT_ITERS = 30
export const MAX_EVENT_ROUNDS = 8

export class TransientSandSim {
  t = 0
  nBounce = 0
  nDetach = 0
  nEventCap = 0
  lastDetachT: number | null = null
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

  step(dtIn?: number): void {
    const { plate: pl, state: st, p } = this
    const dt = dtIn ?? this.dtNow(); const t = this.t; const tNew = t + dt
    pl.beginStep(t, tNew)
    const vStick = p.vStickFrac * pl.vScale()
    const n = st.id.length
    for (let i = 0; i < n; i++) {
      if (!st.alive[i]) continue
      let cur = t
      let rounds = 0
      while (cur < tNew - 1e-15 && rounds < MAX_EVENT_ROUNDS) {
        rounds++
        if (st.sitting[i]) {
          // wtt is linear in tau at fixed (x,y): exact detach root.
          pl.surfAt(st.x[i], st.y[i], cur, S5); const w0 = S5[4]
          pl.surfAt(st.x[i], st.y[i], tNew, S5); const w1 = S5[4]; const h1 = S5[0]
          if (w0 >= -G && w1 >= -G) { st.z[i] = h1; cur = tNew; break }
          const frac = w0 < -G ? 0 : Math.min(1, Math.max(0, (-G - w0) / (w1 - w0)))
          const td = cur + frac * (tNew - cur)
          pl.surfAt(st.x[i], st.y[i], td, S5)
          st.sitting[i] = 0; st.z[i] = S5[0] + 1e-9; st.vz[i] = S5[3]
          cur = td; this.nDetach++; this.lastDetachT = td
          continue
        }
        // Flying: conservative skip if the ballistic minimum stays above
        // the segment surface bound; otherwise sample + bisect first crossing.
        const x0 = st.x[i]; const y0 = st.y[i]; const z0 = st.z[i]
        const vx0 = st.vx[i]; const vy0 = st.vy[i]; const vz0 = st.vz[i]
        const rem = tNew - cur
        {
          const zEnd = z0 + vz0 * rem - 0.5 * G * rem * rem
          const tv = Math.min(Math.max(vz0 / G, 0), rem)
          const zMin = Math.min(z0, zEnd, z0 + vz0 * tv - 0.5 * G * tv * tv)
          if (zMin > pl.hBound) {
            st.x[i] = x0 + vx0 * rem; st.y[i] = y0 + vy0 * rem
            st.z[i] = zEnd; st.vz[i] = vz0 - G * rem
            cur = tNew; break
          }
        }
        let firstBad = -1
        for (let k = 1; k <= CONTACT_SAMPLES; k++) {
          const tau = cur + rem * (k / CONTACT_SAMPLES); const d = tau - cur
          pl.surfAt(x0 + vx0 * d, y0 + vy0 * d, tau, S5)
          if (z0 + vz0 * d - 0.5 * G * d * d - S5[0] <= 0) { firstBad = k; break }
        }
        if (firstBad < 0) {
          const d = rem
          st.x[i] = x0 + vx0 * d; st.y[i] = y0 + vy0 * d
          st.z[i] = z0 + vz0 * d - 0.5 * G * d * d; st.vz[i] = vz0 - G * d
          cur = tNew; break
        }
        let lo = cur + rem * ((firstBad - 1) / CONTACT_SAMPLES)
        let hi = cur + rem * (firstBad / CONTACT_SAMPLES)
        for (let k = 0; k < EVENT_BISECT_ITERS; k++) {
          const tm = 0.5 * (lo + hi); const d = tm - cur
          pl.surfAt(x0 + vx0 * d, y0 + vy0 * d, tm, S5)
          if (z0 + vz0 * d - 0.5 * G * d * d <= S5[0]) hi = tm; else lo = tm
        }
        const tc = 0.5 * (lo + hi); const d = tc - cur
        const xc = x0 + vx0 * d; const yc = y0 + vy0 * d
        const vzc = vz0 - G * d
        pl.surfAt(xc, yc, tc, S5)
        const hc = S5[0]; const hx = S5[1]; const hy = S5[2]; const uz = S5[3]
        const nn = Math.sqrt(hx * hx + hy * hy + 1)
        const nx = -hx / nn; const ny = -hy / nn; const nz = 1 / nn
        const rvz = vzc - uz
        const vn = vx0 * nx + vy0 * ny + rvz * nz
        st.x[i] = xc; st.y[i] = yc; st.z[i] = hc + 1e-9
        if (vn >= 0) { st.vz[i] = vzc; cur = tc; continue }
        const k1 = (1 + p.e) * vn
        const r2x = vx0 - k1 * nx; const r2y = vy0 - k1 * ny; const r2z = rvz - k1 * nz
        const vn2 = r2x * nx + r2y * ny + r2z * nz
        st.vx[i] = vn2 * nx + p.beta * (r2x - vn2 * nx)
        st.vy[i] = vn2 * ny + p.beta * (r2y - vn2 * ny)
        st.vz[i] = vn2 * nz + p.beta * (r2z - vn2 * nz) + uz
        this.nBounce++
        cur = tc
        if (Math.abs(vn2) < vStick) { st.sitting[i] = 1; st.vx[i] = 0; st.vy[i] = 0 }
      }
      if (cur < tNew - 1e-15) {
        // Event-cap fallback: inelastic micro-bounces accumulate in finite
        // time; seat on the surface (counted regularization, as Python).
        this.nEventCap++
        pl.surfAt(st.x[i], st.y[i], tNew, S5)
        st.z[i] = S5[0]; st.vz[i] = S5[3]; st.sitting[i] = 1
        st.vx[i] = 0; st.vy[i] = 0
      }
    }
    for (let i = 0; i < n; i++) {
      if (st.alive[i] && Math.hypot(st.x[i], st.y[i]) > pl.a) st.alive[i] = 0
    }
    this.t = tNew
  }

  advanceTo(tTarget: number): void {
    while (this.t < tTarget - 1e-15) this.step(Math.min(this.dtNow(), tTarget - this.t))
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
