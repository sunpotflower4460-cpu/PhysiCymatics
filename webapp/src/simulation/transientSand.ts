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
  sleepingFrac: number; skippedUpdates: number
}

const S5 = new Float64Array(5)
export const CONTACT_SAMPLES = 8
export const EVENT_BISECT_ITERS = 30
export const MAX_EVENT_ROUNDS = 8
// Sleep (M6-2): computation-management only, not a physical state.
// Bounds are derived in docs/m6_2_event_driven_design.md; margins cover
// integration error of the modal RK4 (validated empirically vs baseline).
export const SLEEP_WINDOW_S = 0.5
// Lookahead ranges are exact for the RK4+Hermite scheme; margin only
// covers floating-point rounding in the range computation itself.
export const SLEEP_ACCEL_MARGIN = 0.01 // m/s^2
export const SLEEP_HEIGHT_MARGIN = 1e-6 // m

export interface SleepOptions { enabled: boolean }

export class TransientSandSim {
  t = 0
  nBounce = 0
  nDetach = 0
  nEventCap = 0
  skippedUpdates = 0
  sleepEvaluations = 0
  wakeCount = 0
  lastDetachT: number | null = null
  readonly sleeping: Uint8Array
  readonly wakeT: Float64Array
  readonly sleepT0: Float64Array
  private nextEvalT: Float64Array
  private lookCache: { end: number; qMin: Float64Array; qMax: Float64Array; aMin: Float64Array; aMax: Float64Array } | null = null
  private phiBuf: Float64Array
  constructor(
    readonly plate: TransientPlate,
    readonly state: ParticleState,
    readonly p: SandParams = DEFAULT_PARAMS,
    readonly sleep: SleepOptions = { enabled: true },
  ) {
    const n = state.id.length
    this.sleeping = new Uint8Array(n)
    this.wakeT = new Float64Array(n)
    this.sleepT0 = new Float64Array(n)
    this.nextEvalT = new Float64Array(n)
    this.phiBuf = new Float64Array(plate.ndof)
  }

  /** Invalidate all sleep schedules (input change). State itself is kept;
   *  flying sleepers are materialized at the current time first. */
  wakeAll(): void {
    const st = this.state
    for (let i = 0; i < st.id.length; i++) {
      if (!this.sleeping[i]) continue
      this.materialize(i, this.t)
      this.sleeping[i] = 0
    }
    this.lookCache = null
  }

  private ensureLookahead(fromT: number): void {
    if (this.lookCache && fromT < this.lookCache.end - SLEEP_WINDOW_S * 0.5) return
    const pl = this.plate
    const dtAt = (t: number): number => {
      const ts = pl.scoreTime(t)
      const f = pl.signal.freqAt(ts) ?? pl.signal.freqAt(ts + 1e-9) ?? 152
      return 1 / (this.p.spp * Math.max(f, 20))
    }
    this.lookCache = pl.lookaheadRanges(fromT, SLEEP_WINDOW_S, dtAt)
  }

  private nominalDtAt(tAbs: number): number {
    const ts = this.plate.scoreTime(tAbs)
    const f = this.plate.signal.freqAt(ts) ?? this.plate.signal.freqAt(ts + 1e-9) ?? 152
    return 1 / (this.p.spp * Math.max(f, 20))
  }

  private materialize(i: number, tAbs: number): void {
    const st = this.state
    if (st.sitting[i]) {
      // z is restored by the caller via surfAt when a segment is active.
      return
    }
    // Replay the ballistic updates in nominal-dt chunks so floating-point
    // accumulation matches the baseline stepwise path.
    while (this.sleepT0[i] < tAbs - 1e-15) {
      const d = Math.min(this.nominalDtAt(this.sleepT0[i]), tAbs - this.sleepT0[i])
      st.x[i] += st.vx[i] * d; st.y[i] += st.vy[i] * d
      st.z[i] += st.vz[i] * d - 0.5 * G * d * d; st.vz[i] -= G * d
      this.sleepT0[i] += d
    }
  }

  get nSleeping(): number { let s = 0; for (let i = 0; i < this.sleeping.length; i++) if (this.state.alive[i]) s += this.sleeping[i]; return s }

  // M6-2R: sleep bounds come from the lookahead cache, not force bounds.
  // Kept for worker API compatibility; invalidates the lookahead cache so
  // changed drive conditions are re-looked-ahead (callers also wakeAll()).
  refreshForceBounds(): void { this.lookCache = null }

  get nAlive(): number { let s = 0; for (const v of this.state.alive) s += v; return s }

  dtNow(): number {
    const ts = this.plate.scoreTime(this.t)
    const f = this.plate.signal.freqAt(ts) ?? this.plate.signal.freqAt(ts + 1e-9) ?? 152
    return 1 / (this.p.spp * Math.max(f, 20))
  }

  private trySleepSitting(i: number, tNew: number): void {
    // Rolling (cRoll>0) moves seated particles horizontally: never sleep
    // them on a vertical-acceleration argument alone (Codex #4410299745).
    if (this.p.cRoll > 0) return
    const pl = this.plate; const st = this.state
    this.sleepEvaluations++
    this.ensureLookahead(tNew)
    const lc = this.lookCache!
    pl.phiAll(st.x[i], st.y[i], this.phiBuf)
    // Exact lower bound of Hermite wtt over the lookahead window.
    let lb = 0
    for (let j = 0; j < pl.ndof; j++) {
      const phi = this.phiBuf[j]
      lb += Math.min(phi * lc.aMin[j], phi * lc.aMax[j])
    }
    if (G + lb < SLEEP_ACCEL_MARGIN) return
    if (!(lc.end > tNew + 8 * this.dtNow())) return
    this.sleeping[i] = 1
    this.wakeT[i] = lc.end
    this.sleepT0[i] = tNew
  }

  private trySleepFlying(i: number, tNew: number): void {
    const pl = this.plate; const st = this.state
    this.sleepEvaluations++
    this.ensureLookahead(tNew)
    const lc = this.lookCache!
    const z0 = st.z[i]; const vz0 = st.vz[i]
    // Exact Hermite height bound over the window: |h| <= sum max|q_j|.
    let hBoundW = 0
    for (let j = 0; j < pl.ndof; j++) hBoundW += Math.max(Math.abs(lc.qMin[j]), Math.abs(lc.qMax[j]))
    const safe = (d: number): boolean => {
      const zEnd = z0 + vz0 * d - 0.5 * G * d * d
      const tv = Math.min(Math.max(vz0 / G, 0), d)
      const zMin = Math.min(z0, zEnd, z0 + vz0 * tv - 0.5 * G * tv * tv)
      return zMin > hBoundW + SLEEP_HEIGHT_MARGIN
    }
    if (!safe(0)) return
    let lo = 0; let hi = Math.min(SLEEP_WINDOW_S, lc.end - tNew)
    if (!safe(hi)) {
      for (let k = 0; k < 40; k++) {
        const mid = 0.5 * (lo + hi)
        if (safe(mid)) lo = mid; else hi = mid
      }
    } else lo = hi
    // Edge-exit time for straight xy motion from inside the disk.
    let wake = tNew + lo
    const px = st.x[i]; const py = st.y[i]; const vx = st.vx[i]; const vy = st.vy[i]
    const A = vx * vx + vy * vy
    if (A > 0) {
      const B = 2 * (px * vx + py * vy)
      const C = px * px + py * py - pl.a * pl.a
      const disc = B * B - 4 * A * C
      if (disc >= 0) {
        const root = (-B + Math.sqrt(disc)) / (2 * A)
        if (root > 0) wake = Math.min(wake, tNew + root)
      }
    }
    if (!(wake > tNew + 8 * this.dtNow())) return
    this.sleeping[i] = 1
    this.wakeT[i] = wake
    this.sleepT0[i] = tNew
  }

  step(dtIn?: number): void {
    const { plate: pl, state: st, p } = this
    // Lookahead ranges assume the nominal dt sequence; a non-nominal dt
    // breaks segment alignment, so wake everything and skip sleep.
    const nominalDt = this.dtNow()
    const customDt = dtIn !== undefined && Math.abs(dtIn - nominalDt) > 1e-12 * nominalDt
    if (customDt && this.sleep.enabled) this.wakeAll()
    const dt = dtIn ?? nominalDt; const t = this.t; const tNew = t + dt
    pl.beginStep(t, tNew)
    // Inelastic-collapse regularization: contacts place at h+1e-9 m; a
    // rebound apex below 1e-8 m above the surface is unresolvable and
    // sustains a numerical limit cycle; stick floor sqrt(2*g*9e-9).
    // Resolution limit, not an event-count rule.
    const vStick = Math.max(p.vStickFrac * pl.vScale(), Math.sqrt(2 * G * 9e-9))
    const n = st.id.length
    for (let i = 0; i < n; i++) {
      if (!st.alive[i]) continue
      let cur = t
      if (this.sleeping[i]) {
        if (this.wakeT[i] > tNew) { this.skippedUpdates++; continue }
        // Wake at the STEP START (<= wakeT, still inside the guaranteed
        // interval) so this step is processed exactly as the baseline.
        this.sleeping[i] = 0; this.wakeCount++
        cur = t
        if (st.sitting[i]) {
          pl.surfAt(st.x[i], st.y[i], t, S5); st.z[i] = S5[0]
        } else {
          this.materialize(i, t)
        }
      }
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
        // Event cap reached (counted): subdivide the remaining interval
        // into 8 legacy single-contact substeps instead of forcing a
        // physical state change. Seating still occurs only via the stick
        // criterion. True inelastic collapse is thereby approximated by
        // ever-smaller bounces, not by fiat.
        this.nEventCap++
        const subN = 8
        for (let sIdx = 0; sIdx < subN && cur < tNew - 1e-15; sIdx++) {
          const subEnd = cur + (tNew - cur) / (subN - sIdx)
          if (st.sitting[i]) {
            pl.surfAt(st.x[i], st.y[i], subEnd, S5)
            st.z[i] = S5[0]
            if (S5[4] < -G) { st.sitting[i] = 0; st.vz[i] = S5[3]; st.z[i] = S5[0] + 1e-9; this.nDetach++ }
            cur = subEnd
            continue
          }
          // Within the substep, keep processing contacts segment by
          // segment: a ballistic endpoint is accepted only at/above the
          // surface, and after each contact the search resumes from tc
          // (never an unchecked remainder). The stick floor bounds
          // collapse iteration; SEG_CAP is a safety valve that leaves
          // the particle at its last contact state, not below surface.
          const SEG_CAP = 256
          for (let seg = 0; seg < SEG_CAP && cur < subEnd - 1e-15 && !st.sitting[i]; seg++) {
            const d = subEnd - cur
            const xf = st.x[i] + st.vx[i] * d; const yf = st.y[i] + st.vy[i] * d
            const zf = st.z[i] + st.vz[i] * d - 0.5 * G * d * d
            pl.surfAt(xf, yf, subEnd, S5)
            if (zf >= S5[0]) {
              st.x[i] = xf; st.y[i] = yf; st.z[i] = zf; st.vz[i] -= G * d; cur = subEnd; break
            }
            const x0 = st.x[i]; const y0 = st.y[i]; const z0 = st.z[i]
            const vx0 = st.vx[i]; const vy0 = st.vy[i]; const vz0 = st.vz[i]
            let lo2 = cur; let hi2 = subEnd
            for (let k = 0; k < EVENT_BISECT_ITERS; k++) {
              const tm = 0.5 * (lo2 + hi2); const dd = tm - cur
              pl.surfAt(x0 + vx0 * dd, y0 + vy0 * dd, tm, S5)
              if (z0 + vz0 * dd - 0.5 * G * dd * dd <= S5[0]) hi2 = tm; else lo2 = tm
            }
            const tc = 0.5 * (lo2 + hi2); const dd = tc - cur
            const xc = x0 + vx0 * dd; const yc = y0 + vy0 * dd; const vzc = vz0 - G * dd
            pl.surfAt(xc, yc, tc, S5)
            const nn2 = Math.sqrt(S5[1] * S5[1] + S5[2] * S5[2] + 1)
            const nx2 = -S5[1] / nn2; const ny2 = -S5[2] / nn2; const nz2 = 1 / nn2
            const rvz2 = vzc - S5[3]
            const vn = vx0 * nx2 + vy0 * ny2 + rvz2 * nz2
            st.x[i] = xc; st.y[i] = yc; st.z[i] = S5[0] + 1e-9
            if (vn < 0) {
              const k1 = (1 + p.e) * vn
              const r2x = vx0 - k1 * nx2; const r2y = vy0 - k1 * ny2; const r2z = rvz2 - k1 * nz2
              const vn2 = r2x * nx2 + r2y * ny2 + r2z * nz2
              st.vx[i] = vn2 * nx2 + p.beta * (r2x - vn2 * nx2)
              st.vy[i] = vn2 * ny2 + p.beta * (r2y - vn2 * ny2)
              st.vz[i] = vn2 * nz2 + p.beta * (r2z - vn2 * nz2) + S5[3]
              this.nBounce++
              if (Math.abs(vn2) < vStick) { st.sitting[i] = 1; st.vx[i] = 0; st.vy[i] = 0 }
            } else st.vz[i] = vzc
            cur = tc
          }
          if (st.sitting[i]) {
            pl.surfAt(st.x[i], st.y[i], subEnd, S5)
            st.z[i] = S5[0]
            if (S5[4] < -G) { st.sitting[i] = 0; st.vz[i] = S5[3]; st.z[i] = S5[0] + 1e-9; this.nDetach++ }
            cur = subEnd
          }
        }
      }
      if (this.sleep.enabled && !customDt && st.alive[i] && tNew >= this.nextEvalT[i]) {
        if (st.sitting[i]) this.trySleepSitting(i, tNew)
        else this.trySleepFlying(i, tNew)
        // Failed evaluations back off 16 steps (management only; the
        // particle simply stays on the baseline path meanwhile).
        if (!this.sleeping[i]) this.nextEvalT[i] = tNew + 16 * dt
      }
    }
    for (let i = 0; i < n; i++) {
      if (st.alive[i] && Math.hypot(st.x[i], st.y[i]) > pl.a) { st.alive[i] = 0; this.sleeping[i] = 0 }
    }
    this.t = tNew
  }

  /** Bring sleeping particles' stored state to the current time
   *  (analytic for flyers; surface height for sitting). Sleep flags stay. */
  materializeAll(): void {
    const st = this.state
    for (let i = 0; i < st.id.length; i++) {
      if (!this.sleeping[i] || !st.alive[i]) continue
      if (st.sitting[i]) {
        this.plate.surfAt(st.x[i], st.y[i], this.t, S5)
        st.z[i] = S5[0]
        this.sleepT0[i] = this.t
      } else {
        this.materialize(i, this.t)
      }
    }
  }

  advanceTo(tTarget: number): void {
    while (this.t < tTarget - 1e-15) this.step(Math.min(this.dtNow(), tTarget - this.t))
    this.materializeAll()
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
      sleepingFrac: alive ? this.nSleeping / alive : 0, skippedUpdates: this.skippedUpdates,
    }
  }

  snapshot(stepsPerSec = 0) {
    const st = this.state
    const xy = new Float32Array(this.nAlive * 2)
    let k = 0
    for (let i = 0; i < st.id.length; i++) {
      if (!st.alive[i]) continue
      if (this.sleeping[i] && !st.sitting[i]) {
        // Sleeping flyer: analytic position at snapshot time (no state change).
        const d = this.t - this.sleepT0[i]
        xy[k++] = st.x[i] + st.vx[i] * d; xy[k++] = st.y[i] + st.vy[i] * d
      } else {
        xy[k++] = st.x[i]; xy[k++] = st.y[i]
      }
    }
    return { t: this.t, xy, aliveCount: k / 2, initialCount: st.id.length, stats: this.stats(stepsPerSec) }
  }
}

export function createTransientSim(
  score: PerformanceScore, count: number, seed: number, driveXY: [number, number] = [0, 0],
  sleep: SleepOptions = { enabled: true },
): TransientSandSim {
  const plate = new TransientPlate(score, driveXY)
  return new TransientSandSim(plate, initialState(plate.a, count, seed), DEFAULT_PARAMS, sleep)
}

/** Force giving ~860 um steady amplitude at the 152 Hz mode (center drive). */
export function calibratedForceN(driveXY: [number, number] = [0, 0]): number {
  const probe = new TransientPlate(
    { schema: 'physicymatics.score.v1', events: [{ start_s: 0, end_s: 1, frequency_hz: 152, forcePeak_N: 1, attack_s: 0, release_s: 0 }], duration_s: 1, initialPhase_rad: 0, audioGain: 0.2 },
    driveXY)
  return probe.forceForAmplitude(probe.modeIndex152(), 860e-6)
}
