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

export interface SleepOptions {
  enabled: boolean
  /** M6-3A event scheduler (active set + wake heap); default on when enabled.
   *  false keeps the M6-2 full-scan loop for comparison. */
  scheduler?: boolean
  /** M6-3B region-shared sitting-sleep evaluation; default on when enabled. */
  regions?: boolean
}

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
  // M6-3A scheduler state: active set (swap-remove list) + wake min-heap.
  // Heap entries are validated against sleeping/wakeT, so stale entries
  // (after wakeAll or rescheduling) can never cause a wrong wake.
  private activeList: Int32Array
  private activePos: Int32Array
  private activeCount = 0
  private heapT = new Float64Array(1024)
  private heapI = new Int32Array(1024)
  private heapSize = 0
  private sleepingNow = 0
  // M6-3B region sharing for SITTING sleep evaluation: exact per-cell Phi
  // ranges (TransientPlate.cellPhiRanges) combined with the lookahead
  // acceleration ranges give a cell-wide lower bound valid at every point
  // of the cell. Inconclusive cells fall back to per-particle evaluation.
  regionHits = 0
  regionRejects = 0
  regionChecks = 0
  private static readonly GRID = 64
  private cellPhiLo: Float64Array | null = null
  private cellPhiHi: Float64Array | null = null
  private cellLb: Float64Array | null = null
  private cellLbHi: Float64Array | null = null
  private cellLbFor: object | null = null
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
    this.activeList = new Int32Array(n)
    this.activePos = new Int32Array(n)
    for (let i = 0; i < n; i++) { this.activeList[i] = i; this.activePos[i] = i }
    this.activeCount = n
    if (this.useRegions()) this.initRegions()
  }

  private useRegions(): boolean { return this.sleep.enabled && this.sleep.regions !== false }

  private initRegions(): void {
    const G = TransientSandSim.GRID; const a = this.plate.a; const nd = this.plate.ndof
    const w = (2 * a) / G
    this.cellPhiLo = new Float64Array(G * G * nd)
    this.cellPhiHi = new Float64Array(G * G * nd)
    this.cellLb = new Float64Array(G * G)
    this.cellLbHi = new Float64Array(G * G)
    const lo = new Float64Array(nd); const hi = new Float64Array(nd)
    for (let gy = 0; gy < G; gy++) {
      for (let gx = 0; gx < G; gx++) {
        this.plate.cellPhiRanges(-a + gx * w, -a + (gx + 1) * w, -a + gy * w, -a + (gy + 1) * w, lo, hi)
        const base = (gy * G + gx) * nd
        this.cellPhiLo.set(lo, base)
        this.cellPhiHi.set(hi, base)
      }
    }
  }

  /** Cell-wide wtt lower bounds for the current lookahead cache (once). */
  private ensureCellLb(): void {
    const lc = this.lookCache!
    if (this.cellLbFor === lc) return
    this.cellLbFor = lc
    const G = TransientSandSim.GRID; const nd = this.plate.ndof
    const loA = this.cellPhiLo!; const hiA = this.cellPhiHi!; const out = this.cellLb!
    const outHi = this.cellLbHi!
    for (let c = 0; c < G * G; c++) {
      let lb = 0; let ub = 0
      const base = c * nd
      for (let j = 0; j < nd; j++) {
        const plo = loA[base + j]; const phi = hiA[base + j]
        const a0 = lc.aMin[j]; const a1 = lc.aMax[j]
        lb += Math.min(plo * a0, plo * a1, phi * a0, phi * a1)
        // Upper bound on the exact per-particle term min(phi*a0, phi*a1)
        // for phi in [plo, phi]: the term is piecewise linear in phi with
        // a kink at 0, so its max is at an endpoint or at 0.
        const g = (v: number) => Math.min(v * a0, v * a1)
        ub += Math.max(g(plo), g(phi), plo <= 0 && 0 <= phi ? 0 : -Infinity)
      }
      out[c] = lb
      outHi[c] = ub
    }
  }

  private cellIndex(x: number, y: number): number {
    const G = TransientSandSim.GRID; const a = this.plate.a
    const gx = Math.min(G - 1, Math.max(0, Math.floor(((x + a) / (2 * a)) * G)))
    const gy = Math.min(G - 1, Math.max(0, Math.floor(((y + a) / (2 * a)) * G)))
    return gy * G + gx
  }

  private useScheduler(): boolean { return this.sleep.enabled && this.sleep.scheduler !== false }

  private setSleeping(i: number, v: 0 | 1): void {
    if (this.sleeping[i] !== v) { this.sleeping[i] = v; this.sleepingNow += v ? 1 : -1 }
  }

  private addActive(i: number): void {
    this.activeList[this.activeCount] = i
    this.activePos[i] = this.activeCount++
  }

  private removeActiveAt(idx: number): void {
    const removed = this.activeList[idx]
    const last = this.activeList[--this.activeCount]
    this.activeList[idx] = last
    this.activePos[last] = idx
    this.activePos[removed] = -1
  }

  private heapPush(t: number, i: number): void {
    if (this.heapSize === this.heapT.length) {
      const nt = new Float64Array(this.heapT.length * 2); nt.set(this.heapT); this.heapT = nt
      const ni = new Int32Array(this.heapI.length * 2); ni.set(this.heapI); this.heapI = ni
    }
    let c = this.heapSize++
    this.heapT[c] = t; this.heapI[c] = i
    while (c > 0) {
      const par = (c - 1) >> 1
      if (this.heapT[par] <= this.heapT[c]) break
      const tt = this.heapT[par]; this.heapT[par] = this.heapT[c]; this.heapT[c] = tt
      const ti = this.heapI[par]; this.heapI[par] = this.heapI[c]; this.heapI[c] = ti
      c = par
    }
  }

  private heapPop(): [number, number] {
    const t = this.heapT[0]; const i = this.heapI[0]
    const last = --this.heapSize
    if (last > 0) {
      this.heapT[0] = this.heapT[last]; this.heapI[0] = this.heapI[last]
      let c = 0
      for (;;) {
        const l = 2 * c + 1; const r = l + 1
        let m = c
        if (l < last && this.heapT[l] < this.heapT[m]) m = l
        if (r < last && this.heapT[r] < this.heapT[m]) m = r
        if (m === c) break
        const tt = this.heapT[m]; this.heapT[m] = this.heapT[c]; this.heapT[c] = tt
        const ti = this.heapI[m]; this.heapI[m] = this.heapI[c]; this.heapI[c] = ti
        c = m
      }
    }
    return [t, i]
  }

  /** Invalidate all sleep schedules (input change). State itself is kept;
   *  flying sleepers are materialized at the current time first. */
  wakeAll(): void {
    const st = this.state
    for (let i = 0; i < st.id.length; i++) {
      if (!this.sleeping[i]) continue
      this.materialize(i, this.t)
      this.setSleeping(i, 0)
    }
    this.lookCache = null
    // Scheduler: every schedule is void; all living particles are active.
    this.heapSize = 0
    this.activeCount = 0
    for (let i = 0; i < st.id.length; i++) {
      this.activePos[i] = -1
      if (st.alive[i]) this.addActive(i)
    }
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
    if (this.cellLb) {
      // Region shortcut: the cell bound holds at every point of the cell,
      // so a safe cell bound proves this particle safe without phiAll.
      this.regionChecks++
      this.ensureCellLb()
      const cell = this.cellIndex(st.x[i], st.y[i])
      // Certain-fail shortcut: even the best exact bound in this cell
      // cannot reach the margin, so the exact evaluation would fail too.
      if (G + this.cellLbHi![cell] < SLEEP_ACCEL_MARGIN) { this.regionRejects++; return }
      const clb = this.cellLb[cell]
      if (G + clb >= SLEEP_ACCEL_MARGIN && lc.end > tNew + 8 * this.dtNow()) {
        this.regionHits++
        this.setSleeping(i, 1)
        this.wakeT[i] = lc.end
        this.sleepT0[i] = tNew
        return
      }
    }
    pl.phiAll(st.x[i], st.y[i], this.phiBuf)
    // Exact lower bound of Hermite wtt over the lookahead window.
    let lb = 0
    for (let j = 0; j < pl.ndof; j++) {
      const phi = this.phiBuf[j]
      lb += Math.min(phi * lc.aMin[j], phi * lc.aMax[j])
    }
    if (G + lb < SLEEP_ACCEL_MARGIN) return
    if (!(lc.end > tNew + 8 * this.dtNow())) return
    this.setSleeping(i, 1)
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
    this.setSleeping(i, 1)
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
    if (this.useScheduler()) {
      this.stepScheduled(t, tNew, dt, vStick, customDt)
    } else {
      // M6-2 legacy path: full scan every step (kept for comparison).
      for (let i = 0; i < n; i++) {
        if (!st.alive[i]) continue
        if (this.sleeping[i]) {
          if (this.wakeT[i] > tNew) { this.skippedUpdates++; continue }
          // Wake at the STEP START (<= wakeT, still inside the guaranteed
          // interval) so this step is processed exactly as the baseline.
          this.setSleeping(i, 0); this.wakeCount++
          if (st.sitting[i]) {
            pl.surfAt(st.x[i], st.y[i], t, S5); st.z[i] = S5[0]
          } else {
            this.materialize(i, t)
          }
        }
        this.advanceParticle(i, t, tNew, dt, vStick, customDt)
      }
      for (let i = 0; i < n; i++) {
        if (st.alive[i] && Math.hypot(st.x[i], st.y[i]) > pl.a) { st.alive[i] = 0; this.setSleeping(i, 0) }
      }
    }
    this.t = tNew
  }

  // M6-3A: only active particles are touched; due sleepers are woken from
  // the heap at the step start with exactly the legacy wake semantics
  // (wake when wakeT <= tNew, materialized at the step-start time).
  private stepScheduled(t: number, tNew: number, dt: number, vStick: number, customDt: boolean): void {
    const { plate: pl, state: st } = this
    while (this.heapSize > 0 && this.heapT[0] <= tNew) {
      const [wt, i] = this.heapPop()
      if (!this.sleeping[i] || this.wakeT[i] !== wt || !st.alive[i]) continue // stale entry
      this.setSleeping(i, 0); this.wakeCount++
      if (st.sitting[i]) {
        pl.surfAt(st.x[i], st.y[i], t, S5); st.z[i] = S5[0]
      } else {
        this.materialize(i, t)
      }
      this.addActive(i)
    }
    this.skippedUpdates += this.sleepingNow
    for (let idx = 0; idx < this.activeCount;) {
      const i = this.activeList[idx]
      if (!st.alive[i]) { this.removeActiveAt(idx); continue }
      this.advanceParticle(i, t, tNew, dt, vStick, customDt)
      if (st.alive[i] && Math.hypot(st.x[i], st.y[i]) > pl.a) { st.alive[i] = 0; this.setSleeping(i, 0) }
      if (!st.alive[i]) { this.removeActiveAt(idx); continue }
      if (this.sleeping[i]) { this.heapPush(this.wakeT[i], i); this.removeActiveAt(idx); continue }
      idx++
    }
  }

  /** Advance one awake particle from t to tNew (events, fallback, sleep
   *  evaluation). Shared verbatim by the legacy scan and the scheduler. */
  private advanceParticle(i: number, t: number, tNew: number, dt: number, vStick: number, customDt: boolean): void {
    const { plate: pl, state: st, p } = this
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

  // M6-3C differential frame for display transfer. Particle identity is
  // the state-array index (stable for the life of the sim). The first
  // call (and any call after resetFrameDelta) sends every living slot;
  // later calls send only slots whose DISPLAY xy changed since the last
  // frame, plus slots that died. Display xy here is exactly snapshot()'s
  // (analytic for sleeping flyers, stored xy otherwise) at the same
  // physics time — display optimization never touches physical state.
  // Sitting particles (sleeping or not) keep constant xy, so formed
  // patterns mostly stop being resent.
  private prevSentXy: Float32Array | null = null
  private prevSentAlive: Uint8Array | null = null
  private frameFirst = true
  private deltaIdx: Int32Array | null = null
  private deltaDead: Int32Array | null = null

  resetFrameDelta(): void {
    this.prevSentXy = null; this.prevSentAlive = null; this.frameFirst = true
  }

  frameDelta(stepsPerSec = 0) {
    const st = this.state
    const n = st.id.length
    if (!this.prevSentXy) {
      this.prevSentXy = new Float32Array(n * 2).fill(NaN)
      this.prevSentAlive = new Uint8Array(n)
      this.deltaIdx = new Int32Array(n)
      this.deltaDead = new Int32Array(n)
    }
    const prev = this.prevSentXy; const prevAlive = this.prevSentAlive!
    const idxBuf = this.deltaIdx!; const deadBuf = this.deltaDead!
    const full = this.frameFirst
    this.frameFirst = false
    let nIdx = 0; let nDead = 0; let aliveCount = 0
    for (let i = 0; i < n; i++) {
      if (st.alive[i]) {
        aliveCount++
        let x: number; let y: number
        if (this.sleeping[i] && !st.sitting[i]) {
          const d = this.t - this.sleepT0[i]
          x = st.x[i] + st.vx[i] * d; y = st.y[i] + st.vy[i] * d
        } else { x = st.x[i]; y = st.y[i] }
        const fx = Math.fround(x); const fy = Math.fround(y)
        if (!prevAlive[i] || prev[2 * i] !== fx || prev[2 * i + 1] !== fy) {
          idxBuf[nIdx++] = i
          prev[2 * i] = fx; prev[2 * i + 1] = fy
        }
        prevAlive[i] = 1
      } else if (prevAlive[i]) {
        deadBuf[nDead++] = i
        prevAlive[i] = 0
      }
    }
    const xy = new Float32Array(nIdx * 2)
    for (let k = 0; k < nIdx; k++) {
      xy[2 * k] = prev[2 * idxBuf[k]]; xy[2 * k + 1] = prev[2 * idxBuf[k] + 1]
    }
    return {
      t: this.t, full, idx: idxBuf.slice(0, nIdx), xy, dead: deadBuf.slice(0, nDead),
      aliveCount, initialCount: n, stats: this.stats(stepsPerSec),
    }
  }
}

export type TransientFrame = ReturnType<TransientSandSim['frameDelta']>

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
