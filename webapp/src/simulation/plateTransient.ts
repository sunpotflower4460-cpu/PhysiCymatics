// plateTransient.ts — multi-mode transient plate (TS port of
// reference/plate_transient.py). See that file for PUT-IN / floors.
// qdd + 2*zeta*omega*qd + omega^2 q = Phi(drive)*F(t)/M, RK4 substeps,
// cubic-Hermite state evaluation inside a sand step.
import basisData from '../data/plate_modes_transient.json'
import { ScoreSignal, type PerformanceScore, type DriveEvent } from './sequence'
import { G } from './surface'

interface ModeRow { m: number; root_idx: number; f_hz: number; M_kg: number; R: number[]; dR_drho: number[] }
interface BasisData { a_m: number; ngrid: number; zeta: number; modes: ModeRow[] }
const BD = basisData as unknown as BasisData
// Reused trig-recurrence scratch (surfAt/phiAll are synchronous).
const CM = new Float64Array(6)
const SM = new Float64Array(6)

export interface Dof { mode: number; kind: 'cos' | 'sin' }

/** Tri-state result of the interval-assisted detach search (M6-12R):
 *  'found' = first crossing located; 'none' = absence proven by
 *  interval bounds; 'unknown' = search budget exhausted (absence NOT
 *  proven). 'unknown' must never be treated as 'none'. */
export type DetachSearchResult =
  | { kind: 'found'; t: number }
  | { kind: 'none' }
  | { kind: 'unknown' }

/** M6-12 step-boundary plate checkpoint (structured-cloneable). */
export interface PlateCheckpoint {
  t: number; epoch: number
  driveXY: [number, number]
  accelMode: 'hermite' | 'ode'
  scoreKey: string
  q: Float64Array; qd: Float64Array
  seg: {
    t0: number; t1: number
    q0: Float64Array; qd0: Float64Array
    q1: Float64Array; qd1: Float64Array
  } | null
  counters: [number, number, number]
}

export class TransientPlate {
  readonly a = BD.a_m
  readonly zeta = BD.zeta
  readonly modes = BD.modes
  readonly dofs: Dof[] = []
  readonly omega: Float64Array
  readonly mass: Float64Array
  phiDrive: Float64Array
  q: Float64Array; qd: Float64Array
  hBound = 0
  t = 0
  driveXY: [number, number] = [0, 0]
  private seg: { t0: number; t1: number; q0: Float64Array; qd0: Float64Array; q1: Float64Array; qd1: Float64Array } | null = null
  private R: Float64Array[]; private dR: Float64Array[]
  private ng1 = BD.ngrid - 1
  // M6-4A hot-path tables: flattened R and dR/a (same float values as the
  // per-call divisions), typed DOF descriptors, and a direct-mapped cache
  // of per-tau Hermite modal coefficients (exact tau key, cleared per
  // segment). Arithmetic and summation order are unchanged.
  private Rflat!: Float64Array
  private dRflat!: Float64Array
  private dofMode!: Int32Array
  private dofM!: Int32Array
  private dofCos!: Uint8Array
  private static readonly TAU_SLOTS = 64
  private tauSlotTau = new Float64Array(TransientPlate.TAU_SLOTS).fill(NaN)
  private tauQ!: Float64Array
  private tauQd!: Float64Array
  private tauQdd!: Float64Array
  // M6-6: q-only coefficient cache for heightAt (same A..D arithmetic
  // as modalCoeffs, so cached q values are bitwise identical).
  private tauSlotTauQ = new Float64Array(TransientPlate.TAU_SLOTS).fill(NaN)
  private tauQQ!: Float64Array
  // M6-5: exact per-DOF Hermite q ranges over the current segment.
  private segQMin!: Float64Array
  private segQMax!: Float64Array
  // M6-6: lazily computed per-slice q ranges for the current segment,
  // cached against the segment object and slice count.
  private sliceCache: { seg: object; k: number; min: Float64Array; max: Float64Array } | null = null
  signal: ScoreSignal
  epoch = 0 // score time = physics t - epoch (live score replacement)
  // Evaluation counters (M6-7 profiling; read-only diagnostics).
  surfCalls = 0
  heightCalls = 0
  heightAccelCalls = 0
  // M6-7C: q+qdd coefficient cache for heightAccelAt/accelAt (same
  // Hermite expressions as modalCoeffs, so values are bitwise identical).
  private tauSlotTauQH = new Float64Array(TransientPlate.TAU_SLOTS).fill(NaN)
  private tauQH!: Float64Array
  private tauQddH!: Float64Array
  // M6-8A: DOFs live in the current segment (any Hermite endpoint
  // nonzero). Dead DOFs have q(tau) identically zero on the segment,
  // so their height/accel terms are exact zeros and can be skipped
  // bitwise (finite state assumed, as everywhere in the physics).
  // With center drive only the m=0 modes are ever live.
  liveDofs!: Int32Array
  liveCount = 0
  /** A/B measurement switch (bench only): false evaluates all DOFs. */
  liveMaskEnabled = true
  /** M6-11 research switch: 'ode' evaluates modal acceleration from the
   *  equation of motion applied to the Hermite q/qd, instead of the
   *  Hermite second derivative. Only wtt changes (h/uz stay bitwise the
   *  Hermite values); results are NOT bitwise vs the default. Default
   *  'hermite' is the production model. See docs/m6_11_results.md.
   *  Switching invalidates the tau coefficient caches so evaluations
   *  never mix definitions within a segment. Sim-level caches (sleep
   *  lookahead) must additionally be invalidated by the caller when
   *  switching mid-run (wakeAll + refreshForceBounds), exactly as for
   *  drive-condition changes. */
  private _accelMode: 'hermite' | 'ode' = 'hermite'
  get accelMode(): 'hermite' | 'ode' { return this._accelMode }
  set accelMode(m: 'hermite' | 'ode') {
    if (m === this._accelMode) return
    this._accelMode = m
    this.invalidateAccelCaches()
  }

  private computeLiveMask(q0: Float64Array, qd0: Float64Array, q1: Float64Array, qd1: Float64Array, t0: number, t1: number): void {
    const n = this.ndof
    let lc = 0
    if (this.liveMaskEnabled) {
      let forcedLive = false
      // Effective acceleration form for THIS segment: tiny segments
      // (h < 1e-7) use the ODE form in every accel-bearing cache
      // (modalCoeffs / tauQH fill), so the mask must keep directly
      // forced DOFs live under the same criterion (M6-12R review).
      if (this.accelMode === 'ode' || t1 - t0 < 1e-7) {
        const ts0 = t0 - this.epoch; const ts1 = t1 - this.epoch
        for (const e of this.signal.score.events) {
          if (e.forcePeak_N !== 0 && e.start_s <= ts1 && e.end_s >= ts0) { forcedLive = true; break }
        }
      }
      for (let j = 0; j < n; j++) {
        if (q0[j] !== 0 || qd0[j] !== 0 || q1[j] !== 0 || qd1[j] !== 0
          || (forcedLive && this.phiDrive[j] !== 0)) this.liveDofs[lc++] = j
      }
    } else {
      for (let j = 0; j < n; j++) this.liveDofs[lc++] = j
    }
    this.liveCount = lc
  }

  constructor(score: PerformanceScore, driveXY: [number, number] = [0, 0]) {
    this.signal = new ScoreSignal(score)
    this.modes.forEach((md, i) => {
      this.dofs.push({ mode: i, kind: 'cos' })
      if (md.m > 0) this.dofs.push({ mode: i, kind: 'sin' })
    })
    const n = this.dofs.length
    this.omega = new Float64Array(n); this.mass = new Float64Array(n)
    this.dofs.forEach((d, j) => {
      this.omega[j] = 2 * Math.PI * this.modes[d.mode].f_hz
      this.mass[j] = this.modes[d.mode].M_kg
    })
    this.q = new Float64Array(n); this.qd = new Float64Array(n)
    this.phiDrive = new Float64Array(n)
    this.R = this.modes.map((m) => Float64Array.from(m.R))
    this.dR = this.modes.map((m) => Float64Array.from(m.dR_drho))
    const stride = BD.ngrid
    this.Rflat = new Float64Array(this.modes.length * stride)
    this.dRflat = new Float64Array(this.modes.length * stride)
    this.modes.forEach((md, mi) => {
      for (let k = 0; k < stride; k++) {
        this.Rflat[mi * stride + k] = md.R[k]
        this.dRflat[mi * stride + k] = md.dR_drho[k] / this.a
      }
    })
    this.dofMode = new Int32Array(n); this.dofM = new Int32Array(n); this.dofCos = new Uint8Array(n)
    this.dofs.forEach((d, j) => {
      this.dofMode[j] = d.mode
      this.dofM[j] = this.modes[d.mode].m
      this.dofCos[j] = d.kind === 'cos' ? 1 : 0
    })
    this.tauQ = new Float64Array(TransientPlate.TAU_SLOTS * n)
    this.tauQd = new Float64Array(TransientPlate.TAU_SLOTS * n)
    this.tauQdd = new Float64Array(TransientPlate.TAU_SLOTS * n)
    this.tauQQ = new Float64Array(TransientPlate.TAU_SLOTS * n)
    this.tauQH = new Float64Array(TransientPlate.TAU_SLOTS * n)
    this.tauQddH = new Float64Array(TransientPlate.TAU_SLOTS * n)
    this.liveDofs = new Int32Array(n)
    for (let j = 0; j < n; j++) this.liveDofs[j] = j
    this.liveCount = n
    this.segQMin = new Float64Array(n)
    this.segQMax = new Float64Array(n)
    this.setDrive(driveXY[0], driveXY[1])
  }

  get ndof(): number { return this.dofs.length }

  /** Phi for one DOF at (x,y) using precomputed radial index + angle terms. */
  private phi(modeIdx: number, kind: string, idx: number, cm: number, sm: number): number {
    const R = this.R[modeIdx][idx]
    if (this.modes[modeIdx].m === 0) return R
    return kind === 'cos' ? R * cm : R * sm
  }

  /** Invalidate tau caches and rebuild the live-DOF mask for the
   *  current segment. Required whenever anything the ODE
   *  acceleration or the mask's forced-DOF criterion depends on
   *  changes mid-segment: drive position (phiDrive), score/epoch
   *  (force events), or accel mode (M6-14R reviews). */
  private invalidateAccelCaches(): void {
    this.tauSlotTau.fill(NaN)
    this.tauSlotTauQH.fill(NaN)
    const seg = this.seg
    if (seg) this.computeLiveMask(seg.q0, seg.qd0, seg.q1, seg.qd1, seg.t0, seg.t1)
  }

  setDrive(x: number, y: number): void {
    const changed = this.driveXY[0] !== x || this.driveXY[1] !== y
    this.driveXY = [x, y]
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const th = Math.atan2(y, x)
    this.dofs.forEach((d, j) => {
      const m = this.modes[d.mode].m
      this.phiDrive[j] = this.phi(d.mode, d.kind, idx, Math.cos(m * th), Math.sin(m * th))
    })
    if (changed) this.invalidateAccelCaches()
  }

  /** Replace the score at physics time atTime; q/qd are preserved. */
  setScore(score: PerformanceScore, atTime: number): void {
    this.signal = new ScoreSignal(score)
    this.epoch = atTime
    this.invalidateAccelCaches()
  }

  /** Switch to a single tone at freqHz, phase-continuous with the current
   *  signal at atTime (same theta, same force peak as current event).
   *  q/qd and all particle state are preserved; no snapping. */
  setFrequency(freqHz: number, atTime: number): void {
    const ts = atTime - this.epoch
    const theta = this.signal.thetaAt(ts)
    const force = this.signal.forcePeakAt(ts) || this.signal.score.events[0]?.forcePeak_N || 0
    const dur = 3600
    this.setScore({
      schema: 'physicymatics.score.v1', duration_s: dur, initialPhase_rad: theta,
      audioGain: this.signal.score.audioGain,
      events: [{ start_s: 0, end_s: dur, frequency_hz: freqHz, forcePeak_N: force, attack_s: 0, release_s: 0 }],
    }, atTime)
  }

  scoreTime(t: number): number { return t - this.epoch }

  /** Per-DOF forcing bound |Phi_j(drive)| * maxForcePeak / M_j [m/s^2 per unit q]. */
  forceBounds(): Float64Array {
    const fMax = Math.max(0, ...this.signal.score.events.map((e) => e.forcePeak_N))
    const out = new Float64Array(this.ndof)
    for (let j = 0; j < this.ndof; j++) out[j] = (Math.abs(this.phiDrive[j]) * fMax) / this.mass[j]
    return out
  }

  /** Per-DOF modal energy amplitude E_j = sqrt(qd^2 + (omega q)^2) at segment end. */
  energyPerDof(out: Float64Array): void {
    const q = this.seg ? this.seg.q1 : this.q
    const qd = this.seg ? this.seg.qd1 : this.qd
    for (let j = 0; j < this.ndof; j++) out[j] = Math.hypot(qd[j], this.omega[j] * q[j])
  }

  /** Exact per-DOF Phi ranges over an axis-aligned rectangle, for the
   *  AS-IMPLEMENTED Phi: Phi_j(x,y) = R_j[bin(r)] * trig_m(theta) with the
   *  same radial binning as phiAll/surfAt. r-range gives a bin range
   *  (R min/max over those bins); theta-range gives the exact trig range
   *  (endpoints + interior extrema). The interval product therefore
   *  contains Phi_j at EVERY point of the rectangle — no representative
   *  point, no sampling. M6-3B region-shared sleep bounds use this. */
  cellPhiRanges(x0: number, x1: number, y0: number, y1: number,
                outLo: Float64Array, outHi: Float64Array): void {
    const containsOrigin = x0 <= 0 && 0 <= x1 && y0 <= 0 && 0 <= y1
    const dx = x0 > 0 ? x0 : x1 < 0 ? -x1 : 0
    const dy = y0 > 0 ? y0 : y1 < 0 ? -y1 : 0
    const rLo = Math.hypot(dx, dy)
    const rHi = Math.max(
      Math.hypot(x0, y0), Math.hypot(x0, y1), Math.hypot(x1, y0), Math.hypot(x1, y1))
    const bLo = Math.min(Math.floor(rLo * (this.ng1 / this.a)), this.ng1)
    const bHi = Math.min(Math.floor(rHi * (this.ng1 / this.a)), this.ng1)
    // Covering angular arc of the rectangle (theta is monotonic along
    // each edge, so corner angles suffice when origin is outside).
    let thLo = 0; let thHi = 0
    if (!containsOrigin) {
      const angs = [Math.atan2(y0, x0), Math.atan2(y0, x1), Math.atan2(y1, x0), Math.atan2(y1, x1)]
        .sort((u, v) => u - v)
      let bestGap = -1; let start = angs[0]
      for (let k = 0; k < 4; k++) {
        const nxt = k < 3 ? angs[k + 1] : angs[0] + 2 * Math.PI
        const gap = nxt - angs[k]
        if (gap > bestGap) { bestGap = gap; start = nxt }
      }
      thLo = start; thHi = start + (2 * Math.PI - bestGap)
    }
    const trigRange = (m: number, kind: 'cos' | 'sin'): [number, number] => {
      if (containsOrigin || m === 0) return kind === 'cos' || m === 0 ? (m === 0 ? [1, 1] : [-1, 1]) : [-1, 1]
      const f = kind === 'cos' ? Math.cos : Math.sin
      let lo = Math.min(f(m * thLo), f(m * thHi)); let hi = Math.max(f(m * thLo), f(m * thHi))
      // Interior extrema: cos(m t) at m t = k*pi; sin(m t) at m t = pi/2 + k*pi.
      const off = kind === 'cos' ? 0 : Math.PI / 2
      const kLo = Math.ceil((m * thLo - off) / Math.PI)
      const kHi = Math.floor((m * thHi - off) / Math.PI)
      for (let k = kLo; k <= kHi; k++) {
        const v = f(off + k * Math.PI)
        if (v < lo) lo = v
        if (v > hi) hi = v
      }
      return [lo, hi]
    }
    for (let j = 0; j < this.ndof; j++) {
      const d = this.dofs[j]; const md = this.modes[d.mode]
      let rMin = Infinity; let rMax = -Infinity
      const R = this.R[d.mode]
      for (let b = bLo; b <= bHi; b++) {
        if (R[b] < rMin) rMin = R[b]
        if (R[b] > rMax) rMax = R[b]
      }
      if (md.m === 0) { outLo[j] = rMin; outHi[j] = rMax; continue }
      const [tLo, tHi] = trigRange(md.m, d.kind)
      const p1 = rMin * tLo; const p2 = rMin * tHi; const p3 = rMax * tLo; const p4 = rMax * tHi
      outLo[j] = Math.min(p1, p2, p3, p4)
      outHi[j] = Math.max(p1, p2, p3, p4)
    }
  }

  /** All DOF shape values Phi_j(x, y) (no gradients). */
  phiAll(x: number, y: number, out: Float64Array): void {
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const rs = r < 1e-12 ? 1e-12 : r
    const c = x / rs; const sn = y / rs
    CM[0] = 1; CM[1] = c; SM[0] = 0; SM[1] = sn
    for (let m = 2; m <= 4; m++) {
      CM[m] = CM[m - 1] * c - SM[m - 1] * sn
      SM[m] = SM[m - 1] * c + CM[m - 1] * sn
    }
    const stride = BD.ngrid
    for (let j = 0; j < this.ndof; j++) {
      const m = this.dofM[j]
      const R = this.Rflat[this.dofMode[j] * stride + idx]
      out[j] = m === 0 ? R : this.dofCos[j] ? R * CM[m] : R * SM[m]
    }
  }

  private accel(q: Float64Array, qd: Float64Array, t: number, out: Float64Array): void {
    const ts = t - this.epoch
    const F = this.signal.forcePeakAt(ts) * this.signal.u(ts)
    for (let j = 0; j < this.ndof; j++) {
      out[j] = (this.phiDrive[j] * F) / this.mass[j]
        - 2 * this.zeta * this.omega[j] * qd[j] - this.omega[j] ** 2 * q[j]
    }
  }

  /** RK4-advance a modal state over [t0, t1] (same substeps as beginStep). */
  advanceState(qIn: Float64Array, qdIn: Float64Array, t0: number, t1: number): [Float64Array, Float64Array] {
    const n = this.ndof
    let q = Float64Array.from(qIn); let qd = Float64Array.from(qdIn)
    const hMax = 0.3 / Math.max(...this.omega)
    const nsub = Math.max(1, Math.ceil((t1 - t0) / hMax))
    const h = (t1 - t0) / nsub
    const k1v = new Float64Array(n); const k2v = new Float64Array(n)
    const k3v = new Float64Array(n); const k4v = new Float64Array(n)
    const qt = new Float64Array(n); const qdt = new Float64Array(n)
    let t = t0
    for (let s = 0; s < nsub; s++) {
      this.accel(q, qd, t, k1v)
      for (let j = 0; j < n; j++) { qt[j] = q[j] + 0.5 * h * qd[j]; qdt[j] = qd[j] + 0.5 * h * k1v[j] }
      this.accel(qt, qdt, t + 0.5 * h, k2v)
      const k2q = Float64Array.from(qdt)
      for (let j = 0; j < n; j++) { qt[j] = q[j] + 0.5 * h * k2q[j]; qdt[j] = qd[j] + 0.5 * h * k2v[j] }
      this.accel(qt, qdt, t + 0.5 * h, k3v)
      const k3q = Float64Array.from(qdt)
      for (let j = 0; j < n; j++) { qt[j] = q[j] + h * k3q[j]; qdt[j] = qd[j] + h * k3v[j] }
      this.accel(qt, qdt, t + h, k4v)
      for (let j = 0; j < n; j++) {
        const k1q = qd[j]; const k4q = qdt[j]
        q[j] += (h * (k1q + 2 * k2q[j] + 2 * k3q[j] + k4q)) / 6
        qd[j] += (h * (k1v[j] + 2 * k2v[j] + 2 * k3v[j] + k4v[j])) / 6
      }
      t += h
    }
    return [q, qd]
  }

  beginStep(t0: number, t1: number): void {
    const q0 = Float64Array.from(this.q); const qd0 = Float64Array.from(this.qd)
    const [q, qd] = this.advanceState(q0, qd0, t0, t1)
    this.q = q; this.qd = qd; this.t = t1
    this.seg = { t0, t1, q0, qd0, q1: Float64Array.from(q), qd1: Float64Array.from(qd) }
    this.finalizeSegment()
  }

  /** Derived per-segment state (height bound, tau-cache invalidation,
   *  live mask, exact q ranges). Shared by beginStep and checkpoint
   *  restore so both build identical derived state from the endpoints. */
  private finalizeSegment(): void {
    const n = this.ndof
    const { t0, t1, q0, qd0, q1: q, qd1: qd } = this.seg!
    // Conservative |h| bound over the segment (see Python plate_transient).
    // Hermite tangent basis max is 4/27 (at s=1/3), not 0.105.
    let hb = 0
    for (let j = 0; j < n; j++) {
      hb += Math.max(Math.abs(q0[j]), Math.abs(q[j])) + (4 / 27) * (t1 - t0) * (Math.abs(qd0[j]) + Math.abs(qd[j]))
    }
    this.hBound = hb
    this.tauSlotTau.fill(NaN)
    this.tauSlotTauQ.fill(NaN)
    this.tauSlotTauQH.fill(NaN)
    // In 'ode' accel mode a DOF with all-zero Hermite endpoints can
    // still have nonzero acceleration inside the segment when the
    // drive force is nonzero there (the RK4 substeps may not sample
    // a short event). computeLiveMask keeps such DOFs live so every
    // evaluation path (surfAt / heightAccelAt / accelAt) sums the same
    // terms; their q/qd contributions are exact zeros, so h/uz are
    // bitwise unaffected. Default 'hermite' masking is unchanged.
    this.computeLiveMask(q0, qd0, q, qd, t0, t1)
    for (let j = 0; j < n; j++) {
      const r = hermiteQRange(q0[j], qd0[j], q[j], qd[j], t1 - t0)
      // Widen by a rounding hair for bound use (M6-5R): modalCoeffs
      // evaluates the same cubic in a different operation order, so the
      // implemented q(tau) can differ from the exact range by ulps.
      // hermiteQRange itself stays exact for its other callers.
      const wden = 1e-12 * (1 + Math.abs(r[0]) + Math.abs(r[1]))
      this.segQMin[j] = r[0] - wden; this.segQMax[j] = r[1] + wden
    }
  }

  /** M6-12: full plate state at a step boundary (between steps, i.e.
   *  after a completed beginStep-driven step). Restoring into a plate
   *  built with the same score/drive reproduces the exact state:
   *  q/qd, epoch, drive, accel mode, and the last segment with all
   *  derived state rebuilt by finalizeSegment (no re-integration, so
   *  a score/frequency switch applied at the boundary is preserved). */
  checkpointPlate(): PlateCheckpoint {
    const seg = this.seg
    return {
      t: this.t, epoch: this.epoch,
      driveXY: [this.driveXY[0], this.driveXY[1]],
      accelMode: this._accelMode,
      scoreKey: JSON.stringify(this.signal.score),
      q: Float64Array.from(this.q), qd: Float64Array.from(this.qd),
      seg: seg ? {
        t0: seg.t0, t1: seg.t1,
        q0: Float64Array.from(seg.q0), qd0: Float64Array.from(seg.qd0),
        q1: Float64Array.from(seg.q1), qd1: Float64Array.from(seg.qd1),
      } : null,
      counters: [this.surfCalls, this.heightCalls, this.heightAccelCalls],
    }
  }

  restorePlate(cp: PlateCheckpoint): void {
    if (cp.q.length !== this.ndof) throw new Error('checkpoint ndof mismatch')
    if (cp.scoreKey !== JSON.stringify(this.signal.score)) {
      throw new Error('checkpoint score mismatch: restore target must use the same score')
    }
    this.epoch = cp.epoch
    this.t = cp.t
    this.setDrive(cp.driveXY[0], cp.driveXY[1])
    this.q = Float64Array.from(cp.q); this.qd = Float64Array.from(cp.qd)
    this.sliceCache = null
    if (cp.seg) {
      this.seg = {
        t0: cp.seg.t0, t1: cp.seg.t1,
        q0: Float64Array.from(cp.seg.q0), qd0: Float64Array.from(cp.seg.qd0),
        q1: Float64Array.from(cp.seg.q1), qd1: Float64Array.from(cp.seg.qd1),
      }
      this.finalizeSegment()
    } else {
      this.seg = null
      this.tauSlotTau.fill(NaN); this.tauSlotTauQ.fill(NaN); this.tauSlotTauQH.fill(NaN)
    }
    // Set the mode last via the public setter: it invalidates tau
    // caches and rebuilds the live mask for the restored segment
    // (finalizeSegment already built it under the previous mode;
    // the setter is a no-op when the mode already matches).
    this.accelMode = cp.accelMode
    this.surfCalls = cp.counters[0]; this.heightCalls = cp.counters[1]; this.heightAccelCalls = cp.counters[2]
  }

  /** Exact per-DOF q range over the current segment (M6-5 local bound). */
  segQRanges(): { min: Float64Array; max: Float64Array } {
    return { min: this.segQMin, max: this.segQMax }
  }

  /** Per-DOF q ranges over each of k equal time slices of the current
   *  segment (M6-6): min/max are k*ndof arrays, slice s at [s*ndof, ...].
   *  The Hermite cubic restricted to a slice is the Hermite cubic through
   *  its endpoint values/derivatives (evaluated with the same modalCoeffs
   *  arithmetic surfAt uses), so hermiteQRange on those endpoints is the
   *  exact range over the slice; widened by the same rounding hair as
   *  segQRanges. Slice s covers [t0 + hseg*s/k, t0 + hseg*(s+1)/k]. */
  segQSliceRanges(k: number): { min: Float64Array; max: Float64Array } {
    const seg = this.seg!
    const c = this.sliceCache
    if (c && c.seg === seg && c.k === k) return { min: c.min, max: c.max }
    const n = this.ndof
    const min = new Float64Array(k * n); const max = new Float64Array(k * n)
    const hseg = seg.t1 - seg.t0
    // Copy endpoint modal states first (modalCoeffs slots are reused).
    const qs = new Float64Array((k + 1) * n); const qds = new Float64Array((k + 1) * n)
    for (let s = 0; s <= k; s++) {
      const base = this.modalCoeffs(seg.t0 + (hseg * s) / k)
      qs.set(this.tauQ.subarray(base, base + n), s * n)
      qds.set(this.tauQd.subarray(base, base + n), s * n)
    }
    for (let s = 0; s < k; s++) {
      const ta = seg.t0 + (hseg * s) / k; const tb = seg.t0 + (hseg * (s + 1)) / k
      for (let j = 0; j < n; j++) {
        const r = hermiteQRange(qs[s * n + j], qds[s * n + j], qs[(s + 1) * n + j], qds[(s + 1) * n + j], tb - ta)
        const wden = 1e-12 * (1 + Math.abs(r[0]) + Math.abs(r[1]))
        min[s * n + j] = r[0] - wden; max[s * n + j] = r[1] + wden
      }
    }
    this.sliceCache = { seg, k, min, max }
    return { min, max }
  }

  /** Per-tau Hermite modal coefficients (q, qd, qdd per DOF) for the
   *  current segment, direct-mapped cache keyed by exact tau. Returns
   *  the slot base offset into tauQ/tauQd/tauQdd. The arithmetic is the
   *  same expressions surfAtReference evaluates per DOF. */
  private modalCoeffs(tau: number): number {
    const seg = this.seg!
    const slots = TransientPlate.TAU_SLOTS
    const hbits = (tau * 2654435761) | 0
    const slot = (hbits ^ (hbits >>> 3)) & (slots - 1)
    const n = this.ndof
    const base = slot * n
    if (this.tauSlotTau[slot] === tau) return base
    this.tauSlotTau[slot] = tau
    const hseg = seg.t1 - seg.t0
    const s = (tau - seg.t0) / hseg
    const s2 = s * s; const s3 = s2 * s
    const A = 2 * s3 - 3 * s2 + 1; const B = (s3 - 2 * s2 + s) * hseg
    const C = -2 * s3 + 3 * s2; const D = (s3 - s2) * hseg
    const E = 6 * s2 - 6 * s; const F = (3 * s2 - 4 * s + 1) * hseg
    const Gc = -6 * s2 + 6 * s; const H = (3 * s2 - 2 * s) * hseg
    const P = 12 * s - 6; const Q = (6 * s - 4) * hseg
    const Rr = -12 * s + 6; const Ss = (6 * s - 2) * hseg
    const hh = hseg * hseg
    const { q0, qd0, q1, qd1 } = seg
    // Tiny segments (h < 1e-7 s) always use the ODE form: the Hermite
    // qdd divides by h^2 and degenerates into cancellation garbage as
    // h -> 0, while the ODE form is its stable limit (M6-11R).
    if (this.accelMode === 'ode' || hseg < 1e-7) {
      const ts = tau - this.epoch
      const Fd = this.signal.forcePeakAt(ts) * this.signal.u(ts)
      for (let j = 0; j < n; j++) {
        const qv = A * q0[j] + B * qd0[j] + C * q1[j] + D * qd1[j]
        const qdv = (E * q0[j] + F * qd0[j] + Gc * q1[j] + H * qd1[j]) / hseg
        this.tauQ[base + j] = qv
        this.tauQd[base + j] = qdv
        this.tauQdd[base + j] = (this.phiDrive[j] * Fd) / this.mass[j]
          - 2 * this.zeta * this.omega[j] * qdv - this.omega[j] ** 2 * qv
      }
      return base
    }
    for (let j = 0; j < n; j++) {
      this.tauQ[base + j] = A * q0[j] + B * qd0[j] + C * q1[j] + D * qd1[j]
      this.tauQd[base + j] = (E * q0[j] + F * qd0[j] + Gc * q1[j] + H * qd1[j]) / hseg
      this.tauQdd[base + j] = (P * q0[j] + Q * qd0[j] + Rr * q1[j] + Ss * qd1[j]) / hh
    }
    return base
  }

  /** (h, hx, hy, uz, wtt) at (x, y, tau) inside the current segment.
   *  M6-4A optimized path: per-tau modal coefficients are shared through
   *  modalCoeffs(), spatial weights use flattened tables and reused
   *  scratch. Arithmetic and summation order match surfAtReference. */
  surfAt(x: number, y: number, tau: number, out: Float64Array): void {
    this.surfCalls++
    const seg = this.seg
    if (!seg) { out.fill(0); return }
    const base = this.modalCoeffs(tau)
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const rs = r < 1e-12 ? 1e-12 : r
    const c = x / rs; const sn = y / rs
    CM[0] = 1; CM[1] = c; SM[0] = 0; SM[1] = sn
    for (let m = 2; m <= 4; m++) {
      CM[m] = CM[m - 1] * c - SM[m - 1] * sn
      SM[m] = SM[m - 1] * c + CM[m - 1] * sn
    }
    const stride = BD.ngrid
    let hh = 0; let hx = 0; let hy = 0; let uz = 0; let wtt = 0
    for (let j = 0; j < this.ndof; j++) {
      const mi = this.dofMode[j]; const m = this.dofM[j]
      const R = this.Rflat[mi * stride + idx]; const dR = this.dRflat[mi * stride + idx]
      let W: number; let Wx: number; let Wy: number
      if (m === 0) {
        W = R; Wx = dR * c; Wy = dR * sn
      } else {
        const cm = CM[m]; const sm = SM[m]
        let Wr: number; let Wt: number
        if (this.dofCos[j]) { W = R * cm; Wr = dR * cm; Wt = (-m * R * sm) / rs }
        else { W = R * sm; Wr = dR * sm; Wt = (m * R * cm) / rs }
        Wx = Wr * c - Wt * sn; Wy = Wr * sn + Wt * c
      }
      const q = this.tauQ[base + j]; const qd = this.tauQd[base + j]; const qdd = this.tauQdd[base + j]
      hh += q * W; hx += q * Wx; hy += q * Wy; uz += qd * W; wtt += qdd * W
    }
    out[0] = hh; out[1] = hx; out[2] = hy; out[3] = uz; out[4] = wtt
  }

  /** q-only modal coefficients (M6-6): identical A..D expressions to
   *  modalCoeffs, so the returned q values are bitwise the same; only
   *  the qd/qdd work is skipped. Separate slot cache. */
  private modalCoeffsQ(tau: number): number {
    const seg = this.seg!
    const slots = TransientPlate.TAU_SLOTS
    const hbits = (tau * 2654435761) | 0
    const slot = (hbits ^ (hbits >>> 3)) & (slots - 1)
    const n = this.ndof
    const base = slot * n
    if (this.tauSlotTauQ[slot] === tau) return base
    this.tauSlotTauQ[slot] = tau
    const hseg = seg.t1 - seg.t0
    const s = (tau - seg.t0) / hseg
    const s2 = s * s; const s3 = s2 * s
    const A = 2 * s3 - 3 * s2 + 1; const B = (s3 - 2 * s2 + s) * hseg
    const C = -2 * s3 + 3 * s2; const D = (s3 - s2) * hseg
    const { q0, qd0, q1, qd1 } = seg
    for (let t = 0; t < this.liveCount; t++) {
      const j = this.liveDofs[t]
      this.tauQQ[base + j] = A * q0[j] + B * qd0[j] + C * q1[j] + D * qd1[j]
    }
    return base
  }

  /** Plate height h(x, y, tau) only (M6-6): replicates surfAt's hh
   *  accumulation term for term (same radial bin, same trig recurrence,
   *  same summation order), so the result is bitwise identical to
   *  surfAt's out[0]. Gradients / velocity / acceleration are skipped.
   *  Used where only the height is compared (contact sampling and
   *  bisection); the contact resolution itself still uses surfAt. */
  heightAt(x: number, y: number, tau: number): number {
    this.heightCalls++
    const seg = this.seg
    if (!seg) return 0
    const base = this.modalCoeffsQ(tau)
    const idx = this.spatialSetup(x, y)
    const stride = BD.ngrid
    let hh = 0
    for (let t = 0; t < this.liveCount; t++) {
      const j = this.liveDofs[t]
      const m = this.dofM[j]
      const R = this.Rflat[this.dofMode[j] * stride + idx]
      const W = m === 0 ? R : this.dofCos[j] ? R * CM[m] : R * SM[m]
      hh += this.tauQQ[base + j] * W
    }
    return hh
  }

  /** q+qdd modal coefficients (M6-7C): identical expressions to
   *  modalCoeffs for those two arrays; qd work is skipped. */
  private modalCoeffsQH(tau: number): number {
    const seg = this.seg!
    const slots = TransientPlate.TAU_SLOTS
    const hbits = (tau * 2654435761) | 0
    const slot = (hbits ^ (hbits >>> 3)) & (slots - 1)
    const n = this.ndof
    const base = slot * n
    if (this.tauSlotTauQH[slot] === tau) return base
    this.tauSlotTauQH[slot] = tau
    const hseg = seg.t1 - seg.t0
    const s = (tau - seg.t0) / hseg
    const s2 = s * s; const s3 = s2 * s
    const A = 2 * s3 - 3 * s2 + 1; const B = (s3 - 2 * s2 + s) * hseg
    const C = -2 * s3 + 3 * s2; const D = (s3 - s2) * hseg
    const P = 12 * s - 6; const Q = (6 * s - 4) * hseg
    const Rr = -12 * s + 6; const Ss = (6 * s - 2) * hseg
    const hh = hseg * hseg
    const { q0, qd0, q1, qd1 } = seg
    // Tiny segments (h < 1e-7 s) always use the ODE form: the Hermite
    // qdd divides by h^2 and degenerates into cancellation garbage as
    // h -> 0, while the ODE form is its stable limit (M6-11R).
    if (this.accelMode === 'ode' || hseg < 1e-7) {
      const E = 6 * s2 - 6 * s; const F = (3 * s2 - 4 * s + 1) * hseg
      const Gc = -6 * s2 + 6 * s; const H = (3 * s2 - 2 * s) * hseg
      const ts = tau - this.epoch
      const Fd = this.signal.forcePeakAt(ts) * this.signal.u(ts)
      for (let t = 0; t < this.liveCount; t++) {
        const j = this.liveDofs[t]
        const qv = A * q0[j] + B * qd0[j] + C * q1[j] + D * qd1[j]
        const qdv = (E * q0[j] + F * qd0[j] + Gc * q1[j] + H * qd1[j]) / hseg
        this.tauQH[base + j] = qv
        this.tauQddH[base + j] = (this.phiDrive[j] * Fd) / this.mass[j]
          - 2 * this.zeta * this.omega[j] * qdv - this.omega[j] ** 2 * qv
      }
      return base
    }
    for (let t = 0; t < this.liveCount; t++) {
      const j = this.liveDofs[t]
      this.tauQH[base + j] = A * q0[j] + B * qd0[j] + C * q1[j] + D * qd1[j]
      this.tauQddH[base + j] = (P * q0[j] + Q * qd0[j] + Rr * q1[j] + Ss * qd1[j]) / hh
    }
    return base
  }

  /** Spatial weights setup shared by heightAt/heightAccelAt/accelAt:
   *  fills CM/SM exactly as surfAt does and returns the radial bin. */
  private spatialSetup(x: number, y: number): number {
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const rs = r < 1e-12 ? 1e-12 : r
    const c = x / rs; const sn = y / rs
    CM[0] = 1; CM[1] = c; SM[0] = 0; SM[1] = sn
    for (let m = 2; m <= 4; m++) {
      CM[m] = CM[m - 1] * c - SM[m - 1] * sn
      SM[m] = SM[m - 1] * c + CM[m - 1] * sn
    }
    return idx
  }

  /** (h, wtt) at (x, y, tau) only (M6-7C): bitwise identical to
   *  surfAt's out[0] and out[4] (same W, same summation order).
   *  Used by the sitting detach check, which needs exactly these two. */
  heightAccelAt(x: number, y: number, tau: number, out2: Float64Array): void {
    this.heightAccelCalls++
    if (!this.seg) { out2[0] = 0; out2[1] = 0; return }
    const base = this.modalCoeffsQH(tau)
    const idx = this.spatialSetup(x, y)
    const stride = BD.ngrid
    let hh = 0; let wtt = 0
    for (let t = 0; t < this.liveCount; t++) {
      const j = this.liveDofs[t]
      const m = this.dofM[j]
      const R = this.Rflat[this.dofMode[j] * stride + idx]
      const W = m === 0 ? R : this.dofCos[j] ? R * CM[m] : R * SM[m]
      hh += this.tauQH[base + j] * W
      wtt += this.tauQddH[base + j] * W
    }
    out2[0] = hh; out2[1] = wtt
  }

  /** Height with fused coefficient evaluation (M6-7A): same Hermite
   *  expressions as modalCoeffsQ but computed inline per DOF without
   *  the tau slot cache. Bitwise identical to heightAt; intended for
   *  callers whose taus are unique (contact bisection), where the
   *  cache cannot hit. Sample evaluations keep using heightAt so
   *  shared sample taus still hit the cache. */
  heightAtFused(x: number, y: number, tau: number): number {
    this.heightCalls++
    const seg = this.seg
    if (!seg) return 0
    const idx = this.spatialSetup(x, y)
    const hseg = seg.t1 - seg.t0
    const s = (tau - seg.t0) / hseg
    const s2 = s * s; const s3 = s2 * s
    const A = 2 * s3 - 3 * s2 + 1; const B = (s3 - 2 * s2 + s) * hseg
    const C = -2 * s3 + 3 * s2; const D = (s3 - s2) * hseg
    const { q0, qd0, q1, qd1 } = seg
    const stride = BD.ngrid
    let hh = 0
    for (let t = 0; t < this.liveCount; t++) {
      const j = this.liveDofs[t]
      const m = this.dofM[j]
      const R = this.Rflat[this.dofMode[j] * stride + idx]
      const W = m === 0 ? R : this.dofCos[j] ? R * CM[m] : R * SM[m]
      hh += (A * q0[j] + B * qd0[j] + C * q1[j] + D * qd1[j]) * W
    }
    return hh
  }

  /** wtt at (x, y, tau) only (M6-7C): bitwise identical to surfAt's out[4]. */
  accelAt(x: number, y: number, tau: number): number {
    this.heightAccelCalls++
    if (!this.seg) return 0
    const base = this.modalCoeffsQH(tau)
    const idx = this.spatialSetup(x, y)
    const stride = BD.ngrid
    let wtt = 0
    for (let t = 0; t < this.liveCount; t++) {
      const j = this.liveDofs[t]
      const m = this.dofM[j]
      const R = this.Rflat[this.dofMode[j] * stride + idx]
      const W = m === 0 ? R : this.dofCos[j] ? R * CM[m] : R * SM[m]
      wtt += this.tauQddH[base + j] * W
    }
    return wtt
  }

  /** Evaluation counter for wttLowerBound (research accounting). */
  wttBoundCalls = 0
  private boundPhi: Float64Array | null = null

  /** Conservative lower bound of the *implemented* wtt at fixed (x, y)
   *  over [ta, tb] inside the current segment. Hermite form is linear
   *  in tau, so its bound is the exact endpoint minimum. The ODE form
   *  (also used for tiny segments) is bounded by interval evaluation
   *  of the equation of motion over the exact Hermite q/qd ranges on
   *  the subinterval and the drive-force bound |F| <= max overlapping
   *  forcePeak, widened by a rounding hair (same construction as
   *  lookaheadRanges). Returns -Infinity when no segment is active. */
  wttLowerBound(x: number, y: number, ta: number, tb: number): number {
    const seg = this.seg
    if (!seg) return -Infinity
    this.wttBoundCalls++
    const hseg = seg.t1 - seg.t0
    if (this.accelMode !== 'ode' && hseg >= 1e-7) {
      return Math.min(this.accelAt(x, y, ta), this.accelAt(x, y, tb))
    }
    const n = this.ndof
    const baseA = this.modalCoeffs(ta)
    const qa = this.tauQ.slice(baseA, baseA + n)
    const qda = this.tauQd.slice(baseA, baseA + n)
    const baseB = this.modalCoeffs(tb)
    const qb = this.tauQ.slice(baseB, baseB + n)
    const qdb = this.tauQd.slice(baseB, baseB + n)
    if (!this.boundPhi) this.boundPhi = new Float64Array(n)
    const phi = this.boundPhi
    this.phiAll(x, y, phi)
    const [fLo, fHi] = this.forceRange(ta, tb)
    const h = tb - ta
    let lb = 0
    for (let j = 0; j < n; j++) {
      const rq = hermiteQRange(qa[j], qda[j], qb[j], qdb[j], h)
      const rd = hermiteQdRange(qa[j], qda[j], qb[j], qdb[j], h)
      const w = this.omega[j]
      const c = this.phiDrive[j] / this.mass[j]
      const dLo = Math.min(c * fLo, c * fHi); const dHi = Math.max(c * fLo, c * fHi)
      let lo = dLo - 2 * this.zeta * w * rd[1] - w * w * rq[1]
      let hi = dHi - 2 * this.zeta * w * rd[0] - w * w * rq[0]
      // Rounding hair: the discrepancy between modalCoeffs' evaluation
      // order and the exact range math is ulp-scale on q (~1e-20 abs)
      // amplified by w^2 (~1e-12 abs per DOF). 1e-12 relative matches
      // the segQRange widening; the former 1e-9 left an unprovable
      // ~1e-5 m/s2 band around the detach threshold (M6-12R).
      const wd = 1e-12 * (1 + Math.abs(lo) + Math.abs(hi))
      lo -= wd; hi += wd
      lb += Math.min(phi[j] * lo, phi[j] * hi)
    }
    return lb
  }

  /** Conservative [min, max] of the drive force F(t) over [ta, tb]
   *  (physics times). Within one score event the envelope is
   *  up-then-flat-then-down and the phase is monotone, so exact ranges
   *  follow from endpoints plus plateau/extremum containment. Across
   *  event boundaries (or outside events) fall back to the global
   *  |F| <= max overlapping forcePeak bound. */
  private forceRange(ta: number, tb: number): [number, number] {
    const ts0 = ta - this.epoch; const ts1 = tb - this.epoch
    const sig = this.signal
    // Split at every event boundary inside the interval and union the
    // per-piece exact ranges (force is 0 outside events). Without the
    // split, a subinterval straddling a boundary would keep the loose
    // global bound at any width and interval refinement could never
    // conclude around short events (M6-12R).
    const cuts = [ts0, ts1]
    for (const e of sig.score.events) {
      if (e.start_s > ts0 && e.start_s < ts1) cuts.push(e.start_s)
      if (e.end_s > ts0 && e.end_s < ts1) cuts.push(e.end_s)
    }
    cuts.sort((a, b) => a - b)
    let lo = Infinity; let hi = -Infinity
    for (let ci = 0; ci + 1 < cuts.length; ci++) {
      const a = cuts[ci]; const b = cuts[ci + 1]
      const hit = sig.eventAt(a) ?? sig.eventAt(0.5 * (a + b))
      if (!hit) { lo = Math.min(lo, 0); hi = Math.max(hi, 0); continue }
      const [plo, phi] = this.forceRangeInEvent(hit.ev, a, b)
      lo = Math.min(lo, plo); hi = Math.max(hi, phi)
    }
    // Include exact endpoint force values (eventAt is start-inclusive).
    const fA = sig.forcePeakAt(ts0) * sig.u(ts0)
    const fB = sig.forcePeakAt(ts1) * sig.u(ts1)
    lo = Math.min(lo, fA, fB); hi = Math.max(hi, fA, fB)
    return [lo, hi]
  }

  /** Exact force range on [ts0, ts1] fully inside one score event. */
  private forceRangeInEvent(ev: DriveEvent, ts0: number, ts1: number): [number, number] {
    const sig = this.signal
    // Envelope replica of ScoreSignal.u (smoothstep attack/release).
    const envAt = (t: number): number => {
      let env = 1
      if (ev.attack_s > 0 && t < ev.start_s + ev.attack_s) {
        const u = (t - ev.start_s) / ev.attack_s
        env = Math.min(env, u * u * (3 - 2 * u))
      }
      if (ev.release_s > 0 && t > ev.end_s - ev.release_s) {
        const u = (ev.end_s - t) / ev.release_s
        env = Math.min(env, u * u * (3 - 2 * u))
      }
      return env
    }
    let eLo = Math.min(envAt(ts0), envAt(ts1))
    let eHi = Math.max(envAt(ts0), envAt(ts1))
    if (ts0 < ev.end_s - ev.release_s && ts1 > ev.start_s + ev.attack_s) eHi = 1
    // Overlapping attack/release (attack+release > duration, or a
    // sub-interval inside the overlap): env = min(increasing attack,
    // decreasing release) peaks where the two ramps cross, below 1.
    // smooth() is monotone, so the crossing solves
    // (t-start)/attack = (end-t)/release (M6-12R review).
    if (ev.attack_s > 0 && ev.release_s > 0) {
      const tc = (ev.release_s * ev.start_s + ev.attack_s * ev.end_s) / (ev.attack_s + ev.release_s)
      if (tc >= ts0 && tc <= ts1) eHi = Math.max(eHi, envAt(tc))
    }
    if (eLo < 0) eLo = 0
    if (eHi > 1) eHi = 1
    const thA = sig.thetaAt(ts0); const thB = sig.thetaAt(ts1)
    let sLo = Math.min(Math.sin(thA), Math.sin(thB))
    let sHi = Math.max(Math.sin(thA), Math.sin(thB))
    if (thB - thA >= 2 * Math.PI) { sLo = -1; sHi = 1 } else {
      const k0 = Math.ceil((thA - Math.PI / 2) / (2 * Math.PI))
      if (Math.PI / 2 + 2 * Math.PI * k0 <= thB) sHi = 1
      const k1 = Math.ceil((thA - 3 * Math.PI / 2) / (2 * Math.PI))
      if (3 * Math.PI / 2 + 2 * Math.PI * k1 <= thB) sLo = -1
    }
    const pLo = sLo >= 0 ? eLo * sLo : eHi * sLo
    const pHi = sHi >= 0 ? eHi * sHi : eLo * sHi
    const F = ev.forcePeak_N
    return F >= 0 ? [F * pLo, F * pHi] : [F * pHi, F * pLo]
  }

  /** First time in [ta, tb] (inside the current segment) at which the
   *  implemented wtt at fixed (x, y) crosses below -G.
   *  Tri-state result (M6-12R): 'found' (first crossing, time order
   *  guaranteed), 'none' (proven: every subinterval either sampled
   *  clear with a wttLowerBound proof), or 'unknown' (bound budget
   *  exhausted before a proof completed — NOT a proof of absence;
   *  callers must retry with more budget or fail explicitly).
   *  Scan subintervals are resolved strictly left-to-right: a later
   *  sampled crossing is never accepted before earlier subintervals
   *  are proven clear or yield an earlier root. Roots inside a
   *  1e-11 s window are located by bisection (time error <= window). */
  /** True when the current segment's acceleration uses the ODE form:
   *  explicit 'ode' mode, or the tiny-segment (h < 1e-7) stable
   *  fallback. Callers that assume wtt linear in time (Hermite form)
   *  must not run when this is true. Canonical steps are never tiny,
   *  so default-mode production stepping always returns false here. */
  odeFormActive(): boolean {
    return this.accelMode === 'ode' || (!!this.seg && this.seg.t1 - this.seg.t0 < 1e-7)
  }

  firstDetachBelow(x: number, y: number, ta: number, tb: number, budgetIn = 2048): DetachSearchResult {
    if (!this.seg) return { kind: 'unknown' }
    const bisect = (a: number, b: number): number => {
      let lo = a; let hi = b
      for (let it = 0; it < 40; it++) {
        const m = 0.5 * (lo + hi)
        if (this.accelAt(x, y, m) >= -G) lo = m; else hi = m
      }
      return 0.5 * (lo + hi)
    }
    let budget = budgetIn
    // Left-first recursive resolution of one interval whose endpoint
    // values are known (w(a) >= -G invariant). 'none' is returned
    // only when every leaf is bound-proven; a bracketed root at
    // width <= 1e-9 s is located by bisection (documented resolution
    // floor: an earlier dip inside that window cannot be told apart).
    const search = (a: number, b: number, wa: number, wb: number, depth: number): DetachSearchResult => {
      if (wb < -G) {
        if (depth <= 0 || b - a < 1e-9) return { kind: 'found', t: bisect(a, b) }
        const m = 0.5 * (a + b)
        const wm = this.accelAt(x, y, m)
        const left = search(a, m, wa, wm, depth - 1)
        if (left.kind === 'found') return left
        if (left.kind === 'unknown') {
          // A root in [a, b] is certain (sampled bracket); only its
          // "firstness" inside the unresolved left part is in doubt.
          // At this width that doubt is below 1e-9 s of model time.
          if (b - a < 1e-8) return { kind: 'found', t: bisect(a, b) }
          return left
        }
        return search(m, b, wm, wb, depth - 1)
      }
      if (budget-- <= 0) return { kind: 'unknown' }
      if (this.wttLowerBound(x, y, a, b) >= -G) return { kind: 'none' }
      if (depth <= 0 || b - a < 1e-12) return { kind: 'unknown' }
      const m = 0.5 * (a + b)
      const wm = this.accelAt(x, y, m)
      const left = search(a, m, wa, wm, depth - 1)
      if (left.kind !== 'none') return left
      return search(m, b, wm, wb, depth - 1)
    }
    const wa0 = this.accelAt(x, y, ta)
    if (wa0 < -G) return { kind: 'found', t: ta }
    // Quick whole-interval proof (the common firmly-seated case).
    const wb0 = this.accelAt(x, y, tb)
    if (wb0 >= -G && budget-- > 0 && this.wttLowerBound(x, y, ta, tb) >= -G) {
      return { kind: 'none' }
    }
    // Resolve scan subintervals strictly left-to-right.
    const K = 8
    let tp = ta; let wp = wa0
    for (let k = 1; k <= K; k++) {
      const t = k === K ? tb : ta + ((tb - ta) * k) / K
      const w = k === K ? wb0 : this.accelAt(x, y, t)
      const r = search(tp, t, wp, w, 40)
      if (r.kind !== 'none') return r
      tp = t; wp = w
    }
    return { kind: 'none' }
  }

  /** Pre-M6-4A surfAt, kept verbatim as the comparison reference. */
  surfAtReference(x: number, y: number, tau: number, out: Float64Array): void {
    const seg = this.seg
    if (!seg) { out.fill(0); return }
    const hseg = seg.t1 - seg.t0
    const s = (tau - seg.t0) / hseg
    const s2 = s * s; const s3 = s2 * s
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const rs = r < 1e-12 ? 1e-12 : r
    const c = x / rs; const sn = y / rs
    // cos(m*theta)/sin(m*theta) by recurrence from (c, sn) — no atan2/trig.
    const cmByM = [1, c, 0, 0, 0, 0]; const smByM = [0, sn, 0, 0, 0, 0]
    for (let m = 2; m <= 4; m++) {
      cmByM[m] = cmByM[m - 1] * c - smByM[m - 1] * sn
      smByM[m] = smByM[m - 1] * c + cmByM[m - 1] * sn
    }
    let hh = 0; let hx = 0; let hy = 0; let uz = 0; let wtt = 0
    for (let j = 0; j < this.ndof; j++) {
      const q = (2 * s3 - 3 * s2 + 1) * seg.q0[j] + (s3 - 2 * s2 + s) * hseg * seg.qd0[j]
        + (-2 * s3 + 3 * s2) * seg.q1[j] + (s3 - s2) * hseg * seg.qd1[j]
      const qd = ((6 * s2 - 6 * s) * seg.q0[j] + (3 * s2 - 4 * s + 1) * hseg * seg.qd0[j]
        + (-6 * s2 + 6 * s) * seg.q1[j] + (3 * s2 - 2 * s) * hseg * seg.qd1[j]) / hseg
      const qdd = ((12 * s - 6) * seg.q0[j] + (6 * s - 4) * hseg * seg.qd0[j]
        + (-12 * s + 6) * seg.q1[j] + (6 * s - 2) * hseg * seg.qd1[j]) / (hseg * hseg)
      const d = this.dofs[j]; const md = this.modes[d.mode]; const m = md.m
      const R = this.R[d.mode][idx]; const dR = this.dR[d.mode][idx] / this.a
      let W: number; let Wx: number; let Wy: number
      if (m === 0) {
        W = R; Wx = dR * c; Wy = dR * sn
      } else {
        const cm = cmByM[m]; const sm = smByM[m]
        let Wr: number; let Wt: number
        if (d.kind === 'cos') { W = R * cm; Wr = dR * cm; Wt = (-m * R * sm) / rs }
        else { W = R * sm; Wr = dR * sm; Wt = (m * R * cm) / rs }
        Wx = Wr * c - Wt * sn; Wy = Wr * sn + Wt * c
      }
      hh += q * W; hx += q * Wx; hy += q * Wy; uz += qd * W; wtt += qdd * W
    }
    out[0] = hh; out[1] = hx; out[2] = hy; out[3] = uz; out[4] = wtt
  }

  vScale(): number {
    const q = this.seg ? this.seg.q1 : this.q
    const qd = this.seg ? this.seg.qd1 : this.qd
    let s = 0
    for (let j = 0; j < this.ndof; j++) s += qd[j] ** 2 + (this.omega[j] * q[j]) ** 2
    return Math.sqrt(s)
  }

  /** Displacement scale sqrt(sum q^2 + (qd/omega)^2) at segment end [m]. */
  ampScale(): number {
    const q = this.seg ? this.seg.q1 : this.q
    const qd = this.seg ? this.seg.qd1 : this.qd
    let s = 0
    for (let j = 0; j < this.ndof; j++) s += q[j] ** 2 + (qd[j] / this.omega[j]) ** 2
    return Math.sqrt(s)
  }

  forceForAmplitude(modeIndex: number, targetA: number): number {
    const j = this.dofs.findIndex((d) => d.mode === modeIndex)
    return (targetA * this.mass[j] * 2 * this.zeta * this.omega[j] ** 2) / Math.abs(this.phiDrive[j])
  }

  /**
   * Exact per-DOF ranges of the Hermite q and implemented qdd over a
   * future window, obtained by advancing a clone with the same dt
   * sequence (dtAt) and the same RK4 substeps the real advance will
   * use. Hermite qdd is linear per segment (endpoints suffice); q is
   * cubic (endpoints + critical points). In 'ode' accel mode the
   * qdd range is instead the interval evaluation of the ODE right
   * hand side over the exact Hermite q/qd ranges and the segment's
   * drive-force bound |F| <= max overlapping forcePeak (envelope and
   * sine are both within [-1, 1]); widened by a rounding hair.
   * Valid only while real steps follow dtAt exactly (no shortened steps).
   */
  lookaheadRanges(fromT: number, windowS: number, dtAt: (t: number) => number) {
    const n = this.ndof
    const qMin = new Float64Array(n).fill(Infinity); const qMax = new Float64Array(n).fill(-Infinity)
    const aMin = new Float64Array(n).fill(Infinity); const aMax = new Float64Array(n).fill(-Infinity)
    const q = new Float64Array(this.q); const qd = new Float64Array(this.qd)
    let t = fromT
    const end = fromT + windowS
    while (t < end - 1e-15) {
      const dt = Math.min(dtAt(t), end - t)
      const [q1, qd1] = this.advanceState(q, qd, t, t + dt)
      const h = dt
      let fpk = 0
      // Match the implemented acceleration definition per segment
      // (tiny segments use the ODE form even in 'hermite' mode).
      const odeSeg = this.accelMode === 'ode' || h < 1e-7
      if (odeSeg) {
        const ts0 = t - this.epoch; const ts1 = t + dt - this.epoch
        for (const e of this.signal.score.events) {
          if (e.start_s <= ts1 && e.end_s >= ts0) fpk = Math.max(fpk, Math.abs(e.forcePeak_N))
        }
      }
      for (let j = 0; j < n; j++) {
        const r = hermiteQRange(q[j], qd[j], q1[j], qd1[j], h)
        if (r[0] < qMin[j]) qMin[j] = r[0]
        if (r[1] > qMax[j]) qMax[j] = r[1]
        let lo: number; let hi: number
        if (odeSeg) {
          const rd = hermiteQdRange(q[j], qd[j], q1[j], qd1[j], h)
          const w = this.omega[j]
          const dAbs = Math.abs(this.phiDrive[j] / this.mass[j]) * fpk
          lo = -dAbs - 2 * this.zeta * w * rd[1] - w * w * r[1]
          hi = dAbs - 2 * this.zeta * w * rd[0] - w * w * r[0]
          const wd = 1e-9 * (1 + Math.abs(lo) + Math.abs(hi))
          lo -= wd; hi += wd
        } else {
          const a0 = (-6 * q[j] - 4 * h * qd[j] + 6 * q1[j] - 2 * h * qd1[j]) / (h * h)
          const a1 = (6 * q[j] + 2 * h * qd[j] - 6 * q1[j] + 4 * h * qd1[j]) / (h * h)
          lo = Math.min(a0, a1); hi = Math.max(a0, a1)
        }
        if (lo < aMin[j]) aMin[j] = lo
        if (hi > aMax[j]) aMax[j] = hi
      }
      q.set(q1); qd.set(qd1); t += dt
    }
    return { end, qMin, qMax, aMin, aMax }
  }

  modeIndex152(): number { return this.modes.findIndex((m) => m.m === 0 && m.root_idx === 0) }
  resonanceFreqs(): number[] { return this.modes.map((m) => m.f_hz) }
}

/** Exact [min, max] of the cubic Hermite q(s), s in [0,1]. */
export function hermiteQRange(q0: number, qd0: number, q1: number, qd1: number, h: number): [number, number] {
  const qAt = (s: number): number => {
    const s2 = s * s; const s3 = s2 * s
    return (2 * s3 - 3 * s2 + 1) * q0 + (s3 - 2 * s2 + s) * h * qd0
      + (-2 * s3 + 3 * s2) * q1 + (s3 - s2) * h * qd1
  }
  let lo = Math.min(q0, q1); let hi = Math.max(q0, q1)
  // qd numerator: A s^2 + B s + C = 0
  const A = 6 * (q0 - q1) + 3 * h * (qd0 + qd1)
  const B = -6 * q0 + 6 * q1 + h * (-4 * qd0 - 2 * qd1)
  const C = h * qd0
  const roots: number[] = []
  if (Math.abs(A) < 1e-300) {
    if (Math.abs(B) > 0) roots.push(-C / B)
  } else {
    const disc = B * B - 4 * A * C
    if (disc >= 0) {
      const sq = Math.sqrt(disc)
      roots.push((-B - sq) / (2 * A), (-B + sq) / (2 * A))
    }
  }
  for (const s of roots) {
    if (s > 0 && s < 1) {
      const v = qAt(s)
      if (v < lo) lo = v
      if (v > hi) hi = v
    }
  }
  return [lo, hi]
}

/** Exact [min, max] of the Hermite velocity qd(s) = (A s^2 + B s + C)/h,
 *  s in [0,1] (same A, B, C as hermiteQRange's qd numerator). */
export function hermiteQdRange(q0: number, qd0: number, q1: number, qd1: number, h: number): [number, number] {
  const A = 6 * (q0 - q1) + 3 * h * (qd0 + qd1)
  const B = -6 * q0 + 6 * q1 + h * (-4 * qd0 - 2 * qd1)
  const C = h * qd0
  let lo = Math.min(qd0, qd1); let hi = Math.max(qd0, qd1)
  if (Math.abs(A) > 0) {
    const sv = -B / (2 * A)
    if (sv > 0 && sv < 1) {
      const v = (A * sv * sv + B * sv + C) / h
      if (v < lo) lo = v
      if (v > hi) hi = v
    }
  }
  return [lo, hi]
}
