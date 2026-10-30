// plateTransient.ts — multi-mode transient plate (TS port of
// reference/plate_transient.py). See that file for PUT-IN / floors.
// qdd + 2*zeta*omega*qd + omega^2 q = Phi(drive)*F(t)/M, RK4 substeps,
// cubic-Hermite state evaluation inside a sand step.
import basisData from '../data/plate_modes_transient.json'
import { ScoreSignal, type PerformanceScore } from './sequence'

interface ModeRow { m: number; root_idx: number; f_hz: number; M_kg: number; R: number[]; dR_drho: number[] }
interface BasisData { a_m: number; ngrid: number; zeta: number; modes: ModeRow[] }
const BD = basisData as unknown as BasisData

export interface Dof { mode: number; kind: 'cos' | 'sin' }

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
  signal: ScoreSignal
  epoch = 0 // score time = physics t - epoch (live score replacement)

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
    this.setDrive(driveXY[0], driveXY[1])
  }

  get ndof(): number { return this.dofs.length }

  /** Phi for one DOF at (x,y) using precomputed radial index + angle terms. */
  private phi(modeIdx: number, kind: string, idx: number, cm: number, sm: number): number {
    const R = this.R[modeIdx][idx]
    if (this.modes[modeIdx].m === 0) return R
    return kind === 'cos' ? R * cm : R * sm
  }

  setDrive(x: number, y: number): void {
    this.driveXY = [x, y]
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const th = Math.atan2(y, x)
    this.dofs.forEach((d, j) => {
      const m = this.modes[d.mode].m
      this.phiDrive[j] = this.phi(d.mode, d.kind, idx, Math.cos(m * th), Math.sin(m * th))
    })
  }

  /** Replace the score at physics time atTime; q/qd are preserved. */
  setScore(score: PerformanceScore, atTime: number): void {
    this.signal = new ScoreSignal(score)
    this.epoch = atTime
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

  /** All DOF shape values Phi_j(x, y) (no gradients). */
  phiAll(x: number, y: number, out: Float64Array): void {
    const r = Math.hypot(x, y)
    const idx = Math.min((r * (this.ng1 / this.a)) | 0, this.ng1)
    const rs = r < 1e-12 ? 1e-12 : r
    const c = x / rs; const sn = y / rs
    const cm = [1, c, 0, 0, 0]; const sm = [0, sn, 0, 0, 0]
    for (let m = 2; m <= 4; m++) {
      cm[m] = cm[m - 1] * c - sm[m - 1] * sn
      sm[m] = sm[m - 1] * c + cm[m - 1] * sn
    }
    for (let j = 0; j < this.ndof; j++) {
      const d = this.dofs[j]; const m = this.modes[d.mode].m
      const R = this.R[d.mode][idx]
      out[j] = m === 0 ? R : d.kind === 'cos' ? R * cm[m] : R * sm[m]
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
    const n = this.ndof
    const q0 = Float64Array.from(this.q); const qd0 = Float64Array.from(this.qd)
    const [q, qd] = this.advanceState(q0, qd0, t0, t1)
    this.q = q; this.qd = qd; this.t = t1
    this.seg = { t0, t1, q0, qd0, q1: Float64Array.from(q), qd1: Float64Array.from(qd) }
    // Conservative |h| bound over the segment (see Python plate_transient).
    // Hermite tangent basis max is 4/27 (at s=1/3), not 0.105.
    let hb = 0
    for (let j = 0; j < n; j++) {
      hb += Math.max(Math.abs(q0[j]), Math.abs(q[j])) + (4 / 27) * (t1 - t0) * (Math.abs(qd0[j]) + Math.abs(qd[j]))
    }
    this.hBound = hb
  }

  /** (h, hx, hy, uz, wtt) at (x, y, tau) inside the current segment. */
  surfAt(x: number, y: number, tau: number, out: Float64Array): void {
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
   * Exact per-DOF ranges of the Hermite q and qdd over a future window,
   * obtained by advancing a clone with the same dt sequence (dtAt) and the
   * same RK4 substeps the real advance will use. qdd is linear per segment
   * (endpoints suffice); q is cubic (endpoints + critical points).
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
      for (let j = 0; j < n; j++) {
        const r = hermiteQRange(q[j], qd[j], q1[j], qd1[j], h)
        if (r[0] < qMin[j]) qMin[j] = r[0]
        if (r[1] > qMax[j]) qMax[j] = r[1]
        const a0 = (-6 * q[j] - 4 * h * qd[j] + 6 * q1[j] - 2 * h * qd1[j]) / (h * h)
        const a1 = (6 * q[j] + 2 * h * qd[j] - 6 * q1[j] + 4 * h * qd1[j]) / (h * h)
        const lo = Math.min(a0, a1); const hi = Math.max(a0, a1)
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
