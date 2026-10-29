// sequence.ts — shared drive score: one timeline for force AND audio.
// u(t) = envelope(t)*sin(theta(t)), theta = 2*pi*integral f dt.
// Force = forcePeak_N * u(t) [N]; audio PCM = audioGain * u(t) [dimensionless].
export interface DriveEvent {
  start_s: number; end_s: number; frequency_hz: number
  forcePeak_N: number; attack_s: number; release_s: number
}

export interface PerformanceScore {
  schema: 'physicymatics.score.v1'
  events: DriveEvent[]
  duration_s: number
  initialPhase_rad: number
  audioGain: number
}

function smooth(u: number): number { return u * u * (3 - 2 * u) }

export class ScoreSignal {
  private theta0: number[] = []
  constructor(readonly score: PerformanceScore) {
    let th = score.initialPhase_rad
    for (const e of score.events) {
      this.theta0.push(th)
      th += 2 * Math.PI * e.frequency_hz * (e.end_s - e.start_s)
    }
  }

  eventAt(t: number): { index: number; ev: DriveEvent } | null {
    for (let i = 0; i < this.score.events.length; i++) {
      const e = this.score.events[i]
      if (t >= e.start_s && t < e.end_s) return { index: i, ev: e }
    }
    return null
  }

  u(t: number): number {
    const hit = this.eventAt(t)
    if (!hit) return 0
    const { index: i, ev: e } = hit
    let env = 1
    if (e.attack_s > 0 && t < e.start_s + e.attack_s) env = Math.min(env, smooth((t - e.start_s) / e.attack_s))
    if (e.release_s > 0 && t > e.end_s - e.release_s) env = Math.min(env, smooth((e.end_s - t) / e.release_s))
    return env * Math.sin(this.theta0[i] + 2 * Math.PI * e.frequency_hz * (t - e.start_s))
  }

  forcePeakAt(t: number): number { return this.eventAt(t)?.ev.forcePeak_N ?? 0 }
  freqAt(t: number): number | null { return this.eventAt(t)?.ev.frequency_hz ?? null }

  /** Integrated phase theta(t) (continuous across events). Outside events,
   *  holds the phase at the nearest event boundary. */
  thetaAt(t: number): number {
    const evs = this.score.events
    if (!evs.length) return this.score.initialPhase_rad
    if (t <= evs[0].start_s) return this.theta0[0]
    for (let i = 0; i < evs.length; i++) {
      const e = evs[i]
      if (t < e.start_s) return this.theta0[i]
      if (t < e.end_s) return this.theta0[i] + 2 * Math.PI * e.frequency_hz * (t - e.start_s)
    }
    const last = evs[evs.length - 1]
    return this.theta0[evs.length - 1] + 2 * Math.PI * last.frequency_hz * (last.end_s - last.start_s)
  }
}

/** PCM from the same score: sample k is at time k/sampleRate, starting at 0. */
export function renderScorePcm(score: PerformanceScore, sampleRate: number): Float32Array {
  const sig = new ScoreSignal(score)
  const n = Math.round(score.duration_s * sampleRate)
  const out = new Float32Array(n)
  for (let k = 0; k < n; k++) out[k] = score.audioGain * sig.u(k / sampleRate)
  return out
}

export function singleToneScore(freqHz: number, forcePeakN: number, durationS = 3600): PerformanceScore {
  return {
    schema: 'physicymatics.score.v1', duration_s: durationS, initialPhase_rad: 0, audioGain: 0.2,
    events: [{ start_s: 0, end_s: durationS, frequency_hz: freqHz, forcePeak_N: forcePeakN, attack_s: 0.05, release_s: 0 }],
  }
}
