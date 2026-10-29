// scoreWorklet.ts — streaming score audio (no full-length PCM allocation).
// The AudioWorklet processor integrates phase sample-by-sample from the same
// score definition the physics uses. Frequency changes keep theta continuous.
// Memory: O(events), independent of performance length (phone-friendly).
import type { PerformanceScore } from '../simulation/sequence'

const PROCESSOR_JS = `
class ScoreProcessor extends AudioWorkletProcessor {
  constructor() {
    super()
    this.events = []; this.theta0 = []; this.gain = 0.15
    this.theta = 0; this.tScore = 0; this.toneHz = null
    this.port.onmessage = (ev) => {
      const m = ev.data
      if (m.type === 'score') {
        this.events = m.events; this.gain = m.gain
        this.theta = m.initialPhase; this.tScore = 0; this.toneHz = null
        this.theta0 = []; let th = m.initialPhase
        for (const e of this.events) { this.theta0.push(th); th += 2*Math.PI*e.frequency_hz*(e.end_s-e.start_s) }
      } else if (m.type === 'freq') {
        this.toneHz = m.hz // keep theta: phase-continuous switch
      } else if (m.type === 'gain') {
        this.gain = m.gain
      }
    }
  }
  envAt(e, t) {
    let env = 1
    const sm = (u) => u*u*(3-2*u)
    if (e.attack_s > 0 && t < e.start_s + e.attack_s) env = Math.min(env, sm((t-e.start_s)/e.attack_s))
    if (e.release_s > 0 && t > e.end_s - e.release_s) env = Math.min(env, sm((e.end_s-t)/e.release_s))
    return env
  }
  process(_inputs, outputs) {
    const out = outputs[0][0]
    const dt = 1 / sampleRate
    for (let i = 0; i < out.length; i++) {
      let f = 0; let env = 0
      if (this.toneHz !== null) { f = this.toneHz; env = 1 }
      else {
        for (let k = 0; k < this.events.length; k++) {
          const e = this.events[k]
          if (this.tScore >= e.start_s && this.tScore < e.end_s) { f = e.frequency_hz; env = this.envAt(e, this.tScore); break }
        }
      }
      out[i] = env > 0 ? this.gain * Math.sin(this.theta) : 0
      this.theta += 2 * Math.PI * f * dt
      this.tScore += dt
    }
    return true
  }
}
registerProcessor('score-processor', ScoreProcessor)
`

export class ScoreStreamAudio {
  private ctx: AudioContext | null = null
  private node: AudioWorkletNode | null = null
  private offset = 0
  private segmentStart = 0
  playing = false

  private async ensure(): Promise<void> {
    if (this.node) return
    if (!this.ctx) this.ctx = new AudioContext()
    const url = URL.createObjectURL(new Blob([PROCESSOR_JS], { type: 'application/javascript' }))
    await this.ctx.audioWorklet.addModule(url)
    URL.revokeObjectURL(url)
    this.node = new AudioWorkletNode(this.ctx, 'score-processor')
    this.node.connect(this.ctx.destination)
  }

  async load(score: PerformanceScore): Promise<void> {
    await this.ensure()
    this.offset = 0; this.playing = false
    this.node?.port.postMessage({
      type: 'score', events: score.events, gain: score.audioGain,
      initialPhase: score.initialPhase_rad,
    })
  }

  setFrequency(hz: number): void {
    this.node?.port.postMessage({ type: 'freq', hz })
  }

  async play(): Promise<void> {
    await this.ensure()
    if (!this.ctx) return
    await this.ctx.resume()
    this.segmentStart = this.ctx.currentTime
    this.playing = true
  }

  async pause(): Promise<void> {
    if (!this.ctx || !this.playing) return
    this.offset = this.position()
    await this.ctx.suspend()
    this.playing = false
  }

  position(): number {
    if (!this.ctx) return this.offset
    return this.playing ? this.offset + (this.ctx.currentTime - this.segmentStart) : this.offset
  }

  async suspendForLag(): Promise<void> { if (this.playing) await this.ctx?.suspend() }
  async resumeFromLag(): Promise<void> {
    if (!this.playing || !this.ctx) return
    const before = this.ctx.currentTime
    await this.ctx.resume()
    this.segmentStart += this.ctx.currentTime - before
  }
}
