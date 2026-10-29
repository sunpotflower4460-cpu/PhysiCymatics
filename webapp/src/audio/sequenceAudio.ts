// sequenceAudio.ts — plays PCM rendered from the same score as the drive.
// Pause keeps an offset; resume starts a new source at that offset.
// Position = offset + ctx.currentTime - segmentStart (guide §8.5).
import { renderScorePcm, type PerformanceScore } from '../simulation/sequence'

export class ScoreAudio {
  private ctx: AudioContext | null = null
  private buffer: AudioBuffer | null = null
  private source: AudioBufferSourceNode | null = null
  private offset = 0
  private segmentStart = 0
  playing = false

  load(score: PerformanceScore): void {
    this.stopSource()
    if (!this.ctx) this.ctx = new AudioContext()
    const pcm = renderScorePcm(score, this.ctx.sampleRate)
    this.buffer = this.ctx.createBuffer(1, pcm.length, this.ctx.sampleRate)
    this.buffer.copyToChannel(pcm as Float32Array<ArrayBuffer>, 0)
    this.offset = 0
    this.playing = false
  }

  private stopSource(): void {
    if (this.source) { try { this.source.stop() } catch { /* already stopped */ } this.source = null }
  }

  async play(): Promise<void> {
    if (!this.ctx || !this.buffer) return
    await this.ctx.resume()
    this.stopSource()
    this.source = this.ctx.createBufferSource()
    this.source.buffer = this.buffer
    this.source.connect(this.ctx.destination)
    this.segmentStart = this.ctx.currentTime
    this.source.start(0, this.offset)
    this.playing = true
  }

  pause(): void {
    if (!this.ctx || !this.playing) return
    this.offset = this.position()
    this.stopSource()
    this.playing = false
  }

  position(): number {
    if (!this.ctx) return this.offset
    return this.playing ? this.offset + (this.ctx.currentTime - this.segmentStart) : this.offset
  }

  async suspendForLag(): Promise<void> { await this.ctx?.suspend() }
  async resumeFromLag(): Promise<void> { if (this.playing) await this.ctx?.resume() }
}
