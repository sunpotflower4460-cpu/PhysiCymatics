// audio.ts — internal drive tone (Web Audio). NOT physics: a real synthesized
// sine at the drive frequency so the ear hears what the field is computed for.
// Mic/file FFT input is future work (shown as MOCK on the Music screen).
export class ToneEngine {
  private ctx: AudioContext | null = null
  private osc: OscillatorNode | null = null
  private gain: GainNode | null = null
  playing = false

  start(freqHz: number) {
    if (typeof window === 'undefined') return
    if (!this.ctx) this.ctx = new (window.AudioContext ||
      (window as unknown as { webkitAudioContext: typeof AudioContext }).webkitAudioContext)()
    if (this.ctx.state === 'suspended') void this.ctx.resume()
    this.stop()
    const osc = this.ctx.createOscillator()
    const gain = this.ctx.createGain()
    osc.type = 'sine'
    osc.frequency.value = freqHz
    gain.gain.setValueAtTime(0, this.ctx.currentTime)
    gain.gain.linearRampToValueAtTime(0.15, this.ctx.currentTime + 0.05) // soft attack
    osc.connect(gain).connect(this.ctx.destination)
    osc.start()
    this.osc = osc
    this.gain = gain
    this.playing = true
  }

  setFreq(freqHz: number) {
    if (this.osc && this.ctx) {
      this.osc.frequency.setTargetAtTime(freqHz, this.ctx.currentTime, 0.02)
    }
  }

  stop() {
    if (this.osc && this.ctx && this.gain) {
      this.gain.gain.setTargetAtTime(0, this.ctx.currentTime, 0.03)
      const osc = this.osc
      try { osc.stop(this.ctx.currentTime + 0.1) } catch { /* already stopped */ }
    }
    this.osc = null
    this.playing = false
  }
}
