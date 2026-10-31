// transientPool.ts — M6-4B coordinator for sharded transient-sand workers.
//
// Each worker replicates the plate (deterministic, particle-independent)
// and owns one contiguous shard of the SAME (count, seed) ensemble.
// The coordinator sends identical step batches to every worker and only
// publishes a frame batch when all shards have reported the same physics
// time, so the display never mixes times. Control messages are broadcast
// between batches and therefore apply at the same step boundary in every
// worker. A generation counter makes stale frames from before a re-init
// harmless. No SharedArrayBuffer is required.
import { shardRanges, type TransientFrame, type TransientStats } from './transientSand'
import type { PerformanceScore } from './sequence'

export interface PoolBatch {
  t: number
  frames: TransientFrame[] // per shard, in shard order
  starts: number[] // global start index per shard
  aliveCount: number
  initialCount: number
  stats: TransientStats
}

interface ShardMsg {
  type: string
  frame?: TransientFrame
  generation?: number
  shardIndex?: number
  dtNow?: number
  message?: string
}

export class TransientPool {
  onBatch: ((b: PoolBatch) => void) | null = null
  onError: ((msg: string) => void) | null = null
  private workers: Worker[] = []
  private starts: number[] = []
  private gen = 0
  private count = 0
  private timer: ReturnType<typeof setInterval> | null = null
  private busy = false
  private collected = new Map<number, TransientFrame>()
  private dtHint = 1 / (44 * 152)
  private lastTickMs = 0
  private debtSteps = 0
  private speedEma = 0
  private lastT = 0
  private lastWall = 0
  playing = false

  get workerCount(): number { return this.workers.length }

  init(score: PerformanceScore, count: number, seed: number,
       driveXY: [number, number], workerCount: number): void {
    this.gen++
    this.count = count
    this.collected.clear()
    this.pendingControls = []
    this.busy = false
    this.speedEma = 0
    this.lastT = 0
    this.debtSteps = 0
    const K = Math.max(1, Math.min(workerCount, count))
    if (this.workers.length !== K) {
      for (const w of this.workers) w.terminate()
      this.workers = []
      for (let k = 0; k < K; k++) {
        const w = new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' })
        w.onmessage = (ev: MessageEvent<ShardMsg>) => this.handle(ev.data)
        w.onerror = () => this.onError?.('Workerでエラーが起きました')
        this.workers.push(w)
      }
    }
    this.starts = shardRanges(count, K).map(([a]) => a)
    this.workers.forEach((w, k) => {
      w.postMessage({
        type: 'initPerformance', score, count, seed, driveXY,
        shardIndex: k, shardCount: K, generation: this.gen,
      })
    })
  }

  play(): void {
    if (this.timer || this.workers.length === 0) return
    this.playing = true
    this.lastTickMs = 0
    this.debtSteps = 0
    this.timer = setInterval(() => this.tick(), 33)
  }

  pause(): void {
    this.playing = false
    if (this.timer) { clearInterval(this.timer); this.timer = null }
  }

  private pendingControls: unknown[] = []

  setScore(score: PerformanceScore): void { this.control({ type: 'setScore', score }) }
  setFrequency(hz: number): void { this.control({ type: 'setFrequency', hz }) }
  setDrive(driveXY: [number, number]): void { this.control({ type: 'setDrive', driveXY }) }

  /** Controls apply at a batch boundary: queue while a batch is in flight
   *  so control-response frames never interleave with batch frames. */
  private control(msg: unknown): void {
    if (this.busy) this.pendingControls.push(msg)
    else this.broadcast(msg)
  }

  dispose(): void {
    this.pause()
    for (const w of this.workers) w.terminate()
    this.workers = []
  }

  private broadcast(msg: unknown): void {
    for (const w of this.workers) w.postMessage(msg)
  }

  private tick(): void {
    if (this.busy || this.workers.length === 0) return
    const now = performance.now()
    const elapsed = this.lastTickMs ? (now - this.lastTickMs) / 1000 : 0
    this.lastTickMs = now
    this.debtSteps = Math.min(
      this.debtSteps + Math.floor(elapsed / this.dtHint),
      Math.floor(0.25 / this.dtHint))
    if (this.debtSteps <= 0) return
    const steps = this.debtSteps
    this.debtSteps = 0
    this.busy = true
    this.collected.clear()
    this.broadcast({ type: 'stepBatch', steps, generation: this.gen, speed: this.speedEma })
  }

  private handle(msg: ShardMsg): void {
    if (msg.type === 'error') { this.onError?.(String(msg.message ?? 'worker error')); return }
    if (msg.type !== 'frame' || !msg.frame) return
    if (msg.generation !== this.gen) return // stale: from before a re-init
    if (msg.dtNow && msg.dtNow > 0) this.dtHint = msg.dtNow
    const k = msg.shardIndex ?? 0
    if (this.collected.has(k)) {
      // Should not happen with boundary-queued controls; if it does, a
      // newer frame starts a fresh collection rather than mixing times.
      const prev = this.collected.get(k)!
      if (msg.frame.t <= prev.t) return
      this.collected.clear()
    }
    this.collected.set(k, msg.frame)
    if (this.collected.size < this.workers.length) return
    const frames = this.workers.map((_, i) => this.collected.get(i)!)
    this.collected.clear()
    this.busy = false
    const t = frames[0].t
    const now = performance.now()
    if (this.lastWall > 0) {
      const wall = (now - this.lastWall) / 1000
      if (wall > 0) {
        const inst = (t - this.lastT) / wall
        this.speedEma = this.speedEma === 0 ? inst : 0.7 * this.speedEma + 0.3 * inst
      }
    }
    this.lastT = t; this.lastWall = now
    const aliveCount = frames.reduce((s, f) => s + f.aliveCount, 0)
    const stats = this.mergeStats(frames, aliveCount)
    this.onBatch?.({ t, frames, starts: [...this.starts], aliveCount, initialCount: this.count, stats })
    if (this.pendingControls.length > 0) {
      const pending = this.pendingControls
      this.pendingControls = []
      for (const m of pending) this.broadcast(m)
    }
  }

  private mergeStats(frames: TransientFrame[], aliveCount: number): TransientStats {
    const s0 = frames[0].stats
    let sleeping = 0
    let sittingWeighted = 0
    let nBounce = 0
    let skipped = 0
    for (const f of frames) {
      sleeping += f.stats.sleepingFrac * f.aliveCount
      if (!Number.isNaN(f.stats.sitting)) sittingWeighted += f.stats.sitting * f.aliveCount
      nBounce += f.stats.nBounce
      skipped += f.stats.skippedUpdates
    }
    return {
      ...s0,
      alive: aliveCount,
      sitting: aliveCount ? sittingWeighted / aliveCount : NaN,
      nBounce,
      stepsPerSec: this.speedEma,
      sleepingFrac: aliveCount ? sleeping / aliveCount : 0,
      skippedUpdates: skipped,
    }
  }
}
