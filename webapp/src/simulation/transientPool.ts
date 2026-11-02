// transientPool.ts — M6-4B/M6-4R coordinator for sharded transient-sand
// workers with an explicit synchronization protocol.
//
// Each worker replicates the plate (deterministic, particle-independent)
// and owns one contiguous shard of the SAME (count, seed) ensemble.
//
// Synchronization: work proceeds in ROUNDS. A round is one of
//   init    (batchId 0, controlRevision 0)
//   control (one control change, its own controlRevision)
//   batch   (N physics steps, its own batchId, current controlRevision)
// The coordinator starts a round only when no round is in flight, and
// delivers a batch only when every shard has reported a frame matching
// (generation, kind, batchId, controlRevision). Therefore frames from
// different generations, batches, or control states can never be mixed,
// and a control change is applied by every worker at the same step
// boundary before any further batch starts. Stale frames are dropped
// (counted in droppedStale). A worker error fails the pool loudly
// instead of waiting forever. No SharedArrayBuffer is required.
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

export type FrameKind = 'init' | 'batch' | 'control'

export interface ShardFrameMsg {
  type: string
  frame?: TransientFrame
  generation?: number
  shardIndex?: number
  batchId?: number
  controlRevision?: number
  frameKind?: FrameKind
  dtNow?: number
  message?: string
}

/** Minimal worker surface (real Worker or test fake). */
export interface WorkerLike {
  postMessage(msg: unknown): void
  terminate(): void
  onmessage: ((ev: { data: ShardFrameMsg }) => void) | null
  onerror: ((ev: unknown) => void) | null
}

interface Round { kind: FrameKind; batchId: number; controlRevision: number }

export class TransientPool {
  onBatch: ((b: PoolBatch) => void) | null = null
  onError: ((msg: string) => void) | null = null
  droppedStale = 0
  private workers: WorkerLike[] = []
  private starts: number[] = []
  private gen = 0
  private count = 0
  private timer: ReturnType<typeof setInterval> | null = null
  private round: Round | null = null
  private collected = new Map<number, TransientFrame>()
  private batchSeq = 0
  private controlRev = 0
  private pendingControls: unknown[] = []
  private failed = false
  private dtHint = 1 / (44 * 152)
  private lastTickMs = 0
  private debtSteps = 0
  private speedEma = 0
  private lastT = 0
  private lastWall = 0
  playing = false

  constructor(private factory: () => WorkerLike = () =>
    new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' }) as unknown as WorkerLike) {}

  get workerCount(): number { return this.workers.length }
  get isFailed(): boolean { return this.failed }
  get roundInFlight(): boolean { return this.round !== null }

  init(score: PerformanceScore, count: number, seed: number,
       driveXY: [number, number], workerCount: number): void {
    this.gen++
    this.count = count
    this.collected.clear()
    this.pendingControls = []
    this.batchSeq = 0
    this.controlRev = 0
    this.failed = false
    this.speedEma = 0
    this.lastT = 0
    this.lastWall = 0
    this.debtSteps = 0
    const K = Math.max(1, Math.min(workerCount, count))
    if (this.workers.length !== K) {
      for (const w of this.workers) w.terminate()
      this.workers = []
      for (let k = 0; k < K; k++) {
        const w = this.factory()
        w.onmessage = (ev) => this.handle(ev.data)
        w.onerror = () => this.fail('Workerでエラーが起きました')
        this.workers.push(w)
      }
    }
    this.starts = shardRanges(count, K).map(([a]) => a)
    this.round = { kind: 'init', batchId: 0, controlRevision: 0 }
    this.workers.forEach((w, k) => {
      w.postMessage({
        type: 'initPerformance', score, count, seed, driveXY,
        shardIndex: k, shardCount: K, generation: this.gen,
      })
    })
  }

  play(): void {
    if (this.timer || this.workers.length === 0 || this.failed) return
    this.playing = true
    this.lastTickMs = 0
    this.debtSteps = 0
    this.timer = setInterval(() => this.tick(), 33)
  }

  pause(): void {
    this.playing = false
    if (this.timer) { clearInterval(this.timer); this.timer = null }
  }

  setScore(score: PerformanceScore): void { this.control({ type: 'setScore', score }) }
  setFrequency(hz: number): void { this.control({ type: 'setFrequency', hz }) }
  setDrive(driveXY: [number, number]): void { this.control({ type: 'setDrive', driveXY }) }

  dispose(): void {
    this.pause()
    for (const w of this.workers) w.terminate()
    this.workers = []
    this.round = null
  }

  private fail(msg: string): void {
    if (this.failed) return
    this.failed = true
    this.pause()
    this.round = null
    this.onError?.(msg)
  }

  private broadcast(msg: unknown): void {
    for (const w of this.workers) w.postMessage(msg)
  }

  /** Controls form their own rounds; queue while any round is in flight. */
  private control(msg: Record<string, unknown>): void {
    this.pendingControls.push(msg)
    this.startNext()
  }

  private tick(): void {
    if (this.failed || this.workers.length === 0) return
    const now = performance.now()
    const elapsed = this.lastTickMs ? (now - this.lastTickMs) / 1000 : 0
    this.lastTickMs = now
    this.debtSteps = Math.min(
      this.debtSteps + Math.floor(elapsed / this.dtHint),
      Math.floor(0.25 / this.dtHint))
    this.startNext()
  }

  private startNext(): void {
    if (this.failed || this.round || this.workers.length === 0) return
    const ctrl = this.pendingControls.shift()
    if (ctrl) {
      const rev = this.controlRev + 1
      this.round = { kind: 'control', batchId: this.batchSeq, controlRevision: rev }
      this.collected.clear()
      this.broadcast({ ...(ctrl as object), controlRevision: rev, generation: this.gen })
      return
    }
    if (!this.playing || this.debtSteps <= 0) return
    const steps = this.debtSteps
    this.debtSteps = 0
    this.batchSeq++
    this.round = { kind: 'batch', batchId: this.batchSeq, controlRevision: this.controlRev }
    this.collected.clear()
    this.broadcast({
      type: 'stepBatch', steps, generation: this.gen,
      batchId: this.batchSeq, controlRevision: this.controlRev, speed: this.speedEma,
    })
  }

  private matchesRound(msg: ShardFrameMsg): boolean {
    const r = this.round
    if (!r) return false
    return msg.generation === this.gen
      && msg.frameKind === r.kind
      && (msg.batchId ?? 0) === r.batchId
      && (msg.controlRevision ?? 0) === r.controlRevision
  }

  private handle(msg: ShardFrameMsg): void {
    if (msg.type === 'error') { this.fail(String(msg.message ?? 'worker error')); return }
    if (msg.type !== 'frame' || !msg.frame) return
    if (this.failed) return
    if (!this.matchesRound(msg)) { this.droppedStale++; return }
    if (msg.dtNow && msg.dtNow > 0) this.dtHint = msg.dtNow
    const k = msg.shardIndex ?? 0
    if (k < 0 || k >= this.workers.length || this.collected.has(k)) { this.droppedStale++; return }
    this.collected.set(k, msg.frame)
    if (this.collected.size < this.workers.length) return
    const frames = this.workers.map((_, i) => this.collected.get(i)!)
    this.collected.clear()
    const round = this.round!
    this.round = null
    if (round.kind === 'control') this.controlRev = round.controlRevision
    const t = frames[0].t
    // All shards step identically, so times must agree; if they ever do
    // not, surface it instead of silently mixing.
    for (const f of frames) {
      if (f.t !== t) { this.fail('Worker間で物理時刻が一致しません'); return }
    }
    const now = performance.now()
    if (this.lastWall > 0 && round.kind === 'batch') {
      const wall = (now - this.lastWall) / 1000
      if (wall > 0) {
        const inst = (t - this.lastT) / wall
        this.speedEma = this.speedEma === 0 ? inst : 0.7 * this.speedEma + 0.3 * inst
      }
    }
    if (round.kind === 'batch' || round.kind === 'init') { this.lastT = t; this.lastWall = now }
    const aliveCount = frames.reduce((s, f) => s + f.aliveCount, 0)
    const stats = this.mergeStats(frames, aliveCount)
    this.onBatch?.({ t, frames, starts: [...this.starts], aliveCount, initialCount: this.count, stats })
    this.startNext()
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
