// worker.ts — sand physics off the UI thread.
// Fixed physics dt (T/44) regardless of display rate; each tick advances as
// many whole steps as fit the time budget and reports the measured speed.
// Snapshots are fresh Float32 copies (transferable); the Float64 state stays here.
import { createSim, type SandSim } from './sand'
import { createTransientSim, createTransientSimShard, type TransientSandSim } from './transientSand'
import type { PerformanceScore } from './sequence'

export type SandWorkerIn =
  | { type: 'init'; gamma: number; count: number; seed: number }
  | { type: 'initPerformance'; score: PerformanceScore; count: number; seed: number; driveXY: [number, number]; shardIndex?: number; shardCount?: number; generation?: number }
  | { type: 'stepBatch'; steps: number; generation: number; speed: number }
  | { type: 'setScore'; score: PerformanceScore }
  | { type: 'setFrequency'; hz: number }
  | { type: 'setDrive'; driveXY: [number, number] }
  | { type: 'play' }
  | { type: 'pause' }

let sim: SandSim | null = null
let tsim: TransientSandSim | null = null
let timer: ReturnType<typeof setInterval> | null = null
let speedEma = 0
let lastTickMs = 0
let debtSteps = 0
// M6-4B pool mode: when generation is set, an external coordinator drives
// step batches; the internal timer stays off.
let generation = -1
let shardIndex = 0

function postSnapshot(): void {
  if (tsim) {
    // M6-3C: differential frame (unchanged display xy is not resent).
    const frame = tsim.frameDelta(speedEma)
    ;(self as unknown as Worker).postMessage(
      { type: 'frame', kind: 'transient', frame, generation, shardIndex, dtNow: tsim.dtNow() },
      [frame.xy.buffer, frame.idx.buffer, frame.dead.buffer] as unknown as Transferable[])
    return
  }
  if (!sim) return
  const snap = sim.snapshot(speedEma)
  ;(self as unknown as Worker).postMessage(
    { type: 'snapshot', kind: 'single', snap }, [snap.xy.buffer] as unknown as Transferable[])
}

function tick(): void {
  const active = tsim ?? sim
  if (!active) return
  // Realtime pacing: physics may lag behind wall time but never runs ahead.
  // Unpayable debt is dropped (physics simply progresses slower; the UI
  // shows the measured speed) instead of spiralling.
  const now = performance.now()
  const elapsed = lastTickMs ? (now - lastTickMs) / 1000 : 0
  lastTickMs = now
  const tBefore = active.t
  const refDt = tsim ? tsim.dtNow() : sim!.dt
  debtSteps = Math.min(debtSteps + Math.floor(elapsed / refDt),
                       Math.floor(0.25 / refDt))
  const budgetMs = 24
  const t0 = performance.now()
  let steps = 0
  while (debtSteps > 0 && performance.now() - t0 < budgetMs) {
    active.step(); steps++; debtSteps--
  }
  const inst = elapsed > 0 ? (active.t - tBefore) / elapsed : 0
  speedEma = speedEma === 0 ? inst : 0.7 * speedEma + 0.3 * inst
  postSnapshot()
}

self.onmessage = (ev: MessageEvent<SandWorkerIn>) => {
  try {
    const msg = ev.data
    if (msg.type === 'init') {
      if (timer) { clearInterval(timer); timer = null }
      tsim = null
      sim = createSim(msg.gamma, msg.count, msg.seed)
      speedEma = 0
      postSnapshot()
    } else if (msg.type === 'initPerformance') {
      if (timer) { clearInterval(timer); timer = null }
      sim = null
      generation = msg.generation ?? -1
      shardIndex = msg.shardIndex ?? 0
      tsim = msg.shardCount && msg.shardCount > 1
        ? createTransientSimShard(msg.score, msg.count, msg.seed, msg.driveXY, { enabled: true }, shardIndex, msg.shardCount)
        : createTransientSim(msg.score, msg.count, msg.seed, msg.driveXY)
      speedEma = 0
      postSnapshot()
    } else if (msg.type === 'stepBatch') {
      if (tsim && msg.generation === generation) {
        speedEma = msg.speed
        for (let s = 0; s < msg.steps; s++) tsim.step()
        postSnapshot()
      }
    } else if (msg.type === 'setScore') {
      if (tsim) { tsim.wakeAll(); tsim.plate.setScore(msg.score, tsim.t); tsim.refreshForceBounds(); postSnapshot() }
    } else if (msg.type === 'setFrequency') {
      if (tsim) { tsim.wakeAll(); tsim.plate.setFrequency(msg.hz, tsim.t); tsim.refreshForceBounds(); postSnapshot() }
    } else if (msg.type === 'setDrive') {
      if (tsim) { tsim.wakeAll(); tsim.plate.setDrive(msg.driveXY[0], msg.driveXY[1]); tsim.refreshForceBounds(); postSnapshot() }
    } else if (msg.type === 'play') {
      if (generation >= 0) { /* pool-driven: coordinator sends stepBatch */ }
      else if ((sim || tsim) && !timer) { lastTickMs = 0; debtSteps = 0; timer = setInterval(tick, 33) }
    } else if (msg.type === 'pause') {
      if (timer) { clearInterval(timer); timer = null }
      postSnapshot()
    }
  } catch (err) {
    ;(self as unknown as Worker).postMessage(
      { type: 'error', message: err instanceof Error ? err.message : String(err) })
  }
}
