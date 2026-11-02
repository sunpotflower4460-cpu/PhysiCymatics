// M6-4R: TransientPool synchronization protocol tests with scriptable
// fake workers. The fakes speak the worker protocol (init/batch/control
// frames echoing generation/batchId/controlRevision) but let the test
// deliver frames in any order, so races are reproducible deterministically.
// Physics itself is covered by the shard-equivalence test in sleep.test.ts.
import { describe, it, expect } from 'vitest'
import { TransientPool, type ShardFrameMsg, type WorkerLike, type PoolBatch } from '../transientPool'
import { singleToneScore } from '../sequence'
import type { TransientFrame } from '../transientSand'

function fakeFrame(t: number, alive = 10): TransientFrame {
  return {
    t, full: false, idx: new Int32Array(0), xy: new Float32Array(0), dead: new Int32Array(0),
    aliveCount: alive, initialCount: alive,
    stats: {
      t, freqHz: 152, ampUm: 0, alive, sitting: 0, nBounce: 0,
      stepsPerSec: 0, sleepingFrac: 0, skippedUpdates: 0,
    },
  } as unknown as TransientFrame
}

class FakeWorker implements WorkerLike {
  onmessage: ((ev: { data: ShardFrameMsg }) => void) | null = null
  onerror: ((ev: unknown) => void) | null = null
  outbox: ShardFrameMsg[] = []
  log: string[] = []
  terminated = false
  private t = 0
  private gen = -1
  private batchId = 0
  private rev = 0
  private shard = 0
  private static dt = 1 / (44 * 152)

  postMessage(msg: unknown): void {
    const m = msg as Record<string, unknown>
    if (m.type === 'initPerformance') {
      this.gen = m.generation as number
      this.shard = m.shardIndex as number
      this.batchId = 0; this.rev = 0; this.t = 0
      this.log.push('init')
      this.emit('init')
    } else if (m.type === 'stepBatch') {
      if (m.generation !== this.gen || m.controlRevision !== this.rev) {
        this.log.push(`batch-REJECTED rev=${String(m.controlRevision)} mine=${this.rev}`)
        return
      }
      this.batchId = m.batchId as number
      this.t += (m.steps as number) * FakeWorker.dt
      this.log.push(`batch${this.batchId}@rev${this.rev}`)
      this.emit('batch')
    } else if (m.type === 'setFrequency' || m.type === 'setScore' || m.type === 'setDrive') {
      this.rev = m.controlRevision as number
      this.log.push(`ctrl${this.rev}`)
      this.emit('control')
    }
  }

  private emit(frameKind: 'init' | 'batch' | 'control'): void {
    this.outbox.push({
      type: 'frame', frame: fakeFrame(this.t), generation: this.gen,
      shardIndex: this.shard, batchId: this.batchId, controlRevision: this.rev,
      frameKind, dtNow: FakeWorker.dt,
    })
  }

  /** Deliver the oldest queued frame to the pool. */
  flushOne(): void {
    const m = this.outbox.shift()
    if (m) this.onmessage?.({ data: m })
  }

  flushAll(): void { while (this.outbox.length) this.flushOne() }
  terminate(): void { this.terminated = true }
}

function makePool(K: number) {
  const fakes: FakeWorker[] = []
  const pool = new TransientPool(() => {
    const f = new FakeWorker()
    fakes.push(f)
    return f
  })
  const batches: PoolBatch[] = []
  const errors: string[] = []
  pool.onBatch = (b) => batches.push(b)
  pool.onError = (m) => errors.push(m)
  pool.init(singleToneScore(152, 0.05, 30), 100, 7, [0, 0], K)
  return { pool, fakes, batches, errors }
}

const flushInit = (fakes: FakeWorker[]) => fakes.forEach((f) => f.flushAll())
const startBatch = (pool: TransientPool, steps = 10) => {
  const p = pool as unknown as { playing: boolean; debtSteps: number; startNext: () => void }
  p.playing = true
  p.debtSteps = steps
  p.startNext()
}

describe('M6-4R pool synchronization', () => {
  it('init round completes only when every shard reports, any order', () => {
    const { pool, fakes, batches } = makePool(3)
    fakes[2].flushOne(); fakes[0].flushOne()
    expect(batches.length).toBe(0)
    fakes[1].flushOne()
    expect(batches.length).toBe(1)
    expect(batches[0].t).toBe(0)
    expect(pool.roundInFlight).toBe(false)
  })

  it('control during a batch is applied at the boundary, before the next batch', () => {
    const { pool, fakes, batches } = makePool(2)
    flushInit(fakes)
    startBatch(pool)
    expect(pool.roundInFlight).toBe(true)
    pool.setFrequency(649)
    // Control must NOT reach workers while the batch is in flight.
    for (const f of fakes) expect(f.log).toEqual(['init', 'batch1@rev0'])
    fakes.forEach((f) => f.flushAll()) // batch frames -> delivered, control sent
    expect(batches.length).toBe(2) // init + batch
    for (const f of fakes) expect(f.log).toEqual(['init', 'batch1@rev0', 'ctrl1'])
    fakes.forEach((f) => f.flushAll()) // control acks
    expect(batches.length).toBe(3)
    startBatch(pool)
    for (const f of fakes) expect(f.log[f.log.length - 1]).toBe('batch2@rev1')
    fakes.forEach((f) => f.flushAll())
    expect(batches.length).toBe(4)
    // All delivered batches have one consistent time across shards.
    for (const b of batches) for (const fr of b.frames) expect(fr.t).toBe(b.t)
  })

  it('rapid successive controls apply as sequential revisions in order', () => {
    const { pool, fakes, batches } = makePool(2)
    flushInit(fakes)
    pool.setFrequency(300)
    pool.setFrequency(400)
    pool.setDrive([0.01, 0])
    for (let round = 0; round < 3; round++) fakes.forEach((f) => f.flushAll())
    for (const f of fakes) expect(f.log).toEqual(['init', 'ctrl1', 'ctrl2', 'ctrl3'])
    expect(batches.length).toBe(4) // init + 3 control rounds
  })

  it('stale-generation frames after re-init are dropped', () => {
    const { pool, fakes, batches } = makePool(2)
    flushInit(fakes)
    startBatch(pool)
    const stale = fakes.map((f) => f.outbox.shift()!) // hold batch frames back
    pool.init(singleToneScore(152, 0.05, 30), 100, 9, [0, 0], 2) // re-init (new seed)
    // Deliver the held-back old-generation batch frames now.
    const before = pool.droppedStale
    for (const m of stale) fakes[0].onmessage?.({ data: m })
    expect(pool.droppedStale).toBeGreaterThan(before)
    expect(batches.length).toBe(1) // only the first init batch
    flushInit(fakes)
    expect(batches.length).toBe(2) // new init batch
  })

  it('mismatched batchId frames are dropped, never mixed', () => {
    const { pool, fakes, batches } = makePool(2)
    flushInit(fakes)
    startBatch(pool)
    const bogus: ShardFrameMsg = {
      type: 'frame', frame: fakeFrame(123), generation: 1, shardIndex: 0,
      batchId: 999, controlRevision: 0, frameKind: 'batch',
    }
    fakes[0].onmessage?.({ data: bogus })
    expect(pool.droppedStale).toBe(1)
    expect(batches.length).toBe(1)
    fakes.forEach((f) => f.flushAll())
    expect(batches.length).toBe(2)
    expect(batches[1].t).toBeGreaterThan(0)
    expect(batches[1].t).toBeLessThan(1)
  })

  it('worker error fails loudly instead of hanging', () => {
    const { pool, fakes, errors } = makePool(2)
    flushInit(fakes)
    startBatch(pool)
    fakes[1].onerror?.(new Error('boom'))
    expect(errors.length).toBe(1)
    expect(pool.isFailed).toBe(true)
    expect(pool.roundInFlight).toBe(false)
    // Further frames are ignored, no throw, no delivery.
    fakes.forEach((f) => f.flushAll())
  })

  it('pause during a batch still delivers that batch, then stops', () => {
    const { pool, fakes, batches } = makePool(2)
    flushInit(fakes)
    startBatch(pool)
    pool.pause()
    fakes.forEach((f) => f.flushAll())
    expect(batches.length).toBe(2)
    expect(pool.roundInFlight).toBe(false)
    const logs = fakes.map((f) => f.log.length)
    const p = pool as unknown as { debtSteps: number; startNext: () => void }
    p.debtSteps = 10
    p.startNext() // playing=false: must not start a batch
    expect(fakes.map((f) => f.log.length)).toEqual(logs)
  })
})
