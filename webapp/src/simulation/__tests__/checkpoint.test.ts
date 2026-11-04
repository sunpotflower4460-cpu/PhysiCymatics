// M6-12: step-boundary checkpoint save/restore must reproduce the
// uninterrupted run bitwise (same construction, continued with step()).
import { describe, expect, it } from 'vitest'
import { singleToneScore } from '../sequence'
import {
  TransientSandSim, createTransientSim, calibratedForceN,
  checkpointToBytes, checkpointFromBytes,
} from '../transientSand'

function makeSim(freq: number, n: number, seed: number, mode: 'hermite' | 'ode', sleep = true): TransientSandSim {
  const sim = createTransientSim(singleToneScore(freq, calibratedForceN(), 30), n, seed, [0, 0], { enabled: sleep })
  if (mode === 'ode') sim.plate.accelMode = 'ode'
  return sim
}

function expectBitwise(a: TransientSandSim, b: TransientSandSim): void {
  expect(b.t).toBe(a.t)
  expect(b.nBounce).toBe(a.nBounce)
  expect(b.nDetach).toBe(a.nDetach)
  expect(b.nEventCap).toBe(a.nEventCap)
  const sa = a.state; const sb = b.state
  for (const k of ['x', 'y', 'z', 'vx', 'vy', 'vz'] as const) {
    for (let i = 0; i < sa.id.length; i++) {
      if (!Object.is(sa[k][i], sb[k][i])) {
        throw new Error(`${k}[${i}] differs: ${sa[k][i]} vs ${sb[k][i]}`)
      }
    }
  }
  for (let i = 0; i < sa.id.length; i++) {
    expect(sb.alive[i]).toBe(sa.alive[i])
    expect(sb.sitting[i]).toBe(sa.sitting[i])
    expect(b.sleeping[i]).toBe(a.sleeping[i])
    if (!Object.is(b.wakeT[i], a.wakeT[i])) throw new Error(`wakeT[${i}] differs`)
  }
  for (let j = 0; j < a.plate.q.length; j++) {
    if (!Object.is(a.plate.q[j], b.plate.q[j])) throw new Error(`plate.q[${j}] differs`)
    if (!Object.is(a.plate.qd[j], b.plate.qd[j])) throw new Error(`plate.qd[${j}] differs`)
  }
}

describe('M6-12 checkpoint save/restore', () => {
  for (const mode of ['hermite', 'ode'] as const) {
    it(`bitwise resume after restore into a fresh sim (${mode}, sleep on)`, () => {
      const stepsBefore = 400; const stepsAfter = 400
      const ref = makeSim(152, 120, 11, mode)
      for (let i = 0; i < stepsBefore + stepsAfter; i++) ref.step()

      const src = makeSim(152, 120, 11, mode)
      for (let i = 0; i < stepsBefore; i++) src.step()
      const cp = src.checkpoint()

      const dst = makeSim(152, 120, 11, mode)
      dst.restore(cp)
      for (let i = 0; i < stepsAfter; i++) dst.step()
      expectBitwise(ref, dst)
    }, 300_000)
  }

  it('bitwise resume with sleep disabled (legacy full-scan path)', () => {
    const ref = makeSim(649, 80, 3, 'hermite', false)
    for (let i = 0; i < 600; i++) ref.step()
    const src = makeSim(649, 80, 3, 'hermite', false)
    for (let i = 0; i < 250; i++) src.step()
    const dst = makeSim(649, 80, 3, 'hermite', false)
    dst.restore(src.checkpoint())
    for (let i = 0; i < 350; i++) dst.step()
    expectBitwise(ref, dst)
  }, 300_000)

  it('bytes round-trip preserves bitwise resume', () => {
    const ref = makeSim(152, 60, 5, 'ode')
    for (let i = 0; i < 500; i++) ref.step()
    const src = makeSim(152, 60, 5, 'ode')
    for (let i = 0; i < 200; i++) src.step()
    const cp = checkpointFromBytes(checkpointToBytes(src.checkpoint()))
    const dst = makeSim(152, 60, 5, 'ode')
    dst.restore(cp)
    for (let i = 0; i < 300; i++) dst.step()
    expectBitwise(ref, dst)
  }, 300_000)

  it('checkpoint before the first step restores the initial state', () => {
    const ref = makeSim(152, 40, 9, 'hermite')
    for (let i = 0; i < 300; i++) ref.step()
    const src = makeSim(152, 40, 9, 'hermite')
    const dst = makeSim(152, 40, 9, 'hermite')
    dst.restore(src.checkpoint())
    for (let i = 0; i < 300; i++) dst.step()
    expectBitwise(ref, dst)
  }, 300_000)

  it('restore rejects a different score and particle count', () => {
    const src = makeSim(152, 40, 9, 'hermite')
    src.step()
    const cp = src.checkpoint()
    const wrongScore = makeSim(649, 40, 9, 'hermite')
    expect(() => wrongScore.restore(cp)).toThrow(/score mismatch/)
    const wrongCount = makeSim(152, 41, 9, 'hermite')
    expect(() => wrongCount.restore(cp)).toThrow(/count mismatch/)
  }, 60_000)
})
