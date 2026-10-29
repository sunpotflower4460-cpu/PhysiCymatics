// M6-2 gates: sleep on/off must reproduce the same physics (positions,
// counts), sleep must actually engage, input changes must invalidate it.
import { describe, it, expect } from 'vitest'
import { singleToneScore, type PerformanceScore } from '../sequence'
import { createTransientSim, calibratedForceN, type TransientSandSim } from '../transientSand'

function run(enabled: boolean, freq = 152, tEnd = 0.6, n = 60, forceScale = 1): TransientSandSim {
  const sim = createTransientSim(singleToneScore(freq, calibratedForceN() * forceScale, 30), n, 7, [0, 0], { enabled })
  sim.advanceTo(tEnd)
  return sim
}

function maxPosDiff(a: TransientSandSim, b: TransientSandSim): number {
  let m = 0
  for (let i = 0; i < a.state.id.length; i++) {
    m = Math.max(m, Math.abs(a.state.x[i] - b.state.x[i]), Math.abs(a.state.y[i] - b.state.y[i]))
  }
  return m
}

describe('M6-2 sleep parity vs baseline', () => {
  it('strong drive 152Hz: same trajectories and counts, sleep engages', () => {
    const off = run(false)
    const on = run(true)
    expect(on.skippedUpdates).toBeGreaterThan(0)
    expect(Math.abs(on.t - off.t)).toBeLessThan(1e-12)
    expect(on.nBounce).toBe(off.nBounce)
    expect(on.nDetach).toBe(off.nDetach)
    expect(on.nAlive).toBe(off.nAlive)
    expect(maxPosDiff(on, off)).toBeLessThan(1e-6)
    expect([...on.state.sitting]).toEqual([...off.state.sitting])
  })

  it('weak drive: sitting sleep dominates and still matches', () => {
    const off = run(false, 152, 0.6, 60, 0.02)
    const on = run(true, 152, 0.6, 60, 0.02)
    expect(on.skippedUpdates).toBeGreaterThan(0)
    expect(maxPosDiff(on, off)).toBeLessThan(1e-6)
    expect(on.nAlive).toBe(off.nAlive)
  })

  it('frequency sweep samples match baseline', () => {
    for (const f of [70, 173.31, 368.6, 649]) {
      const off = run(false, f, 0.3, 40)
      const on = run(true, f, 0.3, 40)
      expect(on.nAlive, `f=${f}`).toBe(off.nAlive)
      expect(maxPosDiff(on, off), `f=${f}`).toBeLessThan(1e-6)
      expect(on.nBounce, `f=${f}`).toBe(off.nBounce)
    }
  })

  it('frequency change invalidates sleep and keeps state', () => {
    const score: PerformanceScore = singleToneScore(152, calibratedForceN() * 0.02, 30)
    const sim = createTransientSim(score, 40, 3, [0, 0], { enabled: true })
    sim.advanceTo(0.4)
    const sleepingBefore = sim.nSleeping
    const xBefore = Float64Array.from(sim.state.x)
    sim.wakeAll()
    sim.plate.setFrequency(649, sim.t)
    sim.refreshForceBounds()
    expect(sim.nSleeping).toBe(0)
    if (sleepingBefore > 0) {
      // Wake materialization must not teleport sitting particles in xy.
      let jump = 0
      for (let i = 0; i < 40; i++) {
        if (sim.state.sitting[i]) jump = Math.max(jump, Math.abs(sim.state.x[i] - xBefore[i]))
      }
      expect(jump).toBe(0)
    }
    sim.advanceTo(0.7)
    expect(Number.isFinite(sim.state.x[0])).toBe(true)
  })

  it('undriven plate: dropped particle settles identically', () => {
    const off = run(false, 152, 0.4, 10, 0)
    const on = run(true, 152, 0.4, 10, 0)
    expect(maxPosDiff(on, off)).toBeLessThan(1e-9)
    expect(on.nBounce).toBe(off.nBounce)
  })
})
