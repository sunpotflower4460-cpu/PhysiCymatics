// M6-12R: observedState must report sleeping flyers at the current
// sim time (not their stale sleepT0 arrays), must not mutate state,
// and must make sleep on/off measurements agree.
import { describe, expect, it } from 'vitest'
import { singleToneScore } from '../sequence'
import { createTransientSim, calibratedForceN, TransientSandSim } from '../transientSand'

function run(sleep: boolean, steps: number): TransientSandSim {
  const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 120, 42, [0, 0], { enabled: sleep })
  for (let s = 0; s < steps; s++) sim.step()
  return sim
}

describe('M6-12R observedState', () => {
  it('sleeping flyers are observed at sim time (replay-consistent, sleep-off-consistent)', () => {
    const steps1s = Math.round(1 / (1 / (44 * 152)))
    const on = run(true, steps1s)
    const o1 = new Float64Array(6)
    let sleepingFlyers = 0; let staleWouldDiffer = 0; let replayDiff = 0
    for (let i = 0; i < on.state.id.length; i++) {
      if (!on.state.alive[i] || !on.sleeping[i] || on.state.sitting[i]) continue
      sleepingFlyers++
      on.observedState(i, o1)
      // Independent replay of the stored ballistic state to sim.t.
      const d = on.t - on.sleepT0[i]
      const rx = on.state.x[i] + on.state.vx[i] * d
      const rz = on.state.z[i] + on.state.vz[i] * d - 0.5 * 9.81 * d * d
      replayDiff = Math.max(replayDiff, Math.abs(o1[0] - rx), Math.abs(o1[2] - rz))
      staleWouldDiffer = Math.max(staleWouldDiffer, Math.abs(on.state.z[i] - rz))
    }
    expect(sleepingFlyers).toBeGreaterThan(50) // reviewer's config: 201 at 300 grains
    expect(staleWouldDiffer).toBeGreaterThan(1e-4) // raw-array read really was wrong
    expect(replayDiff).toBeLessThan(1e-12)
    // Sleep on/off observed positions agree before chaotic amplification
    // separates the trajectories (0.25 s; established sleep-parity scale).
    const stepsQ = Math.round(0.25 / (1 / (44 * 152)))
    const a = run(true, stepsQ); const b = run(false, stepsQ)
    const o2 = new Float64Array(6)
    let maxPosDiff = 0
    for (let i = 0; i < a.state.id.length; i++) {
      if (!a.state.alive[i] || !b.state.alive[i]) continue
      a.observedState(i, o1); b.observedState(i, o2)
      maxPosDiff = Math.max(maxPosDiff, Math.abs(o1[0] - o2[0]), Math.abs(o1[1] - o2[1]), Math.abs(o1[2] - o2[2]))
    }
    expect(maxPosDiff).toBeLessThan(5e-6)
  }, 300_000)

  it('observation does not mutate simulation state', () => {
    const sim = run(true, 2000)
    const before = {
      x: Float64Array.from(sim.state.x), z: Float64Array.from(sim.state.z),
      vz: Float64Array.from(sim.state.vz), sleepT0: Float64Array.from(sim.sleepT0),
      sleeping: Uint8Array.from(sim.sleeping), t: sim.t,
    }
    const o = new Float64Array(6)
    for (let i = 0; i < sim.state.id.length; i++) sim.observedState(i, o)
    expect(sim.t).toBe(before.t)
    for (let i = 0; i < sim.state.id.length; i++) {
      expect(sim.state.x[i]).toBe(before.x[i])
      expect(sim.state.z[i]).toBe(before.z[i])
      expect(sim.state.vz[i]).toBe(before.vz[i])
      expect(sim.sleepT0[i]).toBe(before.sleepT0[i])
      expect(sim.sleeping[i]).toBe(before.sleeping[i])
    }
    // And the run continues identically to an unobserved twin.
    const twin = run(true, 2000)
    for (let s = 0; s < 500; s++) { sim.step(); twin.step() }
    for (let i = 0; i < sim.state.id.length; i++) {
      expect(sim.state.x[i]).toBe(twin.state.x[i])
      expect(sim.state.z[i]).toBe(twin.state.z[i])
    }
    expect(sim.nBounce).toBe(twin.nBounce)
  }, 300_000)

  it('chord path length from observed positions agrees sleep on/off', () => {
    const pathOf = (sleep: boolean): number => {
      const sim = createTransientSim(singleToneScore(152, calibratedForceN(), 30), 60, 42, [0, 0], { enabled: sleep })
      const o = new Float64Array(6)
      const prev = new Float64Array(60 * 3)
      const path = new Float64Array(60)
      for (let i = 0; i < 60; i++) { sim.observedState(i, o); prev.set([o[0], o[1], o[2]], i * 3) }
      const steps = Math.round(0.5 / (1 / (44 * 152)))
      for (let s = 0; s < steps; s++) {
        sim.step()
        for (let i = 0; i < 60; i++) {
          if (!sim.state.alive[i]) continue
          sim.observedState(i, o)
          path[i] += Math.hypot(o[0] - prev[i * 3], o[1] - prev[i * 3 + 1], o[2] - prev[i * 3 + 2])
          prev.set([o[0], o[1], o[2]], i * 3)
        }
      }
      let m = 0; let n = 0
      for (let i = 0; i < 60; i++) if (sim.state.alive[i]) { m += path[i]; n++ }
      return m / n
    }
    const pOn = pathOf(true); const pOff = pathOf(false)
    expect(pOff).toBeGreaterThan(0)
    expect(Math.abs(pOn - pOff) / pOff).toBeLessThan(1e-6)
  }, 300_000)
})
