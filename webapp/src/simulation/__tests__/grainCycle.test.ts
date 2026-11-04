// M6-19: grain cycle harness validation on a resting flat plate.
import { describe, expect, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { G } from '../surface'

const base = { muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0 }
const DT = 1 / (44 * 152)

function run(model: 'roll' | 'point', vx: number, seconds: number) {
  const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
  const sim = new GrainCycleSim(plate, { ...base, model } as GrainCycleParams,
    [{ x: -0.05, y: 0, dropM: 0.002, vx }])
  const n = Math.round(seconds / DT)
  for (let k = 0; k < n; k++) sim.step(DT)
  return sim.stats()[0]
}

describe('M6-19 grain cycle harness', () => {
  it('vertical drop: one impact, capture, no horizontal drift', () => {
    for (const model of ['roll', 'point'] as const) {
      const st = run(model, 0, 0.3)
      expect(st.impacts).toBe(1)
      expect(st.contactTime).toBeGreaterThan(0.25)
      expect(Math.abs(st.x - (-0.05))).toBeLessThan(1e-9)
      expect(st.lost).toBe(false)
    }
  })

  it('rolling persists where sliding stops (flat, horizontal throw)', () => {
    const roll = run('roll', 0.2, 0.6)
    const point = run('point', 0.2, 0.6)
    const dxRoll = roll.x - (-0.05)
    const dxPoint = point.x - (-0.05)
    // Point grain: friction stops it (finite stopping distance).
    // Sphere: transitions to rolling at vx/(1+beta) and keeps going.
    expect(dxRoll).toBeGreaterThan(3 * dxPoint)
    expect(dxPoint).toBeGreaterThan(0)
    expect(dxPoint).toBeLessThan(0.02)
    // Analytic rolling speed after transition: 0.2/1.4 = 0.1429;
    // distance in ~0.55 s of contact ~ 0.075 m.
    expect(dxRoll).toBeGreaterThan(0.05)
    void G
  })
})

describe('M6-20 transport decomposition conservation', () => {
  it('displacement and momentum decompositions telescope', () => {
    for (const model of ['roll', 'point'] as const) {
      const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
      const sim = new GrainCycleSim(plate, { ...base, model } as GrainCycleParams,
        [{ x: -0.03, y: 0.01, dropM: 0.003, vx: 0.15 }])
      for (let k = 0; k < Math.round(0.4 / DT); k++) sim.step(DT)
      const st = sim.stats()[0]
      // Displacement: flight + contact = total (center-based).
      expect(st.netFlightX + st.netContactX).toBeCloseTo(st.x - (-0.03), 9)
      expect(st.netFlightY + st.netContactY).toBeCloseTo(st.y - 0.01, 9)
      // Momentum: initial + impact + contact + conversion = final.
      expect(0.15 + st.impactImpulseX + st.contactImpulseX + st.conversionImpulseX)
        .toBeCloseTo(st.vxNow, 6)
      expect(0 + st.impactImpulseY + st.contactImpulseY + st.conversionImpulseY)
        .toBeCloseTo(st.vyNow, 6)
      expect(st.pathFlightM + st.pathContactM).toBeGreaterThan(0)
    }
  })
})
