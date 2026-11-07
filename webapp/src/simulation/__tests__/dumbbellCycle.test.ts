// M6-27B: dumbbell sustained contact on a resting flat plate
// (force-free TransientPlate, h=0 analytically) and a driven
// plate smoke test. docs/m6_27b_design.md.
import { describe, expect, it } from 'vitest'
import {
  DumbbellCycleSim, type DumbbellCycleParams, type DumbbellInitial,
} from '../dumbbellCycle'
import { kineticEnergyPerMass, type DumbbellState } from '../dumbbell'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'
import { G } from '../surface'

const DT = 1 / (44 * 152)
const R = 2.5e-4
const baseParams = (d: number, e = 0.3): DumbbellCycleParams => ({
  lobeRadius: R, halfSeparation: d, beta: 2 / 5,
  muS: 0.3, muK: 0.2, restitution: e,
})
const restSim = (
  initials: DumbbellInitial[],
  params: DumbbellCycleParams,
): DumbbellCycleSim => new DumbbellCycleSim(
  new TransientPlate(singleToneScore(152, 0, 30), [0, 0]), params, initials)

const run = (sim: DumbbellCycleSim, T: number, dt = DT): void => {
  for (let k = 0; k < Math.round(T / dt); k++) sim.step(dt)
}

describe('M6-27B sustained contact (resting flat plate)', () => {
  it('two-lobe simultaneous contact: horizontal dumbbell rests without drift', () => {
    const sim = restSim(
      [{ x: 0.01, y: -0.02, dropM: 0, axis: [1, 0, 0] }], baseParams(R))
    run(sim, 0.5)
    const s = sim.stats()[0]
    expect(s.lost).toBe(false)
    expect(s.failed).toBe(false)
    expect(s.unconsumedTotal).toBe(0)
    expect(Math.abs(s.x - 0.01)).toBeLessThan(1e-6)
    expect(Math.abs(s.y + 0.02)).toBeLessThan(1e-6)
    expect(Math.abs(s.czNow - R)).toBeLessThan(1e-6)
    expect(Math.abs(s.vxNow)).toBeLessThan(1e-5)
    expect(s.contactImpulses).toBeGreaterThan(0)
    expect(s.contactTime).toBeGreaterThan(0.45)
  })

  it('single-lobe support: exactly vertical dumbbell stays standing (symmetry)', () => {
    const d = R
    const sim = restSim(
      [{ x: 0, y: 0, dropM: 0, axis: [0, 0, 1] }], baseParams(d))
    run(sim, 0.2)
    const s = sim.stats()[0]
    expect(s.ez).toBeGreaterThan(0.999)
    expect(Math.abs(s.czNow - (R + d))).toBeLessThan(2e-6)
    expect(Math.abs(s.x)).toBeLessThan(1e-6)
  })

  it('axis across motion: slip-to-roll transition matches v0/(1+beta)', () => {
    const sim = restSim(
      [{ x: -0.05, y: 0, dropM: 0, vx: 0.2, axis: [0, 1, 0] }], baseParams(R, 0))
    run(sim, 1.0)
    const s = sim.stats()[0]
    const expected = 0.2 / (1 + 2 / 5)
    expect(s.vxNow).toBeGreaterThan(0.9 * expected)
    expect(s.vxNow).toBeLessThan(1.1 * expected)
    expect(Math.abs(s.spinMag - expected / R) / (expected / R)).toBeLessThan(0.1)
    expect(Math.abs(s.ey)).toBeGreaterThan(0.99)
    expect(s.x).toBeGreaterThan(-0.05 + 0.08)
    expect(s.unconsumedTotal).toBe(0)
  })

  it('contactSubsteps: explicit default is identical and invalid counts fail', () => {
    const initial: DumbbellInitial = {
      x: 0, y: 0, dropM: 0.001, axis: [Math.SQRT1_2, 0, Math.SQRT1_2],
    }
    const baseline = restSim([initial], baseParams(R, 0))
    const explicit = restSim([initial], { ...baseParams(R, 0), contactSubsteps: 4 })
    run(baseline, 0.1)
    run(explicit, 0.1)
    expect(explicit.stats()).toEqual(baseline.stats())
    for (const contactSubsteps of [0, -1, 1.5, Number.NaN]) {
      expect(() => restSim([initial], { ...baseParams(R, 0), contactSubsteps }))
        .toThrow(RangeError)
    }
    const refined = restSim([initial], { ...baseParams(R, 0), contactSubsteps: 8 })
    run(refined, 0.1)
    expect(refined.stats()[0].failed).toBe(false)
    expect(refined.stats()[0].contactTime).toBeGreaterThan(0)
  })

  it('single-lobe pivot: 45-degree dumbbell topples to two-lobe rest without energy gain', () => {
    const p = baseParams(R, 0)
    const sim = restSim([{
      x: 0, y: 0, dropM: 0,
      axis: [Math.SQRT1_2, 0, Math.SQRT1_2],
    }], p)
    sim.step(DT)
    const grains = (sim as unknown as { grains: { st: DumbbellState }[] }).grains
    const mech = (): number => kineticEnergyPerMass(grains[0].st, p) + G * grains[0].st.cz
    const E0 = mech()
    run(sim, 1.5)
    const s = sim.stats()[0]
    expect(Math.abs(s.ez)).toBeLessThan(0.25) // lying down
    expect(Math.abs(s.czNow - R)).toBeLessThan(0.25 * R)
    expect(mech()).toBeLessThanOrEqual(E0 + 2e-4)
    expect(s.contactImpulses).toBeGreaterThan(0)
    expect(s.failed).toBe(false)
  })
})

describe('M6-27B driven plate smoke', () => {
  it('152Hz eccentric: impacts occur, decomposition telescopes, no dropped time', () => {
    const plate = new TransientPlate(
      singleToneScore(152, calibratedForceN(), 30), [0.06, 0])
    const sim = new DumbbellCycleSim(plate, {
      ...baseParams(R / 2), diagFreqHz: 152, plateHalfThickness: 5e-4,
    }, [{ x: 0.03, y: 0, dropM: 0.003, axis: [0.6, 0.48, 0.64] }])
    run(sim, 0.3)
    const s = sim.stats()[0]
    expect(Number.isFinite(s.x) && Number.isFinite(s.y)).toBe(true)
    expect(s.unconsumedTotal).toBe(0)
    expect(s.failed).toBe(false)
    expect(s.impacts).toBeGreaterThan(0)
    expect(s.netFlightX + s.netContactX + s.netCorrectionX).toBeCloseTo(s.x - 0.03, 7)
    expect(s.netFlightY + s.netContactY + s.netCorrectionY).toBeCloseTo(s.y - 0, 7)
  })
})
