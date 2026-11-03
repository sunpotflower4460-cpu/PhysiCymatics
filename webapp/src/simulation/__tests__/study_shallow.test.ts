// M6-8C shallow-contact convergence study. Skipped unless PHYSI_BENCH=1.
// Compares the production time step (spp=44) against finer steps
// (spp=88/176): with finer steps the 8 contact samples are denser in
// time, so systematic bounce-count differences indicate shallow
// contacts the production sampling can miss between samples.
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientSandSim, createTransientSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

benchIt('shallow contact convergence (spp scaling)', () => {
  const F = calibratedForceN()
  const mk = (freq: number, spp: number) => {
    const plate = new TransientPlate(singleToneScore(freq, F, 30))
    return new TransientSandSim(plate, initialState(plate.a, 60, 7),
      { ...DEFAULT_PARAMS, spp }, { enabled: false })
  }
  for (const freq of [152, 649, 1241.6]) {
    for (const spp of [44, 88, 176]) {
      const sim = mk(freq, spp)
      sim.advanceTo(0.15)
      console.log(`freq=${freq} spp=${spp}: bounce=${sim.nBounce} detach=${sim.nDetach} alive=${sim.nAlive} eventCap=${sim.nEventCap}`)
    }
  }
}, 1800000)

benchIt('in-sim sample count sensitivity (8 vs 32 samples)', () => {
  const F2 = calibratedForceN()
  for (const freq of [649, 1241.6]) {
    for (const samples of [8, 32]) {
      const sim = createTransientSim(singleToneScore(freq, F2, 30), 60, 7, [0, 0], { enabled: false })
      sim.contactSamples = samples
      sim.advanceTo(0.15)
      console.log(`samples freq=${freq} n=${samples}: bounce=${sim.nBounce} detach=${sim.nDetach} alive=${sim.nAlive}`)
    }
  }
}, 1800000)
