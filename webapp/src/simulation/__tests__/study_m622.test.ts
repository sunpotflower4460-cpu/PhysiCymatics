// M6-22 study (PHYSI_BENCH=1): interparticle collisions ON vs OFF
// for rolling grains on the driven plate (dense packing).
// Docs: docs/m6_22_design.md.
import { describe, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

describe('M6-22 pair collision effect on transport', () => {
  benchIt('24 rolling grains, dense spiral, 152Hz ecc, 3s, collisions on/off', () => {
    // Dense sunflower cluster (spacing ~0.6mm vs diameter 0.5mm)
    // around (0.03, 0.01) so pair encounters actually happen.
    const initials = Array.from({ length: 24 }, (_, k) => {
      const ang = k * 2.399963
      const r = 0.0006 * Math.sqrt(k)
      return { x: 0.03 + r * Math.cos(ang), y: 0.01 + r * Math.sin(ang), dropM: 0.003 }
    })
    for (const on of [false, true]) {
      const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.06, 0])
      const params: GrainCycleParams = {
        muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3,
        model: 'roll', interparticle: on,
      }
      const sim = new GrainCycleSim(plate, params, initials)
      const dt = 1 / (44 * 152)
      const t0 = performance.now()
      const n = Math.round(3 / dt)
      for (let k = 0; k < n; k++) sim.step(dt)
      const stats = sim.stats()
      const live = stats.filter((s) => !s.lost)
      const mean = (f: (s: typeof stats[number]) => number) =>
        live.length ? live.reduce((a, s) => a + f(s), 0) / live.length : 0
      console.log(JSON.stringify({
        interparticle: on,
        wallMs: Math.round(performance.now() - t0),
        meanNetM: Number(mean((s) => Math.hypot(
          s.netFlightX + s.netContactX + s.netCorrectionX, s.netFlightY + s.netContactY + s.netCorrectionY)).toFixed(5)),
        meanPathM: Number(mean((s) => s.pathM).toFixed(5)),
        meanRInit: Number((initials.reduce((a, i2) => a + Math.hypot(i2.x, i2.y), 0)
          / initials.length).toFixed(4)),
        meanRFinal: Number(mean((s) => Math.hypot(s.x, s.y)).toFixed(4)),
        totalPairCollisions: stats.reduce((a, s) => a + s.pairCollisions, 0) / 2,
        totalPlateImpacts: stats.reduce((a, s) => a + s.impacts, 0),
        lost: stats.filter((s) => s.lost).length,
        finals: stats.map((s, i) => ({
          i, r0: Number(Math.hypot(initials[i].x, initials[i].y).toFixed(4)),
          r1: Number(Math.hypot(s.x, s.y).toFixed(4)), lost: s.lost,
        })),
      }))
    }
  }, 1800000)
})

describe('M6-22 ensemble check', () => {
  benchIt('3 cluster variants x on/off, 16 grains, 2s', () => {
    const variants = [
      { cx: 0.03, cy: 0.01, rot: 0 },
      { cx: -0.02, cy: 0.025, rot: 1.1 },
      { cx: 0.01, cy: -0.035, rot: 2.2 },
    ]
    for (const v of variants) {
      const initials = Array.from({ length: 16 }, (_, k) => {
        const ang = k * 2.399963 + v.rot
        const r = 0.0006 * Math.sqrt(k)
        return { x: v.cx + r * Math.cos(ang), y: v.cy + r * Math.sin(ang), dropM: 0.003 }
      })
      for (const on of [false, true]) {
        const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.06, 0])
        const sim = new GrainCycleSim(plate, {
          muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3,
          model: 'roll', interparticle: on,
        }, initials)
        const dt = 1 / (44 * 152)
        const t0 = performance.now()
        for (let k = 0; k < Math.round(2 / dt); k++) sim.step(dt)
        const stats = sim.stats()
        const live = stats.filter((s) => !s.lost)
        const mean = (f: (s: typeof stats[number]) => number) =>
          live.length ? live.reduce((a, s) => a + f(s), 0) / live.length : 0
        console.log(JSON.stringify({
          variant: [v.cx, v.cy, v.rot], interparticle: on,
          wallMs: Math.round(performance.now() - t0),
          meanNetM: Number(mean((s) => Math.hypot(
            s.netFlightX + s.netContactX + s.netCorrectionX, s.netFlightY + s.netContactY + s.netCorrectionY)).toFixed(5)),
          meanRFinal: Number(mean((s) => Math.hypot(s.x, s.y)).toFixed(4)),
          pairs: stats.reduce((a, s) => a + s.pairCollisions, 0) / 2,
          lost: stats.filter((s) => s.lost).length,
        }))
      }
    }
  }, 1800000)
})
