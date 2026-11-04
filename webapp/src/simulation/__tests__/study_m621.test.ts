// M6-21 studies (PHYSI_BENCH=1): still-air linear drag effect on
// transport, and impact-impulse cancellation at 152 vs 649 Hz.
// Docs: docs/m6_21_design.md.
import { describe, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function spiral(n: number, r0: number, dr: number) {
  return Array.from({ length: n }, (_, k) => {
    const ang = k * 2.399963
    const r = r0 + dr * k
    return { x: r * Math.cos(ang), y: r * Math.sin(ang), dropM: 0.003 }
  })
}

interface Opts {
  freq: number; T: number; dragGamma?: number; spinGamma?: number
  initials: { x: number; y: number; dropM: number }[]
}

function run(o: Opts) {
  const plate = new TransientPlate(singleToneScore(o.freq, calibratedForceN(), 30), [0.06, 0])
  const params: GrainCycleParams = {
    muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3,
    model: 'roll', dragGamma: o.dragGamma, spinGamma: o.spinGamma,
  }
  const sim = new GrainCycleSim(plate, params, o.initials)
  const dt = 1 / (44 * o.freq)
  const t0 = performance.now()
  const n = Math.round(o.T / dt)
  for (let k = 0; k < n; k++) sim.step(dt)
  const stats = sim.stats()
  const live = stats.filter((s) => !s.lost)
  const mean = (f: (s: typeof stats[number]) => number) =>
    live.length ? live.reduce((a, s) => a + f(s), 0) / live.length : 0
  const netImpulse = live.map((s) => Math.hypot(s.impactImpulseX, s.impactImpulseY))
  return {
    wallMs: Math.round(performance.now() - t0),
    meanNetM: Number(mean((s) => Math.hypot(
      s.netFlightX + s.netContactX, s.netFlightY + s.netContactY)).toFixed(5)),
    meanPathM: Number(mean((s) => s.pathM).toFixed(5)),
    meanR0: Number(mean((_s) => 0).toFixed(4)),
    meanRFinal: Number((live.reduce((a, s) => a + Math.hypot(s.x, s.y), 0)
      / Math.max(1, live.length)).toFixed(4)),
    meanImpulseMagSum: Number(mean((s) => s.impactImpulseMagSum).toFixed(4)),
    meanNetImpulseMag: Number((netImpulse.reduce((a, v) => a + v, 0)
      / Math.max(1, live.length)).toFixed(4)),
    cancellationRatio: Number((() => {
      const ms = mean((s) => s.impactImpulseMagSum)
      const nm = netImpulse.reduce((a, v) => a + v, 0) / Math.max(1, live.length)
      return ms > 0 ? nm / ms : 0
    })().toFixed(4)),
    totalImpacts: stats.reduce((a, s) => a + s.impacts, 0),
    lost: stats.filter((s) => s.lost).length,
    finals: stats.map((s, i) => ({
      i, r0: Number(Math.hypot(o.initials[i].x, o.initials[i].y).toFixed(4)),
      r1: Number(Math.hypot(s.x, s.y).toFixed(4)), lost: s.lost,
    })),
  }
}

describe('M6-21 air drag and cancellation studies', () => {
  benchIt('drag scan: none / physical estimate / 5x (152Hz ecc, 8 grains, 2s)', () => {
    const initials = spiral(8, 0.03, 0.008)
    const cases = [
      { label: 'no-drag', dragGamma: undefined, spinGamma: undefined },
      { label: 'stokes-estimate', dragGamma: 0.51, spinGamma: 1.7 },
      { label: 'stokes-5x', dragGamma: 2.55, spinGamma: 8.5 },
    ]
    for (const c of cases) {
      const r = run({ freq: 152, T: 2, initials, dragGamma: c.dragGamma, spinGamma: c.spinGamma })
      console.log(JSON.stringify({ study: 'drag', label: c.label, ...r }))
    }
  }, 1800000)

  benchIt('cancellation: 152Hz (2s) vs 649Hz (1s), 8 grains', () => {
    const initials = spiral(8, 0.03, 0.008)
    for (const [label, freq, T] of [['152', 152, 2], ['649', 649, 1]] as const) {
      const r = run({ freq, T, initials })
      console.log(JSON.stringify({ study: 'cancel', label, ...r }))
    }
  }, 1800000)
})
