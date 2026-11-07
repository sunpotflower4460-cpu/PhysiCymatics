// M6-19 study (PHYSI_BENCH=1): horizontal transport of rolling
// spheres vs point grains on the driven consistent-geometry plate.
// Same drive, initial drops, friction and restitution; only the
// contact model differs. Exploratory research measurement — the
// production sim is untouched.
import { describe, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const T = 2.0
const initials = Array.from({ length: 8 }, (_, k) => {
  const ang = k * 2.399963
  const r = 0.03 + 0.008 * k
  return { x: r * Math.cos(ang), y: r * Math.sin(ang), dropM: 0.003 }
})

function runModel(model: 'roll' | 'point', freq: number, drive: [number, number]) {
  const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), drive)
  const params: GrainCycleParams = {
    muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3, model,
  }
  const sim = new GrainCycleSim(plate, params, initials)
  const dt = 1 / (44 * freq)
  const t0 = performance.now()
  const n = Math.round(T / dt)
  for (let k = 0; k < n; k++) sim.step(dt)
  const wallMs = performance.now() - t0
  const stats = sim.stats()
  const rows = stats.map((s, i) => ({
    i,
    netM: Number(Math.hypot(s.x - initials[i].x, s.y - initials[i].y).toFixed(5)),
    pathM: Number(s.pathM.toFixed(5)),
    impacts: s.impacts, detaches: s.detaches,
    contactFrac: Number((s.contactTime / T).toFixed(3)),
    rInitM: Number(Math.hypot(initials[i].x, initials[i].y).toFixed(4)),
    rFinalM: Number(Math.hypot(s.x, s.y).toFixed(4)),
    lost: s.lost,
  }))
  const live = rows.filter((r) => !r.lost)
  const mean = (f: (r: typeof rows[number]) => number) =>
    live.length ? live.reduce((a, r) => a + f(r), 0) / live.length : 0
  return {
    model, freq, drive, wallMs: Math.round(wallMs),
    meanNetM: Number(mean((r) => r.netM).toFixed(5)),
    meanPathM: Number(mean((r) => r.pathM).toFixed(5)),
    meanContactFrac: Number(mean((r) => r.contactFrac).toFixed(3)),
    totalImpacts: rows.reduce((a, r) => a + r.impacts, 0),
    totalDetaches: rows.reduce((a, r) => a + r.detaches, 0),
    lost: rows.filter((r) => r.lost).length,
    rows,
  }
}

describe('M6-19 rolling vs point transport on driven plate', () => {
  benchIt('152Hz eccentric drive, 8 grains, 2s, both models', () => {
    for (const model of ['roll', 'point'] as const) {
      console.log(JSON.stringify(runModel(model, 152, [0.06, 0])))
    }
  }, 1800000)
})
