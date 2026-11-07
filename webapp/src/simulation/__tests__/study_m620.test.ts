// M6-20 studies (PHYSI_BENCH=1): transport mechanism decomposition,
// distribution vs plate activity, sensitivity scan. Exploratory
// research measurements on the grain-cycle harness; production
// sim untouched. Docs: docs/m6_20_design.md.
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

interface RunOpts {
  model: 'roll' | 'point'; freq: number; drive: [number, number]
  T: number; muS?: number; muK?: number; restitution?: number
  initials: { x: number; y: number; dropM: number }[]
}

function runModel(o: RunOpts) {
  const plate = new TransientPlate(singleToneScore(o.freq, calibratedForceN(), 30), o.drive)
  const params: GrainCycleParams = {
    muS: o.muS ?? 0.3, muK: o.muK ?? 0.2, radius: 2.5e-4, beta: 2 / 5,
    restitution: o.restitution ?? 0.3, model: o.model,
  }
  const sim = new GrainCycleSim(plate, params, o.initials)
  const dt = 1 / (44 * o.freq)
  const t0 = performance.now()
  const n = Math.round(o.T / dt)
  for (let k = 0; k < n; k++) sim.step(dt)
  return { stats: sim.stats(), wallMs: Math.round(performance.now() - t0), plate }
}

function agg(stats: ReturnType<GrainCycleSim['stats']>, initials: { x: number; y: number }[]) {
  const live = stats.filter((s) => !s.lost)
  const mean = (f: (s: typeof stats[number]) => number) =>
    live.length ? live.reduce((a, s) => a + f(s), 0) / live.length : 0
  const netFlight = live.map((s) => Math.hypot(s.netFlightX, s.netFlightY))
  const netContact = live.map((s) => Math.hypot(s.netContactX, s.netContactY))
  return {
    meanNetM: Number(mean((s) => Math.hypot(
      s.netFlightX + s.netContactX + s.netCorrectionX, s.netFlightY + s.netContactY + s.netCorrectionY)).toFixed(5)),
    meanPathM: Number(mean((s) => s.pathM).toFixed(5)),
    meanNetFlightM: Number((netFlight.reduce((a, v) => a + v, 0) / Math.max(1, live.length)).toFixed(5)),
    meanNetContactM: Number((netContact.reduce((a, v) => a + v, 0) / Math.max(1, live.length)).toFixed(5)),
    meanPathFlightM: Number(mean((s) => s.pathFlightM).toFixed(5)),
    meanPathContactM: Number(mean((s) => s.pathContactM).toFixed(5)),
    meanImpactImpulse: Number(mean((s) => Math.hypot(s.impactImpulseX, s.impactImpulseY)).toFixed(4)),
    meanContactImpulse: Number(mean((s) => Math.hypot(s.contactImpulseX, s.contactImpulseY)).toFixed(4)),
    meanContactTimeS: Number(mean((s) => s.contactTime).toFixed(4)),
    meanContactFrac: Number(mean((s) => (s.t > 0 ? s.contactTime / s.t : 0)).toFixed(4)),
    totalImpacts: stats.reduce((a, s) => a + s.impacts, 0),
    totalDetaches: stats.reduce((a, s) => a + s.detaches, 0),
    lost: stats.filter((s) => s.lost).length,
    finals: stats.map((s, i) => ({
      i, r0: Number(Math.hypot(initials[i].x, initials[i].y).toFixed(4)),
      r1: Number(Math.hypot(s.x, s.y).toFixed(4)),
      x: Number(s.x.toFixed(4)), y: Number(s.y.toFixed(4)),
    })),
  }
}

describe('M6-20 transport studies', () => {
  benchIt('mechanism decomposition: 152Hz eccentric, 8 grains, 2s, both models', () => {
    const initials = spiral(8, 0.03, 0.008)
    for (const model of ['roll', 'point'] as const) {
      const r = runModel({ model, freq: 152, drive: [0.06, 0], T: 2, initials })
      console.log(JSON.stringify({ study: 'mechanism', model, wallMs: r.wallMs, ...agg(r.stats, initials) }))
    }
  }, 1800000)

  benchIt('distribution vs activity: 16 rolling grains, 5s', () => {
    const initials = spiral(16, 0.02, 0.005)
    const freq = 152; const drive: [number, number] = [0.06, 0]; const T = 5
    const r = runModel({ model: 'roll', freq, drive, T, initials })
    // Activity map: plate-only run, mean |wtt| on a 25x25 grid.
    const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), drive)
    const GN = 25; const a = plate.a
    const acc = new Float64Array(GN * GN)
    const out = new Float64Array(5)
    const dt = 1 / (44 * freq)
    const steps = Math.round(T / dt)
    let samples = 0
    for (let k = 0; k < steps; k++) {
      plate.beginStep(k * dt, (k + 1) * dt)
      if (k % 25 !== 0) continue
      samples++
      for (let gy = 0; gy < GN; gy++) {
        for (let gx = 0; gx < GN; gx++) {
          const x = -a + (2 * a * gx) / (GN - 1)
          const y = -a + (2 * a * gy) / (GN - 1)
          if (Math.hypot(x, y) > a) continue
          plate.surfAtConsistent(x, y, (k + 1) * dt, out)
          acc[gy * GN + gx] += Math.abs(out[4])
        }
      }
    }
    const activityAt = (x: number, y: number): number => {
      const fx = Math.min(GN - 1.001, Math.max(0, ((x + a) / (2 * a)) * (GN - 1)))
      const fy = Math.min(GN - 1.001, Math.max(0, ((y + a) / (2 * a)) * (GN - 1)))
      const ix = fx | 0; const iy = fy | 0
      const tx = fx - ix; const ty = fy - iy
      const v = (gx: number, gy: number) => acc[gy * GN + gx] / samples
      return v(ix, iy) * (1 - tx) * (1 - ty) + v(ix + 1, iy) * tx * (1 - ty)
        + v(ix, iy + 1) * (1 - tx) * ty + v(ix + 1, iy + 1) * tx * ty
    }
    // Activity quartile threshold over in-disk cells (area-weighted by cell count).
    const vals: number[] = []
    for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
      const x = -a + (2 * a * gx) / (GN - 1); const y = -a + (2 * a * gy) / (GN - 1)
      if (Math.hypot(x, y) <= a) vals.push(acc[gy * GN + gx] / samples)
    }
    vals.sort((p, q) => p - q)
    const q25 = vals[(vals.length * 0.25) | 0]
    const rows = r.stats.map((s, i) => ({
      i,
      act0: Number(activityAt(initials[i].x, initials[i].y).toFixed(2)),
      act1: Number(activityAt(s.x, s.y).toFixed(2)),
      r0: Number(Math.hypot(initials[i].x, initials[i].y).toFixed(4)),
      r1: Number(Math.hypot(s.x, s.y).toFixed(4)),
      lost: s.lost,
    }))
    const live = rows.filter((row) => !row.lost)
    console.log(JSON.stringify({
      study: 'activity', wallMs: r.wallMs, activityQ25: Number(q25.toFixed(2)),
      meanAct0: Number((live.reduce((s2, row) => s2 + row.act0, 0) / live.length).toFixed(2)),
      meanAct1: Number((live.reduce((s2, row) => s2 + row.act1, 0) / live.length).toFixed(2)),
      fracFinalInLowQuartile: Number((live.filter((row) => row.act1 <= q25).length / live.length).toFixed(3)),
      fracInitInLowQuartile: Number((live.filter((row) => row.act0 <= q25).length / live.length).toFixed(3)),
      rows,
    }))
  }, 1800000)

  benchIt('sensitivity: center drive, 649Hz, friction variants (roll, 8 grains)', () => {
    const initials = spiral(8, 0.03, 0.008)
    const cases: { label: string; o: Omit<RunOpts, 'model' | 'initials'> }[] = [
      { label: '152-center', o: { freq: 152, drive: [0, 0], T: 2 } },
      { label: '649-eccentric-1s', o: { freq: 649, drive: [0.06, 0], T: 1 } },
      { label: '152-ecc-muLow', o: { freq: 152, drive: [0.06, 0], T: 2, muS: 0.1, muK: 0.08 } },
      { label: '152-ecc-muHigh', o: { freq: 152, drive: [0.06, 0], T: 2, muS: 0.6, muK: 0.5 } },
    ]
    for (const c of cases) {
      const r = runModel({ ...c.o, model: 'roll', initials })
      console.log(JSON.stringify({ study: 'sensitivity', label: c.label, wallMs: r.wallMs, ...agg(r.stats, initials) }))
    }
  }, 1800000)
})
