// M6-24 studies (PHYSI_BENCH=1): drift-field map vs activity
// gradient, impact-phase/slope diagnostics, restitution and
// grain-size sensitivity. Docs: docs/m6_24_design.md.
import { describe, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const FREQ = 152
const DRIVE: [number, number] = [0.06, 0]

function probes() {
  const out: { x: number; y: number; dropM: number }[] = []
  for (const r of [0.01, 0.025, 0.04, 0.055, 0.07, 0.085, 0.10]) {
    for (const a of [0, 2.1, 4.2]) out.push({ x: r * Math.cos(a), y: r * Math.sin(a), dropM: 0.003 })
  }
  return out
}

function runProbes(T: number, extra: Partial<GrainCycleParams>) {
  const plate = new TransientPlate(singleToneScore(FREQ, calibratedForceN(), 30), DRIVE)
  const initials = probes()
  const sim = new GrainCycleSim(plate, {
    muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3,
    model: 'roll', diagFreqHz: FREQ, ...extra,
  }, initials)
  const dt = 1 / (44 * FREQ)
  const t0 = performance.now()
  for (let k = 0; k < Math.round(T / dt); k++) sim.step(dt)
  return { stats: sim.stats(), initials, wallMs: Math.round(performance.now() - t0) }
}

function activityGrid(T: number) {
  const plate = new TransientPlate(singleToneScore(FREQ, calibratedForceN(), 30), DRIVE)
  const GN = 25; const a = plate.a
  const acc = new Float64Array(GN * GN)
  const out = new Float64Array(5)
  const dt = 1 / (44 * FREQ)
  const steps = Math.round(T / dt)
  let samples = 0
  for (let k = 0; k < steps; k++) {
    plate.beginStep(k * dt, (k + 1) * dt)
    if (k % 25 !== 0) continue
    samples++
    for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
      const x = -a + (2 * a * gx) / (GN - 1)
      const y = -a + (2 * a * gy) / (GN - 1)
      if (Math.hypot(x, y) > a) continue
      plate.surfAtConsistent(x, y, (k + 1) * dt, out)
      acc[gy * GN + gx] += Math.abs(out[4])
    }
  }
  const at = (x: number, y: number): number => {
    const fx = Math.min(GN - 1.001, Math.max(0, ((x + a) / (2 * a)) * (GN - 1)))
    const fy = Math.min(GN - 1.001, Math.max(0, ((y + a) / (2 * a)) * (GN - 1)))
    const ix = fx | 0; const iy = fy | 0
    const tx = fx - ix; const ty = fy - iy
    const v = (gx: number, gy: number) => acc[gy * GN + gx] / samples
    return v(ix, iy) * (1 - tx) * (1 - ty) + v(ix + 1, iy) * tx * (1 - ty)
      + v(ix, iy + 1) * (1 - tx) * ty + v(ix + 1, iy + 1) * tx * ty
  }
  return at
}

describe('M6-24 drift diagnosis', () => {
  benchIt('drift field vs activity gradient (21 probes, 1.5s)', () => {
    const T = 1.5
    const r = runProbes(T, {})
    const act = activityGrid(T)
    const eps = 0.004
    let downGrad = 0; let live = 0
    let projSum = 0
    const rows = r.stats.map((s, i) => {
      const x0 = r.initials[i].x; const y0 = r.initials[i].y
      const dx = s.x - x0; const dy = s.y - y0
      const gAx = (act(x0 + eps, y0) - act(x0 - eps, y0)) / (2 * eps)
      const gAy = (act(x0, y0 + eps) - act(x0, y0 - eps)) / (2 * eps)
      const gMag = Math.hypot(gAx, gAy)
      const proj = gMag > 0 ? -(dx * gAx + dy * gAy) / gMag : 0
      if (!s.lost) {
        live++
        projSum += proj
        if (proj > 0) downGrad++
      }
      return {
        i, r0: Number(Math.hypot(x0, y0).toFixed(3)),
        driftM: Number(Math.hypot(dx, dy).toFixed(5)),
        downGradM: Number(proj.toFixed(5)),
        act0: Number(act(x0, y0).toFixed(1)),
        lost: s.lost,
      }
    })
    const totW = r.stats.reduce((a2, s) => a2 + s.impactImpulseMagSum, 0)
    const pc = r.stats.reduce((a2, s) => a2 + s.phaseCosW, 0)
    const ps = r.stats.reduce((a2, s) => a2 + s.phaseSinW, 0)
    const sd = r.stats.reduce((a2, s) => a2 + s.slopeDotSum, 0)
    const sm = r.stats.reduce((a2, s) => a2 + s.slopeMagWSum, 0)
    console.log(JSON.stringify({
      study: 'drift-field', wallMs: r.wallMs,
      fracDownGradient: Number((downGrad / Math.max(1, live)).toFixed(3)),
      meanDownGradM: Number((projSum / Math.max(1, live)).toFixed(5)),
      impactPhaseConcentration: Number((Math.hypot(pc, ps) / Math.max(totW, 1e-30)).toFixed(3)),
      impactPhaseRad: Number(Math.atan2(ps, pc).toFixed(3)),
      slopeAlignment: Number((sm > 0 ? sd / sm : 0).toFixed(4)),
      rows,
    }))
  }, 1800000)

  benchIt('sensitivity: restitution and grain size (21 probes, 1.5s each)', () => {
    const T = 1.5
    const act = activityGrid(T)
    const eps = 0.004
    const summarize = (label: string, r: ReturnType<typeof runProbes>) => {
      let down = 0; let live = 0; let projSum = 0; let netSum = 0
      r.stats.forEach((s, i) => {
        if (s.lost) return
        live++
        const x0 = r.initials[i].x; const y0 = r.initials[i].y
        const dx = s.x - x0; const dy = s.y - y0
        netSum += Math.hypot(dx, dy)
        const gAx = (act(x0 + eps, y0) - act(x0 - eps, y0)) / (2 * eps)
        const gAy = (act(x0, y0 + eps) - act(x0, y0 - eps)) / (2 * eps)
        const gMag = Math.hypot(gAx, gAy)
        const proj = gMag > 0 ? -(dx * gAx + dy * gAy) / gMag : 0
        projSum += proj
        if (proj > 0) down++
      })
      console.log(JSON.stringify({
        study: 'sensitivity', label, wallMs: r.wallMs,
        meanNetM: Number((netSum / Math.max(1, live)).toFixed(5)),
        fracDownGradient: Number((down / Math.max(1, live)).toFixed(3)),
        meanDownGradM: Number((projSum / Math.max(1, live)).toFixed(5)),
        lost: r.stats.filter((s) => s.lost).length,
      }))
    }
    summarize('e=0', runProbes(T, { restitution: 0 }))
    summarize('e=0.6', runProbes(T, { restitution: 0.6 }))
    summarize('R=0.1mm', runProbes(T, { radius: 1e-4 }))
    summarize('R=0.5mm', runProbes(T, { radius: 5e-4 }))
  }, 1800000)
})
