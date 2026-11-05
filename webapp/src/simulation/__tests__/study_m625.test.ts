// M6-25 study (PHYSI_BENCH=1): drift field with Kirchhoff
// top-surface motion (ht2 = H/2). Physical H = 1mm from
// reference/plate_transient.py; an exaggerated H = 5mm case
// shows the trend direction only. Docs: docs/m6_25_design.md.
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

function activityGrid(T: number) {
  const plate = new TransientPlate(singleToneScore(FREQ, calibratedForceN(), 30), DRIVE)
  const GN = 25; const a = plate.a
  const acc = new Float64Array(GN * GN)
  const out = new Float64Array(5)
  const dt = 1 / (44 * FREQ)
  let samples = 0
  for (let k = 0; k < Math.round(T / dt); k++) {
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
  return (x: number, y: number): number => {
    const fx = Math.min(GN - 1.001, Math.max(0, ((x + a) / (2 * a)) * (GN - 1)))
    const fy = Math.min(GN - 1.001, Math.max(0, ((y + a) / (2 * a)) * (GN - 1)))
    const ix = fx | 0; const iy = fy | 0
    const tx = fx - ix; const ty = fy - iy
    const v = (gx: number, gy: number) => acc[gy * GN + gx] / samples
    return v(ix, iy) * (1 - tx) * (1 - ty) + v(ix + 1, iy) * tx * (1 - ty)
      + v(ix, iy + 1) * (1 - tx) * ty + v(ix + 1, iy + 1) * tx * ty
  }
}

describe('M6-25 top-surface drift', () => {
  benchIt('drift field with ht2 = 0.5mm (physical) and 2.5mm (trend)', () => {
    const T = 1.5
    const act = activityGrid(T)
    const eps = 0.004
    for (const ht2 of [5e-4, 2.5e-3]) {
      const plate = new TransientPlate(singleToneScore(FREQ, calibratedForceN(), 30), DRIVE)
      const initials = probes()
      const sim = new GrainCycleSim(plate, {
        muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3,
        model: 'roll', diagFreqHz: FREQ, plateHalfThickness: ht2,
      } as GrainCycleParams, initials)
      const dt = 1 / (44 * FREQ)
      const t0 = performance.now()
      for (let k = 0; k < Math.round(T / dt); k++) sim.step(dt)
      const stats = sim.stats()
      let down = 0; let live = 0; let projSum = 0; let netSum = 0
      stats.forEach((s, i) => {
        if (s.lost) return
        live++
        const x0 = initials[i].x; const y0 = initials[i].y
        const dx = s.x - x0; const dy = s.y - y0
        netSum += Math.hypot(dx, dy)
        const gAx = (act(x0 + eps, y0) - act(x0 - eps, y0)) / (2 * eps)
        const gAy = (act(x0, y0 + eps) - act(x0, y0 - eps)) / (2 * eps)
        const gMag = Math.hypot(gAx, gAy)
        const proj = gMag > 0 ? -(dx * gAx + dy * gAy) / gMag : 0
        projSum += proj
        if (proj > 0) down++
      })
      const totW = stats.reduce((a2, s) => a2 + s.impactImpulseMagSum, 0)
      const sd = stats.reduce((a2, s) => a2 + s.slopeDotSum, 0)
      const sm = stats.reduce((a2, s) => a2 + s.slopeMagWSum, 0)
      void totW
      console.log(JSON.stringify({
        study: 'top-surface', ht2, wallMs: Math.round(performance.now() - t0),
        meanNetM: Number((netSum / Math.max(1, live)).toFixed(5)),
        fracDownGradient: Number((down / Math.max(1, live)).toFixed(3)),
        meanDownGradM: Number((projSum / Math.max(1, live)).toFixed(5)),
        slopeAlignment: Number((sm > 0 ? sd / sm : 0).toFixed(4)),
        lost: stats.filter((s) => s.lost).length,
      }))
    }
  }, 1800000)
})
