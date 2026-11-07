// M6-27B drift re-measurement (PHYSI_BENCH=1), identical
// protocol to study_m626: 21 probes, reproduction window
// [0, 1.5s], warmup window [1.8s, 3.3s] (tWarm measured in
// M6-26), H=1mm top surface, activity-grid projection with
// an indeterminate band, direction correlation / MSD /
// low-activity residence. Grain = dumbbell (this module);
// orientations are a deterministic PUT-IN grid, identical
// across aspect ratios. Docs: docs/m6_27b_results.md.
import { describe, expect, it } from 'vitest'
import {
  DumbbellCycleSim, type DumbbellInitial, type DumbbellStats,
} from '../dumbbellCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const DT = 1 / (44 * 152)
const DRIVE: [number, number] = [0.06, 0]
const R = 2.5e-4
const T_WARM = 1.8
const T_TOTAL = T_WARM + 1.5

function initialsFor(): DumbbellInitial[] {
  const out: DumbbellInitial[] = []
  let i = 0
  for (const r of [0.01, 0.025, 0.04, 0.055, 0.07, 0.085, 0.10]) {
    for (const a of [0, 2.1, 4.2]) {
      const polar = [0.3, 0.8, 1.3][i % 3]
      const azim = i * 2.399963
      out.push({
        x: r * Math.cos(a), y: r * Math.sin(a), dropM: 0.003,
        axis: [
          Math.sin(polar) * Math.cos(azim),
          Math.sin(polar) * Math.sin(azim),
          Math.cos(polar),
        ],
      })
      i++
    }
  }
  return out
}

interface Activity {
  at: (x: number, y: number) => number
  lowQuartile: number
}

function activityField(t0: number, t1: number): Activity {
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), DRIVE)
  const GN = 25; const a = plate.a
  const acc = new Float64Array(GN * GN)
  const out = new Float64Array(5)
  let samples = 0
  const steps = Math.round(t1 / DT)
  for (let k = 0; k < steps; k++) {
    plate.beginStep(k * DT, (k + 1) * DT)
    if ((k + 1) * DT < t0 || k % 25 !== 0) continue
    samples++
    for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
      const x = -a + (2 * a * gx) / (GN - 1)
      const y = -a + (2 * a * gy) / (GN - 1)
      if (Math.hypot(x, y) > a) continue
      plate.surfAtConsistent(x, y, (k + 1) * DT, out)
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
  const vals = Array.from(acc).filter((v) => v > 0).sort((p, q) => p - q)
  return { at, lowQuartile: vals[Math.floor(vals.length / 4)] / samples }
}

interface Snap { t: number; stats: DumbbellStats[] }
const snapAt = (snaps: Snap[], t: number): Snap => {
  let best = snaps[0]
  for (const s of snaps) if (Math.abs(s.t - t) < Math.abs(best.t - t)) best = s
  return best
}

function analyzeWindow(start: Snap, end: Snap, act: Activity, label: string) {
  const band = R
  const rows = end.stats.map((s, i) => {
    const s0 = start.stats[i]
    const dx = s.x - s0.x; const dy = s.y - s0.y
    const eps = 0.004
    const gx = (act.at(s0.x + eps, s0.y) - act.at(s0.x - eps, s0.y)) / (2 * eps)
    const gy = (act.at(s0.x, s0.y + eps) - act.at(s0.x, s0.y - eps)) / (2 * eps)
    const gm = Math.hypot(gx, gy)
    const proj = gm > 0 ? -(dx * gx + dy * gy) / gm : 0
    return {
      i, lost: s.lost, x1: +s.x.toFixed(4), y1: +s.y.toFixed(4),
      r0mm: +(1e3 * Math.hypot(s0.x, s0.y)).toFixed(2),
      r1mm: +(1e3 * Math.hypot(s.x, s.y)).toFixed(2),
      netMm: +(Math.hypot(dx, dy) * 1e3).toFixed(3),
      projMm: +(proj * 1e3).toFixed(3),
      cls: s.lost ? 'lost' : proj > band ? 'down' : proj < -band ? 'up' : 'indeterminate',
      act0: +act.at(s0.x, s0.y).toFixed(1), act1: +act.at(s.x, s.y).toFixed(1),
      impacts: s.impacts - s0.impacts, detaches: s.detaches - s0.detaches,
      contactFracWin: +((s.contactTime - s0.contactTime) / (end.t - start.t)).toFixed(4),
    }
  })
  const live = rows.filter((r) => !r.lost)
  const mean = (f: (r: typeof rows[number]) => number) =>
    live.length ? live.reduce((a, r) => a + f(r), 0) / live.length : 0
  return {
    label, tStart: +start.t.toFixed(2), tEnd: +end.t.toFixed(2),
    nLive: live.length, nLost: rows.length - live.length,
    down: live.filter((r) => r.cls === 'down').length,
    up: live.filter((r) => r.cls === 'up').length,
    indeterminate: live.filter((r) => r.cls === 'indeterminate').length,
    meanNetMm: +mean((r) => r.netMm).toFixed(3),
    meanProjMm: +mean((r) => r.projMm).toFixed(3),
    meanAct0: +mean((r) => r.act0).toFixed(1),
    meanAct1: +mean((r) => r.act1).toFixed(1),
    meanR0mm: +mean((r) => r.r0mm).toFixed(2),
    meanR1mm: +mean((r) => r.r1mm).toFixed(2),
    meanContactFracWin: +mean((r) => r.contactFracWin).toFixed(4),
    rows,
  }
}

function spreadStats(snaps: Snap[], t0: number, act: Activity) {
  const win = snaps.filter((s) => s.t >= t0 - 1e-9)
  const base = win[0]
  const msd = win.map((s) => {
    let sum = 0; let n = 0
    s.stats.forEach((g, i) => {
      if (g.lost || base.stats[i].lost) return
      sum += (g.x - base.stats[i].x) ** 2 + (g.y - base.stats[i].y) ** 2; n++
    })
    return { t: +(s.t - base.t).toFixed(2), msdMm2: +(1e6 * sum / Math.max(1, n)).toFixed(3) }
  })
  let dotSum = 0; let cnt = 0
  for (let k = 2; k < win.length; k++) {
    win[k].stats.forEach((g, i) => {
      if (g.lost) return
      const ax = win[k - 1].stats[i].x - win[k - 2].stats[i].x
      const ay = win[k - 1].stats[i].y - win[k - 2].stats[i].y
      const bx = g.x - win[k - 1].stats[i].x
      const by = g.y - win[k - 1].stats[i].y
      const ma = Math.hypot(ax, ay); const mb = Math.hypot(bx, by)
      if (ma > 1e-6 && mb > 1e-6) { dotSum += (ax * bx + ay * by) / (ma * mb); cnt++ }
    })
  }
  let lowTime = 0; let totTime = 0
  for (const s of win) {
    s.stats.forEach((g) => {
      if (g.lost) return
      totTime++
      if (act.at(g.x, g.y) <= act.lowQuartile) lowTime++
    })
  }
  return {
    msd,
    dirCorr: cnt ? +(dotSum / cnt).toFixed(3) : null,
    lowActivityResidenceFrac: totTime ? +(lowTime / totTime).toFixed(3) : null,
  }
}

for (const aspect of [1.0, 1.5, 2.0]) {
  describe(`M6-27B drift aspect ${aspect}`, () => {
    benchIt(`reproduction + warmup windows, aspect ${aspect}`, () => {
      const initials = initialsFor()
      const actRepro = activityField(0, 1.5)
      const actWarm = activityField(T_WARM, T_TOTAL)
      const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), DRIVE)
      const sim = new DumbbellCycleSim(plate, {
        lobeRadius: R, halfSeparation: (aspect - 1) * R, beta: 2 / 5,
        muS: 0.3, muK: 0.2, restitution: 0.3,
        diagFreqHz: 152, plateHalfThickness: 5e-4,
      }, initials)
      const snaps: Snap[] = [{ t: 0, stats: sim.stats() }]
      const stepsPerSnap = Math.round(0.1 / DT)
      const wall0 = performance.now()
      const total = Math.round(T_TOTAL / DT)
      for (let k = 1; k <= total; k++) {
        sim.step(DT)
        if (k % stepsPerSnap === 0 || k === total) snaps.push({ t: k * DT, stats: sim.stats() })
      }
      const wallMs = Math.round(performance.now() - wall0)
      const fin = snaps[snaps.length - 1].stats
      fin.forEach((s, i) => {
        for (const v of [s.x, s.y, s.vxNow, s.vyNow, s.czNow, s.spinMag]) {
          expect(Number.isFinite(v)).toBe(true)
        }
        expect(s.unconsumedTotal).toBe(0)
        expect(s.failed).toBe(false)
        expect(s.netFlightX + s.netContactX + s.netCorrectionX)
          .toBeCloseTo(s.x - initials[i].x, 7)
        expect(s.netFlightY + s.netContactY + s.netCorrectionY)
          .toBeCloseTo(s.y - initials[i].y, 7)
      })
      const repro = analyzeWindow(snaps[0], snapAt(snaps, 1.5), actRepro, 'reproduction')
      const warm = analyzeWindow(
        snapAt(snaps, T_WARM), snaps[snaps.length - 1], actWarm, 'warmup')
      const spread = spreadStats(snaps, T_WARM, actWarm)
      const live = fin.filter((s) => !s.lost)
      const dsum = live.reduce((a, s) => a + s.slopeDotSum, 0)
      const msum = live.reduce((a, s) => a + s.slopeMagWSum, 0)
      console.log(JSON.stringify({
        study: 'm627b-drift', aspect, wallMs, repro, warm,
        slopeAlign: msum > 0 ? +(dsum / msum).toFixed(4) : null,
        totImpacts: fin.reduce((a, s) => a + s.impacts, 0),
        totContactImpulses: fin.reduce((a, s) => a + s.contactImpulses, 0),
        spreadWarm: spread,
      }))
    }, 3600000)
  })
}

describe('M6-27B dt convergence (aspect 1.5, reproduction window)', () => {
  benchIt('plate dt/2 gives the same qualitative transport', () => {
    const aspect = 1.5
    const dt = DT / 2
    const initials = initialsFor()
    const actRepro = activityField(0, 1.5)
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), DRIVE)
    const sim = new DumbbellCycleSim(plate, {
      lobeRadius: R, halfSeparation: (aspect - 1) * R, beta: 2 / 5,
      muS: 0.3, muK: 0.2, restitution: 0.3,
      diagFreqHz: 152, plateHalfThickness: 5e-4,
    }, initials)
    const snaps: Snap[] = [{ t: 0, stats: sim.stats() }]
    const total = Math.round(1.5 / dt)
    const stepsPerSnap = Math.round(0.1 / dt)
    for (let k = 1; k <= total; k++) {
      sim.step(dt)
      if (k % stepsPerSnap === 0 || k === total) snaps.push({ t: k * dt, stats: sim.stats() })
    }
    const fin = snaps[snaps.length - 1].stats
    fin.forEach((s) => {
      expect(Number.isFinite(s.x) && Number.isFinite(s.y)).toBe(true)
      expect(s.unconsumedTotal).toBe(0)
    })
    const repro = analyzeWindow(snaps[0], snaps[snaps.length - 1], actRepro, 'dt2-reproduction')
    console.log(JSON.stringify({ study: 'm627b-dt2', aspect, repro }))
  }, 3600000)
})
