// M6-26 re-measurement studies (PHYSI_BENCH=1), run after the
// stabilization fixes (docs/m6_26_recovery.md). Same SHA and
// same measurement implementation for every compared case.
// Failure conditions (non-finite state, unconsumed time,
// decomposition mismatch) assert; direction of motion never does.
import { describe, expect, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams, type GrainStats } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'
import { G } from '../surface'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const DT152 = 1 / (44 * 152)
const DRIVE: [number, number] = [0.06, 0]
const R_GRAIN = 2.5e-4

function probes21() {
  const out: { x: number; y: number; dropM: number }[] = []
  for (const r of [0.01, 0.025, 0.04, 0.055, 0.07, 0.085, 0.10]) {
    for (const a of [0, 2.1, 4.2]) out.push({ x: r * Math.cos(a), y: r * Math.sin(a), dropM: 0.003 })
  }
  return out
}

interface Activity {
  at: (x: number, y: number) => number
  directMeanAbsWtt: number[]
  gridVsDirectMaxRelErr: number
  lowQuartile: number
}

/** Time-averaged |wtt| over [t0, t1] on a 25x25 grid, plus direct
 *  evaluation at the probe points for grid validation. */
function activityField(
  t0: number, t1: number, freq: number, phaseRad: number,
  points: { x: number; y: number }[],
): Activity {
  const score = singleToneScore(freq, calibratedForceN(), 30)
  score.initialPhase_rad = phaseRad
  const plate = new TransientPlate(score, DRIVE)
  const GN = 25; const a = plate.a
  const acc = new Float64Array(GN * GN)
  const direct = new Float64Array(points.length)
  const out = new Float64Array(5)
  const dt = 1 / (44 * freq)
  let samples = 0
  const steps = Math.round(t1 / dt)
  for (let k = 0; k < steps; k++) {
    plate.beginStep(k * dt, (k + 1) * dt)
    if ((k + 1) * dt < t0 || k % 25 !== 0) continue
    samples++
    for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
      const x = -a + (2 * a * gx) / (GN - 1)
      const y = -a + (2 * a * gy) / (GN - 1)
      if (Math.hypot(x, y) > a) continue
      plate.surfAtConsistent(x, y, (k + 1) * dt, out)
      acc[gy * GN + gx] += Math.abs(out[4])
    }
    points.forEach((pt, i) => {
      plate.surfAtConsistent(pt.x, pt.y, (k + 1) * dt, out)
      direct[i] += Math.abs(out[4])
    })
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
  const directMean = points.map((_, i) => direct[i] / samples)
  let maxRel = 0
  points.forEach((pt, i) => {
    const g = at(pt.x, pt.y)
    if (directMean[i] > 1) maxRel = Math.max(maxRel, Math.abs(g - directMean[i]) / directMean[i])
  })
  const vals = Array.from(acc).filter((v) => v > 0).sort((p, q) => p - q)
  const lowQuartile = vals[Math.floor(vals.length / 4)] / samples
  return { at, directMeanAbsWtt: directMean, gridVsDirectMaxRelErr: maxRel, lowQuartile }
}

/** Warmup time from plate-only RMS wtt windows: first time after
 *  which three consecutive 0.1 s windows change by < 2%. */
function warmupTime(freq: number): { tWarm: number; rmsSeries: number[] } {
  const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), DRIVE)
  const dt = 1 / (44 * freq)
  const out = new Float64Array(5)
  const winSteps = Math.round(0.1 / dt)
  const rms: number[] = []
  for (let w = 0; w < 25; w++) {
    let s = 0; let n = 0
    for (let k = 0; k < winSteps; k++) {
      const t = (w * winSteps + k) * dt
      plate.beginStep(t, t + dt)
      if (k % 10 === 0) {
        plate.surfAtConsistent(0.03, 0, t + dt, out)
        s += out[4] * out[4]; n++
      }
    }
    rms.push(Math.sqrt(s / n))
    if (w >= 5) {
      const rel = [0, 1, 2].map((j) =>
        Math.abs(rms[w - j] - rms[w - j - 1]) / Math.max(1e-9, rms[w - j - 1]))
      if (rel.every((v) => v < 0.02)) return { tWarm: (w - 2) * 0.1, rmsSeries: rms }
    }
  }
  return { tWarm: 2.5, rmsSeries: rms }
}

interface Snap { t: number; stats: GrainStats[] }

function runDrift(
  initials: { x: number; y: number; dropM: number }[],
  freq: number, ht2: number, phaseRad: number, tTotal: number, dt: number,
): { snaps: Snap[]; wallMs: number } {
  const score = singleToneScore(freq, calibratedForceN(), 30)
  score.initialPhase_rad = phaseRad
  const plate = new TransientPlate(score, DRIVE)
  const sim = new GrainCycleSim(plate, {
    muS: 0.3, muK: 0.2, radius: R_GRAIN, beta: 2 / 5, restitution: 0.3,
    model: 'roll', diagFreqHz: freq, plateHalfThickness: ht2,
  } as GrainCycleParams, initials)
  const snaps: Snap[] = [{ t: 0, stats: sim.stats() }]
  const stepsPerSnap = Math.max(1, Math.round(0.1 / dt))
  const wall0 = performance.now()
  const total = Math.round(tTotal / dt)
  for (let k = 1; k <= total; k++) {
    sim.step(dt)
    if (k % stepsPerSnap === 0 || k === total) snaps.push({ t: k * dt, stats: sim.stats() })
  }
  return { snaps, wallMs: Math.round(performance.now() - wall0) }
}

function snapAt(snaps: Snap[], t: number): Snap {
  let best = snaps[0]
  for (const s of snaps) if (Math.abs(s.t - t) < Math.abs(best.t - t)) best = s
  return best
}

function analyzeWindow(
  start: Snap, end: Snap, act: Activity, label: string,
): Record<string, unknown> {
  const band = R_GRAIN // indeterminate projection band (1 grain radius)
  const rows = end.stats.map((s, i) => {
    const s0 = start.stats[i]
    const dx = s.x - s0.x; const dy = s.y - s0.y
    const gradAt = (eps: number): [number, number] => [
      (act.at(s0.x + eps, s0.y) - act.at(s0.x - eps, s0.y)) / (2 * eps),
      (act.at(s0.x, s0.y + eps) - act.at(s0.x, s0.y - eps)) / (2 * eps),
    ]
    const [gx, gy] = gradAt(0.004)
    const gm = Math.hypot(gx, gy)
    const proj = gm > 0 ? -(dx * gx + dy * gy) / gm : 0
    const cls = s.lost ? 'lost' : proj > band ? 'down' : proj < -band ? 'up' : 'indeterminate'
    return {
      i, x0: +s0.x.toFixed(4), y0: +s0.y.toFixed(4),
      x1: +s.x.toFixed(4), y1: +s.y.toFixed(4), lost: s.lost,
      netMm: +(Math.hypot(dx, dy) * 1e3).toFixed(3),
      projMm: +(proj * 1e3).toFixed(3), cls,
      act0: +act.at(s0.x, s0.y).toFixed(1), act1: +act.at(s.x, s.y).toFixed(1),
      impacts: s.impacts - s0.impacts, detaches: s.detaches - s0.detaches,
      contactFracWin: +((s.contactTime - s0.contactTime) / (end.t - start.t)).toFixed(4),
    }
  })
  const live = rows.filter((r) => !r.lost)
  const mean = (f: (r: typeof rows[number]) => number) =>
    live.length ? live.reduce((a, r) => a + f(r), 0) / live.length : 0
  // Gradient-step sensitivity of the mean projection.
  const projByEps: Record<string, number> = {}
  for (const eps of [0.002, 0.004, 0.008]) {
    let sum = 0; let n = 0
    end.stats.forEach((s, i) => {
      if (s.lost) return
      const s0 = start.stats[i]
      const gx = (act.at(s0.x + eps, s0.y) - act.at(s0.x - eps, s0.y)) / (2 * eps)
      const gy = (act.at(s0.x, s0.y + eps) - act.at(s0.x, s0.y - eps)) / (2 * eps)
      const gm = Math.hypot(gx, gy)
      if (gm > 0) { sum += -((s.x - s0.x) * gx + (s.y - s0.y) * gy) / gm; n++ }
    })
    projByEps[`eps${eps}`] = +(1e3 * sum / Math.max(1, n)).toFixed(3)
  }
  return {
    label, tStart: +start.t.toFixed(2), tEnd: +end.t.toFixed(2),
    nInit: rows.length, nLive: live.length, nLost: rows.length - live.length,
    down: live.filter((r) => r.cls === 'down').length,
    up: live.filter((r) => r.cls === 'up').length,
    indeterminate: live.filter((r) => r.cls === 'indeterminate').length,
    meanNetMm: +mean((r) => r.netMm).toFixed(3),
    meanProjMm: +mean((r) => r.projMm).toFixed(3),
    meanAct0: +mean((r) => r.act0).toFixed(1),
    meanAct1: +mean((r) => r.act1).toFixed(1),
    projByEpsMm: projByEps,
    rows,
  }
}

function gateChecks(snaps: Snap[], initials: { x: number; y: number }[]): void {
  const fin = snaps[snaps.length - 1].stats
  fin.forEach((s, i) => {
    for (const v of [s.x, s.y, s.vxNow, s.vyNow, s.pathM]) {
      expect(Number.isFinite(v)).toBe(true)
    }
    expect(s.unconsumedTotal).toBe(0)
    expect(s.failed).toBe(false)
    expect(s.netFlightX + s.netContactX + s.netCorrectionX)
      .toBeCloseTo(s.x - initials[i].x, 7)
    expect(s.netFlightY + s.netContactY + s.netCorrectionY)
      .toBeCloseTo(s.y - initials[i].y, 7)
    expect(s.impactImpulseX + s.contactImpulseX + s.conversionImpulseX)
      .toBeCloseTo(s.vxNow, 6)
    expect(s.impactImpulseY + s.contactImpulseY + s.conversionImpulseY)
      .toBeCloseTo(s.vyNow, 6)
  })
}

/** MSD / direction correlation / low-activity residence over the
 *  snapshot series inside [t0, t1]. */
function spreadStats(snaps: Snap[], t0: number, act: Activity): Record<string, unknown> {
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
  // Direction correlation between consecutive snapshot increments.
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
  // Diffusion fit on the second half of the MSD series.
  const half = msd.slice(Math.floor(msd.length / 2))
  const mt = half.reduce((a, p) => a + p.t, 0) / half.length
  const mm = half.reduce((a, p) => a + p.msdMm2, 0) / half.length
  let num = 0; let den = 0; let ssRes = 0; let ssTot = 0
  for (const p of half) { num += (p.t - mt) * (p.msdMm2 - mm); den += (p.t - mt) ** 2 }
  const slope = den > 0 ? num / den : 0
  for (const p of half) {
    ssRes += (p.msdMm2 - (mm + slope * (p.t - mt))) ** 2
    ssTot += (p.msdMm2 - mm) ** 2
  }
  const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 0
  return {
    msd,
    dirCorr: cnt ? +(dotSum / cnt).toFixed(3) : null,
    lowActivityResidenceFrac: totTime ? +(lowTime / totTime).toFixed(3) : null,
    diffusionFitSecondHalf: {
      slopeMm2PerS: +slope.toFixed(3), Dmm2PerS: +(slope / 4).toFixed(3),
      r2: +r2.toFixed(3), acceptedAsDiffusive: r2 > 0.9,
    },
  }
}

describe('M6-26 plate warmup and local Gamma (plate only)', () => {
  benchIt('warmup criterion and detach-possible fraction at 152/649Hz', () => {
    for (const freq of [152, 649]) {
      const { tWarm, rmsSeries } = warmupTime(freq)
      // Local Gamma-like stats after warmup: max/mean |wtt|/g and
      // fraction of (grid x time) samples with wtt < -g. This is a
      // diagnostic of the drive condition, not a detach proof for
      // general states (normal force also depends on motion).
      const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), DRIVE)
      const dt = 1 / (44 * freq)
      const out = new Float64Array(5)
      let maxAbs = 0; let sumAbs = 0; let below = 0; let n = 0
      const steps = Math.round((tWarm + 0.5) / dt)
      for (let k = 0; k < steps; k++) {
        plate.beginStep(k * dt, (k + 1) * dt)
        if ((k + 1) * dt < tWarm || k % 10 !== 0) continue
        for (let gx = -5; gx <= 5; gx++) for (let gy = -5; gy <= 5; gy++) {
          const x = gx * 0.012; const y = gy * 0.012
          if (Math.hypot(x, y) > plate.a) continue
          plate.surfAtConsistent(x, y, (k + 1) * dt, out)
          maxAbs = Math.max(maxAbs, Math.abs(out[4]))
          sumAbs += Math.abs(out[4]); n++
          if (out[4] < -G) below++
        }
      }
      console.log(JSON.stringify({
        study: 'm626-plate', freq, tWarm,
        rmsWttSeries: rmsSeries.map((v) => +v.toFixed(1)),
        maxGammaLocal: +(maxAbs / G).toFixed(2),
        meanGammaLocal: +(sumAbs / n / G).toFixed(3),
        fracWttBelowMinusG: +(below / n).toFixed(4),
      }))
      expect(Number.isFinite(maxAbs)).toBe(true)
    }
  }, 600000)
})

describe('M6-26 drift re-measurement', () => {
  benchIt('reproduction [0,1.5s] and warmup windows, H=0/1/5mm, 21 probes, 152Hz ecc', () => {
    const initials = probes21()
    const { tWarm } = warmupTime(152)
    const tTotal = tWarm + 1.5
    const actRepro = activityField(0, 1.5, 152, 0, initials)
    const actWarm = activityField(tWarm, tTotal, 152, 0, initials)
    console.log(JSON.stringify({
      study: 'm626-actgrid', tWarm,
      gridVsDirectMaxRelErrRepro: +actRepro.gridVsDirectMaxRelErr.toFixed(4),
      gridVsDirectMaxRelErrWarm: +actWarm.gridVsDirectMaxRelErr.toFixed(4),
    }))
    for (const ht2 of [0, 5e-4, 2.5e-3]) {
      const { snaps, wallMs } = runDrift(initials, 152, ht2, 0, tTotal, DT152)
      gateChecks(snaps, initials)
      const repro = analyzeWindow(snaps[0], snapAt(snaps, 1.5), actRepro, 'reproduction')
      const warm = analyzeWindow(snapAt(snaps, tWarm), snaps[snaps.length - 1], actWarm, 'warmup')
      const spread = spreadStats(snaps, tWarm, actWarm)
      const fin = snaps[snaps.length - 1].stats
      const live = fin.filter((s) => !s.lost)
      const slopeAlign = (() => {
        const d = live.reduce((a, s) => a + s.slopeDotSum, 0)
        const m = live.reduce((a, s) => a + s.slopeMagWSum, 0)
        return m > 0 ? +(d / m).toFixed(4) : null
      })()
      console.log(JSON.stringify({
        study: 'm626-drift', ht2, Hmm: 2 * ht2 * 1e3, wallMs, repro, warm, slopeAlign,
        totImpacts: fin.reduce((a, s) => a + s.impacts, 0),
        spreadWarm: spread,
      }))
    }
  }, 3600000)

  benchIt('release-phase variants (7 probes x 4 phases, 1.0s, H=1mm)', () => {
    const initials = [0.01, 0.025, 0.04, 0.055, 0.07, 0.085, 0.10]
      .map((r) => ({ x: r, y: 0, dropM: 0.003 }))
    for (const phase of [0, Math.PI / 2, Math.PI, 3 * Math.PI / 2]) {
      const act = activityField(0, 1.0, 152, phase, initials)
      const { snaps, wallMs } = runDrift(initials, 152, 5e-4, phase, 1.0, DT152)
      gateChecks(snaps, initials)
      const res = analyzeWindow(snaps[0], snaps[snaps.length - 1], act, `phase${phase}`)
      console.log(JSON.stringify({ study: 'm626-phase', phase, wallMs, ...res }))
    }
  }, 3600000)

  benchIt('dt convergence (2 probes, 0.75s, H=1mm; flat-drop capture vs dt)', () => {
    const initials = [
      { x: 0.04, y: 0, dropM: 0.003 },
      { x: 0.07 * Math.cos(2.1), y: 0.07 * Math.sin(2.1), dropM: 0.003 },
    ]
    for (const scale of [1, 2, 4]) {
      const dt = DT152 / scale
      const { snaps, wallMs } = runDrift(initials, 152, 5e-4, 0, 0.75, dt)
      gateChecks(snaps, initials)
      const fin = snaps[snaps.length - 1].stats
      console.log(JSON.stringify({
        study: 'm626-dtconv', scale, wallMs,
        grains: fin.map((s, i) => ({
          i, x: +s.x.toFixed(5), y: +s.y.toFixed(5), lost: s.lost,
          netMm: +(1e3 * Math.hypot(s.x - initials[i].x, s.y - initials[i].y)).toFixed(3),
          impacts: s.impacts, detaches: s.detaches,
          contactFrac: +(s.contactTime / 0.75).toFixed(4),
        })),
      }))
    }
    for (const scale of [1, 2, 4]) {
      const dt = DT152 / scale
      const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
      const sim = new GrainCycleSim(plate, {
        muS: 0.3, muK: 0.2, radius: R_GRAIN, beta: 2 / 5, restitution: 0.3, model: 'roll',
      } as GrainCycleParams, [{ x: -0.05, y: 0, dropM: 0.002 }])
      for (let k = 0; k < Math.round(0.3 / dt); k++) sim.step(dt)
      const st = sim.stats()[0]
      console.log(JSON.stringify({
        study: 'm626-capture-dt', scale, impacts: st.impacts,
        contactTime: +st.contactTime.toFixed(4), unconsumed: st.unconsumedTotal,
      }))
      expect(st.unconsumedTotal).toBe(0)
    }
  }, 3600000)

  benchIt('drift at 649Hz (7 probes, 0.5s, H=0/1mm)', () => {
    const initials = [0.01, 0.025, 0.04, 0.055, 0.07, 0.085, 0.10]
      .map((r) => ({ x: r * Math.cos(0.7), y: r * Math.sin(0.7), dropM: 0.003 }))
    const act = activityField(0, 0.5, 649, 0, initials)
    for (const ht2 of [0, 5e-4]) {
      const { snaps, wallMs } = runDrift(initials, 649, ht2, 0, 0.5, 1 / (44 * 649))
      gateChecks(snaps, initials)
      const res = analyzeWindow(snaps[0], snaps[snaps.length - 1], act, '649')
      console.log(JSON.stringify({ study: 'm626-649', ht2, wallMs, ...res }))
    }
  }, 3600000)
})
