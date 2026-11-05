// M6-28 nodal-density diagnostics (PHYSI_BENCH=1).
// Nodal set is derived from the plate alone (lowest-decile
// cells of the time-averaged |wtt| grid) and is measurement
// only — it never enters grain motion. Docs:
// docs/m6_28_design.md, docs/m6_28_results.md.
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

interface NodalField {
  distToNode: (x: number, y: number) => number
  actAt: (x: number, y: number) => number
  nNodalCells: number
  threshold: number
}

function nodalField(t0: number, t1: number, drive: [number, number]): NodalField {
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), drive)
  const GN = 49; const a = plate.a
  const acc = new Float64Array(GN * GN)
  const inside: boolean[] = []
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
      if (Math.hypot(x, y) > a) { inside[gy * GN + gx] = false; continue }
      inside[gy * GN + gx] = true
      plate.surfAtConsistent(x, y, (k + 1) * DT, out)
      acc[gy * GN + gx] += Math.abs(out[4])
    }
  }
  const vals: number[] = []
  for (let i = 0; i < GN * GN; i++) if (inside[i]) vals.push(acc[i] / samples)
  vals.sort((p, q) => p - q)
  const threshold = vals[Math.floor(vals.length * 0.1)]
  const nodalCells: [number, number][] = []
  for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
    if (!inside[gy * GN + gx]) continue
    if (acc[gy * GN + gx] / samples <= threshold) {
      nodalCells.push([-a + (2 * a * gx) / (GN - 1), -a + (2 * a * gy) / (GN - 1)])
    }
  }
  const distToNode = (x: number, y: number): number => {
    let best = Infinity
    for (const [nx, ny] of nodalCells) {
      const d = Math.hypot(x - nx, y - ny)
      if (d < best) best = d
    }
    return best
  }
  const actAt = (x: number, y: number): number => {
    const fx = Math.min(GN - 1.001, Math.max(0, ((x + a) / (2 * a)) * (GN - 1)))
    const fy = Math.min(GN - 1.001, Math.max(0, ((y + a) / (2 * a)) * (GN - 1)))
    const ix = fx | 0; const iy = fy | 0
    const tx = fx - ix; const ty = fy - iy
    const v = (gx: number, gy: number) => acc[gy * GN + gx] / samples
    return v(ix, iy) * (1 - tx) * (1 - ty) + v(ix + 1, iy) * tx * (1 - ty)
      + v(ix, iy + 1) * (1 - tx) * ty + v(ix + 1, iy + 1) * tx * ty
  }
  return { distToNode, actAt, nNodalCells: nodalCells.length, threshold }
}

function nullBaseline(field: NodalField, plateRadius: number) {
  const pts: [number, number][] = []
  for (let gx = -23; gx <= 23; gx++) for (let gy = -23; gy <= 23; gy++) {
    const x = gx * 0.005; const y = gy * 0.005
    if (Math.hypot(x, y) <= plateRadius) pts.push([x, y])
  }
  const ds = pts.map(([x, y]) => field.distToNode(x, y)).sort((p, q) => p - q)
  const frac = (band: number) => ds.filter((d) => d <= band).length / ds.length
  return {
    nNull: pts.length,
    medianDistMm: +(1e3 * ds[Math.floor(ds.length / 2)]).toFixed(2),
    frac2p5: +frac(0.0025).toFixed(4),
    frac5: +frac(0.005).toFixed(4),
  }
}

interface Snap { t: number; stats: DumbbellStats[] }
const snapAt = (snaps: Snap[], t: number): Snap => {
  let best = snaps[0]
  for (const s of snaps) if (Math.abs(s.t - t) < Math.abs(best.t - t)) best = s
  return best
}

function nodalMetrics(stats: DumbbellStats[], field: NodalField, nullBase: ReturnType<typeof nullBaseline>) {
  const live = stats.filter((s) => !s.lost)
  const ds = live.map((s) => field.distToNode(s.x, s.y)).sort((p, q) => p - q)
  const frac = (band: number) => (ds.length ? ds.filter((d) => d <= band).length / ds.length : 0)
  const f25 = frac(0.0025); const f5 = frac(0.005)
  return {
    nLive: live.length,
    medianDistMm: ds.length ? +(1e3 * ds[Math.floor(ds.length / 2)]).toFixed(2) : null,
    frac2p5: +f25.toFixed(4), frac5: +f5.toFixed(4),
    enrichment2p5: nullBase.frac2p5 > 0 ? +(f25 / nullBase.frac2p5).toFixed(2) : null,
    enrichment5: nullBase.frac5 > 0 ? +(f5 / nullBase.frac5).toFixed(2) : null,
    meanAct: live.length
      ? +(live.reduce((a, s) => a + field.actAt(s.x, s.y), 0) / live.length).toFixed(1) : null,
    meanRmm: live.length
      ? +(1e3 * live.reduce((a, s) => a + Math.hypot(s.x, s.y), 0) / live.length).toFixed(2) : null,
  }
}

function runCase(
  aspect: number, drive: [number, number], dt: number, tTotal: number,
): { snaps: Snap[]; wallMs: number; plateA: number } {
  const initials = initialsFor()
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), drive)
  const sim = new DumbbellCycleSim(plate, {
    lobeRadius: R, halfSeparation: (aspect - 1) * R, beta: 2 / 5,
    muS: 0.3, muK: 0.2, restitution: 0.3,
    diagFreqHz: 152, plateHalfThickness: 5e-4,
  }, initials)
  const snaps: Snap[] = [{ t: 0, stats: sim.stats() }]
  const stepsPerSnap = Math.max(1, Math.round(0.1 / dt))
  const wall0 = performance.now()
  const total = Math.round(tTotal / dt)
  for (let k = 1; k <= total; k++) {
    sim.step(dt)
    if (k % stepsPerSnap === 0 || k === total) snaps.push({ t: k * dt, stats: sim.stats() })
  }
  const fin = snaps[snaps.length - 1].stats
  fin.forEach((s, i) => {
    expect(Number.isFinite(s.x) && Number.isFinite(s.y)).toBe(true)
    expect(s.unconsumedTotal).toBe(0)
    expect(s.failed).toBe(false)
    expect(s.netFlightX + s.netContactX + s.netCorrectionX)
      .toBeCloseTo(s.x - initials[i].x, 7)
  })
  return { snaps, wallMs: Math.round(performance.now() - wall0), plateA: plate.a }
}

function dirCorrWarm(snaps: Snap[]): number | null {
  const win = snaps.filter((s) => s.t >= T_WARM - 1e-9)
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
  return cnt ? +(dotSum / cnt).toFixed(3) : null
}

const CASES: { driveLabel: string; drive: [number, number]; aspect: number }[] = [
  { driveLabel: 'eccentric', drive: [0.06, 0], aspect: 1.0 },
  { driveLabel: 'eccentric', drive: [0.06, 0], aspect: 1.5 },
  { driveLabel: 'center', drive: [0, 0], aspect: 1.0 },
  { driveLabel: 'center', drive: [0, 0], aspect: 1.5 },
]

for (const c of CASES) {
  describe(`M6-28 nodal density ${c.driveLabel} aspect ${c.aspect}`, () => {
    benchIt(`nodal metrics ${c.driveLabel} aspect ${c.aspect}`, () => {
      const fieldRepro = nodalField(0, 1.5, c.drive)
      const fieldWarm = nodalField(T_WARM, T_TOTAL, c.drive)
      const anyPlate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), c.drive)
      const nullRepro = nullBaseline(fieldRepro, anyPlate.a)
      const nullWarm = nullBaseline(fieldWarm, anyPlate.a)
      const { snaps, wallMs } = runCase(c.aspect, c.drive, DT, T_TOTAL)
      console.log(JSON.stringify({
        study: 'm628-nodal', drive: c.driveLabel, aspect: c.aspect, wallMs,
        nodalCellsRepro: fieldRepro.nNodalCells,
        nullRepro, nullWarm,
        atStart: nodalMetrics(snaps[0].stats, fieldRepro, nullRepro),
        atReproEnd: nodalMetrics(snapAt(snaps, 1.5).stats, fieldRepro, nullRepro),
        atWarmStart: nodalMetrics(snapAt(snaps, T_WARM).stats, fieldWarm, nullWarm),
        atWarmEnd: nodalMetrics(snaps[snaps.length - 1].stats, fieldWarm, nullWarm),
        dirCorrWarm: dirCorrWarm(snaps),
      }))
    }, 3600000)
  })
}

describe('M6-28 dt convergence of nodal metrics', () => {
  benchIt('eccentric aspect 1.5 reproduction at dt/2', () => {
    const drive: [number, number] = [0.06, 0]
    const fieldRepro = nodalField(0, 1.5, drive)
    const anyPlate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), drive)
    const nullRepro = nullBaseline(fieldRepro, anyPlate.a)
    const { snaps, wallMs } = runCase(1.5, drive, DT / 2, 1.5)
    console.log(JSON.stringify({
      study: 'm628-dt2', aspect: 1.5, wallMs, nullRepro,
      atReproEnd: nodalMetrics(snaps[snaps.length - 1].stats, fieldRepro, nullRepro),
    }))
  }, 3600000)
})
