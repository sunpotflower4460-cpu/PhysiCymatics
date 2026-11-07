// M6-29 residence-vs-escape diagnostics (PHYSI_BENCH=1).
// Zone definitions are derived from the plate alone and are
// measurement only. Every plate step is recorded so zone
// definitions can be varied on the identical trajectories.
// Docs: docs/m6_29_design.md, docs/m6_29_results.md.
import { describe, expect, it } from 'vitest'
import {
  DumbbellCycleSim, type DumbbellInitial,
} from '../dumbbellCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const DT = 1 / (44 * 152)
const R = 2.5e-4
const T_WARM = 1.8
const T_TOTAL = 3.3
const ZONE_NAMES = ['node', 'intermediate', 'antinode'] as const

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

interface PlateField {
  actAt: (x: number, y: number) => number
  distForFraction: (fraction: number, x: number, y: number) => number
  q90: number
  thresholds: Record<string, number>
}

function plateField(t0: number, t1: number, drive: [number, number]): PlateField {
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
  const quantile = (f: number): number => vals[Math.floor(vals.length * f)]
  const thresholds: Record<string, number> = {}
  const cellsByFraction = new Map<number, [number, number][]>()
  for (const f of [0.05, 0.1, 0.2]) {
    const th = quantile(f)
    thresholds[`f${f}`] = th
    const cells: [number, number][] = []
    for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
      if (!inside[gy * GN + gx]) continue
      if (acc[gy * GN + gx] / samples <= th) {
        cells.push([-a + (2 * a * gx) / (GN - 1), -a + (2 * a * gy) / (GN - 1)])
      }
    }
    cellsByFraction.set(f, cells)
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
  const distForFraction = (fraction: number, x: number, y: number): number => {
    const cells = cellsByFraction.get(fraction) as [number, number][]
    let best = Infinity
    for (const [nx, ny] of cells) {
      const dx = x - nx; const dy = y - ny
      const d2 = dx * dx + dy * dy
      if (d2 < best) best = d2
    }
    return Math.sqrt(best)
  }
  return { actAt, distForFraction, q90: quantile(0.9), thresholds }
}

interface ZoneDef {
  label: string
  nodalFraction: number
  nodeBandM: number
}
const PRIMARY: ZoneDef = { label: 'primary-f10-band5mm', nodalFraction: 0.1, nodeBandM: 0.005 }
const SENSITIVITIES: ZoneDef[] = [
  PRIMARY,
  { label: 'band2p5mm', nodalFraction: 0.1, nodeBandM: 0.0025 },
  { label: 'band10mm', nodalFraction: 0.1, nodeBandM: 0.01 },
  { label: 'nodal-set5pct', nodalFraction: 0.05, nodeBandM: 0.005 },
  { label: 'nodal-set20pct', nodalFraction: 0.2, nodeBandM: 0.005 },
]

function zoneAt(field: PlateField, def: ZoneDef, x: number, y: number): 0 | 1 | 2 {
  const d = field.distForFraction(def.nodalFraction, x, y)
  if (d <= def.nodeBandM) return 0
  if (d > 0.01 && field.actAt(x, y) >= field.q90) return 2
  return 1
}

function nullAreas(field: PlateField, def: ZoneDef, plateRadius: number): number[] {
  const counts = [0, 0, 0]
  let n = 0
  for (let gx = -23; gx <= 23; gx++) for (let gy = -23; gy <= 23; gy++) {
    const x = gx * 0.005; const y = gy * 0.005
    if (Math.hypot(x, y) > plateRadius) continue
    counts[zoneAt(field, def, x, y)]++
    n++
  }
  return counts.map((c) => c / n)
}

interface Trajectory {
  nGrains: number
  nSamples: number
  dt: number
  x: Float64Array[]
  y: Float64Array[]
  lost: Uint8Array[]
  mode: Uint8Array[] // 0 fly, 1 contact, 2 lost (sample-end observation)
  contactDt: Float64Array[] // per-sample interval ending at index (index 0 = 0)
  impactsDt: Float64Array[]
  detachesDt: Float64Array[]
  contactImpulsesDt: Float64Array[]
}

function runTrajectory(
  aspect: number, drive: [number, number], dt: number, tTotal: number,
): { traj: Trajectory; wallMs: number; plateA: number } {
  const initials = initialsFor()
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), drive)
  const sim = new DumbbellCycleSim(plate, {
    lobeRadius: R, halfSeparation: (aspect - 1) * R, beta: 2 / 5,
    muS: 0.3, muK: 0.2, restitution: 0.3,
    diagFreqHz: 152, plateHalfThickness: 5e-4,
  }, initials)
  const steps = Math.round(tTotal / dt)
  const ng = initials.length
  const traj: Trajectory = {
    nGrains: ng, nSamples: steps + 1, dt,
    x: [], y: [], lost: [], mode: [],
    contactDt: [], impactsDt: [], detachesDt: [], contactImpulsesDt: [],
  }
  for (let g = 0; g < ng; g++) {
    traj.x.push(new Float64Array(steps + 1))
    traj.y.push(new Float64Array(steps + 1))
    traj.lost.push(new Uint8Array(steps + 1))
    traj.mode.push(new Uint8Array(steps + 1))
    traj.contactDt.push(new Float64Array(steps + 1))
    traj.impactsDt.push(new Float64Array(steps + 1))
    traj.detachesDt.push(new Float64Array(steps + 1))
    traj.contactImpulsesDt.push(new Float64Array(steps + 1))
  }
  let prev = sim.stats()
  const record = (idx: number, stats: typeof prev): void => {
    stats.forEach((s, g) => {
      traj.x[g][idx] = s.x
      traj.y[g][idx] = s.y
      traj.lost[g][idx] = s.lost ? 1 : 0
      traj.mode[g][idx] = s.modeNow === 'contact' ? 1 : s.modeNow === 'lost' ? 2 : 0
      if (idx > 0) {
        traj.contactDt[g][idx] = s.contactTime - prev[g].contactTime
        traj.impactsDt[g][idx] = s.impacts - prev[g].impacts
        traj.detachesDt[g][idx] = s.detaches - prev[g].detaches
        traj.contactImpulsesDt[g][idx] = s.contactImpulses - prev[g].contactImpulses
      }
    })
    prev = stats
  }
  record(0, prev)
  const wall0 = performance.now()
  for (let k = 1; k <= steps; k++) {
    sim.step(dt)
    record(k, sim.stats())
  }
  const fin = prev
  fin.forEach((s, i) => {
    expect(Number.isFinite(s.x) && Number.isFinite(s.y)).toBe(true)
    expect(s.unconsumedTotal).toBe(0)
    expect(s.failed).toBe(false)
    expect(s.netFlightX + s.netContactX + s.netCorrectionX)
      .toBeCloseTo(s.x - initials[i].x, 7)
    expect(s.netFlightY + s.netContactY + s.netCorrectionY)
      .toBeCloseTo(s.y - initials[i].y, 7)
  })
  return { traj, wallMs: Math.round(performance.now() - wall0), plateA: plate.a }
}

interface Episode {
  zone: 0 | 1 | 2
  duration: number
  contactDuration: number
  impacts: number
  detaches: number
  contactImpulses: number
  endType: 'exit' | 'lost' | 'rightCensored' | 'open'
  leftCensored: boolean
  firstContactFrac: number
  lastContactFrac: number
  lastMode: number
}

function median(vals: number[]): number | null {
  if (!vals.length) return null
  const s = [...vals].sort((a, b) => a - b)
  return s[Math.floor(s.length / 2)]
}
const round = (v: number | null, digits = 4): number | null =>
  v === null || !Number.isFinite(v) ? null : +v.toFixed(digits)

function analyzeWindow(
  traj: Trajectory, field: PlateField, def: ZoneDef,
  startS: number, endS: number, areas: number[],
): Record<string, unknown> {
  const { dt } = traj
  const startIdx = Math.round(startS / dt)
  const endIdx = Math.min(traj.nSamples - 1, Math.round(endS / dt))
  const windowDuration = (endIdx - startIdx) * dt
  const zoneExposure = [0, 0, 0]
  const zoneContact = [0, 0, 0]
  const entries = [0, 0, 0]
  const exits = [0, 0, 0]
  const transitions = [
    [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0],
  ] // to: node, intermediate, antinode, lost, rightCensored
  const episodes: Episode[] = []
  let totalLiveExposure = 0
  let lostGrainsInWindow = 0

  for (let g = 0; g < traj.nGrains; g++) {
    if (traj.lost[g][startIdx]) continue
    let zone = zoneAt(field, def, traj.x[g][startIdx], traj.y[g][startIdx])
    let ep: Episode = {
      zone, duration: 0, contactDuration: 0, impacts: 0, detaches: 0,
      contactImpulses: 0, endType: 'open', leftCensored: true,
      firstContactFrac: 0, lastContactFrac: 0, lastMode: traj.mode[g][startIdx],
    }
    let firstInterval = true
    for (let k = startIdx; k < endIdx; k++) {
      const contactInc = Math.min(dt, Math.max(0, traj.contactDt[g][k + 1]))
      ep.duration += dt
      ep.contactDuration += contactInc
      ep.impacts += traj.impactsDt[g][k + 1]
      ep.detaches += traj.detachesDt[g][k + 1]
      ep.contactImpulses += traj.contactImpulsesDt[g][k + 1]
      if (firstInterval) ep.firstContactFrac = contactInc / dt
      ep.lastContactFrac = contactInc / dt
      ep.lastMode = traj.mode[g][k + 1]
      firstInterval = false
      totalLiveExposure += dt
      if (traj.lost[g][k + 1]) {
        ep.endType = 'lost'
        episodes.push(ep)
        transitions[ep.zone][3]++
        lostGrainsInWindow++
        ep = null as unknown as Episode
        break
      }
      const nextZone = zoneAt(field, def, traj.x[g][k + 1], traj.y[g][k + 1])
      if (nextZone !== zone) {
        ep.endType = 'exit'
        episodes.push(ep)
        exits[zone]++
        transitions[zone][nextZone]++
        entries[nextZone]++
        zone = nextZone
        ep = {
          zone, duration: 0, contactDuration: 0, impacts: 0, detaches: 0,
          contactImpulses: 0, endType: 'open', leftCensored: false,
          firstContactFrac: 0, lastContactFrac: 0, lastMode: traj.mode[g][k + 1],
        }
        firstInterval = true
      }
    }
    if (ep) {
      ep.endType = 'rightCensored'
      episodes.push(ep)
      transitions[ep.zone][4]++
    }
  }

  for (const ep of episodes) {
    zoneExposure[ep.zone] += ep.duration
    zoneContact[ep.zone] += ep.contactDuration
  }

  const zoneSummary = ZONE_NAMES.map((name, z) => {
    const eps = episodes.filter((e) => e.zone === z)
    const exitEps = eps.filter((e) => e.endType === 'exit')
    const byDominance = (contactDominant: boolean) => {
      const group = eps.filter(
        (e) => (e.duration > 0 && e.contactDuration / e.duration >= 0.5) === contactDominant)
      const groupExits = group.filter((e) => e.endType === 'exit')
      return {
        episodes: group.length,
        exitEpisodes: groupExits.length,
        durationShareOfZone: zoneExposure[z] > 0
          ? round(group.reduce((a, e) => a + e.duration, 0) / zoneExposure[z]) : null,
        medianExitDurationS: round(median(groupExits.map((e) => e.duration)), 4),
      }
    }
    const exitFinalContact = exitEps.filter((e) => e.lastContactFrac >= 0.5).length
    return {
      zone: name,
      areaFractionNull: round(areas[z]),
      exposureS: round(zoneExposure[z], 3),
      occupancyLive: totalLiveExposure > 0 ? round(zoneExposure[z] / totalLiveExposure) : null,
      occupancyInitialPop: round(zoneExposure[z] / (traj.nGrains * windowDuration)),
      occupancyEnrichmentVsArea: areas[z] > 0 && totalLiveExposure > 0
        ? round((zoneExposure[z] / totalLiveExposure) / areas[z], 3) : null,
      contactTimeS: round(zoneContact[z], 3),
      flightTimeS: round(zoneExposure[z] - zoneContact[z], 3),
      contactFraction: zoneExposure[z] > 0 ? round(zoneContact[z] / zoneExposure[z]) : null,
      episodes: eps.length,
      entries: entries[z],
      exits: exits[z],
      entryRatePerOutsideS: totalLiveExposure - zoneExposure[z] > 0
        ? round(entries[z] / (totalLiveExposure - zoneExposure[z]), 3) : null,
      escapeRatePerS: zoneExposure[z] > 0 ? round(exits[z] / zoneExposure[z], 3) : null,
      leftCensoredEpisodes: eps.filter((e) => e.leftCensored).length,
      rightCensoredEpisodes: eps.filter((e) => e.endType === 'rightCensored').length,
      lostEpisodes: eps.filter((e) => e.endType === 'lost').length,
      medianExitDurationS: round(median(exitEps.map((e) => e.duration)), 4),
      medianAllObservedDurationS: round(median(eps.map((e) => e.duration)), 4),
      exitSurvivalGt0p1: exitEps.length
        ? round(exitEps.filter((e) => e.duration > 0.1).length / exitEps.length, 3) : null,
      exitSurvivalGt0p5: exitEps.length
        ? round(exitEps.filter((e) => e.duration > 0.5).length / exitEps.length, 3) : null,
      exitSurvivalGt1p0: exitEps.length
        ? round(exitEps.filter((e) => e.duration > 1.0).length / exitEps.length, 3) : null,
      medianEpisodeContactFraction: round(median(
        eps.map((e) => (e.duration > 0 ? e.contactDuration / e.duration : 0))), 4),
      exitsFromContactInterval: exitFinalContact,
      exitsFromFlightInterval: exitEps.length - exitFinalContact,
      impacts: eps.reduce((a, e) => a + e.impacts, 0),
      detaches: eps.reduce((a, e) => a + e.detaches, 0),
      contactImpulses: eps.reduce((a, e) => a + e.contactImpulses, 0),
      impactsPerExposureS: zoneExposure[z] > 0
        ? round(eps.reduce((a, e) => a + e.impacts, 0) / zoneExposure[z], 2) : null,
      contactDominant: byDominance(true),
      flightDominant: byDominance(false),
    }
  })

  return {
    window: [startS, endIdx * dt],
    definition: def.label,
    totalLiveExposureS: round(totalLiveExposure, 3),
    initialPopulationExposureS: round(traj.nGrains * windowDuration, 3),
    lostGrainsInWindow,
    transitionsToNodeInterAntiLostCensored: transitions,
    zones: zoneSummary,
  }
}

function condensedSensitivity(
  analysis: Record<string, unknown>,
): Record<string, unknown> {
  const zones = analysis.zones as Record<string, unknown>[]
  const pick = (z: Record<string, unknown>) => ({
    zone: z.zone,
    occupancyLive: z.occupancyLive,
    enrichment: z.occupancyEnrichmentVsArea,
    entryRate: z.entryRatePerOutsideS,
    escapeRate: z.escapeRatePerS,
    medianExitDurationS: z.medianExitDurationS,
    contactFraction: z.contactFraction,
    episodes: z.episodes,
    exits: z.exits,
  })
  return {
    definition: analysis.definition,
    window: analysis.window,
    node: pick(zones[0]),
    antinode: pick(zones[2]),
  }
}

interface CaseDef {
  label: string
  driveLabel: string
  drive: [number, number]
  aspect: number
  dt: number
}

const CASES: CaseDef[] = [
  { label: 'eccentric-aspect1', driveLabel: 'eccentric', drive: [0.06, 0], aspect: 1.0, dt: DT },
  { label: 'eccentric-aspect1p5', driveLabel: 'eccentric', drive: [0.06, 0], aspect: 1.5, dt: DT },
  { label: 'center-aspect1', driveLabel: 'center', drive: [0, 0], aspect: 1.0, dt: DT },
  { label: 'center-aspect1p5', driveLabel: 'center', drive: [0, 0], aspect: 1.5, dt: DT },
  { label: 'eccentric-aspect1p5-dt2', driveLabel: 'eccentric', drive: [0.06, 0], aspect: 1.5, dt: DT / 2 },
]

for (const c of CASES) {
  describe(`M6-29 residence/escape ${c.label}`, () => {
    benchIt(`residence escape ${c.label}`, () => {
      const fieldRepro = plateField(0, 1.5, c.drive)
      const fieldWarm = plateField(T_WARM, T_TOTAL, c.drive)
      const { traj, wallMs, plateA } = runTrajectory(c.aspect, c.drive, c.dt, T_TOTAL)
      const windows = [
        { field: fieldRepro, start: 0, end: 1.5, label: 'reproduction' },
        { field: fieldWarm, start: T_WARM, end: T_TOTAL, label: 'warmup' },
      ]
      const primary: Record<string, unknown> = {}
      const sensitivity: Record<string, unknown>[] = []
      for (const w of windows) {
        const areasPrimary = nullAreas(w.field, PRIMARY, plateA)
        primary[w.label] = analyzeWindow(traj, w.field, PRIMARY, w.start, w.end, areasPrimary)
        for (const def of SENSITIVITIES) {
          const areas = nullAreas(w.field, def, plateA)
          sensitivity.push({
            window: w.label,
            ...condensedSensitivity(analyzeWindow(traj, w.field, def, w.start, w.end, areas)),
          })
        }
      }
      console.log(JSON.stringify({
        study: 'm629-residence',
        case: c.label, drive: c.driveLabel, aspect: c.aspect, dt: c.dt, wallMs,
        fieldThresholds: {
          repro: fieldRepro.thresholds, warm: fieldWarm.thresholds,
          q90Repro: fieldRepro.q90, q90Warm: fieldWarm.q90,
        },
        primaryAreas: {
          repro: nullAreas(fieldRepro, PRIMARY, plateA),
          warm: nullAreas(fieldWarm, PRIMARY, plateA),
        },
        primary,
        sensitivity,
      }))
    }, 3600000)
  })
}
