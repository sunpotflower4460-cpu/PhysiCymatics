// M6-31 initial-ensemble and per-grain contribution diagnostics
// (PHYSI_BENCH=1). Paired P1C4/P2C4 runs on five deterministic
// initial ensembles. Zone definitions are plate-only and
// measurement only. Docs: docs/m6_31_design.md,
// docs/m6_31_results.md.
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
const T_TOTAL = 3.3
const DRIVE: [number, number] = [0.06, 0]
const RADII = [0.01, 0.025, 0.04, 0.055, 0.07, 0.085, 0.10]
const ZONE_NAMES = ['node', 'intermediate', 'antinode'] as const

function mulberry32(seed: number): () => number {
  let a = seed >>> 0
  return () => {
    a |= 0
    a = (a + 0x6D2B79F5) | 0
    let t = Math.imul(a ^ (a >>> 15), 1 | a)
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

function baseInitials(): DumbbellInitial[] {
  const out: DumbbellInitial[] = []
  let i = 0
  for (const r of RADII) {
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

function seededInitials(seed: number): DumbbellInitial[] {
  const rand = mulberry32(seed)
  const out: DumbbellInitial[] = []
  for (const r of RADII) {
    for (let j = 0; j < 3; j++) {
      const azimuth = 2 * Math.PI * rand()
      const ez = 2 * rand() - 1
      const axisAzimuth = 2 * Math.PI * rand()
      const horizontal = Math.sqrt(Math.max(0, 1 - ez * ez))
      out.push({
        x: r * Math.cos(azimuth), y: r * Math.sin(azimuth), dropM: 0.003,
        axis: [
          horizontal * Math.cos(axisAzimuth),
          horizontal * Math.sin(axisAzimuth),
          ez,
        ],
      })
    }
  }
  return out
}

function initialsForSeed(seed: number | null): DumbbellInitial[] {
  const initials = seed === null ? baseInitials() : seededInitials(seed)
  const radii = initials.map((i) => Math.hypot(i.x, i.y)).sort((a, b) => a - b)
  const expected = RADII.flatMap((r) => [r, r, r]).sort((a, b) => a - b)
  radii.forEach((r, i) => expect(r).toBeCloseTo(expected[i], 12))
  initials.forEach((i) => {
    expect(Math.hypot(i.axis[0], i.axis[1], i.axis[2])).toBeCloseTo(1, 12)
  })
  return initials
}

interface PlateField {
  actAt: (x: number, y: number) => number
  nodalDistance: (x: number, y: number) => number
  q90: number
  threshold: number
}

function makePlateField(t0: number, t1: number): PlateField {
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), DRIVE)
  const GN = 49
  const a = plate.a
  const acc = new Float64Array(GN * GN)
  const inside = new Uint8Array(GN * GN)
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
      inside[gy * GN + gx] = 1
      plate.surfAtConsistent(x, y, (k + 1) * DT, out)
      acc[gy * GN + gx] += Math.abs(out[4])
    }
  }
  const vals: number[] = []
  for (let i = 0; i < GN * GN; i++) if (inside[i]) vals.push(acc[i] / samples)
  vals.sort((p, q) => p - q)
  const quantile = (f: number): number => vals[Math.floor(vals.length * f)]
  const threshold = quantile(0.1)
  const nodalCells: [number, number][] = []
  for (let gy = 0; gy < GN; gy++) for (let gx = 0; gx < GN; gx++) {
    if (!inside[gy * GN + gx]) continue
    if (acc[gy * GN + gx] / samples <= threshold) {
      nodalCells.push([-a + (2 * a * gx) / (GN - 1), -a + (2 * a * gy) / (GN - 1)])
    }
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
  const nodalDistance = (x: number, y: number): number => {
    let best = Infinity
    for (const [nx, ny] of nodalCells) {
      const dx = x - nx; const dy = y - ny
      const d2 = dx * dx + dy * dy
      if (d2 < best) best = d2
    }
    return Math.sqrt(best)
  }
  return { actAt, nodalDistance, q90: quantile(0.9), threshold }
}

let fieldsCache: { repro: PlateField; warm: PlateField } | null = null
function plateFields(): { repro: PlateField; warm: PlateField } {
  if (!fieldsCache) {
    fieldsCache = {
      repro: makePlateField(0, 1.5),
      warm: makePlateField(T_WARM, T_TOTAL),
    }
  }
  return fieldsCache
}

function zoneAt(field: PlateField, x: number, y: number): 0 | 1 | 2 {
  const d = field.nodalDistance(x, y)
  if (d <= 0.005) return 0
  if (d > 0.01 && field.actAt(x, y) >= field.q90) return 2
  return 1
}

function nullAreas(field: PlateField, plateRadius: number): number[] {
  const counts = [0, 0, 0]
  let n = 0
  for (let gx = -23; gx <= 23; gx++) for (let gy = -23; gy <= 23; gy++) {
    const x = gx * 0.005; const y = gy * 0.005
    if (Math.hypot(x, y) > plateRadius) continue
    counts[zoneAt(field, x, y)]++
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
  mode: Uint8Array[] // sample-end: 0 fly, 1 contact, 2 lost
  contactDt: Float64Array[] // interval ending at each sample
  impactsDt: Float64Array[]
  detachesDt: Float64Array[]
  contactImpulsesDt: Float64Array[]
}

function emptyTrajectory(nGrains: number, nSamples: number, dt: number): Trajectory {
  const makeF = () => Array.from({ length: nGrains }, () => new Float64Array(nSamples))
  const makeU = () => Array.from({ length: nGrains }, () => new Uint8Array(nSamples))
  return {
    nGrains, nSamples, dt,
    x: makeF(), y: makeF(), lost: makeU(), mode: makeU(),
    contactDt: makeF(), impactsDt: makeF(), detachesDt: makeF(),
    contactImpulsesDt: makeF(),
  }
}

function runTrajectory(
  dt: number, contactSubsteps: number, initials: DumbbellInitial[],
): { traj: Trajectory; finalStats: DumbbellStats[]; wallMs: number; plateA: number } {
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), DRIVE)
  const sim = new DumbbellCycleSim(plate, {
    lobeRadius: R, halfSeparation: 0.5 * R, beta: 2 / 5,
    muS: 0.3, muK: 0.2, restitution: 0.3,
    diagFreqHz: 152, plateHalfThickness: 5e-4, contactSubsteps,
  }, initials)
  const steps = Math.round(T_TOTAL / dt)
  const traj = emptyTrajectory(initials.length, steps + 1, dt)
  let prev = sim.stats()
  const record = (idx: number, stats: DumbbellStats[]): void => {
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
  const finalStats = prev
  finalStats.forEach((s, i) => {
    expect(Number.isFinite(s.x) && Number.isFinite(s.y)).toBe(true)
    expect(s.unconsumedTotal).toBe(0)
    expect(s.failed).toBe(false)
    expect(s.netFlightX + s.netContactX + s.netCorrectionX)
      .toBeCloseTo(s.x - initials[i].x, 7)
    expect(s.netFlightY + s.netContactY + s.netCorrectionY)
      .toBeCloseTo(s.y - initials[i].y, 7)
  })
  return {
    traj, finalStats, wallMs: Math.round(performance.now() - wall0), plateA: plate.a,
  }
}

function downsample(traj: Trajectory, factor: number): Trajectory {
  if (factor === 1) return traj
  const coarseSteps = Math.floor((traj.nSamples - 1) / factor)
  const out = emptyTrajectory(traj.nGrains, coarseSteps + 1, traj.dt * factor)
  for (let g = 0; g < traj.nGrains; g++) {
    for (let j = 0; j <= coarseSteps; j++) {
      const src = j * factor
      out.x[g][j] = traj.x[g][src]
      out.y[g][j] = traj.y[g][src]
      out.lost[g][j] = traj.lost[g][src]
      out.mode[g][j] = traj.mode[g][src]
      if (j === 0) continue
      for (let k = src - factor + 1; k <= src; k++) {
        out.contactDt[g][j] += traj.contactDt[g][k]
        out.impactsDt[g][j] += traj.impactsDt[g][k]
        out.detachesDt[g][j] += traj.detachesDt[g][k]
        out.contactImpulsesDt[g][j] += traj.contactImpulsesDt[g][k]
      }
    }
  }
  return out
}

interface Episode {
  grain: number
  zone: 0 | 1 | 2
  duration: number
  contactDuration: number
  impacts: number
  detaches: number
  contactImpulses: number
  endType: 'exit' | 'lost' | 'rightCensored' | 'open'
  leftCensored: boolean
  lastContactFrac: number
}

function median(vals: number[]): number | null {
  if (!vals.length) return null
  const s = [...vals].sort((a, b) => a - b)
  return s[Math.floor(s.length / 2)]
}
const round = (v: number | null, digits = 4): number | null =>
  v === null || !Number.isFinite(v) ? null : +v.toFixed(digits)

function concentration(values: number[]): Record<string, unknown> {
  const total = values.reduce((a, v) => a + v, 0)
  const desc = [...values].sort((a, b) => b - a)
  let grainsForHalf: number | null = null
  if (total > 0) {
    let acc = 0
    for (let i = 0; i < desc.length; i++) {
      acc += desc[i]
      if (acc >= total / 2) { grainsForHalf = i + 1; break }
    }
  }
  let gini: number | null = null
  if (total > 0 && values.length > 0) {
    const asc = [...values].sort((a, b) => a - b)
    let weighted = 0
    asc.forEach((v, i) => { weighted += (i + 1) * v })
    gini = (2 * weighted) / (values.length * total) - (values.length + 1) / values.length
  }
  return {
    population: values.length,
    total: round(total, 6),
    nonzeroGrains: values.filter((v) => v > 0).length,
    top1Share: total > 0 ? round(desc[0] / total) : null,
    top3Share: total > 0
      ? round(desc.slice(0, 3).reduce((a, v) => a + v, 0) / total) : null,
    grainsForHalf,
    gini: round(gini),
  }
}

function analyzeWindow(
  traj: Trajectory, field: PlateField, areas: number[],
  startS: number, endS: number,
  initials: DumbbellInitial[], initialField: PlateField,
): Record<string, unknown> {
  const { dt } = traj
  const startIdx = Math.round(startS / dt)
  const endIdx = Math.min(traj.nSamples - 1, Math.round(endS / dt))
  const windowDuration = (endIdx - startIdx) * dt
  const entries = [0, 0, 0]
  const exits = [0, 0, 0]
  const episodes: Episode[] = []
  const lostAtInWindow: (number | null)[] = Array(traj.nGrains).fill(null)
  let totalLiveExposure = 0
  let lostGrainsInWindow = 0

  for (let g = 0; g < traj.nGrains; g++) {
    if (traj.lost[g][startIdx]) continue
    let zone = zoneAt(field, traj.x[g][startIdx], traj.y[g][startIdx])
    let ep: Episode | null = {
      grain: g, zone, duration: 0, contactDuration: 0, impacts: 0, detaches: 0,
      contactImpulses: 0, endType: 'open', leftCensored: true, lastContactFrac: 0,
    }
    for (let k = startIdx; k < endIdx; k++) {
      const current = ep as Episode
      const contactInc = Math.min(dt, Math.max(0, traj.contactDt[g][k + 1]))
      current.duration += dt
      current.contactDuration += contactInc
      current.impacts += traj.impactsDt[g][k + 1]
      current.detaches += traj.detachesDt[g][k + 1]
      current.contactImpulses += traj.contactImpulsesDt[g][k + 1]
      current.lastContactFrac = contactInc / dt
      totalLiveExposure += dt
      if (traj.lost[g][k + 1]) {
        current.endType = 'lost'
        episodes.push(current)
        lostGrainsInWindow++
        lostAtInWindow[g] = (k + 1) * dt
        ep = null
        break
      }
      const nextZone = zoneAt(field, traj.x[g][k + 1], traj.y[g][k + 1])
      if (nextZone !== zone) {
        current.endType = 'exit'
        episodes.push(current)
        exits[zone]++
        entries[nextZone]++
        zone = nextZone
        ep = {
          grain: g, zone, duration: 0, contactDuration: 0, impacts: 0, detaches: 0,
          contactImpulses: 0, endType: 'open', leftCensored: false, lastContactFrac: 0,
        }
      }
    }
    if (ep) {
      ep.endType = 'rightCensored'
      episodes.push(ep)
    }
  }

  const exposure = [0, 0, 0]
  const contact = [0, 0, 0]
  for (const ep of episodes) {
    exposure[ep.zone] += ep.duration
    contact[ep.zone] += ep.contactDuration
  }

  const endpointDistances: number[] = []
  let aliveAtEnd = 0
  let within5AtEnd = 0
  for (let g = 0; g < traj.nGrains; g++) {
    if (traj.lost[g][endIdx]) continue
    aliveAtEnd++
    const d = field.nodalDistance(traj.x[g][endIdx], traj.y[g][endIdx])
    endpointDistances.push(d)
    if (d <= 0.005) within5AtEnd++
  }

  const zones = ZONE_NAMES.map((name, z) => {
    const eps = episodes.filter((e) => e.zone === z)
    const exitEps = eps.filter((e) => e.endType === 'exit')
    const impacts = eps.reduce((a, e) => a + e.impacts, 0)
    const detaches = eps.reduce((a, e) => a + e.detaches, 0)
    const contactImpulses = eps.reduce((a, e) => a + e.contactImpulses, 0)
    const contactDominant = eps.filter(
      (e) => e.duration > 0 && e.contactDuration / e.duration >= 0.5)
    return {
      zone: name,
      areaFractionNull: round(areas[z]),
      exposureS: round(exposure[z], 3),
      occupancyLive: totalLiveExposure > 0 ? round(exposure[z] / totalLiveExposure) : null,
      occupancyInitialPop: round(exposure[z] / (traj.nGrains * windowDuration)),
      occupancyEnrichmentVsArea: areas[z] > 0 && totalLiveExposure > 0
        ? round((exposure[z] / totalLiveExposure) / areas[z], 3) : null,
      contactFraction: exposure[z] > 0 ? round(contact[z] / exposure[z]) : null,
      episodes: eps.length,
      entries: entries[z],
      exits: exits[z],
      entryRatePerOutsideS: totalLiveExposure - exposure[z] > 0
        ? round(entries[z] / (totalLiveExposure - exposure[z]), 3) : null,
      escapeRatePerS: exposure[z] > 0 ? round(exits[z] / exposure[z], 3) : null,
      leftCensoredEpisodes: eps.filter((e) => e.leftCensored).length,
      rightCensoredEpisodes: eps.filter((e) => e.endType === 'rightCensored').length,
      lostEpisodes: eps.filter((e) => e.endType === 'lost').length,
      medianExitDurationS: round(median(exitEps.map((e) => e.duration)), 4),
      medianAllObservedDurationS: round(median(eps.map((e) => e.duration)), 4),
      medianEpisodeContactFraction: round(median(
        eps.map((e) => (e.duration > 0 ? e.contactDuration / e.duration : 0))), 4),
      contactDominantEpisodes: contactDominant.length,
      contactDominantDurationShare: exposure[z] > 0
        ? round(contactDominant.reduce((a, e) => a + e.duration, 0) / exposure[z]) : null,
      exitsFromContactInterval: exitEps.filter((e) => e.lastContactFrac >= 0.5).length,
      exitsFromFlightInterval: exitEps.filter((e) => e.lastContactFrac < 0.5).length,
      impactsPerExposureS: exposure[z] > 0 ? round(impacts / exposure[z], 2) : null,
      detachesPerExposureS: exposure[z] > 0 ? round(detaches / exposure[z], 2) : null,
      contactImpulsesPerExposureS: exposure[z] > 0
        ? round(contactImpulses / exposure[z], 2) : null,
    }
  })

  interface GrainRaw {
    grain: number
    aliveAtStart: boolean
    aliveAtEnd: boolean
    liveExposureS: number
    exposureByZone: number[]
    contactByZone: number[]
    episodesByZone: number[]
    entriesByZone: number[]
    exitsByZone: number[]
    impactsByZone: number[]
    detachesByZone: number[]
    endpointNodalDistanceM: number | null
  }
  const grainRaw: GrainRaw[] = []
  for (let g = 0; g < traj.nGrains; g++) {
    const eps = episodes.filter((e) => e.grain === g)
    const byZone = (pick: (e: Episode) => number): number[] =>
      [0, 1, 2].map((z) => eps.filter((e) => e.zone === z)
        .reduce((a, e) => a + pick(e), 0))
    grainRaw.push({
      grain: g,
      aliveAtStart: !traj.lost[g][startIdx],
      aliveAtEnd: !traj.lost[g][endIdx],
      liveExposureS: eps.reduce((a, e) => a + e.duration, 0),
      exposureByZone: byZone((e) => e.duration),
      contactByZone: byZone((e) => e.contactDuration),
      episodesByZone: byZone(() => 1),
      entriesByZone: byZone((e) => (e.leftCensored ? 0 : 1)),
      exitsByZone: byZone((e) => (e.endType === 'exit' ? 1 : 0)),
      impactsByZone: byZone((e) => e.impacts),
      detachesByZone: byZone((e) => e.detaches),
      endpointNodalDistanceM: traj.lost[g][endIdx]
        ? null : field.nodalDistance(traj.x[g][endIdx], traj.y[g][endIdx]),
    })
  }
  const aliveGrains = grainRaw.filter((g) => g.aliveAtStart)
  const nodeConcentration = {
    exposureS: concentration(aliveGrains.map((g) => g.exposureByZone[0])),
    episodes: concentration(aliveGrains.map((g) => g.episodesByZone[0])),
    entries: concentration(aliveGrains.map((g) => g.entriesByZone[0])),
    exits: concentration(aliveGrains.map((g) => g.exitsByZone[0])),
    impacts: concentration(aliveGrains.map((g) => g.impactsByZone[0])),
    detaches: concentration(aliveGrains.map((g) => g.detachesByZone[0])),
  }
  const fullNode = {
    occupancy: totalLiveExposure > 0 ? exposure[0] / totalLiveExposure : null,
    entryRate: totalLiveExposure - exposure[0] > 0
      ? entries[0] / (totalLiveExposure - exposure[0]) : null,
    escapeRate: exposure[0] > 0 ? exits[0] / exposure[0] : null,
  }
  const loo = aliveGrains.map((g) => {
    const live = totalLiveExposure - g.liveExposureS
    const nodeExp = exposure[0] - g.exposureByZone[0]
    return {
      occupancy: live > 0 ? nodeExp / live : null,
      entryRate: live - nodeExp > 0
        ? (entries[0] - g.entriesByZone[0]) / (live - nodeExp) : null,
      escapeRate: nodeExp > 0
        ? (exits[0] - g.exitsByZone[0]) / nodeExp : null,
    }
  })
  const range = (vals: (number | null)[]): (number | null)[] => {
    const nums = vals.filter((v): v is number => v !== null)
    return nums.length ? [round(Math.min(...nums)), round(Math.max(...nums))] : [null, null]
  }
  const maxChange = (vals: (number | null)[], full: number | null): number | null => {
    if (full === null) return null
    const nums = vals.filter((v): v is number => v !== null)
    return nums.length ? round(Math.max(...nums.map((v) => Math.abs(v - full)))) : null
  }
  const grains = grainRaw.map((g) => {
    const ini = initials[g.grain]
    const initialDistance = initialField.nodalDistance(ini.x, ini.y)
    return {
      grain: g.grain,
      initial: {
        x: round(ini.x, 6), y: round(ini.y, 6),
        radiusM: round(Math.hypot(ini.x, ini.y), 6),
        axis: ini.axis.map((v) => round(v, 6)),
        nodalDistanceM: round(initialDistance, 6),
        zone: ZONE_NAMES[zoneAt(initialField, ini.x, ini.y)],
      },
      aliveAtStart: g.aliveAtStart,
      aliveAtEnd: g.aliveAtEnd,
      lostAtInWindowS: round(lostAtInWindow[g.grain], 6),
      liveExposureS: round(g.liveExposureS, 6),
      exposureByZoneS: g.exposureByZone.map((v) => round(v, 6)),
      contactByZoneS: g.contactByZone.map((v) => round(v, 6)),
      episodesByZone: g.episodesByZone,
      entriesByZone: g.entriesByZone,
      exitsByZone: g.exitsByZone,
      impactsByZone: g.impactsByZone,
      detachesByZone: g.detachesByZone,
      node: {
        exposureShareOfLive: g.liveExposureS > 0
          ? round(g.exposureByZone[0] / g.liveExposureS) : null,
        exposureShareOfNodeTotal: exposure[0] > 0
          ? round(g.exposureByZone[0] / exposure[0]) : null,
        contactFraction: g.exposureByZone[0] > 0
          ? round(g.contactByZone[0] / g.exposureByZone[0]) : null,
        entryShare: entries[0] > 0 ? round(g.entriesByZone[0] / entries[0]) : null,
        exitShare: exits[0] > 0 ? round(g.exitsByZone[0] / exits[0]) : null,
      },
      endpointNodalDistanceM: round(g.endpointNodalDistanceM, 6),
      endpointWithin5mm: g.endpointNodalDistanceM === null
        ? null : g.endpointNodalDistanceM <= 0.005,
    }
  })

  return {
    window: [startIdx * dt, endIdx * dt],
    measurementDt: dt,
    totalLiveExposureS: round(totalLiveExposure, 3),
    lostGrainsInWindow,
    endpoint: {
      alive: aliveAtEnd,
      medianNodalDistanceM: round(median(endpointDistances), 6),
      within5mmAliveFraction: aliveAtEnd > 0 ? round(within5AtEnd / aliveAtEnd) : null,
      within5mmInitialFraction: round(within5AtEnd / traj.nGrains),
    },
    zones,
    nodeConcentration,
    leaveOneOutNode: {
      occupancyLiveRange: range(loo.map((v) => v.occupancy)),
      entryRateRange: range(loo.map((v) => v.entryRate)),
      escapeRateRange: range(loo.map((v) => v.escapeRate)),
      maxAbsChange: {
        occupancyLive: maxChange(loo.map((v) => v.occupancy), fullNode.occupancy),
        entryRate: maxChange(loo.map((v) => v.entryRate), fullNode.entryRate),
        escapeRate: maxChange(loo.map((v) => v.escapeRate), fullNode.escapeRate),
      },
    },
    grains,
  }
}

interface EnsembleDef {
  label: string
  seed: number | null
}

interface SettingDef {
  label: string
  plateDt: number
  contactSubsteps: number
}

const ENSEMBLES: EnsembleDef[] = [
  { label: 'E0', seed: null },
  { label: 'E1', seed: 110501 },
  { label: 'E2', seed: 110502 },
  { label: 'E3', seed: 110503 },
  { label: 'E4', seed: 110504 },
]

const SETTINGS: SettingDef[] = [
  { label: 'P1C4', plateDt: DT, contactSubsteps: 4 },
  { label: 'P2C4', plateDt: DT / 2, contactSubsteps: 4 },
]

function initialDiagnostics(
  initials: DumbbellInitial[], field: PlateField,
): Record<string, unknown> {
  const distances = initials.map((i) => field.nodalDistance(i.x, i.y))
  const counts = [0, 0, 0]
  initials.forEach((i) => { counts[zoneAt(field, i.x, i.y)]++ })
  return {
    medianNodalDistanceM: round(median(distances), 6),
    within5mmCount: distances.filter((d) => d <= 0.005).length,
    within5mmFraction: round(distances.filter((d) => d <= 0.005).length / distances.length),
    zoneCounts: {
      node: counts[0], intermediate: counts[1], antinode: counts[2],
    },
  }
}

function analyzeTrajectory(
  traj: Trajectory, plateA: number, initials: DumbbellInitial[],
): { reproduction: Record<string, unknown>; warmup: Record<string, unknown> } {
  const fields = plateFields()
  return {
    reproduction: analyzeWindow(
      traj, fields.repro, nullAreas(fields.repro, plateA),
      0, 1.5, initials, fields.repro),
    warmup: analyzeWindow(
      traj, fields.warm, nullAreas(fields.warm, plateA),
      T_WARM, T_TOTAL, initials, fields.repro),
  }
}

for (const ensemble of ENSEMBLES) {
  for (const setting of SETTINGS) {
    const caseLabel = `${ensemble.label}-${setting.label}`
    describe(`M6-31 ensemble contributions ${caseLabel}`, () => {
      benchIt(`ensemble ${caseLabel}`, () => {
        const initials = initialsForSeed(ensemble.seed)
        const { traj, finalStats, wallMs, plateA } = runTrajectory(
          setting.plateDt, setting.contactSubsteps, initials)
        const commonTraj = downsample(traj, Math.round(DT / setting.plateDt))
        const fields = plateFields()
        console.log(JSON.stringify({
          study: 'm631-ensemble',
          case: caseLabel,
          ensemble: ensemble.label,
          seed: ensemble.seed,
          setting: setting.label,
          drive: 'eccentric',
          aspect: 1.5,
          plateDt: setting.plateDt,
          contactSubsteps: setting.contactSubsteps,
          contactSubstepDt: setting.plateDt / setting.contactSubsteps,
          measurementDt: commonTraj.dt,
          wallMs,
          survivors: finalStats.filter((s) => !s.lost).length,
          initialDiagnostics: initialDiagnostics(initials, fields.repro),
          fieldThresholds: {
            repro: fields.repro.threshold, warm: fields.warm.threshold,
            q90Repro: fields.repro.q90, q90Warm: fields.warm.q90,
          },
          analyses: analyzeTrajectory(commonTraj, plateA, initials),
        }))
      }, 3600000)
    })
  }
}
