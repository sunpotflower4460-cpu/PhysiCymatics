// M6-19 RESEARCH ONLY — not wired into TransientSandSim.
// Flight -> impact -> contact (rolling or seated friction) ->
// detach cycle for a single grain on the consistent plate field,
// used to measure horizontal transport of rolling vs point
// grains (docs/m6_19_design.md). Contact search: 8 samples +
// bisection per plate segment (same shallow-contact limitation
// as production). Capture after impact when the rebound normal
// speed <= g*dt (exact when restitution = 0).
import { G } from './surface'
import type { TransientPlate } from './plateTransient'
import {
  stepDynRolling, dynCenterPartials, dynCenterVelocity,
  type SurfaceFnDyn,
} from './rollingDynamic'
import {
  sphereImpact, flightToRolling, rollingToFlight, solveContactParams,
  type FlightState, type ImpactParams,
} from './rollingImpact'
import {
  stepSeated, seatedVz,
  type SurfaceFn, type SurfaceField, type SeatedState,
} from './seatedFriction'
import type { RollingState } from './rollingSphere'

export type GrainModel = 'roll' | 'point'

export interface GrainCycleParams extends ImpactParams {
  model: GrainModel
  /** M6-21: linear air drag in flight (1/s), still air. 0 = off. */
  dragGamma?: number
  /** M6-21: viscous spin decay in flight (1/s). 0 = off. */
  spinGamma?: number
}

export interface GrainStats {
  t: number
  /** Center horizontal position (roll model) / point position. */
  x: number; y: number
  pathM: number
  impacts: number
  detaches: number
  contactTime: number
  lost: boolean
  // M6-20 transport decomposition (docs/m6_20_design.md).
  netFlightX: number; netFlightY: number
  netContactX: number; netContactY: number
  pathFlightM: number; pathContactM: number
  impactImpulseX: number; impactImpulseY: number
  /** Sum of per-impact |dVh| magnitudes (cancellation analysis). */
  impactImpulseMagSum: number
  contactImpulseX: number; contactImpulseY: number
  conversionImpulseX: number; conversionImpulseY: number
  vxNow: number; vyNow: number
}

interface Grain {
  mode: 'fly' | 'contact'
  fly: FlightState
  roll: RollingState
  seated: SeatedState
  stats: GrainStats
  lastHx: number; lastHy: number
  lastQ: [number, number] | null
}

const CONTACT_SAMPLES = 8

export class GrainCycleSim {
  private grains: Grain[] = []
  private pendingDrops: number[] | null = null
  private currentDt = 0
  t = 0
  private dynSurf: SurfaceFnDyn
  private surfFn: SurfaceFn
  private out5 = new Float64Array(5)
  private field: SurfaceField = {
    h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0,
  }

  constructor(
    readonly plate: TransientPlate,
    readonly params: GrainCycleParams,
    initials: { x: number; y: number; dropM: number; vx?: number; vy?: number }[],
  ) {
    this.dynSurf = (x, y, t, o) => plate.surfaceFieldDynamicConsistent(x, y, t, o)
    this.surfFn = (x, y, t, o) => plate.surfaceFieldConsistent(x, y, t, o)
    // Plate is assumed at rest at t=0 (h = 0), as in the driven
    // scores used by the studies; initial heights are set on the
    // first step()'s segment in initDrop().
    this.pendingDrops = initials.map((i) => i.dropM)
    for (const ini of initials) {
      const cz = params.model === 'roll' ? params.radius + ini.dropM : ini.dropM
      this.grains.push({
        mode: 'fly',
        fly: {
          cx: ini.x, cy: ini.y, cz,
          vx: ini.vx ?? 0, vy: ini.vy ?? 0, vz: 0, wx: 0, wy: 0, wz: 0,
        },
        roll: { x: ini.x, y: ini.y, vx: 0, vy: 0, wx: 0, wy: 0, wz: 0 },
        seated: { x: ini.x, y: ini.y, vx: 0, vy: 0 },
        stats: {
          t: 0, x: ini.x, y: ini.y, pathM: 0,
          impacts: 0, detaches: 0, contactTime: 0, lost: false,
          netFlightX: 0, netFlightY: 0, netContactX: 0, netContactY: 0,
          pathFlightM: 0, pathContactM: 0,
          impactImpulseX: 0, impactImpulseY: 0, impactImpulseMagSum: 0,
          contactImpulseX: 0, contactImpulseY: 0,
          conversionImpulseX: 0, conversionImpulseY: 0,
          vxNow: ini.vx ?? 0, vyNow: ini.vy ?? 0,
        },
        lastHx: ini.x, lastHy: ini.y, lastQ: null,
      })
    }
  }

  stats(): GrainStats[] {
    return this.grains.map((g) => ({ ...g.stats }))
  }

  /** Advance one plate segment [t, t+dt] (plate.beginStep inside). */
  step(dt: number): void {
    const t0 = this.t
    this.currentDt = dt
    this.plate.beginStep(t0, t0 + dt)
    if (this.pendingDrops) {
      this.grains.forEach((g, i) => {
        this.plate.surfAtConsistent(g.fly.cx, g.fly.cy, t0, this.out5)
        g.fly.cz = this.out5[0]
          + (this.params.model === 'roll' ? this.params.radius : 0)
          + (this.pendingDrops as number[])[i]
      })
      this.pendingDrops = null
    }
    for (const g of this.grains) {
      if (!g.stats.lost) this.advanceGrain(g, t0, dt)
      const [hx, hy] = this.horizontalOf(g, t0 + dt)
      g.stats.pathM += Math.hypot(hx - g.lastHx, hy - g.lastHy)
      g.lastHx = hx; g.lastHy = hy
      g.stats.x = hx; g.stats.y = hy; g.stats.t = t0 + dt
      const [vx, vy] = this.centerVelH(g, t0 + dt)
      g.stats.vxNow = vx; g.stats.vyNow = vy
      if (Math.hypot(hx, hy) > this.plate.a) g.stats.lost = true
    }
    this.t = t0 + dt
  }

  /** Center horizontal position (consistent across transitions). */
  private horizontalOf(g: Grain, tAbs: number): [number, number] {
    if (this.params.model === 'roll') {
      if (g.mode === 'fly') return [g.fly.cx, g.fly.cy]
      const cp = dynCenterPartials(g.roll.x, g.roll.y, tAbs, this.dynSurf, this.params.radius)
      return [cp.C[0], cp.C[1]]
    }
    return g.mode === 'fly' ? [g.fly.cx, g.fly.cy] : [g.seated.x, g.seated.y]
  }

  /** Center horizontal velocity in the current mode. */
  private centerVelH(g: Grain, tAbs: number): [number, number] {
    if (g.mode === 'fly') return [g.fly.vx, g.fly.vy]
    if (this.params.model === 'roll') {
      const v: [number, number, number] = [0, 0, 0]
      dynCenterVelocity(g.roll, this.dynSurf, this.params, tAbs, v)
      return [v[0], v[1]]
    }
    return [g.seated.vx, g.seated.vy]
  }

  private sphereGap(g: Grain, fly: FlightState, tAbs: number): number {
    const [qx, qy] = solveContactParams(
      fly.cx, fly.cy, this.dynSurf, this.params.radius, tAbs, g.lastQ ?? undefined)
    g.lastQ = [qx, qy]
    const cp = dynCenterPartials(qx, qy, tAbs, this.dynSurf, this.params.radius)
    return (fly.cx - cp.C[0]) * cp.n[0] + (fly.cy - cp.C[1]) * cp.n[1]
      + (fly.cz - cp.C[2]) * cp.n[2]
  }

  private ballisticAt(fly: FlightState, tau: number, out: FlightState): void {
    const gam = this.params.dragGamma ?? 0
    const gs = this.params.spinGamma ?? 0
    if (gam === 0) {
      out.cx = fly.cx + fly.vx * tau
      out.cy = fly.cy + fly.vy * tau
      out.cz = fly.cz + fly.vz * tau - 0.5 * G * tau * tau
      out.vx = fly.vx; out.vy = fly.vy; out.vz = fly.vz - G * tau
    } else {
      // Linear drag closed form (still air): v' = -gam v - g zhat.
      const e = Math.exp(-gam * tau)
      const k = (1 - e) / gam
      out.cx = fly.cx + fly.vx * k
      out.cy = fly.cy + fly.vy * k
      out.cz = fly.cz + (fly.vz + G / gam) * k - (G / gam) * tau
      out.vx = fly.vx * e; out.vy = fly.vy * e
      out.vz = (fly.vz + G / gam) * e - G / gam
    }
    const es = gs === 0 ? 1 : Math.exp(-gs * tau)
    out.wx = fly.wx * es; out.wy = fly.wy * es; out.wz = fly.wz * es
  }

  private advanceGrain(g: Grain, t0: number, dt: number): void {
    let rem = dt
    let tc = t0
    for (let iter = 0; iter < 32 && rem > 0; iter++) {
      if (g.mode === 'fly') {
        const [hx0, hy0] = this.horizontalOf(g, tc)
        const used = this.advanceFlight(g, tc, rem)
        rem -= used
        tc = t0 + (dt - rem)
        const [hx1, hy1] = this.horizontalOf(g, tc)
        g.stats.netFlightX += hx1 - hx0
        g.stats.netFlightY += hy1 - hy0
        g.stats.pathFlightM += Math.hypot(hx1 - hx0, hy1 - hy0)
      } else {
        const [hx0, hy0] = this.horizontalOf(g, tc)
        const [vx0, vy0] = this.centerVelH(g, tc)
        const before = rem
        rem -= this.advanceContact(g, tc, rem)
        g.stats.contactTime += before - rem
        tc = t0 + (dt - rem)
        const [hx1, hy1] = this.horizontalOf(g, tc)
        const [vx1, vy1] = this.centerVelH(g, tc)
        g.stats.netContactX += hx1 - hx0
        g.stats.netContactY += hy1 - hy0
        g.stats.pathContactM += Math.hypot(hx1 - hx0, hy1 - hy0)
        g.stats.contactImpulseX += vx1 - vx0
        g.stats.contactImpulseY += vy1 - vy0
      }
    }
  }

  /** Returns elapsed time consumed. */
  private advanceFlight(g: Grain, tc: number, rem: number): number {
    const p = this.params
    const fly = g.fly
    const probe: FlightState = { ...fly }
    const gapAt = (tau: number): number => {
      this.ballisticAt(fly, tau, probe)
      if (p.model === 'roll') return this.sphereGap(g, probe, tc + tau)
      this.plate.surfAtConsistent(probe.cx, probe.cy, tc + tau, this.out5)
      return probe.cz - this.out5[0]
    }
    let prevTau = 0
    let prevGap = gapAt(0)
    if (prevGap <= 0) {
      this.doImpact(g, tc)
      if (g.mode === 'contact') return 0
      // Separating despite the nonpositive gap (numeric overlap):
      // advance ballistically rather than loop forever.
      this.ballisticAt(fly, rem, fly)
      return rem
    }
    for (let k = 1; k <= CONTACT_SAMPLES; k++) {
      const tau = (rem * k) / CONTACT_SAMPLES
      const gap = gapAt(tau)
      if (gap <= 0) {
        let lo = prevTau; let hi = tau
        for (let it = 0; it < 50; it++) {
          const mid = 0.5 * (lo + hi)
          if (gapAt(mid) > 0) lo = mid; else hi = mid
        }
        this.ballisticAt(fly, hi, fly)
        this.doImpact(g, tc + hi)
        return hi
      }
      prevTau = tau; prevGap = gap
    }
    void prevGap
    this.ballisticAt(fly, rem, fly)
    return rem
  }

  private doImpact(g: Grain, tc: number): void {
    const p = this.params
    g.stats.impacts++
    if (p.model === 'roll') {
      const vxPre = g.fly.vx; const vyPre = g.fly.vy
      const r = sphereImpact(g.fly, this.dynSurf, p, tc)
      if (r.mode === 'separating') return
      g.stats.impactImpulseX += g.fly.vx - vxPre
      g.stats.impactImpulseY += g.fly.vy - vyPre
      g.stats.impactImpulseMagSum += Math.hypot(g.fly.vx - vxPre, g.fly.vy - vyPre)
      // Capture if rebound normal speed is negligible: recompute
      // normal relative speed after impact.
      const [qx, qy] = r.contact
      const cp = dynCenterPartials(qx, qy, tc, this.dynSurf, p.radius)
      this.plate.surfAtConsistent(qx, qy, tc, this.out5)
      const unPost = g.fly.vx * cp.n[0] + g.fly.vy * cp.n[1]
        + (g.fly.vz - this.out5[3]) * cp.n[2]
      // Capture threshold g*dt (docs/m6_19_design.md; exact at e=0).
      if (unPost <= G * this.currentDt) {
        const vxPre = g.fly.vx; const vyPre = g.fly.vy
        g.roll = flightToRolling(g.fly, this.dynSurf, p, tc)
        g.mode = 'contact'
        const vPost: [number, number, number] = [0, 0, 0]
        dynCenterVelocity(g.roll, this.dynSurf, p, tc, vPost)
        g.stats.conversionImpulseX += vPost[0] - vxPre
        g.stats.conversionImpulseY += vPost[1] - vyPre
      }
      return
    }
    // Point impact on the consistent surface.
    const fly = g.fly
    this.plate.surfaceFieldConsistent(fly.cx, fly.cy, tc, this.field)
    const s = this.field
    const sn = Math.hypot(s.hx, s.hy, 1)
    const nx = -s.hx / sn; const ny = -s.hy / sn; const nz = 1 / sn
    const ux = fly.vx; const uy = fly.vy; const uz = fly.vz - s.ht
    const un = ux * nx + uy * ny + uz * nz
    if (un >= 0) return
    const Jn = -(1 + p.restitution) * un
    const utx = ux - un * nx; const uty = uy - un * ny; const utz = uz - un * nz
    const uts = Math.hypot(utx, uty, utz)
    let jfx = 0; let jfy = 0; let jfz = 0
    if (uts > 0) {
      const mag = Math.min(uts, (uts <= p.muS * Jn ? uts : p.muK * Jn))
      jfx = -mag * utx / uts; jfy = -mag * uty / uts; jfz = -mag * utz / uts
    }
    fly.vx += Jn * nx + jfx
    fly.vy += Jn * ny + jfy
    fly.vz += Jn * nz + jfz
    g.stats.impactImpulseX += Jn * nx + jfx
    g.stats.impactImpulseY += Jn * ny + jfy
    g.stats.impactImpulseMagSum += Math.hypot(Jn * nx + jfx, Jn * ny + jfy)
    const unPost = fly.vx * nx + fly.vy * ny + (fly.vz - s.ht) * nz
    if (unPost <= G * this.currentDt) {
      g.seated = { x: fly.cx, y: fly.cy, vx: fly.vx, vy: fly.vy }
      g.mode = 'contact'
    }
  }

  /** Returns elapsed time consumed in contact. */
  private advanceContact(g: Grain, tc: number, rem: number): number {
    const p = this.params
    if (p.model === 'roll') {
      const r = stepDynRolling(g.roll, this.dynSurf, p, tc, rem)
      if (r.mode === 'detach') {
        g.stats.detaches++
        g.fly = rollingToFlight(g.roll, this.dynSurf, p, tc + r.advanced)
        g.mode = 'fly'
      }
      return r.advanced
    }
    const r = stepSeated(g.seated, this.surfFn, p, tc, rem)
    if (r.mode === 'detach') {
      g.stats.detaches++
      const at = tc + r.advanced
      this.plate.surfaceFieldConsistent(g.seated.x, g.seated.y, at, this.field)
      g.fly = {
        cx: g.seated.x, cy: g.seated.y, cz: this.field.h,
        vx: g.seated.vx, vy: g.seated.vy,
        vz: seatedVz(g.seated, this.field),
        wx: 0, wy: 0, wz: 0,
      }
      g.mode = 'fly'
    }
    return r.advanced
  }
}
