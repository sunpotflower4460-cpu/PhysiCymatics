// M6-27B RESEARCH ONLY — not wired into the production sim.
// Dumbbell flight -> impact -> sustained contact -> detach
// cycle on the consistent plate surface.
// docs/m6_27b_design.md (impulse time-stepping floors are
// stated there; this is NOT an exact constrained ODE).
import type { TransientPlate } from './plateTransient'
import { G } from './surface'
import {
  dumbbellFlightStep, dumbbellImpact, lobeCenters, contactVelocity,
  type DumbbellParams, type DumbbellState,
} from './dumbbell'
import type { SurfaceField } from './seatedFriction'

export interface DumbbellCycleParams extends DumbbellParams {
  plateHalfThickness?: number
  diagFreqHz?: number
  /** Research only: sustained-contact substeps per plate step.
   *  Omit for the established default (4). */
  contactSubsteps?: number
}

export interface DumbbellInitial {
  x: number; y: number; dropM: number
  vx?: number; vy?: number
  axis: [number, number, number]
  wx?: number; wy?: number; wz?: number
}

export interface DumbbellStats {
  x: number; y: number; t: number
  vxNow: number; vyNow: number
  czNow: number; vzNow: number
  ex: number; ey: number; ez: number
  spinMag: number
  /** Grain mode at the observation step end (measurement only). */
  modeNow: 'fly' | 'contact' | 'lost'
  pathM: number
  impacts: number
  contactImpulses: number
  detaches: number
  contactTime: number
  lost: boolean
  netFlightX: number; netFlightY: number
  netContactX: number; netContactY: number
  netCorrectionX: number; netCorrectionY: number
  impactImpulseMagSum: number
  phaseCosW: number; phaseSinW: number
  slopeDotSum: number; slopeMagWSum: number
  unconsumedTotal: number
  failed: boolean
}

type V3 = [number, number, number]
const CONTACT_SAMPLES = 8
const CONTACT_SUBSTEPS = 4
const GAP_TOL = 1e-7

interface LobeContact {
  gap: number
  point: V3
  n: V3
  vs: V3
  q: [number, number]
  center: V3
}

interface Grain {
  st: DumbbellState
  mode: 'fly' | 'contact' | 'lost'
  lastQ: [([number, number] | undefined), ([number, number] | undefined)]
  stats: DumbbellStats
  lastX: number; lastY: number
}

const fc: SurfaceField = {
  h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0,
}

export class DumbbellCycleSim {
  t = 0
  private grains: Grain[] = []
  private currentDt = 0
  private pendingDrops: number[] | null
  private readonly contactSubsteps: number

  constructor(
    private plate: TransientPlate,
    private params: DumbbellCycleParams,
    initials: DumbbellInitial[],
  ) {
    const contactSubsteps = params.contactSubsteps ?? CONTACT_SUBSTEPS
    if (!Number.isInteger(contactSubsteps) || contactSubsteps <= 0) {
      throw new RangeError('DumbbellCycleSim: contactSubsteps must be a positive integer')
    }
    this.contactSubsteps = contactSubsteps
    this.pendingDrops = initials.map((i) => i.dropM)
    this.grains = initials.map((ini) => {
      const em = Math.hypot(...ini.axis) || 1
      const st: DumbbellState = {
        cx: ini.x, cy: ini.y, cz: params.lobeRadius + ini.dropM,
        vx: ini.vx ?? 0, vy: ini.vy ?? 0, vz: 0,
        ex: ini.axis[0] / em, ey: ini.axis[1] / em, ez: ini.axis[2] / em,
        wx: ini.wx ?? 0, wy: ini.wy ?? 0, wz: ini.wz ?? 0,
      }
      return {
        st, mode: 'fly' as const, lastQ: [undefined, undefined],
        lastX: ini.x, lastY: ini.y,
        stats: {
          x: ini.x, y: ini.y, t: 0, vxNow: st.vx, vyNow: st.vy, czNow: st.cz, vzNow: st.vz,
          ex: st.ex, ey: st.ey, ez: st.ez,
          spinMag: Math.hypot(st.wx, st.wy, st.wz), modeNow: 'fly',
          pathM: 0, impacts: 0, contactImpulses: 0, detaches: 0,
          contactTime: 0, lost: false,
          netFlightX: 0, netFlightY: 0, netContactX: 0, netContactY: 0,
          netCorrectionX: 0, netCorrectionY: 0,
          impactImpulseMagSum: 0, phaseCosW: 0, phaseSinW: 0,
          slopeDotSum: 0, slopeMagWSum: 0,
          unconsumedTotal: 0, failed: false,
        },
      }
    })
  }

  stats(): DumbbellStats[] {
    return this.grains.map((g) => ({ ...g.stats }))
  }

  step(dt: number): void {
    const t0 = this.t
    this.currentDt = dt
    this.plate.beginStep(t0, t0 + dt)
    if (this.pendingDrops) {
      const drops = this.pendingDrops
      this.grains.forEach((g, i) => {
        // Place so the minimum lobe gap equals dropM at t0
        // (vertical shifts, iterated against the solved gap).
        g.st.cz = drops[i] + this.params.lobeRadius
          + this.params.halfSeparation * Math.abs(g.st.ez)
        for (let it = 0; it < 3; it++) {
          const infos = this.lobeContacts(g, g.st, t0)
          const min = infos[0].gap <= infos[1].gap ? infos[0] : infos[1]
          g.st.cz += (drops[i] - min.gap) / Math.max(0.2, min.n[2])
        }
      })
      this.pendingDrops = null
    }
    for (const g of this.grains) {
      if (g.mode !== 'lost') this.advanceGrain(g, t0, dt, 0)
    }
    this.t = t0 + dt
    for (const g of this.grains) {
      const s = g.stats
      s.x = g.st.cx; s.y = g.st.cy; s.t = this.t
      s.vxNow = g.st.vx; s.vyNow = g.st.vy; s.czNow = g.st.cz; s.vzNow = g.st.vz
      s.ex = g.st.ex; s.ey = g.st.ey; s.ez = g.st.ez
      s.spinMag = Math.hypot(g.st.wx, g.st.wy, g.st.wz)
      s.modeNow = g.mode
      s.pathM += Math.hypot(g.st.cx - g.lastX, g.st.cy - g.lastY)
      g.lastX = g.st.cx; g.lastY = g.st.cy
      if (Math.hypot(g.st.cx, g.st.cy) > this.plate.a) {
        g.mode = 'lost'
        s.lost = true
      }
    }
  }

  /** Lightweight offset-center geometry at contact params q:
   *  only first/second spatial derivatives and first time
   *  derivatives are needed for C(q), its Jacobian, n and the
   *  material velocity. Same definitions as dynCenterPartials
   *  (docs/m6_25_design.md) without the third-derivative field,
   *  which contact detection does not use. */
  private lightGeom(qx: number, qy: number, t: number): {
    C: V3; n: V3; vs: V3; a: number; b: number; c: number; d: number
  } {
    const ht2 = this.params.plateHalfThickness ?? 0
    const R = this.params.lobeRadius
    this.plate.surfaceFieldConsistent(qx, qy, t, fc)
    const m: V3 = [-fc.hx, -fc.hy, 1]
    const sn = Math.hypot(m[0], m[1], m[2])
    const n: V3 = [m[0] / sn, m[1] / sn, m[2] / sn]
    const nDeriv = (md: V3): V3 => {
      const snd = (m[0] * md[0] + m[1] * md[1] + m[2] * md[2]) / sn
      return [
        md[0] / sn - n[0] * snd / sn,
        md[1] / sn - n[1] * snd / sn,
        md[2] / sn - n[2] * snd / sn,
      ]
    }
    const nx = nDeriv([-fc.hxx, -fc.hxy, 0])
    const ny = nDeriv([-fc.hxy, -fc.hyy, 0])
    const Sx: V3 = [1 - ht2 * fc.hxx, -ht2 * fc.hxy, fc.hx]
    const Sy: V3 = [-ht2 * fc.hxy, 1 - ht2 * fc.hyy, fc.hy]
    const Cx: V3 = [Sx[0] + R * nx[0], Sx[1] + R * nx[1], Sx[2] + R * nx[2]]
    const Cy: V3 = [Sy[0] + R * ny[0], Sy[1] + R * ny[1], Sy[2] + R * ny[2]]
    const C: V3 = [
      qx - ht2 * fc.hx + R * n[0],
      qy - ht2 * fc.hy + R * n[1],
      fc.h + ht2 + R * n[2],
    ]
    return {
      C, n, a: Cx[0], b: Cy[0], c: Cx[1], d: Cy[1],
      vs: [-ht2 * fc.htx, -ht2 * fc.hty, fc.ht],
    }
  }

  private lobeContacts(g: Grain, st: DumbbellState, t: number): [LobeContact, LobeContact] {
    const centers = lobeCenters(st, this.params)
    const R = this.params.lobeRadius
    const one = (idx: 0 | 1): LobeContact => {
      const c = centers[idx]
      let qx = g.lastQ[idx] ? g.lastQ[idx][0] : c[0]
      let qy = g.lastQ[idx] ? g.lastQ[idx][1] : c[1]
      let geom = this.lightGeom(qx, qy, t)
      for (let it = 0; it < 12; it++) {
        const ex = geom.C[0] - c[0]; const ey = geom.C[1] - c[1]
        if (Math.hypot(ex, ey) < 1e-12) break
        const det = geom.a * geom.d - geom.b * geom.c
        qx -= (geom.d * ex - geom.b * ey) / det
        qy -= (-geom.c * ex + geom.a * ey) / det
        geom = this.lightGeom(qx, qy, t)
      }
      const q: [number, number] = [qx, qy]
      g.lastQ[idx] = q
      const gap = (c[0] - geom.C[0]) * geom.n[0]
        + (c[1] - geom.C[1]) * geom.n[1]
        + (c[2] - geom.C[2]) * geom.n[2]
      return {
        gap, q, center: c, n: geom.n,
        point: [geom.C[0] - R * geom.n[0], geom.C[1] - R * geom.n[1], geom.C[2] - R * geom.n[2]],
        vs: geom.vs,
      }
    }
    return [one(0), one(1)]
  }

  private minGap(g: Grain, st: DumbbellState, t: number): number {
    const infos = this.lobeContacts(g, st, t)
    return Math.min(infos[0].gap, infos[1].gap)
  }

  private flightStateAt(start: DumbbellState, tau: number): DumbbellState {
    const out: DumbbellState = { ...start }
    if (tau > 0) dumbbellFlightStep(out, this.params, tau)
    return out
  }

  private recordDiag(g: Grain, tc: number, dvx: number, dvy: number, q: [number, number]): void {
    const mag = Math.hypot(dvx, dvy)
    g.stats.impactImpulseMagSum += mag
    if (this.params.diagFreqHz) {
      const ph = 2 * Math.PI * this.params.diagFreqHz * tc
      g.stats.phaseCosW += mag * Math.cos(ph)
      g.stats.phaseSinW += mag * Math.sin(ph)
    }
    this.plate.surfaceFieldConsistent(q[0], q[1], tc, fc)
    g.stats.slopeDotSum += fc.hx * dvx + fc.hy * dvy
    g.stats.slopeMagWSum += Math.hypot(fc.hx, fc.hy) * mag
  }

  /** Impact the deepest lobe at the grain's current state.
   *  Returns 'captured' | 'rebounded' | 'separating'. */
  private doFlightImpact(g: Grain, tc: number): 'captured' | 'rebounded' | 'separating' {
    const infos = this.lobeContacts(g, g.st, tc)
    const info = infos[0].gap <= infos[1].gap ? infos[0] : infos[1]
    const vxPre = g.st.vx; const vyPre = g.st.vy
    const r = dumbbellImpact(g.st, info.point, info.n, info.vs, this.params)
    if (r.mode === 'separating') return 'separating'
    g.stats.impacts++
    this.recordDiag(g, tc, g.st.vx - vxPre, g.st.vy - vyPre, info.q)
    const uPost = contactVelocity(g.st, info.point, info.vs)
    const unPost = uPost[0] * info.n[0] + uPost[1] * info.n[1] + uPost[2] * info.n[2]
    if (unPost <= G * this.currentDt) {
      g.mode = 'contact'
      return 'captured'
    }
    return 'rebounded'
  }

  /** Returns elapsed time consumed in flight. */
  private advanceFlight(g: Grain, tc: number, rem: number): number {
    const start: DumbbellState = { ...g.st }
    const gapAt = (tau: number): number =>
      this.minGap(g, this.flightStateAt(start, tau), tc + tau)
    const setFlight = (tau: number): void => {
      Object.assign(g.st, this.flightStateAt(start, tau))
    }
    let prevTau = 0
    let prevGap = gapAt(0)
    if (prevGap <= 0) {
      const res = this.doFlightImpact(g, tc)
      if (res === 'captured') return 0
      // Separating (or rebounded) from zero gap: scan the
      // remainder for a fresh crossing (M6-26 discipline).
      const start2: DumbbellState = { ...g.st }
      const gapAt2 = (tau: number): number =>
        this.minGap(g, this.flightStateAt(start2, tau), tc + tau)
      let lastPos = -1
      for (let k = 1; k <= CONTACT_SAMPLES; k++) {
        const tau = (rem * k) / CONTACT_SAMPLES
        const gap = gapAt2(tau)
        if (gap > 0) {
          lastPos = tau
        } else if (lastPos >= 0) {
          let lo = lastPos; let hi = tau
          for (let it = 0; it < 50; it++) {
            const mid = 0.5 * (lo + hi)
            if (gapAt2(mid) > 0) lo = mid; else hi = mid
          }
          Object.assign(g.st, this.flightStateAt(start2, hi))
          this.doFlightImpact(g, tc + hi)
          return hi
        }
      }
      Object.assign(g.st, this.flightStateAt(start2, rem))
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
        setFlight(hi)
        this.doFlightImpact(g, tc + hi)
        return hi
      }
      prevTau = tau
      prevGap = gap
    }
    void prevGap
    setFlight(rem)
    return rem
  }

  /** Sustained-contact advance by impulse time stepping
   *  (docs/m6_27b_design.md). Returns elapsed consumed;
   *  a detach transitions to flight at a substep boundary. */
  private advanceContact(g: Grain, tc: number, rem: number): number {
    const h = rem / this.contactSubsteps
    const e0Params: DumbbellCycleParams = { ...this.params, restitution: 0 }
    let elapsed = 0
    for (let sub = 0; sub < this.contactSubsteps; sub++) {
      const tSub = tc + elapsed
      const infos0 = this.lobeContacts(g, g.st, tSub)
      const activePrev: [boolean, boolean] = [
        infos0[0].gap <= GAP_TOL, infos0[1].gap <= GAP_TOL,
      ]
      dumbbellFlightStep(g.st, this.params, h)
      elapsed += h
      const tEnd = tc + elapsed
      let infos = this.lobeContacts(g, g.st, tEnd)
      // Detach: both lobes clear and separating.
      const separating = infos.every((info) => {
        if (info.gap <= GAP_TOL) return false
        const u = contactVelocity(g.st, info.point, info.vs)
        return u[0] * info.n[0] + u[1] * info.n[1] + u[2] * info.n[2] >= -1e-12
      })
      if (separating) {
        g.mode = 'fly'
        g.stats.detaches++
        return elapsed
      }
      // Positional correction for the deepest penetration
      // (accounted separately from physical transport).
      const deep = infos[0].gap <= infos[1].gap ? infos[0] : infos[1]
      if (deep.gap < 0) {
        const pen = -deep.gap
        g.st.cx += deep.n[0] * pen
        g.st.cy += deep.n[1] * pen
        g.st.cz += deep.n[2] * pen
        g.stats.netCorrectionX += deep.n[0] * pen
        g.stats.netCorrectionY += deep.n[1] * pen
        infos = this.lobeContacts(g, g.st, tEnd)
      }
      // Sequential contact impulses (deepest first, iterated).
      for (let round = 0; round < 8; round++) {
        let applied = false
        const order = infos[0].gap <= infos[1].gap ? [0, 1] : [1, 0]
        for (const idx of order) {
          const info = infos[idx]
          if (info.gap > GAP_TOL) continue
          const u = contactVelocity(g.st, info.point, info.vs)
          const un = u[0] * info.n[0] + u[1] * info.n[1] + u[2] * info.n[2]
          if (un >= -1e-12) continue
          const isStrike = !activePrev[idx] && un < -2 * G * h
          const vxPre = g.st.vx; const vyPre = g.st.vy
          const r = dumbbellImpact(
            g.st, info.point, info.n, info.vs, isStrike ? this.params : e0Params)
          if (r.mode === 'separating') continue
          applied = true
          g.stats.contactImpulses++
          if (isStrike) {
            g.stats.impacts++
            this.recordDiag(g, tEnd, g.st.vx - vxPre, g.st.vy - vyPre, info.q)
          }
        }
        if (!applied) break
      }
    }
    return rem
  }

  private advanceGrain(g: Grain, t0: number, dt: number, depth: number): void {
    let rem = dt
    let tc = t0
    for (let iter = 0; iter < 32 && rem > 0; iter++) {
      if (g.mode === 'fly') {
        const x0 = g.st.cx; const y0 = g.st.cy
        const used = this.advanceFlight(g, tc, rem)
        rem -= used
        tc = t0 + (dt - rem)
        g.stats.netFlightX += g.st.cx - x0
        g.stats.netFlightY += g.st.cy - y0
      } else if (g.mode === 'contact') {
        const x0 = g.st.cx; const y0 = g.st.cy
        const corrX0 = g.stats.netCorrectionX; const corrY0 = g.stats.netCorrectionY
        const before = rem
        rem -= this.advanceContact(g, tc, rem)
        g.stats.contactTime += before - rem
        tc = t0 + (dt - rem)
        // Positional corrections inside advanceContact are
        // accounted in netCorrection, not as contact transport.
        g.stats.netContactX += (g.st.cx - x0) - (g.stats.netCorrectionX - corrX0)
        g.stats.netContactY += (g.st.cy - y0) - (g.stats.netCorrectionY - corrY0)
      } else {
        return
      }
    }
    if (rem > 0) {
      if (depth < 8) {
        const half = rem / 2
        this.advanceGrain(g, tc, half, depth + 1)
        this.advanceGrain(g, tc + half, rem - half, depth + 1)
        return
      }
      g.stats.unconsumedTotal += rem
      g.stats.failed = true
      throw new Error(
        `DumbbellCycleSim: event budget exhausted at t=${tc} with ${rem}s unconsumed`)
    }
  }
}
