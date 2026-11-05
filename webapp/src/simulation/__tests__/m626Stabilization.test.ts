// M6-26 stabilization regression tests. Each case reproduces a
// review counterexample from PRs #16/#18/#19/#22/#24/#25/#28 with
// expectations from analytic formulas or fine-step reference runs,
// not from the functions under test. Docs: docs/m6_26_recovery.md.
import { describe, expect, it } from 'vitest'
import { GrainCycleSim, type GrainCycleParams } from '../grainCycle'
import { TransientPlate } from '../plateTransient'
import { singleToneScore } from '../sequence'
import { calibratedForceN } from '../transientSand'
import { G } from '../surface'
import {
  stepSeated, type SurfaceFn, type SeatedState,
} from '../seatedFriction'
import {
  stepRolling, type RollingState, type RollingParams,
} from '../rollingSphere'
import {
  stepCurvedRolling, curvedProjectToRolling, type SurfaceFn3,
} from '../rollingSphereCurved'
import {
  dynCenterPartials, dynRollingAccel, dynProjectToRolling,
  stepDynRolling, type SurfaceFnDyn,
  type DynRollingParams,
} from '../rollingDynamic'

const DT = 1 / (44 * 152)
const base = { muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0 }

describe('M6-26 3-A top-surface geometry', () => {
  it('flat drop from 3mm top gap contacts at sqrt(2h/g) for H=0/1/5mm', () => {
    const expected = Math.sqrt(2 * 0.003 / G)
    for (const ht2 of [0, 5e-4, 2.5e-3]) {
      const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
      const sim = new GrainCycleSim(plate, {
        ...base, model: 'roll', plateHalfThickness: ht2,
      } as GrainCycleParams, [{ x: -0.05, y: 0, dropM: 0.003 }])
      let contactStep = -1
      for (let k = 0; k < Math.round(0.06 / DT); k++) {
        sim.step(DT)
        if (sim.stats()[0].impacts > 0) { contactStep = k; break }
      }
      expect(contactStep).toBeGreaterThan(0)
      // First impact must lie in the step interval containing the
      // analytic free-fall time (bisection locates it inside).
      expect(contactStep * DT).toBeLessThanOrEqual(expected + 1e-12)
      expect((contactStep + 1) * DT).toBeGreaterThanOrEqual(expected - 1e-12)
      const st = sim.stats()[0]
      expect(Math.abs(st.x - (-0.05))).toBeLessThan(1e-9)
    }
  })

  it('rolling solution satisfies differentiated no-slip with top surface (K term)', () => {
    const A0 = 0.002; const K = 2 * Math.PI / 0.24; const OM = 60
    const surf: SurfaceFnDyn = (x, y, t, o) => {
      const s = Math.sin(OM * t); const c = Math.cos(OM * t)
      const px = Math.cos(K * x); const py = Math.cos(K * y)
      const sx = Math.sin(K * x); const sy = Math.sin(K * y)
      o.h = A0 * s * (px + py)
      o.hx = -A0 * s * K * sx; o.hy = -A0 * s * K * sy
      o.hxx = -A0 * s * K * K * px; o.hyy = -A0 * s * K * K * py; o.hxy = 0
      o.hxxx = A0 * s * K ** 3 * sx; o.hyyy = A0 * s * K ** 3 * sy; o.hxxy = 0; o.hxyy = 0
      o.ht = A0 * OM * c * (px + py); o.htt = -A0 * OM * OM * s * (px + py)
      o.htx = -A0 * OM * c * K * sx; o.hty = -A0 * OM * c * K * sy
      o.htxx = -A0 * OM * c * K * K * px; o.htyy = -A0 * OM * c * K * K * py; o.htxy = 0
      o.httx = A0 * OM * OM * s * K * sx; o.htty = A0 * OM * OM * s * K * sy
    }
    const ht2 = 5e-3 // exaggerated: analytic surface, math check only
    const p: DynRollingParams = {
      muS: 10, muK: 5, radius: 1e-3, beta: 2 / 5, plateHalfThickness: ht2,
    }
    const st: RollingState = { x: 0.0317, y: -0.0213, vx: 0.08, vy: -0.05, wx: 0, wy: 0, wz: 0 }
    const t = 0.37
    dynProjectToRolling(st, surf, p, t)
    // Slip from its definition u = J qdot + R n_t - R (w x n).
    const slipOf = (s2: RollingState, tt: number): [number, number, number] => {
      const cp = dynCenterPartials(s2.x, s2.y, tt, surf, p.radius, ht2)
      const w = [s2.wx, s2.wy, s2.wz]
      const wxn = [
        w[1] * cp.n[2] - w[2] * cp.n[1],
        w[2] * cp.n[0] - w[0] * cp.n[2],
        w[0] * cp.n[1] - w[1] * cp.n[0],
      ]
      return [0, 1, 2].map((i) =>
        cp.D[0][i] * s2.vx + cp.D[1][i] * s2.vy
        + p.radius * cp.nD[2][i] - p.radius * wxn[i]) as [number, number, number]
    }
    const u0 = slipOf(st, t)
    // The linearized top surface (mid-surface normal) leaves a
    // small NORMAL component in J qdot of order ht2*|grad h||H||v|;
    // the no-slip constraint is on the tangential part. The normal
    // residual is measured here and must stay at that scale.
    const cp0 = dynCenterPartials(st.x, st.y, t, surf, p.radius, ht2)
    const un0 = u0[0] * cp0.n[0] + u0[1] * cp0.n[1] + u0[2] * cp0.n[2]
    const ut0 = u0.map((v, i) => v - un0 * cp0.n[i])
    expect(Math.hypot(...ut0)).toBeLessThan(1e-9)
    expect(Math.abs(un0)).toBeLessThan(1e-5)
    const a = dynRollingAccel(st, surf, p, t)
    expect(a.mode).toBe('roll')
    const eps = 1e-6
    const s1: RollingState = {
      x: st.x + st.vx * eps, y: st.y + st.vy * eps,
      vx: st.vx + a.ax * eps, vy: st.vy + a.ay * eps,
      wx: st.wx + a.dwx * eps, wy: st.wy + a.dwy * eps, wz: st.wz + a.dwz * eps,
    }
    const u1 = slipOf(s1, t + eps)
    const cp1 = dynCenterPartials(s1.x, s1.y, t + eps, surf, p.radius, ht2)
    const un1 = u1[0] * cp1.n[0] + u1[1] * cp1.n[1] + u1[2] * cp1.n[2]
    const ut1 = u1.map((v, i) => v - un1 * cp1.n[i])
    const udot = Math.hypot(...ut1.map((v, i) => (v - ut0[i]) / eps))
    // Missing K term would give ~ht2*|htxx*v| ~ 0.03 m/s^2 here.
    expect(udot).toBeLessThan(5e-3)
  })
})

describe('M6-26 3-B event handling', () => {
  it('seated: transverse velocity does not reverse on h=0.5x (fine-step reference)', () => {
    const surf: SurfaceFn = (x, _y, _t, o) => {
      o.h = 0.5 * x; o.hx = 0.5; o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0; o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
    }
    const p = { muS: 0.3, muK: 0.2 }
    const coarse: SeatedState = { x: 0, y: 0, vx: 0.01, vy: 0.01 }
    stepSeated(coarse, surf, p, 0, 0.01)
    const fine: SeatedState = { x: 0, y: 0, vx: 0.01, vy: 0.01 }
    for (let k = 0; k < 2000; k++) stepSeated(fine, surf, p, k * 5e-6, 5e-6)
    expect(fine.vy).toBeGreaterThan(0)
    expect(coarse.vy).toBeGreaterThan(0)
    expect(Math.abs(coarse.vy - fine.vy)).toBeLessThan(0.3 * Math.abs(fine.vy) + 1e-4)
    expect(Math.abs(coarse.vx - fine.vx)).toBeLessThan(2e-3)
  })

  it('seated slip: mid-step detach on h=A sin(wt) at analytic time', () => {
    const om = 2 * Math.PI * 100
    const A = 1.5 * G / (om * om)
    const surf: SurfaceFn = (_x, _y, t, o) => {
      o.h = A * Math.sin(om * t); o.hx = 0; o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0
      o.ht = A * om * Math.cos(om * t); o.htt = -A * om * om * Math.sin(om * t)
      o.htx = 0; o.hty = 0
    }
    const st: SeatedState = { x: 0, y: 0, vx: 0.1, vy: 0 }
    const r = stepSeated(st, surf, { muS: 0, muK: 0 }, 0, 0.01)
    expect(r.mode).toBe('detach')
    expect(Math.abs(r.advanced - Math.asin(2 / 3) / om)).toBeLessThan(2e-4)
  })

  it('plane rolling: pending detach survives false slip-zero rejection', () => {
    const om = 2 * Math.PI * 50
    const A = 2 * G / (om * om)
    const surf: SurfaceFn = (x, _y, t, o) => {
      o.h = 2 * x + A * Math.sin(om * t); o.hx = 2; o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0
      o.ht = A * om * Math.cos(om * t); o.htt = -A * om * om * Math.sin(om * t)
      o.htx = 0; o.hty = 0
    }
    const p: RollingParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5 }
    const st: RollingState = { x: 0, y: 0, vx: 0.001, vy: 0.001, wx: 0, wy: 0, wz: 0 }
    const r = stepRolling(st, surf, p, 0, 0.01)
    expect(r.mode).toBe('detach')
    expect(Math.abs(r.advanced - (Math.PI / 6) / om)).toBeLessThan(3e-4)
  })

  it('dynamic rolling: pending detach survives false slip-zero rejection', () => {
    const om = 2 * Math.PI * 50
    const A = 2 * G / (om * om)
    const surf: SurfaceFnDyn = (x, _y, t, o) => {
      o.h = 2 * x + A * Math.sin(om * t); o.hx = 2; o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0; o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
      o.ht = A * om * Math.cos(om * t); o.htt = -A * om * om * Math.sin(om * t)
      o.htx = 0; o.hty = 0; o.htxx = 0; o.htxy = 0; o.htyy = 0; o.httx = 0; o.htty = 0
    }
    const p: RollingParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5 }
    const st: RollingState = { x: 0, y: 0, vx: 0.001, vy: 0.001, wx: 0, wy: 0, wz: 0 }
    const r = stepDynRolling(st, surf, p, 0, 0.01)
    expect(r.mode).toBe('detach')
    expect(Math.abs(r.advanced - (Math.PI / 6) / om)).toBeLessThan(3e-4)
  })

  it('curved rolling: intermediate friction-margin violation is not accepted as rolling', () => {
    const surf: SurfaceFn3 = (x, _y, o) => {
      o.h = 0.001 * Math.sin(100 * x); o.hx = 0.1 * Math.cos(100 * x); o.hy = 0
      o.hxx = -10 * Math.sin(100 * x); o.hxy = 0; o.hyy = 0
      o.hxxx = -1000 * Math.cos(100 * x); o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
    }
    const p: RollingParams = { muS: 0.01, muK: 0.005, radius: 0.001, beta: 0.4 }
    const mk = (): RollingState => {
      const s: RollingState = { x: -Math.PI / 200, y: 0, vx: 0.3, vy: 0, wx: 0, wy: 0, wz: 0 }
      curvedProjectToRolling(s, surf, p)
      return s
    }
    const coarse = mk()
    const rc = stepCurvedRolling(coarse, surf, p, 0.11)
    const fine = mk()
    let fineSlipped = false
    for (let k = 0; k < 1100; k++) {
      const r = stepCurvedRolling(fine, surf, p, 1e-4)
      if (r.mode === 'slip') fineSlipped = true
    }
    expect(fineSlipped).toBe(true)
    // One coarse step must agree with the subdivided reference
    // (pure-rolling acceptance disagrees strongly).
    expect(rc.advanced).toBeCloseTo(0.11, 12)
    expect(Math.abs(coarse.x - fine.x)).toBeLessThan(5e-4)
    expect(Math.abs(coarse.vx - fine.vx)).toBeLessThan(0.01)
    // Spin is the most sensitive observable of the missed slip
    // episode (pure rolling gives wy~242 vs slipped ~283 rad/s).
    expect(Math.abs(coarse.wy - fine.wy)).toBeLessThan(15)
  })

  it('grainCycle: event-limit case consumes the full step (no dropped time)', () => {
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.06, 0])
    const sim = new GrainCycleSim(plate, {
      muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3, model: 'roll',
    } as GrainCycleParams, [{ x: 0.03, y: 0, dropM: 0.003 }])
    for (let k = 0; k < Math.round(0.06 / DT); k++) sim.step(DT)
    const st = sim.stats()[0]
    expect(st.unconsumedTotal).toBe(0)
    expect(st.failed).toBe(false)
    expect(st.t).toBeCloseTo(Math.round(0.06 / DT) * DT, 12)
    expect(Number.isFinite(st.x)).toBe(true)
  })
})

describe('M6-26 3-C numerics, pair collisions, accounting', () => {
  it('impact counter counts only impulses: flat drop e=0.3 gives 5 impacts', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
    const sim = new GrainCycleSim(plate, {
      ...base, model: 'roll', restitution: 0.3,
    } as GrainCycleParams, [{ x: -0.05, y: 0, dropM: 0.002 }])
    for (let k = 0; k < Math.round(0.3 / DT); k++) sim.step(DT)
    const st = sim.stats()[0]
    expect(st.impacts).toBe(5)
    expect(st.contactTime).toBeGreaterThan(0.2)
  })

  it('dragGamma -> 0 flight is continuous with the no-drag flight', () => {
    const firstImpactStep = (gam: number): number => {
      const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
      const sim = new GrainCycleSim(plate, {
        ...base, model: 'roll', dragGamma: gam,
      } as GrainCycleParams, [{ x: -0.05, y: 0, dropM: 0.003 }])
      for (let k = 0; k < Math.round(0.06 / DT); k++) {
        sim.step(DT)
        if (sim.stats()[0].impacts > 0) return k
      }
      return -1
    }
    const s0 = firstImpactStep(0)
    expect(s0).toBeGreaterThan(0)
    expect(firstImpactStep(1e-8)).toBe(s0)
    expect(firstImpactStep(1e-12)).toBe(s0)
  })

  it('pair separation displacement is accounted: flight+contact+correction = net', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
    const initials = [
      { x: -0.004, y: 0, dropM: 0.005, vx: 0.2 },
      { x: 0.004, y: 0, dropM: 0.005, vx: -0.2 },
    ]
    const sim = new GrainCycleSim(plate, {
      ...base, model: 'roll', restitution: 1, interparticle: true,
    } as GrainCycleParams, initials)
    for (let k = 0; k < Math.round(0.025 / DT); k++) sim.step(DT)
    const stats = sim.stats()
    expect(stats[0].pairCollisions).toBeGreaterThanOrEqual(1)
    stats.forEach((s, i) => {
      expect(
        s.netFlightX + s.netContactX + s.netCorrectionX,
      ).toBeCloseTo(s.x - initials[i].x, 9)
      expect(
        s.netFlightY + s.netContactY + s.netCorrectionY,
      ).toBeCloseTo(s.y - initials[i].y, 9)
    })
  })

  it('pair impulse into the plate resolves a plate impact before recapture', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
    const sim = new GrainCycleSim(plate, {
      ...base, model: 'roll', restitution: 1, interparticle: true,
    } as GrainCycleParams, [
      { x: 0, y: 0, dropM: 0.001 }, { x: 0.001, y: 0, dropM: 0.001 },
    ])
    sim.step(DT) // consume pending-drop initialization
    const grains = (sim as unknown as { grains: Record<string, unknown>[] }).grains
    grains[0].mode = 'contact'
    grains[0].roll = { x: 0, y: 0, vx: 0, vy: 0, wx: 0, wy: 0, wz: 0 }
    grains[1].mode = 'fly'
    grains[1].fly = {
      cx: 0.0002, cy: 0, cz: base.radius + 0.00035,
      vx: -0.05, vy: 0, vz: -0.25, wx: 0, wy: 0, wz: 0,
    }
    const impactsBefore = sim.stats()[0].impacts
    sim.step(DT)
    expect(sim.stats()[0].pairCollisions).toBeGreaterThanOrEqual(1)
    expect(sim.stats()[0].impacts).toBe(impactsBefore + 1)
    expect(grains[0].mode).toBe('fly') // rebounded, not silently recaptured
  })

  it('third derivatives near the origin approach the exact center limit', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
    const j = plate.dofs.findIndex(
      (d) => plate.modes[d.mode].m === 1 && d.kind === 'cos')
    expect(j).toBeGreaterThanOrEqual(0)
    plate.q[j] = 1e-4
    plate.beginStep(0, DT)
    const t3 = { hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0 }
    plate.surfaceThirdConsistent(0, 0, 0, t3)
    const center = t3.hxxx
    expect(Math.abs(center)).toBeGreaterThan(0.1)
    for (const r of [1e-11, 1e-10, 1e-8, 1e-6, 1e-5]) {
      const o = { hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0 }
      plate.surfaceThirdConsistent(r * Math.cos(0.9), r * Math.sin(0.9), 0, o)
      expect(Math.abs(o.hxxx - center) / Math.abs(center)).toBeLessThan(0.05)
    }
  })

  it('center coordinates are consistent across capture/detach (decomposition telescopes)', () => {
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.06, 0])
    const sim = new GrainCycleSim(plate, {
      muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3, model: 'roll',
    } as GrainCycleParams, [{ x: 0.03, y: 0, dropM: 0.003 }])
    for (let k = 0; k < Math.round(0.3 / DT); k++) sim.step(DT)
    const st = sim.stats()[0]
    expect(st.impacts).toBeGreaterThan(0)
    expect(
      st.netFlightX + st.netContactX + st.netCorrectionX,
    ).toBeCloseTo(st.x - 0.03, 9)
    expect(
      st.netFlightY + st.netContactY + st.netCorrectionY,
    ).toBeCloseTo(st.y - 0, 9)
  })
})

describe('M6-26 gate: center polynomial and zero-gap recollision', () => {
  it('center-interval derivatives match FD of the height function itself', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 30), [0, 0])
    let seed = 777
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
    plate.beginStep(0, DT)
    const tau = DT / 2
    const o5 = new Float64Array(5)
    const hAt = (x: number, y: number): number => {
      plate.surfAtConsistent(x, y, tau, o5); return o5[0]
    }
    const fld = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
    // Points inside the center radial interval (a/511 ~ 2.35e-4 m).
    for (const [x, y] of [[3e-5, -2e-5], [-1e-4, 5e-5], [8e-5, 9e-5]] as const) {
      plate.surfaceFieldConsistent(x, y, tau, fld)
      const e = 2e-6
      const hxFD = (hAt(x + e, y) - hAt(x - e, y)) / (2 * e)
      const hyFD = (hAt(x, y + e) - hAt(x, y - e)) / (2 * e)
      const hxxFD = (hAt(x + e, y) - 2 * hAt(x, y) + hAt(x - e, y)) / (e * e)
      const hyyFD = (hAt(x, y + e) - 2 * hAt(x, y) + hAt(x, y - e)) / (e * e)
      const hxyFD = (hAt(x + e, y + e) - hAt(x + e, y - e)
        - hAt(x - e, y + e) + hAt(x - e, y - e)) / (4 * e * e)
      expect(Math.abs(fld.hx - hxFD)).toBeLessThan(1e-6 + 1e-3 * Math.abs(hxFD))
      expect(Math.abs(fld.hy - hyFD)).toBeLessThan(1e-6 + 1e-3 * Math.abs(hyFD))
      expect(Math.abs(fld.hxx - hxxFD)).toBeLessThan(1e-3 + 2e-3 * Math.abs(hxxFD))
      expect(Math.abs(fld.hyy - hyyFD)).toBeLessThan(1e-3 + 2e-3 * Math.abs(hyyFD))
      expect(Math.abs(fld.hxy - hxyFD)).toBeLessThan(1e-3 + 2e-3 * Math.abs(hxyFD))
      // Third derivatives via FD of the (now polynomial) Hessian.
      const t3 = { hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0 }
      plate.surfaceThirdConsistent(x, y, tau, t3)
      const e3 = 5e-6
      const fp = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
      const fm = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
      plate.surfaceFieldConsistent(x + e3, y, tau, fp)
      plate.surfaceFieldConsistent(x - e3, y, tau, fm)
      const hxxxFD = (fp.hxx - fm.hxx) / (2 * e3)
      plate.surfaceFieldConsistent(x, y + e3, tau, fp)
      plate.surfaceFieldConsistent(x, y - e3, tau, fm)
      const hyyyFD = (fp.hyy - fm.hyy) / (2 * e3)
      expect(Math.abs(t3.hxxx - hxxxFD)).toBeLessThan(0.05 + 0.02 * Math.abs(hxxxFD))
      expect(Math.abs(t3.hyyy - hyyyFD)).toBeLessThan(0.05 + 0.02 * Math.abs(hyyyFD))
    }
  })

  it('separating from zero gap still recollides within the same step on an accelerating plate', () => {
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.06, 0])
    const sim = new GrainCycleSim(plate, {
      muS: 0.3, muK: 0.2, radius: 2.5e-4, beta: 2 / 5, restitution: 0.3, model: 'roll',
    } as GrainCycleParams, [{ x: 0.08, y: 0.05, dropM: 0.05 }])
    const o5 = new Float64Array(5)
    let placed = false
    for (let k = 0; k < 1500 && !placed; k++) {
      sim.step(DT)
      // Search for a point with large positive plate acceleration.
      let best = { x: 0, y: 0, wtt: 0, ht: 0, hx: 0, hy: 0 }
      for (let gx = -4; gx <= 4; gx++) for (let gy = -4; gy <= 4; gy++) {
        const x = gx * 0.015; const y = gy * 0.015
        plate.surfAtConsistent(x, y, sim.t, o5)
        if (o5[4] > best.wtt) best = { x, y, wtt: o5[4], ht: o5[3], hx: o5[1], hy: o5[2] }
      }
      if (k < 200 || best.wtt < 80) continue
      // Place the grain exactly at the offset center, separating
      // slowly (2 mm/s along the normal). With wtt > 80 the gap
      // must close again in ~2*u/(g+wtt) << DT.
      const surfDyn: SurfaceFnDyn = (x, y, t, o) => plate.surfaceFieldDynamicConsistent(x, y, t, o)
      const cp = dynCenterPartials(best.x, best.y, sim.t, surfDyn, 2.5e-4, 0)
      const u = 0.002
      const grains = (sim as unknown as { grains: Record<string, unknown>[] }).grains
      grains[0].mode = 'fly'
      grains[0].fly = {
        cx: cp.C[0], cy: cp.C[1], cz: cp.C[2],
        vx: u * cp.n[0], vy: u * cp.n[1], vz: best.ht + u * cp.n[2],
        wx: 0, wy: 0, wz: 0,
      }
      const grainsT = grains[0] as { stats?: unknown }
      void grainsT
      const before = sim.stats()[0]
      sim.step(DT)
      const after = sim.stats()[0]
      expect(after.impacts).toBeGreaterThan(before.impacts)
      expect(after.impactImpulseMagSum).toBeGreaterThan(before.impactImpulseMagSum)
      placed = true
    }
    expect(placed).toBe(true)
  })
})
