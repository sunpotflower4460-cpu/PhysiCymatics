// M6-13: standalone validation of the seated-grain Coulomb friction
// formulation (research module, not wired into the sand sim).
// Expectations are analytic (incline statics/kinematics, detach
// condition), not fitted to the implementation.
import { describe, expect, it } from 'vitest'
import {
  seatedAccel, stepSeated, seatedVz, tangentialSpeed,
  type SurfaceFn, type SeatedState,
} from '../seatedFriction'
import { G } from '../surface'

const flat: SurfaceFn = (_x, _y, _t, o) => {
  o.h = 0; o.hx = 0; o.hy = 0; o.hxx = 0; o.hxy = 0; o.hyy = 0
  o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
}

const incline = (theta: number): SurfaceFn => (x, _y, _t, o) => {
  const tn = Math.tan(theta)
  o.h = x * tn; o.hx = tn; o.hy = 0; o.hxx = 0; o.hxy = 0; o.hyy = 0
  o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
}

const inclineVibrating = (theta: number, amp: number, om: number): SurfaceFn => (x, _y, t, o) => {
  const tn = Math.tan(theta)
  o.h = x * tn + amp * Math.sin(om * t); o.hx = tn; o.hy = 0
  o.hxx = 0; o.hxy = 0; o.hyy = 0
  o.ht = amp * om * Math.cos(om * t); o.htt = -amp * om * om * Math.sin(om * t)
  o.htx = 0; o.hty = 0
}

describe('M6-13 seated friction formulation', () => {
  it('flat static plate: no spontaneous drift', () => {
    const st: SeatedState = { x: 0.01, y: -0.02, vx: 0, vy: 0 }
    for (let k = 0; k < 1000; k++) {
      expect(stepSeated(st, flat, { muS: 0.3, muK: 0.2 }, k * 1e-3, 1e-3)).toBe('stick')
    }
    expect(st.x).toBe(0.01); expect(st.y).toBe(-0.02)
  })

  it('static incline: stick boundary at tan(theta) = muS', () => {
    const th = 20 * Math.PI / 180 // tan = 0.364
    const stStick: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    expect(seatedAccel(stStick, incline(th), { muS: 0.4, muK: 0.2 }, 0).mode).toBe('stick')
    const stSlip: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    expect(seatedAccel(stSlip, incline(th), { muS: 0.3, muK: 0.2 }, 0).mode).toBe('slip')
    // Stuck grain stays exactly put.
    for (let k = 0; k < 500; k++) stepSeated(stStick, incline(th), { muS: 0.4, muK: 0.2 }, k * 1e-3, 1e-3)
    expect(stStick.x).toBe(0)
  })

  it('frictionless incline: ax = -g sin(theta)cos(theta) (down-slope)', () => {
    const th = 0.5
    const st: SeatedState = { x: 0, y: 0, vx: 0.01, vy: 0 }
    const a = seatedAccel(st, incline(th), { muS: 0, muK: 0 }, 0)
    expect(a.mode).toBe('slip')
    expect(a.ax).toBeCloseTo(-G * Math.sin(th) * Math.cos(th), 12)
    expect(a.ay).toBeCloseTo(0, 12)
  })

  it('kinetic friction incline: tangential accel = g(sin - muK cos)', () => {
    const th = 25 * Math.PI / 180; const muK = 0.2
    const st: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    const dt = 2e-4; const T = 0.4
    for (let t = 0; t < T; t += dt) stepSeated(st, incline(th), { muS: 0.1, muK }, t, dt)
    // Tangential speed along the slope after time T.
    const o = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
    incline(th)(st.x, st.y, 0, o)
    const vt = tangentialSpeed(st, o)
    const want = G * (Math.sin(th) - muK * Math.cos(th)) * T
    expect(Math.abs(vt - want) / want).toBeLessThan(2e-3)
    expect(st.x).toBeLessThan(0) // slid down-slope (-x for h = +x tan)
  })

  it('vibrating flat plate: detach exactly when htt < -g at rest', () => {
    const om = 2 * Math.PI * 50
    // amp*om^2 = 0.8 g: never detaches
    const aSmall = (0.8 * G) / (om * om)
    const st1: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    let detached = false
    for (let t = 0; t < 0.1; t += 1e-4) {
      if (stepSeated(st1, inclineVibrating(0, aSmall, om), { muS: 0.5, muK: 0.4 }, t, 1e-4) === 'detach') detached = true
    }
    expect(detached).toBe(false)
    // amp*om^2 = 1.5 g: detaches near sin(om t) = g/(1.5g) = 2/3 rising
    const aBig = (1.5 * G) / (om * om)
    const st2: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    let tDetach = -1
    for (let t = 0; t < 0.02; t += 2e-5) {
      if (stepSeated(st2, inclineVibrating(0, aBig, om), { muS: 0.5, muK: 0.4 }, t, 2e-5) === 'detach') { tDetach = t; break }
    }
    const want = Math.asin(2 / 3) / om
    expect(tDetach).toBeGreaterThan(0)
    expect(Math.abs(tDetach - want)).toBeLessThan(3e-5)
  })

  it('pure vertical vibration cannot break stick on a fixed slope', () => {
    // Required static friction and the normal force both scale with
    // (g + htt), so their ratio tan(theta)/muS is invariant: vertical
    // shaking alone never causes slip while attached (only detach).
    const th = 20 * Math.PI / 180
    const om = 2 * Math.PI * 30
    const amp = (0.9 * G) / (om * om)
    const st: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    for (let t = 0; t < 0.2; t += 2e-5) {
      expect(stepSeated(st, inclineVibrating(th, amp, om), { muS: 0.4, muK: 0.35 }, t, 2e-5)).toBe('stick')
    }
    expect(st.x).toBe(0)
  })

  it('oscillating slope causes stick/slip cycles and net creep', () => {
    // h = x * tan(theta0 + dth*sin(om t)): the local slope itself
    // oscillates across the static boundary atan(muS), which is the
    // mechanism a vibrating plate actually provides (grad h at a
    // fixed point is time-dependent).
    const th0 = 18 * Math.PI / 180; const dth = 8 * Math.PI / 180
    const om = 2 * Math.PI * 5 // slow enough that x*htt << g (no detach)
    const muS = Math.tan(22 * Math.PI / 180)
    const surf: SurfaceFn = (x, _y, t, o) => {
      const th = th0 + dth * Math.sin(om * t)
      const tn = Math.tan(th)
      const dT = (1 + tn * tn) * dth * om * Math.cos(om * t)
      const ddT = 2 * tn * (1 + tn * tn) * (dth * om * Math.cos(om * t)) ** 2
        - (1 + tn * tn) * dth * om * om * Math.sin(om * t)
      o.h = x * tn; o.hx = tn; o.hy = 0; o.hxx = 0; o.hxy = 0; o.hyy = 0
      o.ht = x * dT; o.htt = x * ddT; o.htx = dT; o.hty = 0
    }
    const st: SeatedState = { x: 0.01, y: 0, vx: 0, vy: 0 }
    const dt = 2e-5
    let sawSlip = false; let sawStick = false
    for (let t = 0; t < 0.6; t += dt) {
      const m = stepSeated(st, surf, { muS, muK: 0.35 }, t, dt)
      if (m === 'slip') sawSlip = true
      if (m === 'stick') sawStick = true
      if (m === 'detach') break
    }
    expect(sawStick).toBe(true)
    expect(sawSlip).toBe(true)
    expect(st.x).toBeLessThan(0.01 - 1e-4) // net down-slope creep
  })

  it('seatedVz matches the constraint velocity', () => {
    const th = 0.3
    const o = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
    incline(th)(0.01, 0, 0, o)
    const st: SeatedState = { x: 0.01, y: 0, vx: 0.2, vy: -0.1 }
    expect(seatedVz(st, o)).toBeCloseTo(Math.tan(th) * 0.2, 15)
  })
})
