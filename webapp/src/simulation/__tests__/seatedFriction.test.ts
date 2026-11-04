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
      expect(stepSeated(st, flat, { muS: 0.3, muK: 0.2 }, k * 1e-3, 1e-3).mode).toBe('stick')
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
      if (stepSeated(st1, inclineVibrating(0, aSmall, om), { muS: 0.5, muK: 0.4 }, t, 1e-4).mode === 'detach') detached = true
    }
    expect(detached).toBe(false)
    // amp*om^2 = 1.5 g: detaches near sin(om t) = g/(1.5g) = 2/3 rising
    const aBig = (1.5 * G) / (om * om)
    const st2: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    let tDetach = -1
    for (let t = 0; t < 0.02; t += 2e-5) {
      const r2 = stepSeated(st2, inclineVibrating(0, aBig, om), { muS: 0.5, muK: 0.4 }, t, 2e-5); if (r2.mode === 'detach') { tDetach = t + r2.advanced; break }
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
      expect(stepSeated(st, inclineVibrating(th, amp, om), { muS: 0.4, muK: 0.35 }, t, 2e-5).mode).toBe('stick')
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
      const m = stepSeated(st, surf, { muS, muK: 0.35 }, t, dt).mode
      if (m === 'slip') sawSlip = true
      if (m === 'stick') sawStick = true
      if (m === 'detach') break
    }
    expect(sawStick).toBe(true)
    expect(sawSlip).toBe(true)
    expect(st.x).toBeLessThan(0.01 - 1e-4) // net down-slope creep
  })

  it('M6-13R dynamic surface: h=x(1+t), vx=1 gives ax=-5.905 (mixed derivative)', () => {
    // h = x(1+t): hx=1+t, ht=x, htt=0, htx=1, H=0. At t=0, vx=1:
    // B = 2*htx*vx = 2; M a = -grad h (g+B) -> 2 ax = -(g+2).
    const surf: SurfaceFn = (x, _y, t, o) => {
      o.h = x * (1 + t); o.hx = 1 + t; o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0
      o.ht = x; o.htt = 0; o.htx = 1; o.hty = 0
    }
    const st: SeatedState = { x: 0.3, y: 0, vx: 1, vy: 0 }
    const a = seatedAccel(st, surf, { muS: 0, muK: 0 }, 0)
    expect(a.ax).toBeCloseTo(-(G + 2) / 2, 9) // -5.905
    expect(a.ax).toBeCloseTo(-5.905, 3)
    // Normal force from the same constraint B: (B+g)/sqrt(1+hx^2).
    expect(a.nPerMass).toBeCloseTo((G + 2) / Math.SQRT2, 9)
  })

  it('M6-13R dynamic curved surface: acceleration satisfies the constraint identity', () => {
    // h = 0.5*k*(x^2+y^2) * (1 + 0.3 sin(w t)): check the computed
    // acceleration against an independent finite-difference of the
    // constrained trajectory at small dt (frictionless), verifying
    // curvature + mixed-derivative terms together.
    const k = 40; const om = 30
    const surf: SurfaceFn = (x, y, t, o) => {
      const c = 1 + 0.3 * Math.sin(om * t)
      const cd = 0.3 * om * Math.cos(om * t)
      const cdd = -0.3 * om * om * Math.sin(om * t)
      o.h = 0.5 * k * c * (x * x + y * y)
      o.hx = k * c * x; o.hy = k * c * y
      o.hxx = k * c; o.hxy = 0; o.hyy = k * c
      o.ht = 0.5 * k * cd * (x * x + y * y)
      o.htt = 0.5 * k * cdd * (x * x + y * y)
      o.htx = k * cd * x; o.hty = k * cd * y
    }
    const st: SeatedState = { x: 0.02, y: -0.015, vx: 0.1, vy: 0.07 }
    const a = seatedAccel(st, surf, { muS: 0, muK: 0 }, 0.01)
    // Integrate 2 ms with tiny steps as the reference trajectory and
    // compare the initial acceleration by finite difference.
    const ref: SeatedState = { ...st }
    const h = 1e-6
    let prevVx = st.vx
    stepSeated(ref, surf, { muS: 0, muK: 0 }, 0.01, h)
    const fdAx = (ref.vx - prevVx) / h
    prevVx = ref.vx
    stepSeated(ref, surf, { muS: 0, muK: 0 }, 0.01 + h, h)
    const fdAx2 = (ref.vx - prevVx) / h
    // RK4 self-consistency is weak evidence alone, so also verify the
    // algebraic identity M a + grad h (g+B) = Qf (=0) directly.
    const o = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
    surf(st.x, st.y, 0.01, o)
    const B = st.vx * (o.hxx * st.vx + o.hxy * st.vy) + st.vy * (o.hxy * st.vx + o.hyy * st.vy)
      + 2 * (o.htx * st.vx + o.hty * st.vy) + o.htt
    const lhsX = a.ax + o.hx * (o.hx * a.ax + o.hy * a.ay) + o.hx * (G + B)
    const lhsY = a.ay + o.hy * (o.hx * a.ax + o.hy * a.ay) + o.hy * (G + B)
    expect(Math.abs(lhsX)).toBeLessThan(1e-9)
    expect(Math.abs(lhsY)).toBeLessThan(1e-9)
    expect(Math.abs(fdAx - a.ax)).toBeLessThan(Math.abs(fdAx2 - fdAx) + 0.05)
  })

  it('M6-13R stopping event: exact stop distance, dt-convergent, no drift', () => {
    // Flat plate, vx=0.01, muK=0.2: analytic stop distance
    // v^2/(2 muK g) = 2.548e-5 m, stop time v/(muK g) = 5.097e-3 s.
    const want = (0.01 * 0.01) / (2 * 0.2 * G)
    for (const dt of [1e-2, 1e-3, 1e-4]) {
      const st: SeatedState = { x: 0, y: 0, vx: 0.01, vy: 0 }
      let mode = ''
      for (let t = 0; t < 0.05; t += dt) {
        mode = stepSeated(st, flat, { muS: 0.3, muK: 0.2 }, t, dt).mode
      }
      expect(mode).toBe('stick')
      expect(st.vx).toBe(0)
      expect(Math.abs(st.x - want)).toBeLessThan(1e-9)
      expect(st.y).toBe(0)
    }
    // Negative initial velocity is symmetric.
    const st: SeatedState = { x: 0, y: 0, vx: -0.01, vy: 0.004 }
    for (let t = 0; t < 0.05; t += 1e-3) stepSeated(st, flat, { muS: 0.3, muK: 0.2 }, t, 1e-3)
    const v0 = Math.hypot(0.01, 0.004)
    const wantD = (v0 * v0) / (2 * 0.2 * G)
    expect(Math.abs(Math.hypot(st.x, st.y) - wantD)).toBeLessThan(1e-9)
    expect(st.vx).toBe(0); expect(st.vy).toBe(0)
  })

  it('M6-13R stop on incline: sticks if tan<muS, re-slides if not', () => {
    const th = 20 * Math.PI / 180
    const slide = (muS: number): SeatedState => {
      const st: SeatedState = { x: 0, y: 0, vx: 0.05, vy: 0 } // up-slope (+x)
      for (let t = 0; t < 0.6; t += 1e-3) stepSeated(st, incline(th), { muS, muK: 0.2 }, t, 1e-3)
      return st
    }
    const stuck = slide(0.5)
    expect(stuck.vx).toBe(0)
    expect(stuck.x).toBeGreaterThan(0) // stopped partway up, stays
    const xStop = stuck.x
    for (let t = 0; t < 0.2; t += 1e-3) stepSeated(stuck, incline(th), { muS: 0.5, muK: 0.2 }, t, 1e-3)
    expect(stuck.x).toBe(xStop)
    const reslid = slide(0.3) // tan20=0.364 > 0.3: stops, then slides back down
    expect(reslid.vx).toBeLessThan(0)
    expect(reslid.x).toBeLessThan(0)
  })

  it('M6-13R mid-step detach is located inside a coarse step', () => {
    const om = 2 * Math.PI * 50
    const aBig = (1.5 * G) / (om * om)
    const surf = inclineVibrating(0, aBig, om)
    const st: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    const want = Math.asin(2 / 3) / om
    // One coarse 5 ms step straddling the detach instant (~2.33 ms).
    const r = stepSeated(st, surf, { muS: 0.5, muK: 0.4 }, 0, 5e-3)
    expect(r.mode).toBe('detach')
    expect(Math.abs(r.advanced - want)).toBeLessThan(1e-6)
  })

  it('seatedVz matches the constraint velocity', () => {
    const th = 0.3
    const o = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
    incline(th)(0.01, 0, 0, o)
    const st: SeatedState = { x: 0.01, y: 0, vx: 0.2, vy: -0.1 }
    expect(seatedVz(st, o)).toBeCloseTo(Math.tan(th) * 0.2, 15)
  })
})

describe('M6-14R slip onset and event ordering', () => {
  it('slip from rest feels kinetic friction from the first step', () => {
    // h = 0.5 x, rest, muS=0.1 (breaks), muK=0.2, one dt=0.1 step.
    // Analytic: a_slope = g(sin - muK cos), vx = -a_slope cos t.
    const th = Math.atan(0.5)
    const st: SeatedState = { x: 0, y: 0, vx: 0, vy: 0 }
    const r = stepSeated(st, incline(th), { muS: 0.1, muK: 0.2 }, 0, 0.1)
    expect(r.mode).toBe('slip')
    const vxExact = -G * (Math.sin(th) - 0.2 * Math.cos(th)) * Math.cos(th) * 0.1
    expect(st.vx).toBeCloseTo(vxExact, 6) // ~ -0.23544
    expect(st.vx).toBeCloseTo(-0.23544, 4)
  })

  it('stop precedes detach inside one step (event ordering)', () => {
    // Flat surface falling away: h(t) with htt = -g t/0.02, so the
    // normal reaches zero at t where g+htt=0 -> t=0.02 (starting at
    // htt=0). Grain sliding at 0.01 m/s stops at ~v/(muK g)=5.1ms,
    // BEFORE the detach. Correct order: stop (stick) then detach;
    // the step must end stuck at the analytic stop point when the
    // detach is reported, not mid-slide past it.
    const surf: SurfaceFn = (_x, _y, t, o) => {
      o.h = -G * t * t * t / (6 * 0.02); o.hx = 0; o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0
      o.ht = -G * t * t / (2 * 0.02); o.htt = -G * t / 0.02
      o.htx = 0; o.hty = 0
    }
    const st: SeatedState = { x: 0, y: 0, vx: 0.01, vy: 0 }
    const r = stepSeated(st, surf, { muS: 0.3, muK: 0.2 }, 0, 0.03)
    expect(r.mode).toBe('detach')
    // Detach at t=0.02; grain stopped before that (vx ~ 0, position
    // near the flat-plate stop distance, slightly less friction as
    // normal decays, so within [stop_flat, 2*stop_flat]).
    expect(r.advanced).toBeCloseTo(0.02, 3)
    expect(Math.abs(st.vx)).toBeLessThan(1e-6)
    const stopFlat = 0.01 * 0.01 / (2 * 0.2 * G)
    expect(st.x).toBeGreaterThan(stopFlat * 0.9)
    expect(st.x).toBeLessThan(stopFlat * 2.2)
  })
})
