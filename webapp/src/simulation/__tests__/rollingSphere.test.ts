// M6-15: rolling sphere on plane surfaces — analytic validation.
// Expectations are classical rigid-body results (solid sphere
// beta=2/5), not fitted to the implementation.
import { describe, expect, it } from 'vitest'
import {
  rollingAccel, stepRolling, kineticEnergyPerMass, centerHeight,
  slipVelocity, type RollingState, type RollingParams,
} from '../rollingSphere'
import type { SurfaceFn } from '../seatedFriction'
import { G } from '../surface'

const P: RollingParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5 }

const flat: SurfaceFn = (_x, _y, _t, o) => {
  o.h = 0; o.hx = 0; o.hy = 0; o.hxx = 0; o.hxy = 0; o.hyy = 0
  o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
}
const incline = (th: number): SurfaceFn => (x, _y, _t, o) => {
  o.h = x * Math.tan(th); o.hx = Math.tan(th); o.hy = 0
  o.hxx = 0; o.hxy = 0; o.hyy = 0; o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
}
const flatVibe = (amp: number, om: number): SurfaceFn => (_x, _y, t, o) => {
  o.h = amp * Math.sin(om * t); o.hx = 0; o.hy = 0
  o.hxx = 0; o.hxy = 0; o.hyy = 0
  o.ht = amp * om * Math.cos(om * t); o.htt = -amp * om * om * Math.sin(om * t)
  o.htx = 0; o.hty = 0
}

const rest = (): RollingState => ({ x: 0, y: 0, vx: 0, vy: 0, wx: 0, wy: 0, wz: 0 })

describe('M6-15 rolling sphere (plane)', () => {
  it('incline: rolling acceleration is g*sin(theta)/(1+beta)', () => {
    const th = 20 * Math.PI / 180
    const a = rollingAccel(rest(), incline(th), P, 0)
    expect(a.mode).toBe('roll')
    // Tangential acceleration magnitude from (ax, ay): |v_t'| with
    // v_t' = ax*e1 + ay*e2.
    const vtAcc = Math.hypot(a.ax, a.ay, Math.tan(th) * a.ax)
    expect(vtAcc).toBeCloseTo(G * Math.sin(th) / (1 + P.beta), 12)
    // Rolling downhill (negative x for upward slope +x).
    expect(a.ax).toBeLessThan(0)
  })

  it('incline: stick boundary at tan(theta) = muS (1+beta)/beta', () => {
    const tanCrit = P.muS * (1 + P.beta) / P.beta // 1.05
    const below = Math.atan(tanCrit * 0.99)
    const above = Math.atan(tanCrit * 1.01)
    expect(rollingAccel(rest(), incline(below), P, 0).mode).toBe('roll')
    expect(rollingAccel(rest(), incline(above), P, 0).mode).toBe('slip')
  })

  it('incline: energy conserved while rolling (KE + g*z_C)', () => {
    const th = 15 * Math.PI / 180
    const st = rest()
    const surf = incline(th)
    const e0 = kineticEnergyPerMass(st, surf, P, 0) + G * centerHeight(st, surf, P, 0)
    let t = 0
    for (let k = 0; k < 2000; k++) {
      const r = stepRolling(st, surf, P, t, 1e-3)
      expect(r.mode).toBe('roll')
      t += 1e-3
    }
    const e1 = kineticEnergyPerMass(st, surf, P, t) + G * centerHeight(st, surf, P, t)
    expect(Math.abs(e1 - e0)).toBeLessThan(1e-6 * Math.abs(e0) + 1e-9)
    // Speed check: |v_t| = a_t * t
    const vt = Math.hypot(st.vx, st.vy, Math.tan(th) * st.vx)
    expect(vt).toBeCloseTo(G * Math.sin(th) / (1 + P.beta) * t, 4)
    // No slip remains.
    const u: [number, number, number] = [0, 0, 0]
    slipVelocity(st, surf, P, t, u)
    expect(Math.hypot(...u)).toBeLessThan(1e-9)
  })

  it('flat: pure sliding transitions to rolling at v0/(1+beta)', () => {
    const st: RollingState = { ...rest(), vx: 0.7 }
    let mode = ''
    let t = 0
    for (let k = 0; k < 20000 && mode !== 'roll'; k++) {
      mode = stepRolling(st, flat, P, t, 1e-3).mode
      t += 1e-3
    }
    expect(mode).toBe('roll')
    expect(st.vx).toBeCloseTo(0.7 / (1 + P.beta), 6) // 0.5
    // Rolling without slipping: u = vx - R*wy = 0 (n=+z).
    expect(st.wy).toBeCloseTo(st.vx / P.radius, 3)
    // Transition time: u decays at muK*g*(1+1/beta).
    const tExpect = 0.7 / (P.muK * G * (1 + 1 / P.beta))
    expect(Math.abs(t - tExpect)).toBeLessThan(2e-3)
  })

  it('flat: excessive backspin ends rolling backward', () => {
    // During slip on flat, C = vx + beta*R*wy is conserved
    // (dvx = f, d(R*wy) = -f/beta). Rolling needs wy = vx/R, so
    // vx_final = C/(1+beta). Backspin wy<0 with |beta*R*wy| > vx
    // makes C < 0: the sphere ends rolling in -x.
    const st: RollingState = { ...rest(), vx: 0.2, wy: -0.6 / P.radius }
    const C = 0.2 + P.beta * (P.radius * st.wy)
    let mode = ''
    let t = 0
    for (let k = 0; k < 20000 && mode !== 'roll'; k++) {
      mode = stepRolling(st, flat, P, t, 1e-3).mode
      t += 1e-3
    }
    expect(mode).toBe('roll')
    expect(st.vx).toBeCloseTo(C / (1 + P.beta), 6)
    expect(st.vx).toBeLessThan(0)
    expect(st.wy).toBeCloseTo(st.vx / P.radius, 3)
  })

  it('flat vibrating: detach when htt < -g', () => {
    const om = 2 * Math.PI * 50
    const amp = (2 * G) / (om * om) // htt min = -2g
    const st: RollingState = { ...rest(), vx: 0.1 }
    const r = stepRolling(st, flatVibe(amp, om), P, 0, 0.01)
    expect(r.mode).toBe('detach')
    // htt = -amp*om^2*sin(om t) = -g at sin = 1/2 -> t = (pi/6)/om.
    expect(r.advanced).toBeCloseTo((Math.PI / 6) / om, 4)
  })

  it('curved surface is rejected, not approximated', () => {
    const curved: SurfaceFn = (x, _y, _t, o) => {
      o.h = x * x; o.hx = 2 * x; o.hy = 0; o.hxx = 2; o.hxy = 0; o.hyy = 0
      o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
    }
    expect(() => rollingAccel(rest(), curved, P, 0)).toThrow()
  })

  it('step-size convergence on incline rolling', () => {
    const th = 10 * Math.PI / 180
    const surf = incline(th)
    const run = (dt: number): number => {
      const st = rest(); let t = 0
      while (t < 0.5 - 1e-12) { stepRolling(st, surf, P, t, dt); t += dt }
      return st.x
    }
    const x1 = run(1e-2); const x2 = run(1e-3); const x3 = run(1e-4)
    const exact = -0.5 * (G * Math.sin(th) / (1 + P.beta)) * Math.cos(th) * 0.25
    expect(Math.abs(x2 - exact)).toBeLessThan(Math.abs(x1 - exact) + 1e-12)
    expect(Math.abs(x3 - exact)).toBeLessThan(1e-6)
  })
})
