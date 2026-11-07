// M6-18: frictional sphere impact — analytic stick/slide impulses,
// angular-momentum conservation about the contact, energy
// dissipation, and a full drop -> impact -> rolling cycle.
import { describe, expect, it } from 'vitest'
import {
  sphereImpact, flightToRolling, rollingToFlight,
  contactAngularMomentum, flightEnergyPerMass,
  type FlightState, type ImpactParams,
} from '../rollingImpact'
import { stepDynRolling, type SurfaceFnDyn } from '../rollingDynamic'
import { G } from '../surface'

const P: ImpactParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5, restitution: 0.5 }

const zeroDyn = (o: Parameters<SurfaceFnDyn>[3]): void => {
  o.h = 0; o.hx = 0; o.hy = 0; o.hxx = 0; o.hxy = 0; o.hyy = 0
  o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
  o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
  o.htxx = 0; o.htxy = 0; o.htyy = 0; o.httx = 0; o.htty = 0
}
const flatDyn: SurfaceFnDyn = (_x, _y, _t, o) => zeroDyn(o)
const inclineDyn = (th: number): SurfaceFnDyn => (x, _y, _t, o) => {
  zeroDyn(o)
  o.h = x * Math.tan(th); o.hx = Math.tan(th)
}

describe('M6-18 frictional impact', () => {
  it('stick impact on flat: analytic post-impact state', () => {
    // Drop with vx=0.7, vz=-1, no spin. Jn = 1.5, stick needs
    // beta/(1+beta)*0.7 = 0.2 <= muS*Jn = 0.45 -> stick.
    const st: FlightState = {
      cx: 0, cy: 0, cz: P.radius, vx: 0.7, vy: 0, vz: -1, wx: 0, wy: 0, wz: 0,
    }
    const r = sphereImpact(st, flatDyn, P, 0)
    expect(r.mode).toBe('stick')
    expect(st.vz).toBeCloseTo(0.5, 12) // e * 1
    expect(st.vx).toBeCloseTo(0.7 / (1 + P.beta), 12) // 0.5
    expect(st.wy).toBeCloseTo(st.vx / P.radius, 9) // contact slip zero
  })

  it('slide impact on flat: analytic impulses', () => {
    // Grazing impact: vz=-0.1 -> Jn=0.15; stick needs 0.2 > muS*Jn.
    const st: FlightState = {
      cx: 0, cy: 0, cz: P.radius, vx: 0.7, vy: 0, vz: -0.1, wx: 0, wy: 0, wz: 0,
    }
    const r = sphereImpact(st, flatDyn, P, 0)
    expect(r.mode).toBe('slide')
    expect(r.normalImpulse).toBeCloseTo(0.15, 12)
    expect(st.vx).toBeCloseTo(0.7 - P.muK * 0.15, 12)
    expect(st.vz).toBeCloseTo(0.05, 12)
    // wy = muK*Jn/(beta*R)
    expect(st.wy).toBeCloseTo((P.muK * 0.15) / (P.beta * P.radius), 9)
    // Residual slip: u' = u - muK*Jn*(1+1/beta) > 0 (no reversal).
    const uAfter = st.vx - P.radius * st.wy
    expect(uAfter).toBeCloseTo(0.7 - P.muK * 0.15 * (1 + 1 / P.beta), 12)
    expect(uAfter).toBeGreaterThan(0)
  })

  it('separating contact applies no impulse', () => {
    const st: FlightState = {
      cx: 0, cy: 0, cz: P.radius, vx: 0.3, vy: 0, vz: 0.2, wx: 0, wy: 0, wz: 0,
    }
    const r = sphereImpact(st, flatDyn, P, 0)
    expect(r.mode).toBe('separating')
    expect(st.vz).toBe(0.2)
  })

  it('angular momentum about contact conserved; energy never gains (random)', () => {
    let seed = 12345
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    const surf = inclineDyn(0.31)
    for (let k = 0; k < 200; k++) {
      const st: FlightState = {
        cx: 0.01 * rnd(), cy: 0.01 * rnd(), cz: 0,
        vx: rnd(), vy: rnd(), vz: -Math.abs(rnd()) - 0.05,
        wx: 300 * rnd(), wy: 300 * rnd(), wz: 300 * rnd(),
      }
      // Place center on the offset surface above (cx, cy):
      // contact params solve uses (cx,cy); set cz from the surface.
      const th = 0.31
      st.cz = st.cx * Math.tan(th) + P.radius / Math.cos(th)
      const L0: [number, number, number] = [0, 0, 0]
      const L1: [number, number, number] = [0, 0, 0]
      contactAngularMomentum(st, surf, P, 0, L0)
      const e0 = flightEnergyPerMass(st, P)
      const r = sphereImpact(st, surf, P, 0)
      if (r.mode === 'separating') continue
      contactAngularMomentum(st, surf, P, 0, L1)
      for (let i = 0; i < 3; i++) {
        expect(Math.abs(L1[i] - L0[i])).toBeLessThan(1e-9 * (1 + Math.abs(L0[i])) + 1e-12)
      }
      expect(flightEnergyPerMass(st, P)).toBeLessThanOrEqual(e0 * (1 + 1e-12))
    }
  })

  it('full cycle: drop with e=0 -> stick impact -> rolling at v0/(1+beta)', () => {
    const p0: ImpactParams = { ...P, restitution: 0 }
    // Ballistic drop from 1 cm above contact, vx = 0.7.
    const drop = 0.01
    const tFall = Math.sqrt(2 * drop / G)
    const st: FlightState = {
      cx: 0, cy: 0, cz: P.radius + drop, vx: 0.7, vy: 0, vz: 0, wx: 0, wy: 0, wz: 0,
    }
    st.cx += st.vx * tFall
    st.cz = P.radius
    st.vz = -G * tFall
    // Jn = sqrt(2 g drop) = 0.443; stick needs 0.2 > muS*Jn=0.133,
    // so the impact slides; the seated slip phase then finishes the
    // transition. The conserved C = vx + beta*R*wy still gives the
    // same final rolling speed 0.7/(1+beta).
    const r = sphereImpact(st, flatDyn, p0, tFall)
    expect(r.mode).toBe('slide')
    expect(st.vz).toBeCloseTo(0, 12)
    const roll = flightToRolling(st, flatDyn, p0, tFall)
    let t = tFall
    let mode = ''
    for (let k = 0; k < 2000 && mode !== 'roll'; k++) {
      mode = stepDynRolling(roll, flatDyn, p0, t, 1e-3).mode
      t += 1e-3
    }
    expect(mode).toBe('roll')
    for (let k = 0; k < 500; k++) {
      expect(stepDynRolling(roll, flatDyn, p0, t, 1e-3).mode).toBe('roll')
      t += 1e-3
    }
    expect(roll.vx).toBeCloseTo(0.7 / (1 + P.beta), 9)
    expect(roll.x).toBeGreaterThan(0.2)
  })

  it('conversion round-trip on an incline', () => {
    const surf = inclineDyn(0.2)
    const roll0 = { x: 0.015, y: -0.007, vx: 0.11, vy: -0.05, wx: 12, wy: -34, wz: 5 }
    const fl = rollingToFlight(roll0, surf, P, 0)
    const roll1 = flightToRolling(fl, surf, P, 0)
    expect(roll1.x).toBeCloseTo(roll0.x, 9)
    expect(roll1.y).toBeCloseTo(roll0.y, 9)
    expect(roll1.vx).toBeCloseTo(roll0.vx, 9)
    expect(roll1.vy).toBeCloseTo(roll0.vy, 9)
    expect(roll1.wx).toBe(roll0.wx)
  })

  it('second bounce after stick impact needs no friction impulse', () => {
    const st: FlightState = {
      cx: 0, cy: 0, cz: P.radius, vx: 0.7, vy: 0, vz: -1, wx: 0, wy: 0, wz: 0,
    }
    sphereImpact(st, flatDyn, P, 0) // stick, leaves with vz=+0.5
    // Flight up and down: returns with vz=-0.5, same vx, wy.
    st.vz = -0.5
    const vxBefore = st.vx; const wyBefore = st.wy
    const r2 = sphereImpact(st, flatDyn, P, 1)
    expect(r2.mode).toBe('stick')
    expect(Math.hypot(...r2.frictionImpulse)).toBeLessThan(1e-12)
    expect(st.vx).toBeCloseTo(vxBefore, 12)
    expect(st.wy).toBeCloseTo(wyBefore, 9)
    expect(st.vz).toBeCloseTo(0.25, 12)
  })
})
