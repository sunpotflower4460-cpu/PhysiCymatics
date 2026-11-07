// M6-16: rolling sphere on static curved surfaces — validation
// against finite differences of the center constraint, the
// cylindrical-trough oscillation period, energy conservation,
// convex-crest detachment, and the plane limit.
import { describe, expect, it } from 'vitest'
import {
  centerPartials, curvedRollingAccel, stepCurvedRolling,
  curvedKineticEnergy, curvedCenterHeight, curvedProjectToRolling,
  type SurfaceFn3, type SurfaceField3,
} from '../rollingSphereCurved'
import type { RollingState, RollingParams } from '../rollingSphere'
import { G } from '../surface'

const P: RollingParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5 }
const rest = (): RollingState => ({ x: 0, y: 0, vx: 0, vy: 0, wx: 0, wy: 0, wz: 0 })

// Sinusoidal surface h = A (cos kx + cos ky), all derivatives closed form.
const A0 = 0.002; const K = 2 * Math.PI / 0.24
const sinusoid: SurfaceFn3 = (x, y, o) => {
  o.h = A0 * (Math.cos(K * x) + Math.cos(K * y))
  o.hx = -A0 * K * Math.sin(K * x); o.hy = -A0 * K * Math.sin(K * y)
  o.hxx = -A0 * K * K * Math.cos(K * x); o.hyy = -A0 * K * K * Math.cos(K * y); o.hxy = 0
  o.hxxx = A0 * K ** 3 * Math.sin(K * x); o.hyyy = A0 * K ** 3 * Math.sin(K * y)
  o.hxxy = 0; o.hxyy = 0
}

// Cylindrical trough radius Rc (axis along x): h(y) = Rc - sqrt(Rc^2-y^2).
const trough = (Rc: number): SurfaceFn3 => (_x, y, o) => {
  const q = Math.sqrt(Rc * Rc - y * y)
  o.h = Rc - q; o.hx = 0; o.hy = y / q
  o.hxx = 0; o.hxy = 0; o.hyy = (Rc * Rc) / (q ** 3)
  o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = (3 * Rc * Rc * y) / (q ** 5)
}

// Paraboloid bowl h = (x^2+y^2)/(2 Rc) (third derivatives vanish).
const paraboloid = (Rc: number): SurfaceFn3 => (x, y, o) => {
  o.h = (x * x + y * y) / (2 * Rc); o.hx = x / Rc; o.hy = y / Rc
  o.hxx = 1 / Rc; o.hyy = 1 / Rc; o.hxy = 0
  o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
}

// Convex parabolic crest h = -y^2/(2 Rc).
const crest = (Rc: number): SurfaceFn3 => (_x, y, o) => {
  o.h = -(y * y) / (2 * Rc); o.hx = 0; o.hy = -y / Rc
  o.hxx = 0; o.hxy = 0; o.hyy = -1 / Rc
  o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
}

describe('M6-16 curved rolling', () => {
  it('center partials match finite differences of C(q)', () => {
    const fld = (): SurfaceField3 => ({
      h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0,
    })
    const cOf = (x: number, y: number): [number, number, number] => {
      const o = fld(); sinusoid(x, y, o)
      const sn = Math.hypot(o.hx, o.hy, 1)
      return [x - P.radius * o.hx / sn, y - P.radius * o.hy / sn, o.h + P.radius / sn]
    }
    for (const [x, y] of [[0.03, -0.02], [-0.05, 0.04], [0.001, 0.002]] as const) {
      const cp = centerPartials(x, y, sinusoid, P.radius)
      const e1 = 1e-6
      for (let i = 0; i < 3; i++) {
        const fdX = (cOf(x + e1, y)[i] - cOf(x - e1, y)[i]) / (2 * e1)
        const fdY = (cOf(x, y + e1)[i] - cOf(x, y - e1)[i]) / (2 * e1)
        expect(Math.abs(fdX - cp.J[0][i])).toBeLessThan(1e-5)
        expect(Math.abs(fdY - cp.J[1][i])).toBeLessThan(1e-5)
      }
      const e2 = 2e-4
      for (let i = 0; i < 3; i++) {
        const fxx = (cOf(x + e2, y)[i] - 2 * cOf(x, y)[i] + cOf(x - e2, y)[i]) / (e2 * e2)
        const fyy = (cOf(x, y + e2)[i] - 2 * cOf(x, y)[i] + cOf(x, y - e2)[i]) / (e2 * e2)
        const fxy = (cOf(x + e2, y + e2)[i] - cOf(x + e2, y - e2)[i]
          - cOf(x - e2, y + e2)[i] + cOf(x - e2, y - e2)[i]) / (4 * e2 * e2)
        expect(Math.abs(fxx - cp.H[0][i])).toBeLessThan(2e-3)
        expect(Math.abs(fxy - cp.H[1][i])).toBeLessThan(2e-3)
        expect(Math.abs(fyy - cp.H[2][i])).toBeLessThan(2e-3)
      }
    }
  })

  it('cylindrical trough: small-oscillation period', () => {
    const Rc = 0.05
    const surf = trough(Rc)
    const st: RollingState = { ...rest(), y: 0.004 }
    const T = 2 * Math.PI * Math.sqrt((1 + P.beta) * (Rc - P.radius) / G)
    let t = 0; let prevY = st.y; const crossings: number[] = []
    const dt = 2e-4
    while (t < 3.2 * T) {
      const r = stepCurvedRolling(st, surf, P, dt)
      expect(r.mode).toBe('roll')
      t += dt
      if (prevY > 0 && st.y <= 0) crossings.push(t)
      prevY = st.y
    }
    expect(crossings.length).toBeGreaterThanOrEqual(3)
    const period = (crossings[2] - crossings[0]) / 2 // one downward crossing per period
    expect(Math.abs(period - T) / T).toBeLessThan(0.01)
  })

  it('paraboloid: energy conserved while rolling', () => {
    const surf = paraboloid(0.05)
    const st: RollingState = { ...rest(), x: 0.02 }
    const e0 = curvedKineticEnergy(st, surf, P) + G * curvedCenterHeight(st, surf, P)
    for (let k = 0; k < 3000; k++) {
      const r = stepCurvedRolling(st, surf, P, 5e-4)
      expect(r.mode).toBe('roll')
    }
    const e1 = curvedKineticEnergy(st, surf, P) + G * curvedCenterHeight(st, surf, P)
    expect(Math.abs(e1 - e0)).toBeLessThan(1e-6 * Math.abs(e0) + 1e-9)
  })

  it('sinusoid: energy conserved and slip stays zero in roll mode', () => {
    const st: RollingState = { ...rest(), x: 0.03, y: -0.01 }
    const e0 = curvedKineticEnergy(st, sinusoid, P) + G * curvedCenterHeight(st, sinusoid, P)
    for (let k = 0; k < 2000; k++) {
      const r = stepCurvedRolling(st, sinusoid, { ...P, muS: 5 }, 5e-4)
      expect(r.mode).toBe('roll')
    }
    const e1 = curvedKineticEnergy(st, sinusoid, P) + G * curvedCenterHeight(st, sinusoid, P)
    expect(Math.abs(e1 - e0)).toBeLessThan(1e-5 * Math.abs(e0) + 1e-9)
  })

  it('convex crest: fast grain detaches, slow grain stays', () => {
    const Rc = 0.05
    const surf = crest(Rc)
    // At crest, N per mass ~= g - V^2/Rc (center follows crest).
    const fast: RollingState = { ...rest(), vy: 1.0 }
    curvedProjectToRolling(fast, surf, P)
    expect(curvedRollingAccel(fast, surf, P).mode).toBe('detach')
    const slow: RollingState = { ...rest(), vy: 0.3 }
    curvedProjectToRolling(slow, surf, P)
    const a = curvedRollingAccel(slow, surf, P)
    expect(a.mode).toBe('roll')
    expect(a.nPerMass).toBeCloseTo(G - 0.09 / Rc, 1)
  })

  it('plane limit matches M6-15 incline acceleration', () => {
    const th = 20 * Math.PI / 180
    const plane: SurfaceFn3 = (x, _y, o) => {
      o.h = x * Math.tan(th); o.hx = Math.tan(th); o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0; o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
    }
    const a = curvedRollingAccel(rest(), plane, P)
    expect(a.mode).toBe('roll')
    const vtAcc = Math.hypot(a.ax, a.ay, Math.tan(th) * a.ax)
    expect(vtAcc).toBeCloseTo(G * Math.sin(th) / (1 + P.beta), 9)
  })
})
