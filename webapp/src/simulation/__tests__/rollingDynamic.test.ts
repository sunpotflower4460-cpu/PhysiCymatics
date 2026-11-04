// M6-17: rolling on time-dependent curved surfaces — FD partials,
// cross-check against the independent M6-15 plane module, static
// limit against the trough period, vibrating-trough detachment.
import { describe, expect, it } from 'vitest'
import {
  dynCenterPartials, dynRollingAccel, stepDynRolling,
  type SurfaceFnDyn, type SurfaceFieldDyn,
} from '../rollingDynamic'
import { stepRolling, type RollingState, type RollingParams } from '../rollingSphere'
import type { SurfaceFn } from '../seatedFriction'
import { G } from '../surface'

const P: RollingParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5 }
const rest = (): RollingState => ({ x: 0, y: 0, vx: 0, vy: 0, wx: 0, wy: 0, wz: 0 })

const A0 = 0.002; const K = 2 * Math.PI / 0.24; const OM = 60
const dynSinusoid: SurfaceFnDyn = (x, y, t, o) => {
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

describe('M6-17 dynamic curved rolling', () => {
  it('center partials (x,y,t) match finite differences', () => {
    const fld = (): SurfaceFieldDyn => ({
      h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, hxxx: 0, hxxy: 0, hxyy: 0, hyyy: 0,
      ht: 0, htt: 0, htx: 0, hty: 0, htxx: 0, htxy: 0, htyy: 0, httx: 0, htty: 0,
    })
    const cOf = (x: number, y: number, t: number): [number, number, number] => {
      const o = fld(); dynSinusoid(x, y, t, o)
      const sn = Math.hypot(o.hx, o.hy, 1)
      return [x - P.radius * o.hx / sn, y - P.radius * o.hy / sn, o.h + P.radius / sn]
    }
    const x = 0.0317; const y = -0.0213; const t = 0.37
    const cp = dynCenterPartials(x, y, t, dynSinusoid, P.radius)
    const es = 1e-6; const et = 1e-7
    for (let i = 0; i < 3; i++) {
      expect(Math.abs((cOf(x + es, y, t)[i] - cOf(x - es, y, t)[i]) / (2 * es) - cp.D[0][i])).toBeLessThan(1e-5)
      expect(Math.abs((cOf(x, y + es, t)[i] - cOf(x, y - es, t)[i]) / (2 * es) - cp.D[1][i])).toBeLessThan(1e-5)
      expect(Math.abs((cOf(x, y, t + et)[i] - cOf(x, y, t - et)[i]) / (2 * et) - cp.D[2][i])).toBeLessThan(1e-4)
    }
    const e2 = 2e-4; const e2t = 2e-5
    for (let i = 0; i < 3; i++) {
      const fxx = (cOf(x + e2, y, t)[i] - 2 * cOf(x, y, t)[i] + cOf(x - e2, y, t)[i]) / (e2 * e2)
      const fxt = (cOf(x + e2, y, t + e2t)[i] - cOf(x + e2, y, t - e2t)[i]
        - cOf(x - e2, y, t + e2t)[i] + cOf(x - e2, y, t - e2t)[i]) / (4 * e2 * e2t)
      const ftt = (cOf(x, y, t + e2t)[i] - 2 * cOf(x, y, t)[i] + cOf(x, y, t - e2t)[i]) / (e2t * e2t)
      expect(Math.abs(fxx - cp.H[0][i])).toBeLessThan(5e-3)
      expect(Math.abs(fxt - cp.H[3][i])).toBeLessThan(5e-2)
      expect(Math.abs(ftt - cp.H[5][i])).toBeLessThan(5e-2)
    }
  })

  it('vibrating plane incline matches M6-15 plane module trajectory', () => {
    const th = 15 * Math.PI / 180
    const om = 2 * Math.PI * 20
    const amp = 0.3 * G / (om * om)
    const planeFn: SurfaceFn = (x, _y, t, o) => {
      o.h = x * Math.tan(th) + amp * Math.sin(om * t); o.hx = Math.tan(th); o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0
      o.ht = amp * om * Math.cos(om * t); o.htt = -amp * om * om * Math.sin(om * t)
      o.htx = 0; o.hty = 0
    }
    const dynFn: SurfaceFnDyn = (x, _y, t, o) => {
      o.h = x * Math.tan(th) + amp * Math.sin(om * t); o.hx = Math.tan(th); o.hy = 0
      o.hxx = 0; o.hxy = 0; o.hyy = 0; o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = 0
      o.ht = amp * om * Math.cos(om * t); o.htt = -amp * om * om * Math.sin(om * t)
      o.htx = 0; o.hty = 0; o.htxx = 0; o.htxy = 0; o.htyy = 0; o.httx = 0; o.htty = 0
    }
    const sa = rest(); const sb = rest()
    let t = 0
    for (let k = 0; k < 1000; k++) {
      expect(stepRolling(sa, planeFn, P, t, 5e-4).mode).toBe('roll')
      expect(stepDynRolling(sb, dynFn, P, t, 5e-4).mode).toBe('roll')
      t += 5e-4
    }
    expect(Math.abs(sa.x - sb.x)).toBeLessThan(1e-9)
    expect(Math.abs(sa.vx - sb.vx)).toBeLessThan(1e-9)
    expect(Math.abs(sa.wy - sb.wy)).toBeLessThan(1e-6)
  })

  it('static limit: trough period matches M6-16 analytic value', () => {
    const Rc = 0.05
    const surf: SurfaceFnDyn = (_x, y, _t, o) => {
      const q = Math.sqrt(Rc * Rc - y * y)
      o.h = Rc - q; o.hx = 0; o.hy = y / q
      o.hxx = 0; o.hxy = 0; o.hyy = (Rc * Rc) / (q ** 3)
      o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = (3 * Rc * Rc * y) / (q ** 5)
      o.ht = 0; o.htt = 0; o.htx = 0; o.hty = 0
      o.htxx = 0; o.htxy = 0; o.htyy = 0; o.httx = 0; o.htty = 0
    }
    const st: RollingState = { ...rest(), y: 0.004 }
    const T = 2 * Math.PI * Math.sqrt((1 + P.beta) * (Rc - P.radius) / G)
    let t = 0; let prevY = st.y; const crossings: number[] = []
    while (t < 2.2 * T) {
      expect(stepDynRolling(st, surf, P, t, 2e-4).mode).toBe('roll')
      t += 2e-4
      if (prevY > 0 && st.y <= 0) crossings.push(t)
      prevY = st.y
    }
    expect(crossings.length).toBeGreaterThanOrEqual(2)
    expect(Math.abs((crossings[1] - crossings[0]) - T) / T).toBeLessThan(0.01)
  })

  it('vibrating trough: detach at analytic phase, contact at small amp', () => {
    const Rc = 0.05; const om = 2 * Math.PI * 30
    const mk = (amp: number): SurfaceFnDyn => (_x, y, t, o) => {
      const q = Math.sqrt(Rc * Rc - y * y)
      o.h = Rc - q + amp * Math.sin(om * t); o.hx = 0; o.hy = y / q
      o.hxx = 0; o.hxy = 0; o.hyy = (Rc * Rc) / (q ** 3)
      o.hxxx = 0; o.hxxy = 0; o.hxyy = 0; o.hyyy = (3 * Rc * Rc * y) / (q ** 5)
      o.ht = amp * om * Math.cos(om * t); o.htt = -amp * om * om * Math.sin(om * t)
      o.htx = 0; o.hty = 0; o.htxx = 0; o.htxy = 0; o.htyy = 0; o.httx = 0; o.htty = 0
    }
    // Small amp: htt min = -0.5 g -> stays seated at the bottom.
    const small = rest()
    let t = 0
    for (let k = 0; k < 400; k++) {
      expect(stepDynRolling(small, mk(0.5 * G / (om * om)), P, t, 1e-4).mode).toBe('roll')
      t += 1e-4
    }
    expect(Math.abs(small.y)).toBeLessThan(1e-9)
    // Large amp: htt min = -2g -> detach at sin(om t) = 1/2.
    const big = rest()
    t = 0
    let detachT = -1
    for (let k = 0; k < 400 && detachT < 0; k++) {
      const r = stepDynRolling(big, mk(2 * G / (om * om)), P, t, 1e-4)
      if (r.mode === 'detach') detachT = t + r.advanced
      t += 1e-4
    }
    expect(detachT).toBeGreaterThan(0)
    expect(Math.abs(detachT - (Math.PI / 6) / om)).toBeLessThan(3e-4)
    expect(dynRollingAccel(rest(), mk(2 * G / (om * om)), P, 0).mode).toBe('roll')
  })
})

describe('M6-17 integration smoke: rolling on the consistent plate field', () => {
  it('single grain steps on an excited plate without NaN', async () => {
    const { TransientPlate } = await import('../plateTransient')
    const { singleToneScore } = await import('../sequence')
    const { calibratedForceN } = await import('../transientSand')
    const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.05, 0.03])
    // Driven from rest (no artificial modal excitation): the grain
    // should stay in contact and roll/slip rather than detach at once.
    const surf: SurfaceFnDyn = (x, y, t, o) => plate.surfaceFieldDynamicConsistent(x, y, t, o)
    const dt = 1 / (44 * 152)
    const st: RollingState = { ...rest(), x: 0.02, y: -0.01 }
    const modes = new Set<string>()
    let t = 0
    for (let k = 0; k < 200; k++) {
      plate.beginStep(t, t + dt)
      const r = stepDynRolling(st, surf, P, t, dt)
      modes.add(r.mode)
      t += dt
      if (r.mode === 'detach') break
    }
    for (const v of [st.x, st.y, st.vx, st.vy, st.wx, st.wy, st.wz]) {
      expect(Number.isFinite(v)).toBe(true)
    }
    expect(modes.size).toBeGreaterThan(0)
    console.log('plate rolling smoke modes:', [...modes].join(','), 'final t=', t)
  })
})
