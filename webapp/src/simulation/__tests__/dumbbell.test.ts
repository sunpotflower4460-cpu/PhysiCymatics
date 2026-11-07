// M6-27 stage A: rigid dumbbell flight + impact, validated
// against analytic rigid-body results (docs/m6_27_design.md):
// sphere limit, contact-point angular-momentum conservation,
// energy non-gain, torque-free precession period, orientation
// scattering. Expectations are derived independently of the
// module's K-matrix code.
import { describe, expect, it } from 'vitest'
import {
  dumbbellFlightStep, dumbbellImpact, angularMomentum,
  angularMomentumAbout, kineticEnergyPerMass, inertiaPerMass,
  lobeCenters, contactVelocity,
  type DumbbellParams, type DumbbellState,
} from '../dumbbell'

const R = 2.5e-4
const mk = (d: number): DumbbellParams => ({
  lobeRadius: R, halfSeparation: d, beta: 2 / 5,
  muS: 0.3, muK: 0.2, restitution: 0.3,
})
const st0 = (over: Partial<DumbbellState>): DumbbellState => ({
  cx: 0, cy: 0, cz: R, vx: 0, vy: 0, vz: 0,
  ex: 0, ey: 0, ez: 1, wx: 0, wy: 0, wz: 0, ...over,
})

let seed = 12345
const rnd = (): number => {
  seed = (seed * 1103515245 + 12345) & 0x7fffffff
  return seed / 0x7fffffff * 2 - 1
}

describe('M6-27 dumbbell: sphere limit and analytic impacts', () => {
  it('d=0 reproduces the sphere stick impact formulas', () => {
    const p = mk(0)
    const beta = p.beta
    const st = st0({ vx: 0.2, vz: -0.5 })
    const r = dumbbellImpact(st, [0, 0, 0], [0, 0, 1], [0, 0, 0], p)
    expect(r.mode).toBe('stick')
    expect(st.vz).toBeCloseTo(0.15, 12) // -e * vz
    expect(st.vx).toBeCloseTo(0.2 / (1 + beta), 12)
    expect(st.wy).toBeCloseTo(0.2 / ((1 + beta) * R), 6)
    const u = contactVelocity(st, [0, 0, 0], [0, 0, 0])
    expect(Math.hypot(u[0], u[1])).toBeLessThan(1e-12) // rolling
  })

  it('vertical dumbbell: stick result from Kxx = 1 + (d+R)^2/Iperp', () => {
    const d = R
    const p = mk(d)
    const { iPerp } = inertiaPerMass(p)
    const st = st0({ cz: d + R, vx: 0.2, vz: -0.5 })
    const r = dumbbellImpact(st, [0, 0, 0], [0, 0, 1], [0, 0, 0], p)
    expect(r.mode).toBe('stick')
    const Kxx = 1 + ((d + R) ** 2) / iPerp
    expect(st.vx).toBeCloseTo(0.2 * (1 - 1 / Kxx), 10)
    expect(st.vz).toBeCloseTo(0.15, 12)
    const u = contactVelocity(st, [0, 0, 0], [0, 0, 0])
    expect(Math.hypot(u[0], u[1])).toBeLessThan(1e-10)
  })

  it('separating contact applies no impulse', () => {
    const p = mk(R)
    const st = st0({ cz: 2 * R, vx: 0.1, vz: 0.4 })
    const r = dumbbellImpact(st, [0, 0, 0], [0, 0, 1], [0, 0, 0], p)
    expect(r.mode).toBe('separating')
    expect(st.vx).toBe(0.1)
    expect(st.vz).toBe(0.4)
  })
})

describe('M6-27 dumbbell: conservation laws (deterministic random cases)', () => {
  it('contact-point angular momentum conserved; energy never gains (e<=1, Vs=0)', () => {
    let checked = 0
    for (let k = 0; k < 200 && checked < 60; k++) {
      const d = R * (0.5 + Math.abs(rnd()))
      const p = { ...mk(d), restitution: Math.abs(rnd()) }
      let e: [number, number, number] = [rnd(), rnd(), rnd()]
      const em = Math.hypot(...e)
      if (em < 0.2) continue
      e = [e[0] / em, e[1] / em, e[2] / em]
      const st = st0({
        cx: 0.01 * rnd(), cy: 0.01 * rnd(), cz: 0.01,
        vx: rnd(), vy: rnd(), vz: rnd(),
        ex: e[0], ey: e[1], ez: e[2],
        wx: 300 * rnd(), wy: 300 * rnd(), wz: 300 * rnd(),
      })
      const lobes = lobeCenters(st, p)
      const lower = lobes[0][2] <= lobes[1][2] ? lobes[0] : lobes[1]
      const nRaw: [number, number, number] = [0.3 * rnd(), 0.3 * rnd(), 1]
      const nm = Math.hypot(...nRaw)
      const n: [number, number, number] = [nRaw[0] / nm, nRaw[1] / nm, nRaw[2] / nm]
      const point: [number, number, number] = [
        lower[0] - R * n[0], lower[1] - R * n[1], lower[2] - R * n[2],
      ]
      const H0: [number, number, number] = [0, 0, 0]
      angularMomentumAbout(st, p, point, H0)
      const E0 = kineticEnergyPerMass(st, p)
      const r = dumbbellImpact(st, point, n, [0, 0, 0], p)
      if (r.mode === 'separating') continue
      checked++
      const H1: [number, number, number] = [0, 0, 0]
      angularMomentumAbout(st, p, point, H1)
      for (let i = 0; i < 3; i++) {
        expect(Math.abs(H1[i] - H0[i])).toBeLessThan(1e-12 + 1e-9 * Math.abs(H0[i]))
      }
      expect(kineticEnergyPerMass(st, p)).toBeLessThanOrEqual(E0 + 1e-12)
    }
    expect(checked).toBeGreaterThan(30)
  })

  it('torque-free flight: L conserved, axis precesses with period 2pi*Iperp/|L|', () => {
    const p = mk(R)
    const { iPerp } = inertiaPerMass(p)
    const th = 0.6
    const st = st0({
      cz: 0.1, ex: Math.sin(th), ey: 0, ez: Math.cos(th), wy: 100,
    })
    const L0: [number, number, number] = [0, 0, 0]
    angularMomentum(st, p, L0)
    const lmag = Math.hypot(...L0)
    const period = 2 * Math.PI * iPerp / lmag
    const rotE = (): number => {
      const E = kineticEnergyPerMass(st, p)
      return E - 0.5 * (st.vx ** 2 + st.vy ** 2 + st.vz ** 2)
    }
    const Er0 = rotE()
    const eStart = [st.ex, st.ey, st.ez]
    const N = 628
    const dt = period / N
    for (let k = 0; k < N; k++) dumbbellFlightStep(st, p, dt)
    expect(Math.hypot(st.ex - eStart[0], st.ey - eStart[1], st.ez - eStart[2]))
      .toBeLessThan(2e-3)
    const L1: [number, number, number] = [0, 0, 0]
    angularMomentum(st, p, L1)
    for (let i = 0; i < 3; i++) {
      expect(Math.abs(L1[i] - L0[i])).toBeLessThan(1e-6 * lmag)
    }
    expect(Math.abs(rotE() - Er0)).toBeLessThan(1e-6 * Er0)
    expect(Math.abs(Math.hypot(st.ex, st.ey, st.ez) - 1)).toBeLessThan(1e-12)
  })

  it('spin about the axis keeps the axis fixed', () => {
    const p = mk(R)
    const st = st0({ cz: 0.1, ex: 0.6, ey: 0.48, ez: 0.64, wx: 60, wy: 48, wz: 64 })
    // w parallel to e (100 * e): e must not move.
    for (let k = 0; k < 100; k++) dumbbellFlightStep(st, p, 1e-4)
    expect(Math.hypot(st.ex - 0.6, st.ey - 0.48, st.ez - 0.64)).toBeLessThan(1e-9)
  })
})

describe('M6-27 dumbbell: orientation scattering at a flat impact', () => {
  const scatter = (d: number): { angles: number[]; speeds: number[] } => {
    const p = mk(d)
    const angles: number[] = []
    const speeds: number[] = []
    for (const polar of [0.3, 0.6, 0.9, 1.2]) {
      for (const azim of [0, 2.1, 4.2]) {
        const e: [number, number, number] = [
          Math.sin(polar) * Math.cos(azim),
          Math.sin(polar) * Math.sin(azim),
          Math.cos(polar),
        ]
        const cz = R + d * Math.abs(e[2])
        const st = st0({
          cz, vx: 0.25, vz: -0.6, ex: e[0], ey: e[1], ez: e[2],
        })
        const lobes = lobeCenters(st, p)
        const lower = lobes[0][2] <= lobes[1][2] ? lobes[0] : lobes[1]
        const point: [number, number, number] = [lower[0], lower[1], lower[2] - R]
        dumbbellImpact(st, point, [0, 0, 1], [0, 0, 0], p)
        angles.push(Math.atan2(st.vy, st.vx))
        speeds.push(Math.hypot(st.vx, st.vy))
      }
    }
    return { angles, speeds }
  }

  it('sphere control (d=0) shows zero angular scatter', () => {
    const { angles } = scatter(0)
    for (const a of angles) expect(Math.abs(a)).toBeLessThan(1e-12)
  })

  it('aspect-2 dumbbell scatters the outgoing direction by orientation', () => {
    const { angles, speeds } = scatter(R)
    const mean = angles.reduce((s, v) => s + v, 0) / angles.length
    const std = Math.sqrt(angles.reduce((s, v) => s + (v - mean) ** 2, 0) / angles.length)
    expect(std).toBeGreaterThan(0.05)
    expect(Math.max(...speeds)).toBeGreaterThan(0)
    // Symmetry: axis in the x-z plane (azimuth 0) gives vy = 0.
    const p = mk(R)
    const e: [number, number, number] = [Math.sin(0.9), 0, Math.cos(0.9)]
    const st = st0({ cz: R + R * Math.abs(e[2]), vx: 0.25, vz: -0.6, ex: e[0], ey: 0, ez: e[2] })
    const lobes = lobeCenters(st, p)
    const lower = lobes[0][2] <= lobes[1][2] ? lobes[0] : lobes[1]
    dumbbellImpact(st, [lower[0], lower[1], lower[2] - R], [0, 0, 1], [0, 0, 0], p)
    expect(Math.abs(st.vy)).toBeLessThan(1e-12)
    expect(Math.abs(st.wx)).toBeLessThan(1e-9)
    expect(Math.abs(st.wz)).toBeLessThan(1e-9)
  })
})
