// M6-22: sphere-sphere frictional impact — analytic validation.
import { describe, expect, it } from 'vitest'
import { pairImpact } from '../pairImpact'
import { flightEnergyPerMass, type FlightState, type ImpactParams } from '../rollingImpact'

const P: ImpactParams = { muS: 0.3, muK: 0.2, radius: 1e-3, beta: 2 / 5, restitution: 1 }
const mk = (cx: number, vx: number, vy = 0, vz = 0, wx = 0, wy = 0, wz = 0): FlightState => ({
  cx, cy: 0, cz: 0, vx, vy, vz, wx, wy, wz,
})

function totalMomentum(a: FlightState, b: FlightState): [number, number, number] {
  return [a.vx + b.vx, a.vy + b.vy, a.vz + b.vz]
}
function totalAngMom(a: FlightState, b: FlightState, p: ImpactParams): [number, number, number] {
  const one = (s: FlightState): [number, number, number] => {
    const L: [number, number, number] = [
      s.cy * s.vz - s.cz * s.vy + p.beta * p.radius * p.radius * s.wx,
      s.cz * s.vx - s.cx * s.vz + p.beta * p.radius * p.radius * s.wy,
      s.cx * s.vy - s.cy * s.vx + p.beta * p.radius * p.radius * s.wz,
    ]
    return L
  }
  const la = one(a); const lb = one(b)
  return [la[0] + lb[0], la[1] + lb[1], la[2] + lb[2]]
}

describe('M6-22 pair impact', () => {
  it('head-on elastic: equal spheres exchange velocities', () => {
    const a = mk(-P.radius, 1)
    const b = mk(P.radius, -0.5)
    expect(pairImpact(a, b, P)).toBe('stick')
    expect(a.vx).toBeCloseTo(-0.5, 12)
    expect(b.vx).toBeCloseTo(1, 12)
  })

  it('oblique stick: tangential contact velocity vanishes', () => {
    const a = mk(-P.radius, 1, 0.05)
    const b = mk(P.radius, 0, 0)
    // Stick needs beta/(2(1+beta))*0.05 = 0.00714 <= muS*Jn(=1*1)=0.3
    expect(pairImpact(a, b, { ...P, restitution: 0.5 })).toBe('stick')
    // Contact relative tangential velocity after impact:
    const uTy = (a.vy - b.vy) + P.radius * (a.wz + b.wz) * 0
      - P.radius * ((a.wx + b.wx) * 0 - 0) // y-comp of R(w_sum x n), n=x: (w x n)_y = wz*1 - 0 = wz_sum... computed below directly
    void uTy
    const wsx = 0; void wsx
    // u = (Va-Vb) + R (wa+wb) x n, n=(1,0,0): (w x n) = (0, wz, -wy)
    const uy = (a.vy - b.vy) + P.radius * (a.wz + b.wz)
    const uz = (a.vz - b.vz) - P.radius * (a.wy + b.wy)
    expect(Math.abs(uy)).toBeLessThan(1e-12)
    expect(Math.abs(uz)).toBeLessThan(1e-12)
  })

  it('momentum/angular momentum conserved and energy never gains (random)', () => {
    let seed = 999
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let k = 0; k < 300; k++) {
      const th = rnd() * Math.PI; const ph = rnd() * 2 * Math.PI
      const nx = Math.sin(th) * Math.cos(ph); const ny = Math.sin(th) * Math.sin(ph); const nz = Math.cos(th)
      const a: FlightState = {
        cx: -P.radius * nx, cy: -P.radius * ny, cz: -P.radius * nz,
        vx: rnd(), vy: rnd(), vz: rnd(),
        wx: 400 * rnd(), wy: 400 * rnd(), wz: 400 * rnd(),
      }
      const b: FlightState = {
        cx: P.radius * nx, cy: P.radius * ny, cz: P.radius * nz,
        vx: rnd(), vy: rnd(), vz: rnd(),
        wx: 400 * rnd(), wy: 400 * rnd(), wz: 400 * rnd(),
      }
      // Ensure closing along n.
      const closing = (a.vx - b.vx) * nx + (a.vy - b.vy) * ny + (a.vz - b.vz) * nz
      if (closing <= 0.05) continue
      const p: ImpactParams = { ...P, restitution: Math.abs(rnd()) }
      const m0 = totalMomentum(a, b)
      const L0 = totalAngMom(a, b, p)
      const e0 = flightEnergyPerMass(a, p) + flightEnergyPerMass(b, p)
      const mode = pairImpact(a, b, p)
      expect(mode).not.toBe('separating')
      const m1 = totalMomentum(a, b)
      const L1 = totalAngMom(a, b, p)
      for (let i = 0; i < 3; i++) {
        expect(Math.abs(m1[i] - m0[i])).toBeLessThan(1e-12)
        expect(Math.abs(L1[i] - L0[i])).toBeLessThan(1e-12 * (1 + Math.abs(L0[i])) + 1e-15)
      }
      expect(flightEnergyPerMass(a, p) + flightEnergyPerMass(b, p))
        .toBeLessThanOrEqual(e0 * (1 + 1e-12))
    }
  })

  it('separating pairs are untouched', () => {
    const a = mk(-P.radius, -1)
    const b = mk(P.radius, 1)
    expect(pairImpact(a, b, P)).toBe('separating')
    expect(a.vx).toBe(-1)
  })
})
