// M2 gate: TS sand core vs the shared M1 fixture (Python oracle output),
// plus single-body mechanics. Tolerance pre-declared: same arithmetic in
// float64; only libm sin/cos ulp differences are expected over 0.3 s.
import { describe, it, expect } from 'vitest'
import { SandSim, initialState } from '../sand'
import { G, SingleModeSurface } from '../surface'
import truth from '../../data/truth/sand_m1_short.json'

interface Truth {
  init_x: number[]; init_y: number[]; x: number[]; y: number[]
  sitting: boolean[]; n_bounce: number
}
const T = truth as unknown as Truth

describe('sand TS port vs M1 fixture (c002 oracle)', () => {
  it('matches final positions to < 1e-9 m over 0.3 s', () => {
    const surf = new SingleModeSurface(80)
    const st = initialState(surf.a, 0, 0, { x: T.init_x, y: T.init_y })
    const sim = new SandSim(surf, st)
    // Same step count as the Python fixture: nsteps = round(t_max / dt).
    const nsteps = Math.round(0.3 / sim.dt)
    for (let i = 0; i < nsteps; i++) sim.step()
    const xs: number[] = []; const ys: number[] = []; const sit: boolean[] = []
    for (let i = 0; i < st.id.length; i++) {
      if (!st.alive[i]) continue
      xs.push(st.x[i]); ys.push(st.y[i]); sit.push(!!st.sitting[i])
    }
    expect(xs.length).toBe(T.x.length)
    expect(sim.nBounce).toBe(T.n_bounce)
    const dx = Math.max(...xs.map((v, i) => Math.abs(v - T.x[i])))
    const dy = Math.max(...ys.map((v, i) => Math.abs(v - T.y[i])))
    expect(dx).toBeLessThan(1e-9)
    expect(dy).toBeLessThan(1e-9)
    expect(sit).toEqual(T.sitting)
  })
})

describe('single-body mechanics', () => {
  it('free flight follows the ballistic parabola', () => {
    const surf = new SingleModeSurface(80)
    const st = initialState(surf.a, 0, 0, { x: [0.01], y: [0] })
    st.z[0] = 0.05
    const sim = new SandSim(surf, st)
    for (let i = 0; i < 10; i++) sim.step()
    expect(Math.abs(st.z[0] - (0.05 - 0.5 * G * sim.t * sim.t))).toBeLessThan(1e-12)
  })

  it('surface sample obeys wtt = -omega^2 w and keeps sign', () => {
    const surf = new SingleModeSurface(80)
    const s = surf.sample(0.03, 0.01, 0.001)
    expect(Math.abs(s.wtt_ms2 + surf.omega ** 2 * s.w_m)).toBeLessThan(1e-9)
    expect(surf.sample(0.2, 0, 0).inside).toBe(false)
    expect(Math.abs(surf.fHz - 152)).toBeLessThan(0.5)
  })

  it('gamma 0.8 control stays trapped with no losses', () => {
    const surf = new SingleModeSurface(0.8)
    const sim = new SandSim(surf, initialState(surf.a, 200, 3))
    sim.advanceTo(0.5)
    expect(sim.nAlive).toBe(200)
    expect(sim.stats().trapped).toBeGreaterThan(0.9)
  })
})
