// M6-27 stage-A scattering diagnostics (PHYSI_BENCH=1):
// single-impact angular scatter vs aspect ratio, and
// multi-bounce direction persistence on a resting flat plate
// (flight + impact only; no sustained contact model yet).
// Docs: docs/m6_27_design.md, docs/m6_27_results.md.
import { describe, expect, it } from 'vitest'
import {
  dumbbellFlightStep, dumbbellImpact, lobeCenters,
  type DumbbellParams, type DumbbellState,
} from '../dumbbell'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const R = 2.5e-4
const mk = (dOverR: number, e = 0.3): DumbbellParams => ({
  lobeRadius: R, halfSeparation: dOverR * R, beta: 2 / 5,
  muS: 0.3, muK: 0.2, restitution: e,
})

const orientations: [number, number, number][] = []
for (const polar of [0.2, 0.5, 0.8, 1.1]) {
  for (const azim of [0, 1.05, 2.1, 3.14, 4.19, 5.24]) {
    orientations.push([
      Math.sin(polar) * Math.cos(azim),
      Math.sin(polar) * Math.sin(azim),
      Math.cos(polar),
    ])
  }
}

function circStats(angles: number[]): { meanResultant: number; circStd: number } {
  const mx = angles.reduce((s, a) => s + Math.cos(a), 0) / angles.length
  const my = angles.reduce((s, a) => s + Math.sin(a), 0) / angles.length
  const mr = Math.hypot(mx, my)
  return { meanResultant: mr, circStd: Math.sqrt(-2 * Math.log(Math.max(1e-12, mr))) }
}

describe('M6-27 scattering diagnostics', () => {
  benchIt('single flat impact: outgoing direction scatter vs aspect', () => {
    for (const dOverR of [0, 0.5, 1.0]) {
      const p = mk(dOverR)
      const angles: number[] = []
      const speeds: number[] = []
      const spins: number[] = []
      for (const e of orientations) {
        const st: DumbbellState = {
          cx: 0, cy: 0, cz: R + p.halfSeparation * Math.abs(e[2]),
          vx: 0.25, vy: 0, vz: -0.6,
          ex: e[0], ey: e[1], ez: e[2], wx: 0, wy: 0, wz: 0,
        }
        const lobes = lobeCenters(st, p)
        const lower = lobes[0][2] <= lobes[1][2] ? lobes[0] : lobes[1]
        dumbbellImpact(st, [lower[0], lower[1], lower[2] - R], [0, 0, 1], [0, 0, 0], p)
        angles.push(Math.atan2(st.vy, st.vx))
        speeds.push(Math.hypot(st.vx, st.vy))
        spins.push(Math.hypot(st.wx, st.wy, st.wz))
      }
      const cs = circStats(angles)
      console.log(JSON.stringify({
        study: 'm627-scatter1', dOverR, aspect: 1 + dOverR,
        meanResultant: +cs.meanResultant.toFixed(4),
        circStdRad: +cs.circStd.toFixed(4),
        meanSpeedOut: +(speeds.reduce((s, v) => s + v, 0) / speeds.length).toFixed(4),
        meanSpinOut: +(spins.reduce((s, v) => s + v, 0) / spins.length).toFixed(1),
      }))
      expect(speeds.every((v) => Number.isFinite(v))).toBe(true)
    }
  }, 600000)

  benchIt('multi-bounce on resting flat plate: direction persistence vs aspect', () => {
    const dt = 2e-5
    const T = 1.0
    for (const dOverR of [0, 0.5, 1.0]) {
      const p = mk(dOverR, 0.8)
      const cosByBounce: number[][] = []
      let totalBounces = 0
      const finalAngles: number[] = []
      for (const e of orientations) {
        const st: DumbbellState = {
          cx: 0, cy: 0, cz: 0.005 + R + p.halfSeparation * Math.abs(e[2]),
          vx: 0.25, vy: 0, vz: 0,
          ex: e[0], ey: e[1], ez: e[2], wx: 0, wy: 0, wz: 0,
        }
        const gapOf = (s: DumbbellState): number => {
          const lobes = lobeCenters(s, p)
          return Math.min(lobes[0][2], lobes[1][2]) - R
        }
        let bounces = 0
        const steps = Math.round(T / dt)
        for (let k = 0; k < steps && bounces < 300; k++) {
          const saved: DumbbellState = { ...st }
          dumbbellFlightStep(st, p, dt)
          if (gapOf(st) > 0) continue
          // Bisect this step from the saved state to first touch.
          let lo = 0; let hi = 1
          for (let it = 0; it < 40; it++) {
            const mid = 0.5 * (lo + hi)
            const probe: DumbbellState = { ...saved }
            dumbbellFlightStep(probe, p, dt * mid)
            if (gapOf(probe) > 0) lo = mid; else hi = mid
          }
          Object.assign(st, saved)
          dumbbellFlightStep(st, p, dt * hi)
          const lobes = lobeCenters(st, p)
          const li = lobes[0][2] <= lobes[1][2] ? 0 : 1
          const lower = lobes[li]
          const r = dumbbellImpact(
            st, [lower[0], lower[1], lower[2] - R], [0, 0, 1], [0, 0, 0], p)
          if (r.mode === 'separating') {
            // Numeric overlap while separating: step over it.
            dumbbellFlightStep(st, p, dt * (1 - hi))
            continue
          }
          bounces++
          const sp = Math.hypot(st.vx, st.vy)
          const cosA = sp > 1e-9 ? st.vx / sp : 1
          if (!cosByBounce[bounces - 1]) cosByBounce[bounces - 1] = []
          cosByBounce[bounces - 1].push(cosA)
          dumbbellFlightStep(st, p, dt * (1 - hi))
        }
        totalBounces += bounces
        finalAngles.push(Math.atan2(st.vy, st.vx))
      }
      const meanCos = cosByBounce.slice(0, 10).map((arr) =>
        +(arr.reduce((s, v) => s + v, 0) / arr.length).toFixed(3))
      const cs = circStats(finalAngles)
      console.log(JSON.stringify({
        study: 'm627-multibounce', dOverR, aspect: 1 + dOverR,
        meanBounces: +(totalBounces / orientations.length).toFixed(1),
        meanCosVsInitialByBounce: meanCos,
        finalMeanResultant: +cs.meanResultant.toFixed(4),
        finalCircStdRad: +cs.circStd.toFixed(4),
      }))
      expect(totalBounces).toBeGreaterThan(0)
    }
  }, 1800000)
})
