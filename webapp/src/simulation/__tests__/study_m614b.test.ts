// M6-14B Phase D: sand A/B — production bin geometry (A) vs the
// geometrically consistent surface (B), identical initial states.
// Sleep and the local non-contact proof are DISABLED on both sides:
// their guarantees are derived for the bin geometry and are not yet
// re-derived for B (docs/m6_14_results.md). This is an exploratory
// research comparison, not a production change.
// Run: PHYSI_BENCH=1 npx vitest run .../study_m614b.test.ts
import { describe, it } from 'vitest'
import { singleToneScore } from '../sequence'
import { createTransientSim, calibratedForceN, type TransientSandSim } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip
const N = 300
const T = 1.0

function run(geo: 'bin' | 'consistent', freq: number, drive: [number, number], seed: number) {
  const sim = createTransientSim(
    singleToneScore(freq, calibratedForceN(), 30), N, seed, drive,
    { enabled: false, localBound: false })
  sim.plate.geometryMode = geo
  return sim
}

function hist(sim: TransientSandSim): { radial: number[]; angular: number[]; path: number[] } {
  const radial = new Array(12).fill(0) as number[]
  const angular = new Array(12).fill(0) as number[]
  const out = new Float64Array(8)
  const path: number[] = []
  const alive = (sim as unknown as { state: { alive: Uint8Array } }).state.alive
  for (let i = 0; i < N; i++) {
    if (!alive[i]) continue
    sim.observedState(i, out)
    const r = Math.hypot(out[0], out[1]) / sim.plate.a
    radial[Math.min(11, (r * 12) | 0)]++
    const th = Math.atan2(out[1], out[0])
    angular[Math.min(11, (((th + Math.PI) / (2 * Math.PI)) * 12) | 0)]++
  }
  return { radial, angular, path }
}

function corr(a: number[], b: number[]): number {
  const ma = a.reduce((s, v) => s + v, 0) / a.length
  const mb = b.reduce((s, v) => s + v, 0) / b.length
  let cov = 0; let va = 0; let vb = 0
  for (let i = 0; i < a.length; i++) {
    cov += (a[i] - ma) * (b[i] - mb); va += (a[i] - ma) ** 2; vb += (b[i] - mb) ** 2
  }
  return va > 0 && vb > 0 ? cov / Math.sqrt(va * vb) : 0
}

describe('M6-14B sand A/B geometry comparison', () => {
  benchIt('152/649Hz x center/eccentric x 2 seeds, 300 grains, 1s', () => {
    for (const freq of [152, 649]) {
      for (const drive of [[0, 0], [0.06, 0]] as [number, number][]) {
        for (const seed of [42, 7]) {
          const t0 = performance.now()
          const A = run('bin', freq, drive, seed)
          const B = run('consistent', freq, drive, seed)
          // Step both in lockstep; find first divergence > 1 nm.
          let firstDiv = -1
          const oa = new Float64Array(8); const ob = new Float64Array(8)
          const steps = 100 // check every 10 ms
          for (let s = 1; s <= steps; s++) {
            A.advanceTo((T * s) / steps)
            B.advanceTo((T * s) / steps)
            if (firstDiv < 0) {
              for (let i = 0; i < N; i++) {
                A.observedState(i, oa); B.observedState(i, ob)
                if (Math.abs(oa[0] - ob[0]) > 1e-9 || Math.abs(oa[1] - ob[1]) > 1e-9
                  || Math.abs(oa[2] - ob[2]) > 1e-9 ) {
                  firstDiv = (T * s) / steps; break
                }
              }
            }
          }
          const sa = A.stats(); const sb = B.stats()
          const ha = hist(A); const hb = hist(B)
          // Mean horizontal displacement from start is not tracked;
          // report distribution correlations and event counts.
          console.log(JSON.stringify({
            freq, drive, seed,
            firstDivergence_s: firstDiv,
            bounceA: sa.nBounce, bounceB: sb.nBounce,
            aliveA: sa.alive, aliveB: sb.alive,
            sittingA: sa.sitting, sittingB: sb.sitting,
            radialCorr: Number(corr(ha.radial, hb.radial).toFixed(3)),
            angularCorr: Number(corr(ha.angular, hb.angular).toFixed(3)),
            radialA: ha.radial, radialB: hb.radial,
            wallMsPair: Math.round(performance.now() - t0),
          }))
        }
      }
    }
  }, 1800000)
})
