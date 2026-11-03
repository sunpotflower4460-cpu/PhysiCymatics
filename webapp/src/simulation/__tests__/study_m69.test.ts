// M6-9 cause-isolation studies. Skipped unless PHYSI_BENCH=1.
//  A: Hermite segment-length error of the plate surface (plate only)
//  B: 8-sample contact detection vs dense reference within ONE segment
//  C: first-bounce time vs time-step resolution (single particles)
//  D: radial-bin discontinuity sizes of the implemented surface
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientSandSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'
import { G } from '../surface'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function excitedPlate(freq: number, seedIn: number): TransientPlate {
  const plate = new TransientPlate(singleToneScore(freq, 0, 30))
  let seed = seedIn
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  return plate
}

benchIt('M6-9 experiments A/B/C/D', () => {
  const F = calibratedForceN()
  void F

  // ---- D: radial bin jumps of the implemented surface ----
  {
    const plate = excitedPlate(649, 111)
    const outA = new Float64Array(5); const outB = new Float64Array(5)
    plate.beginStep(0, 1.5e-4)
    let worst = 0; let worstInfo = ''
    const ng1 = 511
    for (let k = 1; k < ng1; k++) {
      const rEdge = (k * plate.a) / ng1
      for (const th of [0, 0.7, 1.9, 3.3]) {
        const x0 = (rEdge - 1e-12) * Math.cos(th); const y0 = (rEdge - 1e-12) * Math.sin(th)
        const x1 = (rEdge + 1e-12) * Math.cos(th); const y1 = (rEdge + 1e-12) * Math.sin(th)
        plate.surfAt(x0, y0, 0.7e-4, outA); plate.surfAt(x1, y1, 0.7e-4, outB)
        const jump = Math.abs(outB[0] - outA[0])
        if (jump > worst) { worst = jump; worstInfo = `bin=${k} th=${th}` }
      }
    }
    console.log(`D bin-jump: worst height discontinuity=${worst.toExponential(3)} m (${worstInfo}) at |q|~2e-4 all-mode excitation`)
    // Per-mode max |dR| between adjacent bins from the mode table.
    const modes = plate.modes
    const per = modes.map((m) => {
      let mx = 0
      for (let k = 0; k + 1 < m.R.length; k++) mx = Math.max(mx, Math.abs(m.R[k + 1] - m.R[k]))
      return `m${m.m}r${m.root_idx}:${mx.toFixed(4)}`
    })
    console.log(`D per-mode max |dR|: ${per.join(' ')}`)
  }

  // ---- A: Hermite segment-length error (plate only, free decay) ----
  for (const freq of [152, 1241.6]) {
    const mkPlate = () => excitedPlate(freq, 222)
    const sampleHeights = (spp: number): Float64Array => {
      const plate = mkPlate()
      const dt = 1 / (spp * freq)
      const out = new Float64Array(5)
      const res: number[] = []
      let seed = 333
      const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
      const pts: [number, number][] = []
      for (let i = 0; i < 40; i++) pts.push([(rnd() * 2 - 1) * plate.a * 0.95, (rnd() * 2 - 1) * plate.a * 0.95])
      const segs = Math.round(0.01 / dt)
      for (let sIdx = 0; sIdx < segs; sIdx++) {
        plate.beginStep(sIdx * dt, (sIdx + 1) * dt)
        if (sIdx % Math.max(1, Math.round(0.002 / dt)) === 0) {
          for (const fr of [0.23, 0.51, 0.87]) {
            for (const [x, y] of pts) { plate.surfAt(x, y, sIdx * dt + fr * dt, out); res.push(out[0], out[4]) }
          }
        }
      }
      return Float64Array.from(res)
    }
    const ref = sampleHeights(352)
    for (const spp of [44, 88, 176]) {
      const got = sampleHeights(spp)
      let maxH = 0; let maxW = 0
      for (let i = 0; i < Math.min(ref.length, got.length); i += 2) {
        maxH = Math.max(maxH, Math.abs(ref[i] - got[i]))
        maxW = Math.max(maxW, Math.abs(ref[i + 1] - got[i + 1]))
      }
      console.log(`A freq=${freq} spp=${spp}: max|dh|=${maxH.toExponential(3)} m max|dwtt|=${maxW.toExponential(3)} m/s^2 (vs spp=352, sampled at common nominal times)`)
    }
  }

  // ---- B: 8-sample detection vs dense reference in ONE segment ----
  for (const freq of [152, 649, 1241.6]) {
    const plate = excitedPlate(freq, 444)
    const dt = 1 / (44 * freq)
    plate.beginStep(0, dt)
    let seed = 555
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    let trials = 0; let trueContacts = 0; let missed = 0; let jumpsSeen = 0
    for (let n = 0; n < 4000; n++) {
      const x0 = rnd() * 0.1; const y0 = rnd() * 0.1
      const vx = rnd() * 0.6; const vy = rnd() * 0.6; const vz = rnd() * 0.4
      const h0 = plate.heightAt(x0, y0, 0)
      const z0 = h0 + Math.abs(rnd()) * 4e-4 // start above surface, near it
      const gap = (d: number): number => (z0 + vz * d - 0.5 * G * d * d) - plate.heightAt(x0 + vx * d, y0 + vy * d, d)
      if (gap(0) <= 0) continue
      trials++
      let denseMin = Infinity; let prev = gap(0)
      for (let k = 1; k <= 400; k++) {
        const g = gap((k / 400) * dt)
        if (g < denseMin) denseMin = g
        if (Math.abs(g - prev) > 2e-5) jumpsSeen++ // crude jump detector
        prev = g
      }
      let sampledMin = Infinity
      for (let k = 1; k <= 8; k++) sampledMin = Math.min(sampledMin, gap((k / 8) * dt))
      if (denseMin <= 0) { trueContacts++; if (sampledMin > 0) missed++ }
    }
    console.log(`B freq=${freq}: trials=${trials} denseContacts=${trueContacts} missedBy8Samples=${missed} jumpCandidates=${jumpsSeen}`)
  }

  // ---- C: first-bounce time vs spp (single particles, binary search) ----
  for (const freq of [152, 649, 1241.6]) {
    for (const spp of [44, 88, 176]) {
      const times: number[] = []
      for (let seedI = 1; seedI <= 6; seedI++) {
        const bouncedBy = (tEnd: number): boolean => {
          const plate = new TransientPlate(singleToneScore(freq, F, 30))
          const sim = new TransientSandSim(plate, initialState(plate.a, 1, seedI),
            { ...DEFAULT_PARAMS, spp }, { enabled: false })
          sim.advanceTo(tEnd)
          return sim.nBounce > 0
        }
        let lo = 0; let hi = 0.06
        if (!bouncedBy(hi)) continue
        for (let k = 0; k < 22; k++) {
          const mid = 0.5 * (lo + hi)
          if (bouncedBy(mid)) hi = mid; else lo = mid
        }
        times.push(hi)
      }
      const med = times.length ? [...times].sort((a, b) => a - b)[Math.floor(times.length / 2)] : NaN
      console.log(`C freq=${freq} spp=${spp}: first-bounce median=${(med * 1e3).toFixed(3)}ms times=[${times.map((t) => (t * 1e3).toFixed(2)).join(',')}]`)
    }
  }
}, 3600000)
