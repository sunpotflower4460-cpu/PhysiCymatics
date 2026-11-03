// M6-9 follow-up studies. Skipped unless PHYSI_BENCH=1.
//  WARNING (M6-9R): A2 below is FLAWED — it samples at
//  tEnd - 0.37*dt_spp, i.e. different absolute times per spp (up to
//  ~48us apart at 152Hz). Its ~140um figure is a time-offset artifact.
//  Corrected version: study_m69c.test.ts (A3, shared absolute times).
//  A2: Hermite segment error (flawed sampling, kept for the record)
//  B2: 8-sample detection vs dense reference, resonant large amplitude
//  C2: first-bounce time on a pre-excited plate vs spp
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientSandSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'
import { G } from '../surface'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

benchIt('M6-9 experiments A2/B2/C2', () => {
  const F = calibratedForceN()
  void F
  const excited = (freq: number, seedIn: number) => {
    const plate = new TransientPlate(singleToneScore(freq, 0, 30))
    let seed = seedIn
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
    return plate
  }

  // ---- A2: common sample times are multiples of dt44; fractions of
  // the spp=44 segment that are also multiples of finer dts.
  for (const freq of [152, 1241.6]) {
    const sampleAtCommon = (spp: number): Float64Array => {
      const plate = excited(freq, 222)
      const dt = 1 / (spp * freq)
      const dt44 = 1 / (44 * freq)
      const out = new Float64Array(5)
      const res: number[] = []
      let seed = 333
      const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
      const pts: [number, number][] = []
      for (let i = 0; i < 40; i++) pts.push([(rnd() * 2 - 1) * plate.a * 0.95, (rnd() * 2 - 1) * plate.a * 0.95])
      const ratio = Math.round(dt44 / dt) // 1,2,4,8
      const totalSegs44 = Math.round(0.01 / dt44)
      let segIdx = 0
      for (let s44 = 0; s44 < totalSegs44; s44++) {
        for (let sub = 0; sub < ratio; sub++) {
          plate.beginStep(segIdx * dt, (segIdx + 1) * dt)
          segIdx++
        }
        // sample inside the just-completed spp44 segment at its quarter points:
        // those times are segment boundaries of the finer runs only for
        // fr in {0, 0.25, ...} when ratio divides accordingly; use the
        // plate's CURRENT segment (last sub-segment) interior instead:
        // simpler and still aligned: sample at segment ENDS of dt44 grid
        // plus the midpoint of the last sub-segment.
        const tEnd = (s44 + 1) * dt44
        for (const [x, y] of pts) { plate.surfAt(x, y, tEnd - 0.37 * dt, out); res.push(out[0], out[4]) }
      }
      return Float64Array.from(res)
    }
    const ref = sampleAtCommon(352)
    for (const spp of [44, 88, 176]) {
      const got = sampleAtCommon(spp)
      let maxH = 0; let maxW = 0
      for (let i = 0; i < Math.min(ref.length, got.length); i += 2) {
        maxH = Math.max(maxH, Math.abs(ref[i] - got[i]))
        maxW = Math.max(maxW, Math.abs(ref[i + 1] - got[i + 1]))
      }
      console.log(`A2 freq=${freq} spp=${spp}: max|dh|=${maxH.toExponential(3)} max|dwtt|=${maxW.toExponential(3)} (vs spp=352)`)
    }
  }

  // ---- B2: resonant single-mode large amplitude, near-grazing starts ----
  for (const [freq, dofPick] of [[649, 'm0r1'], [1241.6, 'm4r1']] as const) {
    const plate = new TransientPlate(singleToneScore(freq, 0, 30))
    // Excite only the target mode's cos DOF strongly.
    const mi = plate.modes.findIndex((m) => `m${m.m}r${m.root_idx}` === dofPick)
    let dof = -1
    for (let j = 0; j < plate.ndof; j++) if ((plate as unknown as { dofMode: Int32Array; dofCos: boolean[] }).dofMode[j] === mi && (plate as unknown as { dofMode: Int32Array; dofCos: boolean[] }).dofCos[j]) { dof = j; break }
    plate.q[dof] = 8e-4
    const dt = 1 / (44 * freq)
    plate.beginStep(0, dt)
    let seed = 777
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    let trials = 0; let denseContacts = 0; let missed = 0; let minMissedDepth = Infinity
    for (let n = 0; n < 6000; n++) {
      const x0 = rnd() * 0.1; const y0 = rnd() * 0.1
      const vx = rnd() * 0.3; const vy = rnd() * 0.3; const vz = rnd() * 0.2
      const h0 = plate.heightAt(x0, y0, 0)
      const z0 = h0 + Math.abs(rnd()) * 1.5e-4
      const gap = (d: number): number => (z0 + vz * d - 0.5 * G * d * d) - plate.heightAt(x0 + vx * d, y0 + vy * d, d)
      if (gap(0) <= 0) continue
      trials++
      let denseMin = Infinity
      for (let k = 1; k <= 800; k++) denseMin = Math.min(denseMin, gap((k / 800) * dt))
      let sampledMin = Infinity
      for (let k = 1; k <= 8; k++) sampledMin = Math.min(sampledMin, gap((k / 8) * dt))
      if (denseMin <= 0) {
        denseContacts++
        if (sampledMin > 0) { missed++; if (-denseMin < minMissedDepth) minMissedDepth = -denseMin }
      }
    }
    console.log(`B2 freq=${freq} mode=${dofPick}: trials=${trials} denseContacts=${denseContacts} missed=${missed} shallowestMissedDepth=${missed ? minMissedDepth.toExponential(2) : 'n/a'}`)
  }

  // ---- C2: first bounce on pre-excited plate (force off, decaying) ----
  for (const [freq, dofPick] of [[649, 'm0r1'], [1241.6, 'm4r1']] as const) {
    for (const spp of [44, 176]) {
      const times: number[] = []
      for (let seedI = 1; seedI <= 6; seedI++) {
        const bouncedBy = (tEnd: number): boolean => {
          const plate = new TransientPlate(singleToneScore(freq, 0, 30))
          const mi = plate.modes.findIndex((m) => `m${m.m}r${m.root_idx}` === dofPick)
          for (let j = 0; j < plate.ndof; j++) if ((plate as unknown as { dofMode: Int32Array; dofCos: boolean[] }).dofMode[j] === mi && (plate as unknown as { dofMode: Int32Array; dofCos: boolean[] }).dofCos[j]) plate.q[j] = 8e-4
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
      console.log(`C2 freq=${freq} spp=${spp}: firstBounce=[${times.map((t) => (t * 1e3).toFixed(3)).join(',')}] ms`)
    }
  }
}, 3600000)
