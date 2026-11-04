// M6-14 studies (skipped unless PHYSI_BENCH=1):
//  B1: contact-input differences between bin geometry (production)
//      and consistent Hermite geometry (research): height, surface
//      normal angle, and wtt over driven plates; |R| overshoot check
//      for the |Phi|<=1 height-bound assumption.
//  D2: angular structure of emergent sand distribution, center vs
//      off-center drive, multiple seeds (measurement only).
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientPlate } from '../plateTransient'
import { createTransientSim, calibratedForceN } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

benchIt('M6-14 B1 contact-input differences: bin vs consistent geometry', () => {
  const ri = new Float64Array(3)
  for (const freq of [152, 649]) {
    for (const drive of [[0, 0], [0.05, 0.03]] as [number, number][]) {
      const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(drive), 30), drive)
      const dt = 1 / (44 * freq)
      const warm = Math.round(0.2 / dt)
      for (let k = 0; k < warm; k++) plate.beginStep(k * dt, (k + 1) * dt)
      const oA = new Float64Array(5); const oB = new Float64Array(5)
      let maxDh = 0; let sumDh = 0; let maxAng = 0; let sumAng = 0; let maxDwtt = 0; let n = 0
      for (let k = 0; k < 22; k++) {
        const t0 = (warm + k) * dt
        plate.beginStep(t0, t0 + dt)
        const tau = t0 + dt / 2
        for (let gy = 0; gy < 48; gy++) {
          for (let gx = 0; gx < 48; gx++) {
            const x = ((gx + 0.5) * 2 / 48 - 1) * plate.a
            const y = ((gy + 0.5) * 2 / 48 - 1) * plate.a
            if (Math.hypot(x, y) > plate.a) continue
            plate.surfAt(x, y, tau, oA)
            plate.surfAtConsistent(x, y, tau, oB)
            const dh = Math.abs(oA[0] - oB[0])
            maxDh = Math.max(maxDh, dh); sumDh += dh
            // Angle between surface normals (-hx, -hy, 1).
            const na = Math.hypot(oA[1], oA[2], 1); const nb = Math.hypot(oB[1], oB[2], 1)
            const cosA = Math.min(1, Math.max(-1, (oA[1] * oB[1] + oA[2] * oB[2] + 1) / (na * nb)))
            const ang = Math.acos(cosA) * 180 / Math.PI
            maxAng = Math.max(maxAng, ang); sumAng += ang
            maxDwtt = Math.max(maxDwtt, Math.abs(oA[4] - oB[4]))
            n++
          }
        }
      }
      console.log(`B1 f=${freq} drive=(${drive[0]},${drive[1]}): n=${n} |dh| mean=${(sumDh / n * 1e6).toFixed(2)}um max=${(maxDh * 1e6).toFixed(2)}um, normal angle mean=${(sumAng / n).toFixed(3)}deg max=${maxAng.toFixed(2)}deg, |dwtt| max=${maxDwtt.toFixed(2)} m/s2`)
    }
  }
  // |Phi| <= 1 assumption check for the interpolant (hBound uses it).
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0, 0])
  let maxOver = 0
  for (let mi = 0; mi < plate.modes.length; mi++) {
    for (let k = 0; k <= 20000; k++) {
      plate.radialInterp(mi, k / 20000, ri)
      maxOver = Math.max(maxOver, Math.abs(ri[0]))
    }
  }
  console.log(`B1 max |R_interp| over all modes = ${maxOver.toFixed(6)} (production hBound assumes |Phi| <= 1 via max|R| normalization; overshoot beyond 1 would need re-derivation)`)
}, 900_000)

benchIt('M6-14 D2 angular structure: center vs off-center drive, 3 seeds', () => {
  const o = new Float64Array(6)
  for (const drive of [[0, 0], [0.05, 0.03]] as [number, number][]) {
    for (const freq of [152, 649]) {
      for (const seed of [42, 43, 44]) {
        const sim = createTransientSim(singleToneScore(freq, calibratedForceN(drive), 30), 300, seed, drive)
        const dt = sim.dtNow()
        const steps = Math.min(Math.round(3 / dt), 80_000)
        for (let s = 0; s < steps; s++) sim.step()
        // Angular Fourier amplitudes of final particle angles and
        // radially-binned angular amplitudes (inner/outer halves).
        const amps = [0, 0, 0, 0, 0]
        const ampsIn = [0, 0, 0, 0, 0]; const ampsOut = [0, 0, 0, 0, 0]
        const re = [0, 0, 0, 0, 0]; const im = [0, 0, 0, 0, 0]
        const reI = [0, 0, 0, 0, 0]; const imI = [0, 0, 0, 0, 0]
        const reO = [0, 0, 0, 0, 0]; const imO = [0, 0, 0, 0, 0]
        let alive = 0; let nIn = 0; let nOut = 0
        for (let i = 0; i < sim.state.id.length; i++) {
          if (!sim.state.alive[i]) continue
          sim.observedState(i, o)
          const th = Math.atan2(o[1], o[0])
          const r = Math.hypot(o[0], o[1]) / sim.plate.a
          alive++
          for (let m = 1; m <= 4; m++) {
            re[m] += Math.cos(m * th); im[m] += Math.sin(m * th)
            if (r < 0.5) { reI[m] += Math.cos(m * th); imI[m] += Math.sin(m * th) }
            else { reO[m] += Math.cos(m * th); imO[m] += Math.sin(m * th) }
          }
          if (r < 0.5) nIn++; else nOut++
        }
        for (let m = 1; m <= 4; m++) {
          amps[m] = Math.hypot(re[m], im[m]) / alive
          ampsIn[m] = nIn ? Math.hypot(reI[m], imI[m]) / nIn : 0
          ampsOut[m] = nOut ? Math.hypot(reO[m], imO[m]) / nOut : 0
        }
        const fmt = (a: number[]) => a.slice(1).map((v) => v.toFixed(3)).join(',')
        console.log(`D2 drive=(${drive[0]},${drive[1]}) f=${freq} seed=${seed} t=${sim.t.toFixed(2)} alive=${alive} |A_m| m=1..4 all=[${fmt(amps)}] inner=[${fmt(ampsIn)}] outer=[${fmt(ampsOut)}] (uniform-noise scale ~1/sqrt(N)=${(1 / Math.sqrt(alive)).toFixed(3)})`)
      }
    }
  }
}, 3_600_000)
