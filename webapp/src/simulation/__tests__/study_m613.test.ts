// M6-13 studies (skipped unless PHYSI_BENCH=1):
//  G1: contact-geometry audit — implemented height discontinuities at
//      radial bin edges, and implemented gradient (dR samples) vs the
//      finite-difference gradient of the implemented height itself.
//  G2: slope scale on the driven plate — whether |grad h| ever comes
//      near static-friction angles (decides if seated sliding can
//      matter at calibrated drive at all).
//  D1: off-center drive long run (current model, no friction):
//      non-axisymmetric excitation and emergent sand distribution.
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientPlate } from '../plateTransient'
import { createTransientSim, calibratedForceN, TransientSandSim } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function excited(freq: number, drive: [number, number]): TransientPlate {
  const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), drive)
  let seed = 222
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  return plate
}

benchIt('M6-13 G1 geometry audit: bin-edge jumps and gradient consistency', () => {
  const plate = excited(152, [0.05, 0.03])
  const dt = 1 / (44 * 152)
  plate.beginStep(0, dt)
  const tau = dt / 2
  const out = new Float64Array(5)
  const ng1 = 511
  const th = 0.3; const cx = Math.cos(th); const cy = Math.sin(th)
  // Height jump across every radial bin edge along one ray.
  let maxJump = 0; let sumJump = 0; let nJump = 0
  for (let k = 1; k < ng1; k++) {
    const re = (k / ng1) * plate.a
    plate.surfAt((re - 1e-9) * cx, (re - 1e-9) * cy, tau, out); const hL = out[0]
    plate.surfAt((re + 1e-9) * cx, (re + 1e-9) * cy, tau, out); const hR = out[0]
    const j = Math.abs(hR - hL)
    maxJump = Math.max(maxJump, j); sumJump += j; nJump++
  }
  console.log(`G1 bin-edge height jumps (excited plate): max=${(maxJump * 1e6).toFixed(2)} um mean=${((sumJump / nJump) * 1e6).toFixed(3)} um over ${nJump} edges`)
  // Implemented gradient vs finite-difference gradient of implemented
  // height, at random interior points (delta = 2 um).
  let seed = 99
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  const d = 2e-6
  let maxAng = 0; let sumRel = 0; let n = 0; let nearZeroRadialFD = 0
  for (let k = 0; k < 3000; k++) {
    const x = (rnd() * 2 - 1) * plate.a * 0.9
    const y = (rnd() * 2 - 1) * plate.a * 0.9
    if (Math.hypot(x, y) > plate.a * 0.9) continue
    plate.surfAt(x, y, tau, out)
    const gx = out[1]; const gy = out[2]
    const hxFD = (plate.heightAt(x + d, y, tau) - plate.heightAt(x - d, y, tau)) / (2 * d)
    const hyFD = (plate.heightAt(x, y + d, tau) - plate.heightAt(x, y - d, tau)) / (2 * d)
    const gi = Math.hypot(gx, gy); const gf = Math.hypot(hxFD, hyFD)
    if (gi < 1e-9) continue
    const cosA = Math.min(1, Math.max(-1, (gx * hxFD + gy * hyFD) / (gi * (gf || 1e-30))))
    maxAng = Math.max(maxAng, Math.acos(cosA) * 180 / Math.PI)
    sumRel += Math.abs(gi - gf) / gi
    // Radial FD component ~0 while implemented radial gradient is not:
    const rr = Math.hypot(x, y) || 1
    const radFD = (hxFD * x + hyFD * y) / rr
    const radImpl = (gx * x + gy * y) / rr
    if (Math.abs(radFD) < 0.05 * Math.abs(radImpl)) nearZeroRadialFD++
    n++
  }
  console.log(`G1 gradient consistency at ${n} points: max angle(impl, FD)=${maxAng.toFixed(1)} deg, mean | |g_impl|-|g_FD| |/|g_impl|=${(sumRel / n).toFixed(3)}, radial-FD~0 fraction=${(nearZeroRadialFD / n).toFixed(3)}`)
}, 600_000)

benchIt('M6-13 G2 slope scale on the driven plate', () => {
  for (const freq of [152, 649]) {
    const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), [0.05, 0.03])
    const dt = 1 / (44 * freq)
    const warm = Math.round(0.2 / dt)
    for (let k = 0; k < warm; k++) plate.beginStep(k * dt, (k + 1) * dt)
    const out = new Float64Array(5)
    const mags: number[] = []
    for (let k = 0; k < 44; k++) {
      const t0 = (warm + k) * dt
      plate.beginStep(t0, t0 + dt)
      for (let gy = 0; gy < 32; gy++) {
        for (let gx = 0; gx < 32; gx++) {
          const x = ((gx + 0.5) * 2 / 32 - 1) * plate.a
          const y = ((gy + 0.5) * 2 / 32 - 1) * plate.a
          if (Math.hypot(x, y) > plate.a) continue
          plate.surfAt(x, y, t0 + dt / 2, out)
          mags.push(Math.hypot(out[1], out[2]))
        }
      }
    }
    mags.sort((a, b) => a - b)
    const pct = (p: number) => mags[Math.min(mags.length - 1, Math.floor(p * mags.length))]
    console.log(`G2 f=${freq} |grad h| over driven plate: median=${pct(0.5).toExponential(2)} p99=${pct(0.99).toExponential(2)} max=${mags[mags.length - 1].toExponential(2)} (slope angle max=${(Math.atan(mags[mags.length - 1]) * 180 / Math.PI).toFixed(3)} deg; tan(theta) vs typical muS 0.3-0.6)`)
  }
}, 600_000)

benchIt('M6-13 D1 off-center drive long run (current model, 300 grains)', () => {
  const GRID = 32
  const o = new Float64Array(6)
  const densityOf = (sim: TransientSandSim): Float64Array => {
    const den = new Float64Array(GRID * GRID)
    for (let i = 0; i < sim.state.id.length; i++) {
      if (!sim.state.alive[i]) continue
      sim.observedState(i, o)
      const gx = Math.min(GRID - 1, Math.max(0, ((o[0] / sim.plate.a + 1) * GRID / 2) | 0))
      const gy = Math.min(GRID - 1, Math.max(0, ((o[1] / sim.plate.a + 1) * GRID / 2) | 0))
      den[gy * GRID + gx]++
    }
    return den
  }
  for (const drive of [[0, 0], [0.05, 0.03]] as [number, number][]) {
    for (const freq of [152, 649]) {
      const sim = createTransientSim(singleToneScore(freq, calibratedForceN(drive), 30), 300, 42, drive)
      const dt = sim.dtNow()
      const steps = Math.min(Math.round(3 / dt), 80_000)
      const w0 = Date.now()
      for (let s = 0; s < steps; s++) sim.step()
      const den = densityOf(sim)
      // Angular anisotropy of the final density: dipole moment.
      let mx = 0; let my = 0; let alive = 0
      for (let i = 0; i < sim.state.id.length; i++) {
        if (!sim.state.alive[i]) continue
        sim.observedState(i, o); mx += o[0]; my += o[1]; alive++
      }
      // Activity map anisotropy for reference (plate-only, same drive).
      const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(drive), 30), drive)
      const pdt = 1 / (44 * freq)
      const warm = Math.round(0.2 / pdt)
      for (let k = 0; k < warm; k++) plate.beginStep(k * pdt, (k + 1) * pdt)
      const out = new Float64Array(5)
      let actDipX = 0; let actDipY = 0; let actSum = 0
      for (let k = 0; k < 44; k++) {
        const t0 = (warm + k) * pdt
        plate.beginStep(t0, t0 + pdt)
        for (let gy = 0; gy < GRID; gy++) {
          for (let gx = 0; gx < GRID; gx++) {
            const x = ((gx + 0.5) * 2 / GRID - 1) * plate.a
            const y = ((gy + 0.5) * 2 / GRID - 1) * plate.a
            if (Math.hypot(x, y) > plate.a) continue
            plate.surfAt(x, y, t0 + pdt / 2, out)
            const w = Math.abs(out[4]) / 44
            actDipX += w * x; actDipY += w * y; actSum += w
          }
        }
      }
      void den
      console.log(`D1 drive=(${drive[0]},${drive[1]}) f=${freq} t=${sim.t.toFixed(2)}s wall=${((Date.now() - w0) / 1000).toFixed(1)}s alive=${alive} bounce=${sim.nBounce} detach=${sim.nDetach} sandCentroid=(${(mx / alive).toExponential(2)},${(my / alive).toExponential(2)}) activityCentroid=(${(actDipX / actSum).toExponential(2)},${(actDipY / actSum).toExponential(2)})`)
    }
  }
}, 1_800_000)
