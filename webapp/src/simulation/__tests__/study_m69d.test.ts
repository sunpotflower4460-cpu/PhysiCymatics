// M6-9R follow-ups (skipped unless PHYSI_BENCH=1):
//  D1: where does the spp=44 vs spp=176 particle evolution first
//      diverge (locked-increment comparison)?
//  D2: radial-bin audit — height jumps at bin boundaries and whether
//      the gradient used for contact normals matches the implemented
//      (piecewise-constant-in-r) height field.
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientSandSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

benchIt('M6-9R D1 first divergence spp44 vs spp176 + chaos control', () => {
  const F = calibratedForceN()
  const N = 60
  const cmp = (a: TransientSandSim, b: TransientSandSim) => {
    let stateDiv = false; let maxD = 0
    for (let i = 0; i < N; i++) {
      if (a.state.sitting[i] !== b.state.sitting[i] || a.state.alive[i] !== b.state.alive[i]) stateDiv = true
      maxD = Math.max(maxD, Math.abs(a.state.x[i] - b.state.x[i]), Math.abs(a.state.y[i] - b.state.y[i]), Math.abs(a.state.z[i] - b.state.z[i]))
    }
    return { stateDiv, maxD }
  }
  for (const freq of [649, 1241.6]) {
    const mk = (spp: number, dz = 0) => {
      const plate = new TransientPlate(singleToneScore(freq, F, 30))
      const st = initialState(plate.a, N, 12345)
      if (dz) for (let i = 0; i < N; i++) st.z[i] += dz
      return new TransientSandSim(plate, st, { ...DEFAULT_PARAMS, spp }, { enabled: true })
    }
    // Three pairs advanced in 2ms locked increments:
    //  P44: spp44 vs spp176 (surface representation differs)
    //  P88: spp88 vs spp176 (smaller surface difference)
    //  PCH: spp176 vs spp176 with all z perturbed by +6e-7 m
    //       (chaos control: same-size perturbation, identical plate)
    const pairs: [string, TransientSandSim, TransientSandSim][] = [
      ['44v176', mk(44), mk(176)],
      ['88v176', mk(88), mk(176)],
      ['chaos176', mk(176), mk(176, 6e-7)],
    ]
    // Stage 1: 6ms locked increments (incremental advanceTo has a
    // per-call overhead of ~3-5x, measured separately, so keep calls few).
    const dT = 6e-3; const T = 0.12
    const firstPos: Record<string, number> = {}; const firstState: Record<string, number> = {}
    for (let k = 1; k <= Math.round(T / dT); k++) {
      const t = k * dT
      for (const [name, a, b] of pairs) {
        a.advanceTo(t); b.advanceTo(t)
        const { stateDiv, maxD } = cmp(a, b)
        if (!(name in firstPos) && maxD > 1e-6) firstPos[name] = t
        if (!(name in firstState) && stateDiv) firstState[name] = t
        if (k % 5 === 0) console.log(`D1 f=${freq} ${name} t=${t.toFixed(3)} maxPosDiff=${maxD.toExponential(2)} bounce=${a.nBounce}/${b.nBounce} detach=${a.nDetach}/${b.nDetach}`)
      }
    }
    for (const [name, a, b] of pairs) {
      console.log(`D1 f=${freq} ${name}: first|dpos|>1um window<=${firstPos[name] ?? -1}s firstStateDiv window<=${firstState[name] ?? -1}s final bounce=${a.nBounce}/${b.nBounce} detach=${a.nDetach}/${b.nDetach} caps=${a.nEventCap}/${b.nEventCap}`)
    }
    // Stage 2: refine onset (state divergence) to +-0.5ms with fresh
    // pairs: one single advanceTo to the window start, then fine steps.
    const specs: Record<string, [number, number, number, number]> = {
      '44v176': [44, 0, 176, 0], '88v176': [88, 0, 176, 0], 'chaos176': [176, 0, 176, 6e-7],
    }
    for (const name of Object.keys(specs)) {
      const wEnd = firstState[name]
      if (wEnd === undefined) { console.log(`D1-refine f=${freq} ${name}: no state divergence up to ${T}s`); continue }
      const wStart = Math.max(0, wEnd - dT)
      const [s1, z1, s2, z2] = specs[name]
      const a = mk(s1, z1); const b = mk(s2, z2)
      a.advanceTo(wStart); b.advanceTo(wStart)
      let onset = -1
      for (let t = wStart + 5e-4; t <= wEnd + 1e-9; t += 5e-4) {
        a.advanceTo(t); b.advanceTo(t)
        if (cmp(a, b).stateDiv) { onset = t; break }
      }
      console.log(`D1-refine f=${freq} ${name}: state divergence onset t~${onset.toFixed(4)}s (window ${wStart.toFixed(3)}-${wEnd.toFixed(3)}) bounceAtOnset=${a.nBounce}/${b.nBounce} detachAtOnset=${a.nDetach}/${b.nDetach}`)
    }
  }
}, 900000)

benchIt('M6-9R D2 radial bin + gradient consistency audit', () => {
  const plate = new TransientPlate(singleToneScore(649, calibratedForceN(), 30))
  let seed = 222
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  const dt = 1 / (44 * 649)
  plate.beginStep(0, dt)
  const tau = dt * 0.5
  const out = new Float64Array(5)
  const ng1 = (plate as unknown as { ng1: number }).ng1 // boundaries at r = k*a/ng1
  const ngrid = ng1
  const a = plate.a
  // Height jump across each radial bin boundary, several angles.
  let maxJump = 0; let jumpAt = ''
  for (const th of [0, 0.3, 1.0, 2.2]) {
    for (let bIdx = 1; bIdx < ngrid; bIdx++) {
      const rB = (bIdx / ngrid) * a
      const eps = a * 1e-7
      plate.surfAt((rB - eps) * Math.cos(th), (rB - eps) * Math.sin(th), tau, out); const hIn = out[0]
      plate.surfAt((rB + eps) * Math.cos(th), (rB + eps) * Math.sin(th), tau, out); const hOut = out[0]
      if (Math.abs(hOut - hIn) > maxJump) { maxJump = Math.abs(hOut - hIn); jumpAt = `r=${rB.toFixed(4)} th=${th}` }
    }
  }
  console.log(`D2 bin-boundary height jump: max=${maxJump.toExponential(3)} m at ${jumpAt} (ngrid=${ngrid}, bin width=${(a / ngrid * 1000).toFixed(3)}mm)`)
  // Gradient consistency: implemented h is exactly constant in r inside
  // a bin (FD radial slope = 0), while surfAt's gradient uses the dR
  // table. Tangential (theta) dependence is analytic and should match FD.
  let maxRadialGrad = 0; let maxTangErr = 0; let maxGrad = 0
  for (const th of [0.2, 0.9, 1.7, 2.8]) {
    for (let bIdx = 2; bIdx < ngrid - 1; bIdx += 7) {
      const r = ((bIdx + 0.5) / ngrid) * a
      const x = r * Math.cos(th); const y = r * Math.sin(th)
      plate.surfAt(x, y, tau, out)
      const hx = out[1]; const hy = out[2]
      const gRad = hx * Math.cos(th) + hy * Math.sin(th)
      const gTan = -hx * Math.sin(th) + hy * Math.cos(th)
      maxRadialGrad = Math.max(maxRadialGrad, Math.abs(gRad))
      maxGrad = Math.max(maxGrad, Math.hypot(hx, hy))
      // FD tangential slope of implemented h (arc-length derivative)
      const dth = 1e-6
      plate.surfAt(r * Math.cos(th + dth), r * Math.sin(th + dth), tau, out); const hP = out[0]
      plate.surfAt(r * Math.cos(th - dth), r * Math.sin(th - dth), tau, out); const hM = out[0]
      const fdTan = (hP - hM) / (2 * dth * r)
      maxTangErr = Math.max(maxTangErr, Math.abs(fdTan - gTan))
      // FD radial slope inside the bin (must be exactly 0 if piecewise const)
      const dr = (a / ngrid) * 0.1
      plate.surfAt((r + dr) * Math.cos(th), (r + dr) * Math.sin(th), tau, out); const hR1 = out[0]
      plate.surfAt((r - dr) * Math.cos(th), (r - dr) * Math.sin(th), tau, out); const hR0 = out[0]
      if (hR1 !== hR0) console.log(`D2 NOTE radial FD nonzero inside bin at r=${r}: ${hR1 - hR0}`)
    }
  }
  console.log(`D2 gradient: max|grad|=${maxGrad.toExponential(3)}, max|radial component|=${maxRadialGrad.toExponential(3)} vs implemented-h intra-bin radial slope = 0 exactly; tangential FD vs surfAt max err=${maxTangErr.toExponential(3)}`)
  console.log('D2 interpretation: contact normals use the analytic eigenfunction derivative (dR table), NOT the derivative of the implemented height (0 inside bins, undefined at edges).')
}, 300000)
