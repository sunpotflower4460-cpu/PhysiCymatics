// M6-11R audit studies (skipped unless PHYSI_BENCH=1):
//  R1: 500Hz detach forensics — plate-only radial scan for genuine
//      threshold crossings. Finding: none exist in 0.1 s; the 51
//      detaches reported in M6-11 were an advanceTo artifact (R1b).
//  R1b: tiny-final-step artifact: advanceTo(0.1) at exactly 2200
//      nominal steps used to append a ~1e-15 s segment whose Hermite
//      acceleration is cancellation garbage. Fixed by remnant
//      snapping in advanceTo + stable tiny-segment acceleration.
//  R2: multi-seed reproduction of the 500Hz detach counts.
//  R3: adversarial detach-root search — grazing threshold cases where
//      a scan+bisect finder can miss root pairs; compares K-scan vs
//      the interval-assisted plate finder against a dense reference.
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientSandSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'
import { G } from '../surface'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function hermite3(q0: number, qd0: number, q1: number, qd1: number, h: number, s: number): [number, number, number] {
  const s2 = s * s; const s3 = s2 * s
  const q = (2 * s3 - 3 * s2 + 1) * q0 + (s3 - 2 * s2 + s) * h * qd0 + (-2 * s3 + 3 * s2) * q1 + (s3 - s2) * h * qd1
  const dq = (6 * s2 - 6 * s) * q0 + (3 * s2 - 4 * s + 1) * h * qd0 + (-6 * s2 + 6 * s) * q1 + (3 * s2 - 2 * s) * h * qd1
  const ddq = (12 * s - 6) * q0 + (6 * s - 4) * h * qd0 + (-12 * s + 6) * q1 + (6 * s - 2) * h * qd1
  return [q, dq / h, ddq / (h * h)]
}

benchIt('M6-11R R1 spurious detach zones on the driven plate (500Hz)', () => {
  // Hermite wtt is linear in tau within a segment, so its minimum at
  // fixed (x,y) is at an endpoint: scan endpoints over fine radial
  // points (Phi is radial-bin-wise, so coarse xy grids can miss the
  // extreme radii). Zones where A dips below -G are then re-examined
  // densely under B and a 16x-subdivided hermite representation.
  const plate = new TransientPlate(singleToneScore(500, calibratedForceN(), 30), [0, 0])
  const dt = 1 / (44 * 500)
  const zeta = plate.zeta
  const signal = plate.signal
  const forceAt = (t: number) => signal.forcePeakAt(t - plate.epoch) * signal.u(t - plate.epoch)
  const odeQdd = (j: number, q: number, qd: number, t: number) =>
    (plate.phiDrive[j] * forceAt(t)) / plate.mass[j] - 2 * zeta * plate.omega[j] * qd - plate.omega[j] ** 2 * q
  type Seg = { t0: number; t1: number; q0: Float64Array; qd0: Float64Array; q1: Float64Array; qd1: Float64Array }
  const pts: [number, number][] = []
  for (const th of [0, 0.7]) for (let k = 0; k < 200; k++) {
    const r = (k / 200) * plate.a * 1.5
    pts.push([r * Math.cos(th), r * Math.sin(th)])
  }
  const phis: Float64Array[] = pts.map(([x, y]) => {
    const f = new Float64Array(plate.ndof)
    plate.phiAll(x, y, f)
    return f
  })
  const nSeg = Math.round(0.1 / dt)
  let aBelow = 0; let aBelowSegs = 0
  let first: { seg: Seg; x: number; y: number; end: number } | null = null
  const candSegs: Seg[] = []
  for (let k = 0; k < nSeg; k++) {
    plate.beginStep(k * dt, (k + 1) * dt)
    const sg = (plate as unknown as { seg: Seg }).seg
    const h = sg.t1 - sg.t0
    const qddE: [Float64Array, Float64Array] = [new Float64Array(plate.ndof), new Float64Array(plate.ndof)]
    for (let j = 0; j < plate.ndof; j++) {
      qddE[0][j] = hermite3(sg.q0[j], sg.qd0[j], sg.q1[j], sg.qd1[j], h, 0)[2]
      qddE[1][j] = hermite3(sg.q0[j], sg.qd0[j], sg.q1[j], sg.qd1[j], h, 1)[2]
    }
    let segHas = false
    for (let pi = 0; pi < pts.length; pi++) {
      const phi = phis[pi]
      for (let e = 0; e < 2; e++) {
        let w = 0
        for (let j = 0; j < plate.ndof; j++) w += qddE[e][j] * phi[j]
        if (w + G < 0) {
          aBelow++; segHas = true
          if (!first) first = { seg: { t0: sg.t0, t1: sg.t1, q0: Float64Array.from(sg.q0), qd0: Float64Array.from(sg.qd0), q1: Float64Array.from(sg.q1), qd1: Float64Array.from(sg.qd1) }, x: pts[pi][0], y: pts[pi][1], end: e }
        }
      }
    }
    if (segHas) { aBelowSegs++; if (candSegs.length < 40) candSegs.push({ t0: sg.t0, t1: sg.t1, q0: Float64Array.from(sg.q0), qd0: Float64Array.from(sg.qd0), q1: Float64Array.from(sg.q1), qd1: Float64Array.from(sg.qd1) }) }
  }
  console.log(`R1 radial scan: segments=${nSeg} points=${pts.length} A-below(endpoint,point) pairs=${aBelow} segmentsWithAnyBelow=${aBelowSegs}`)
  const phi1 = new Float64Array(plate.ndof)
  const wttSeg = (sg: Seg, t: number, mode: 'A' | 'B'): number => {
    const s = (t - sg.t0) / (sg.t1 - sg.t0)
    let w = 0
    for (let j = 0; j < plate.ndof; j++) {
      const [q, qd, qddH] = hermite3(sg.q0[j], sg.qd0[j], sg.q1[j], sg.qd1[j], sg.t1 - sg.t0, s)
      w += (mode === 'B' ? odeQdd(j, q, qd, t) : qddH) * phi1[j]
    }
    return w
  }
  // Dense check on candidate segments: does B (and fine16) also dip?
  let bothBelow = 0; let aOnlyBelow = 0; let checked = 0
  for (const seg of candSegs) {
    const fineSegs: Seg[] = []
    let q: Float64Array = Float64Array.from(seg.q0); let qd: Float64Array = Float64Array.from(seg.qd0)
    const hf = (seg.t1 - seg.t0) / 16
    for (let k = 0; k < 16; k++) {
      const [q1, qd1] = plate.advanceState(q, qd, seg.t0 + k * hf, seg.t0 + (k + 1) * hf)
      fineSegs.push({ t0: seg.t0 + k * hf, t1: seg.t0 + (k + 1) * hf, q0: q, qd0: qd, q1, qd1 })
      q = q1; qd = qd1
    }
    for (let pi = 0; pi < pts.length; pi++) {
      phi1.set(phis[pi])
      let mnA = Infinity; let mnB = Infinity; let mnF = Infinity
      for (let sN = 0; sN <= 64; sN++) {
        const t = seg.t0 + ((seg.t1 - seg.t0) * sN) / 64
        mnA = Math.min(mnA, wttSeg(seg, t, 'A') + G)
        mnB = Math.min(mnB, wttSeg(seg, t, 'B') + G)
        const idx = Math.max(0, Math.min(15, Math.ceil((t - seg.t0) / hf) - 1))
        mnF = Math.min(mnF, wttSeg(fineSegs[idx], t, 'A') + G)
      }
      if (mnA < 0) {
        checked++
        if (mnB < 0 || mnF < 0) bothBelow++
        else aOnlyBelow++
      }
    }
  }
  console.log(`R1 candidate dense check: A-below points=${checked} alsoBelow(B or fine)=${bothBelow} A-only(spurious)=${aOnlyBelow}`)
  if (first) {
    const { seg, x, y, end } = first
    phi1.set(phis[pts.findIndex((p) => p[0] === x && p[1] === y)])
    const prof = (name: string, mode: 'A' | 'B') => {
      let mn = Infinity; let tm = 0
      for (let k2 = 0; k2 <= 2000; k2++) {
        const t = seg.t0 + ((seg.t1 - seg.t0) * k2) / 2000
        const v = wttSeg(seg, t, mode) + G
        if (v < mn) { mn = v; tm = t }
      }
      console.log(`R1 first zone seg=[${seg.t0.toPrecision(17)}, ${seg.t1.toPrecision(17)}] end=${end} at (${x}, ${y}) ${name}: min(wtt+G)=${mn.toFixed(4)} at tau=${tm.toPrecision(17)}`)
    }
    prof('A hermite', 'A')
    prof('B ode', 'B')
    let maxEndDiff = 0
    const h = seg.t1 - seg.t0
    for (let j = 0; j < plate.ndof; j++) {
      const aH0 = hermite3(seg.q0[j], seg.qd0[j], seg.q1[j], seg.qd1[j], h, 0)[2]
      const aH1 = hermite3(seg.q0[j], seg.qd0[j], seg.q1[j], seg.qd1[j], h, 1)[2]
      maxEndDiff = Math.max(maxEndDiff,
        Math.abs(aH0 - odeQdd(j, seg.q0[j], seg.qd0[j], seg.t0)) * Math.abs(phi1[j]),
        Math.abs(aH1 - odeQdd(j, seg.q1[j], seg.qd1[j], seg.t1)) * Math.abs(phi1[j]))
    }
    console.log(`R1 endpoint |hermiteQdd - odeQdd|*|Phi| max over DOFs = ${maxEndDiff.toExponential(3)}`)
  }
}, 900000)

benchIt('M6-11R R2 500Hz detach counts across seeds', () => {
  const mk = (mode: 'hermite' | 'ode', spp: number, seed: number) => {
    const plate = new TransientPlate(singleToneScore(500, calibratedForceN(), 30), [0, 0])
    plate.accelMode = mode
    const sim = new TransientSandSim(plate, initialState(plate.a, 60, seed), { ...DEFAULT_PARAMS, spp }, { enabled: true })
    sim.advanceTo(0.1)
    return sim
  }
  for (const seed of [1, 7, 99, 4242]) {
    const a44 = mk('hermite', 44, seed); const a176 = mk('hermite', 176, seed); const b44 = mk('ode', 44, seed)
    console.log(`R2 seed=${seed}: detach A44=${a44.nDetach} A176=${a176.nDetach} B44=${b44.nDetach} | bounce A44=${a44.nBounce} A176=${a176.nBounce} B44=${b44.nBounce}`)
  }
}, 900000)

benchIt('M6-11R R1b tiny final step artifact (500Hz, pre/post fix mechanism)', () => {
  const mk = (mode: 'hermite' | 'ode', spp: number) => {
    const plate = new TransientPlate(singleToneScore(500, calibratedForceN(), 30), [0, 0])
    plate.accelMode = mode
    return new TransientSandSim(plate, initialState(plate.a, 60, 12345), { ...DEFAULT_PARAMS, spp }, { enabled: true })
  }
  // Pure nominal step() driving (production worker pattern): no
  // shortened final step ever occurs.
  for (const mode of ['hermite', 'ode'] as const) {
    const sim = mk(mode, 44)
    const dt = sim.dtNow()
    const nSteps = Math.round(0.1 / dt)
    let firstDetachStep = -1
    for (let k = 0; k < nSteps; k++) {
      sim.step()
      if (firstDetachStep < 0 && sim.nDetach > 0) firstDetachStep = k
    }
    console.log(`R1b step-loop mode=${mode}: steps=${nSteps} t=${sim.t.toPrecision(17)} bounce=${sim.nBounce} detach=${sim.nDetach} firstDetachStep=${firstDetachStep}`)
  }
  // advanceTo with targets that do / do not land on the step grid.
  for (const target of [0.1, 0.100001, 0.0999, 0.05]) {
    const sim = mk('hermite', 44)
    sim.advanceTo(target)
    console.log(`R1b advanceTo(${target}) hermite: bounce=${sim.nBounce} detach=${sim.nDetach} t=${sim.t.toPrecision(17)} lastDetachT=${sim.lastDetachT}`)
  }
  // When do detaches occur under advanceTo(0.1)? Step manually via
  // advanceTo in grid-aligned chunks and watch nDetach.
  {
    const sim = mk('hermite', 44)
    const dt = sim.dtNow()
    let prev = 0
    for (let k = 1; k <= 10; k++) {
      sim.advanceTo(Math.min(0.1, k * Math.round(0.01 / dt) * dt))
      console.log(`R1b aligned-chunk k=${k} t=${sim.t.toPrecision(6)} detach=${sim.nDetach} (+${sim.nDetach - prev})`)
      prev = sim.nDetach
    }
    sim.advanceTo(0.1)
    console.log(`R1b aligned then final advanceTo(0.1): detach=${sim.nDetach} (+${sim.nDetach - prev}) t=${sim.t.toPrecision(17)}`)
  }
}, 600000)

function scanOnlyDetach(plate: TransientPlate, x: number, y: number, ta: number, tb: number): number | null {
  const out = new Float64Array(5)
  const w = (t: number) => { plate.surfAt(x, y, t, out); return out[4] }
  if (w(ta) < -G) return ta
  let tp = ta
  for (let k = 1; k <= 8; k++) {
    const tk = ta + (tb - ta) * (k / 8)
    if (w(tk) < -G) {
      let lo = tp, hi = tk
      for (let i = 0; i < 40; i++) { const m = (lo + hi) / 2; if (w(m) >= -G) lo = m; else hi = m }
      return (lo + hi) / 2
    }
    tp = tk
  }
  return null
}

benchIt('M6-11R R3 adversarial detach-finder search (scan vs interval vs dense)', () => {
  const out = new Float64Array(5)
  const positions: [number, number][] = [[0.01, 0.005], [-0.04, 0.03], [0.07, -0.02], [0.002, 0.09]]
  for (const freq of [152, 649, 1241.6]) {
    const baseF = calibratedForceN()
    for (const fac of [0.02, 0.05, 0.1, 0.2, 0.4]) {
      const plate = new TransientPlate(singleToneScore(freq, baseF * fac, 30), [0, 0])
      plate.accelMode = 'ode'
      const dt = 1 / (44 * freq)
      // Warm up well into the driven transient.
      for (let k = 0; k < 600; k++) plate.beginStep(k * dt, (k + 1) * dt)
      let dense = 0, scanFound = 0, scanMiss = 0, intervalFound = 0, intervalMiss = 0, timeMismatch = 0
      let narrowest = Infinity
      for (let s = 600; s < 720; s++) {
        plate.beginStep(s * dt, (s + 1) * dt)
        const t0 = s * dt, t1 = (s + 1) * dt
        for (const [x, y] of positions) {
          let first: number | null = null, last: number | null = null
          const N = 4000
          for (let k = 0; k <= N; k++) {
            plate.surfAt(x, y, t0 + (t1 - t0) * (k / N), out)
            if (out[4] < -G) { if (first === null) first = t0 + (t1 - t0) * (k / N); last = t0 + (t1 - t0) * (k / N) }
          }
          if (first === null) continue
          dense++; narrowest = Math.min(narrowest, last! - first)
          const sOnly = scanOnlyDetach(plate, x, y, t0, t1)
          const iGot = plate.firstDetachBelow(x, y, t0, t1)
          if (sOnly !== null) scanFound++; else scanMiss++
          if (iGot !== null) intervalFound++; else intervalMiss++
          if (iGot !== null && Math.abs(iGot - first) > (t1 - t0) / N + 2e-6) timeMismatch++
        }
      }
      console.log(`R3 f=${freq} force×${fac}: denseRoots=${dense} scanFound=${scanFound} scanMiss=${scanMiss} intervalFound=${intervalFound} intervalMiss=${intervalMiss} timeMismatch=${timeMismatch} narrowestEpisode=${narrowest === Infinity ? 'none' : narrowest.toExponential(2)}s dt=${dt.toExponential(2)}s`)
    }
  }
}, 900000)
