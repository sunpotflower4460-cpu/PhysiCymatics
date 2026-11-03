// M6-11 sand-level study (skipped unless PHYSI_BENCH=1):
//  S1: bounce/detach counts, hermite (A) vs ode (B), many frequencies
//  S2: divergence onset A vs B (+ chaos control) and final radial
//      distribution distance at 649Hz; eccentric drive; freq switch
//  S3: advanceTo chunking effect (single vs split calls)
//  S4: performance A vs B at 5k/20k grains; scale timing at 100k/600k
// Counts over long horizons are NOT convergent observables (M6-9R D1);
// these numbers measure sensitivity, not correctness by themselves.
import { it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientSandSim, createTransientSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'
import { initialState } from '../sand'
import { DEFAULT_PARAMS } from '../types'

const proc = (globalThis as { process?: { env: Record<string, string | undefined>; memoryUsage: () => { rss: number } } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

function mk(freq: number, mode: 'hermite' | 'ode', n = 60, driveXY: [number, number] = [0, 0], dz = 0, sleep = true): TransientSandSim {
  const sim = createTransientSim(singleToneScore(freq, calibratedForceN(), 30), n, 12345, driveXY, { enabled: sleep })
  sim.plate.accelMode = mode
  if (dz) for (let i = 0; i < sim.state.id.length; i++) sim.state.z[i] += dz
  return sim
}
function maxPosDiff(a: TransientSandSim, b: TransientSandSim): number {
  let m = 0
  for (let i = 0; i < a.state.id.length; i++) {
    m = Math.max(m, Math.abs(a.state.x[i] - b.state.x[i]), Math.abs(a.state.y[i] - b.state.y[i]), Math.abs(a.state.z[i] - b.state.z[i]))
  }
  return m
}
function radialHist(sim: TransientSandSim, bins = 10): number[] {
  const h = new Array(bins).fill(0)
  let alive = 0
  for (let i = 0; i < sim.state.id.length; i++) {
    if (!sim.state.alive[i]) continue
    alive++
    const r = Math.hypot(sim.state.x[i], sim.state.y[i]) / sim.plate.a
    h[Math.min(bins - 1, Math.floor(r * bins))]++
  }
  return h.map((v) => (alive ? v / alive : 0))
}
function histL1(a: number[], b: number[]): number {
  let s = 0
  for (let i = 0; i < a.length; i++) s += Math.abs(a[i] - b[i])
  return s
}

benchIt('M6-11 S1 counts A vs B across frequencies', () => {
  for (const freq of [70, 152, 173.31, 368.6, 500, 649, 1241.6]) {
    const a = mk(freq, 'hermite'); const b = mk(freq, 'ode')
    a.advanceTo(0.1); b.advanceTo(0.1)
    console.log(`S1 f=${freq}: bounce A=${a.nBounce} B=${b.nBounce} detach A=${a.nDetach} B=${b.nDetach} alive A=${a.nAlive} B=${b.nAlive} caps A=${a.nEventCap} B=${b.nEventCap} maxPosDiff=${maxPosDiff(a, b).toExponential(2)} lastDetachT A=${a.lastDetachT} B=${b.lastDetachT}`)
  }
  // Eccentric drive + frequency switch with residual vibration.
  for (const mode of ['hermite', 'ode'] as const) {
    const s = mk(649, mode, 60, [0.05, 0.03])
    s.advanceTo(0.1)
    console.log(`S1 ecc649 mode=${mode}: bounce=${s.nBounce} detach=${s.nDetach} alive=${s.nAlive}`)
  }
  const sw: Record<string, TransientSandSim> = {}
  for (const mode of ['hermite', 'ode'] as const) {
    const s = mk(152, mode)
    s.advanceTo(0.05)
    s.wakeAll(); s.plate.setFrequency(649, s.t); s.refreshForceBounds()
    s.advanceTo(0.1)
    sw[mode] = s
    console.log(`S1 switch152to649 mode=${mode}: bounce=${s.nBounce} detach=${s.nDetach} alive=${s.nAlive}`)
  }
  console.log(`S1 switch maxPosDiff A-vs-B=${maxPosDiff(sw.hermite, sw.ode).toExponential(2)}`)
}, 900000)

benchIt('M6-11 S2 divergence onset and distribution at 649Hz', () => {
  const pairs: [string, TransientSandSim, TransientSandSim][] = [
    ['AvsB', mk(649, 'hermite'), mk(649, 'ode')],
    ['AvsChaos', mk(649, 'hermite'), mk(649, 'hermite', 60, [0, 0], 6e-7)],
  ]
  const dT = 6e-3; const T = 0.12
  const firstState: Record<string, number> = {}
  for (let k = 1; k <= Math.round(T / dT); k++) {
    const t = k * dT
    for (const [name, a, b] of pairs) {
      a.advanceTo(t); b.advanceTo(t)
      let stateDiv = false
      for (let i = 0; i < a.state.id.length; i++) {
        if (a.state.sitting[i] !== b.state.sitting[i] || a.state.alive[i] !== b.state.alive[i]) stateDiv = true
      }
      if (!(name in firstState) && stateDiv) firstState[name] = t
      if (k % 5 === 0) console.log(`S2 ${name} t=${t.toFixed(3)} maxPosDiff=${maxPosDiff(a, b).toExponential(2)} bounce=${a.nBounce}/${b.nBounce} detach=${a.nDetach}/${b.nDetach}`)
    }
  }
  for (const [name, a, b] of pairs) {
    console.log(`S2 ${name}: firstStateDiv<=${firstState[name] ?? -1}s radialHistL1=${histL1(radialHist(a), radialHist(b)).toFixed(3)} final bounce=${a.nBounce}/${b.nBounce} detach=${a.nDetach}/${b.nDetach}`)
  }
}, 900000)

benchIt('M6-11 S3 advanceTo chunking effect', () => {
  const single = mk(649, 'hermite')
  let t0 = performance.now()
  single.advanceTo(0.06)
  const tSingle = performance.now() - t0
  const chunked = mk(649, 'hermite')
  t0 = performance.now()
  for (let k = 1; k <= 10; k++) chunked.advanceTo(k * 6e-3)
  const tChunk = performance.now() - t0
  const fine = mk(649, 'hermite')
  t0 = performance.now()
  for (let k = 1; k <= 60; k++) fine.advanceTo(k * 1e-3)
  const tFine = performance.now() - t0
  console.log(`S3 single: bounce=${single.nBounce} detach=${single.nDetach} wall=${tSingle.toFixed(0)}ms`)
  console.log(`S3 6ms-chunks: bounce=${chunked.nBounce} detach=${chunked.nDetach} wall=${tChunk.toFixed(0)}ms maxPosDiffVsSingle=${maxPosDiff(single, chunked).toExponential(2)}`)
  console.log(`S3 1ms-chunks: bounce=${fine.nBounce} detach=${fine.nDetach} wall=${tFine.toFixed(0)}ms maxPosDiffVsSingle=${maxPosDiff(single, fine).toExponential(2)}`)
  // Aligned chunks (multiples of nominal dt) should reproduce single call.
  const aligned = mk(649, 'hermite')
  const dtN = aligned.dtNow()
  const stepsPerChunk = Math.round(6e-3 / dtN)
  for (let k = 1; k <= Math.floor(0.06 / (stepsPerChunk * dtN)); k++) aligned.advanceTo(k * stepsPerChunk * dtN)
  aligned.advanceTo(0.06)
  console.log(`S3 dt-aligned chunks (dt=${dtN.toExponential(3)}): bounce=${aligned.nBounce} detach=${aligned.nDetach} maxPosDiffVsSingle=${maxPosDiff(single, aligned).toExponential(2)}`)
}, 900000)

benchIt('M6-11 S4 performance A vs B and scale', () => {
  for (const [n, dur] of [[5000, 0.05], [20000, 0.05]] as [number, number][]) {
    for (const mode of ['hermite', 'ode'] as const) {
      const sim = mk(649, mode, n)
      const t0 = performance.now()
      sim.advanceTo(dur)
      const ms = performance.now() - t0
      console.log(`S4 n=${n} mode=${mode} dur=${dur}s: wall=${(ms / 1000).toFixed(2)}s ratioToRealtime=${(ms / 1000 / dur).toFixed(1)} bounce=${sim.nBounce} detach=${sim.nDetach} surfCalls=${sim.plate.surfCalls} heightCalls=${sim.plate.heightCalls} heightAccelCalls=${sim.plate.heightAccelCalls}`)
    }
  }
  for (const [n, dur] of [[100000, 0.02], [600000, 0.02]] as [number, number][]) {
    const rss0 = (proc?.memoryUsage().rss ?? 0) / 1048576
    const sim = mk(649, 'hermite', n)
    const rss1 = (proc?.memoryUsage().rss ?? 0) / 1048576
    const t0 = performance.now()
    sim.advanceTo(dur)
    const ms = performance.now() - t0
    console.log(`S4 scale n=${n} mode=hermite dur=${dur}s: wall=${(ms / 1000).toFixed(2)}s ratioToRealtime=${(ms / 1000 / dur).toFixed(0)} bounce=${sim.nBounce} rssMB~${rss0.toFixed(0)}->${rss1.toFixed(0)} (single-thread, single run)`)
  }
}, 1800000)

benchIt('M6-11 S1b f=500 spurious detach diagnosis', () => {
  const mk = (mode: 'hermite' | 'ode', spp: number, sleep: boolean) => {
    const plate = new TransientPlate(singleToneScore(500, calibratedForceN(), 30), [0, 0])
    plate.accelMode = mode
    const sim = new TransientSandSim(plate, initialState(plate.a, 60, 12345), { ...DEFAULT_PARAMS, spp }, { enabled: sleep })
    sim.advanceTo(0.1)
    return sim
  }
  for (const mode of ['hermite', 'ode'] as const) {
    for (const spp of [44, 88, 176]) {
      for (const sleep of [true, false]) {
        const s = mk(mode, spp, sleep)
        console.log(`S1b mode=${mode} spp=${spp} sleep=${sleep}: bounce=${s.nBounce} detach=${s.nDetach} sleeping=${s.nSleeping}`)
      }
    }
  }
  // What does implemented wtt look like vs analytic at a seated point?
  // Min wtt over time at plate center region under each mode (plate only).
  for (const mode of ['hermite', 'ode'] as const) {
    const plate = new TransientPlate(singleToneScore(500, calibratedForceN(), 30), [0, 0])
    plate.accelMode = mode
    const dt = 1 / (44 * 500)
    const out = new Float64Array(5)
    let minW = Infinity
    for (let k = 0; k < Math.round(0.1 / dt); k++) {
      plate.beginStep(k * dt, (k + 1) * dt)
      for (let sN = 0; sN <= 8; sN++) {
        plate.surfAt(0.02, 0.01, k * dt + (dt * sN) / 8, out)
        minW = Math.min(minW, out[4])
      }
    }
    console.log(`S1b plate-only mode=${mode}: min wtt at (0.02,0.01) over 0.1s = ${minW.toFixed(3)} (detach threshold -9.81)`)
  }
}, 600000)
