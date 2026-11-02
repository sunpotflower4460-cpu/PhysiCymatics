// M6-7 reproducible profile/bench. Skipped unless PHYSI_BENCH=1:
//   PHYSI_BENCH=1 npx vitest run src/simulation/__tests__/bench_m67.test.ts
// Results are appended to webapp/bench_m67_out.txt (not committed).
// PHYSI_BENCH_CASES=initial|formed|all selects the case group.
import { it } from 'vitest'
// @ts-expect-error bench-only node import (runs under vitest node env)
import { appendFileSync } from 'node:fs'
import { singleToneScore } from '../sequence'
import { createTransientSim, calibratedForceN, type SleepOptions } from '../transientSand'
import { TransientPlate } from '../plateTransient'

const proc = (globalThis as { process?: { env: Record<string, string | undefined>; memoryUsage: () => { rss: number } } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip
const OUT = new URL('../../../bench_m67_out.txt', import.meta.url).pathname
const log = (s: string) => { console.log(s); appendFileSync(OUT, s + '\n') }

benchIt('M6-7 profile + A/B', () => {
  const F = calibratedForceN()
  const GROUP = proc?.env.PHYSI_BENCH_CASES ?? 'all'

  // Unit costs: surfAt vs heightAt per call in one live segment.
  {
    const plate = new TransientPlate(singleToneScore(649, F, 30))
    let seed = 987654
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
    plate.beginStep(0, 1.7e-4)
    const N = 200000
    const pts = new Float64Array(N * 3)
    for (let i = 0; i < N; i++) { pts[i * 3] = rnd() * 0.12; pts[i * 3 + 1] = rnd() * 0.12; pts[i * 3 + 2] = Math.abs(rnd()) * 1.7e-4 }
    const out = new Float64Array(5)
    let acc = 0
    let t0 = performance.now()
    for (let i = 0; i < N; i++) { plate.surfAt(pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2], out); acc += out[0] }
    const surfNs = ((performance.now() - t0) * 1e6) / N
    t0 = performance.now()
    for (let i = 0; i < N; i++) acc += plate.heightAt(pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2])
    const hNs = ((performance.now() - t0) * 1e6) / N
    t0 = performance.now()
    for (let i = 0; i < N; i++) acc += plate.heightAtFused(pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2])
    const fNs = ((performance.now() - t0) * 1e6) / N
    const out2 = new Float64Array(2)
    t0 = performance.now()
    for (let i = 0; i < N; i++) { plate.heightAccelAt(pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2], out2); acc += out2[0] + out2[1] }
    const haNs = ((performance.now() - t0) * 1e6) / N
    // Same-tau repeats (cache-hit path, like contact sample taus).
    t0 = performance.now()
    for (let i = 0; i < N; i++) acc += plate.heightAt(pts[i * 3], pts[i * 3 + 1], 8.5e-5)
    const hHitNs = ((performance.now() - t0) * 1e6) / N
    log(`unit: surfAt=${surfNs.toFixed(0)}ns heightAt=${hNs.toFixed(0)}ns heightAtFused=${fNs.toFixed(0)}ns heightAccelAt=${haNs.toFixed(0)}ns heightAtCacheHit=${hHitNs.toFixed(0)}ns acc=${acc.toFixed(3)}`)
  }

  type Mode = { name: string; opts: SleepOptions }
  const MODES: Mode[] = (proc?.env.PHYSI_BENCH_MODES ?? 'base')
    .split(',')
    .map((m): Mode => {
      switch (m) {
        case 'off': return { name: 'off', opts: { enabled: true, localBound: false } }
        case 'k1': return { name: 'k1', opts: { enabled: true, localBound: true } }
        case 'nosleep': return { name: 'nosleep', opts: { enabled: false } }
        case 'locsleep': return { name: 'locsleep', opts: { enabled: true, localSleep: true } }
        case 'base': return { name: 'base', opts: { enabled: true } }
        default: return { name: m, opts: { enabled: true, ...(JSON.parse(m) as object) } }
      }
    })

  const snap = (sim: ReturnType<typeof createTransientSim>) => ({
    surf: sim.plate.surfCalls, height: sim.plate.heightCalls,
    fly: sim.flySegments, glob: sim.globalSkips, loc: sim.localSkips,
    sampled: sim.sampledSegments, bisect: sim.bisectIters, sittingSteps: sim.sittingSteps,
    bounce: sim.nBounce, detach: sim.nDetach, sleepEval: sim.sleepEvaluations,
    wakes: sim.wakeCount, skipped: sim.skippedUpdates, bound: sim.cellBoundComputations,
    regionChecks: sim.regionChecks,
  })
  const run = (opts: SleepOptions, n: number, force: number, dur: number, pre = 0) => {
    const sim = createTransientSim(singleToneScore(152, force, 30), n, 7, [0, 0], opts)
    if (pre > 0) sim.advanceTo(pre)
    const c0 = snap(sim)
    const t0 = performance.now()
    sim.advanceTo(pre + dur)
    const ms = performance.now() - t0
    const c1 = snap(sim)
    const d = Object.fromEntries(Object.keys(c0).map((k) => [k, (c1 as Record<string, number>)[k] - (c0 as Record<string, number>)[k]]))
    return { ms, ...d, bounceTotal: sim.nBounce, detachTotal: sim.nDetach, sleepingFrac: sim.stats().sleepingFrac, rss: Math.round((proc?.memoryUsage().rss ?? 0) / 1048576) }
  }
  run({ enabled: true }, 300, F, 0.2) // warmup

  const cases: [string, string, number, number, number, number, number][] = [
    // group, label, n, force, dur, pre, repeats
    ['initial', '600000-weak-init', 600000, F * 0.05, 0.02, 0, 1],
    ['initial', '100000-weak-init', 100000, F * 0.05, 0.05, 0, 1],
    ['initial', '20000-strong-init', 20000, F, 0.2, 0, 3],
    ['initial', '5000-strong', 5000, F, 0.5, 0, 3],
    ['formed', '5000-weak-formed', 5000, F * 0.05, 0.5, 2, 3],
    ['formed', '20000-weak-formed', 20000, F * 0.05, 0.2, 2, 1],
    ['formed', '500-weak-formed', 500, F * 0.05, 0.5, 2, 3],
  ]
  for (const [group, label, n, force, dur, pre, reps] of cases) {
    if (GROUP !== 'all' && GROUP !== group) continue
    for (const m of MODES) {
      const times: number[] = []
      let last: ReturnType<typeof run> | null = null
      for (let r = 0; r < reps; r++) { last = run(m.opts, n, force, dur, pre); times.push(last.ms) }
      const med = [...times].sort((a, b) => a - b)[Math.floor(times.length / 2)]
      const L = last! as unknown as Record<string, number>
      log(`${label} mode=${m.name}: median=${Math.round(med)}ms surf=${L.surf} height=${L.height} flySeg=${L.fly} globSkip=${L.glob} locSkip=${L.loc} sampled=${L.sampled} bisect=${L.bisect} sittingSteps=${L.sittingSteps} bounce=${L.bounce} detach=${L.detach} sleepEval=${L.sleepEval} wakes=${L.wakes} skippedUpd=${L.skipped} boundComp=${L.bound} regionChecks=${L.regionChecks} bounceTot=${L.bounceTotal} sleepingFrac=${Number(L.sleepingFrac).toFixed(3)} rssMB=${L.rss}`)
    }
  }
}, 3600000)
