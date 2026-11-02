// M6-6 reproducible bench. Skipped unless PHYSI_BENCH=1:
//   PHYSI_BENCH=1 npx vitest run src/simulation/__tests__/bench_m66.test.ts
// Results are appended to webapp/bench_m66_out.txt (not committed).
import { it } from 'vitest'
// @ts-expect-error bench-only node import (runs under vitest node env)
import { appendFileSync } from 'node:fs'
import { singleToneScore } from '../sequence'
import { createTransientSim, calibratedForceN, type SleepOptions } from '../transientSand'

const proc = (globalThis as { process?: { env: Record<string, string | undefined>; memoryUsage: () => { rss: number } } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip
const OUT = new URL('../../../bench_m66_out.txt', import.meta.url).pathname
const log = (s: string) => { console.log(s); appendFileSync(OUT, s + '\n') }

type Mode = { name: string; opts: SleepOptions }
const MODES: Mode[] = [
  { name: 'off', opts: { enabled: true, localBound: false } },
  { name: 'k1', opts: { enabled: true, localBound: true } },
  { name: 'k4a', opts: { enabled: true, localBound: true, localSplits: 4 } },
  { name: 'k8a', opts: { enabled: true, localBound: true, localSplits: 8 } },
]

benchIt('M6-6 sliced proof bench', () => {
  const F = calibratedForceN()
  const run = (opts: SleepOptions, n: number, force: number, dur: number, pre = 0) => {
    const sim = createTransientSim(singleToneScore(152, force, 30), n, 7, [0, 0], opts)
    if (pre > 0) sim.advanceTo(pre)
    const t0 = performance.now()
    sim.advanceTo(pre + dur)
    const ms = performance.now() - t0
    return {
      ms, fly: sim.flySegments, sampled: sim.sampledSegments,
      local: sim.localSkips, sliceTry: sim.localSliceAttempts, sliceSkip: sim.localSliceSkips,
      pieces: sim.localPieceChecks, boundComp: sim.cellBoundComputations,
      bounce: sim.nBounce, detach: sim.nDetach, sleeping: sim.nSleeping,
      rss: Math.round((proc?.memoryUsage().rss ?? 0) / 1048576),
    }
  }
  run({ enabled: true }, 300, F, 0.2) // warmup
  const cases: [string, number, number, number, number, number][] = [
    // label, n, force, dur, pre, repeats
    ['500-strong', 500, F, 0.5, 0, 3],
    ['500-weak', 500, F * 0.05, 0.5, 0, 3],
    ['5000-strong', 5000, F, 0.5, 0, 3],
    ['5000-weak', 5000, F * 0.05, 0.5, 0, 3],
    ['5000-weak-formed', 5000, F * 0.05, 0.5, 2, 3],
    ['20000-strong', 20000, F, 0.2, 0, 3],
    ['20000-weak', 20000, F * 0.05, 0.2, 0, 3],
    ['100000-weak', 100000, F * 0.05, 0.05, 0, 1],
    ['600000-weak', 600000, F * 0.05, 0.02, 0, 1],
  ]
  for (const [label, n, force, dur, pre, reps] of cases) {
    const results = new Map<string, { times: number[]; last: ReturnType<typeof run> }>()
    for (let r = 0; r < reps; r++) {
      const modes = r % 2 === 0 ? MODES : [...MODES].reverse()
      for (const m of modes) {
        const res = run(m.opts, n, force, dur, pre)
        const e = results.get(m.name) ?? { times: [], last: res }
        e.times.push(res.ms); e.last = res
        results.set(m.name, e)
      }
    }
    for (const m of MODES) {
      const e = results.get(m.name)!
      const sorted = [...e.times].sort((a, b) => a - b)
      const med = sorted[Math.floor(sorted.length / 2)]
      const L = e.last
      log(`${label} mode=${m.name}: median=${Math.round(med)}ms min=${Math.round(sorted[0])} max=${Math.round(sorted[sorted.length - 1])} sampled=${L.sampled} localSkip=${L.local} sliceTry=${L.sliceTry} sliceSkip=${L.sliceSkip} pieces=${L.pieces} boundComp=${L.boundComp} bounce=${L.bounce} detach=${L.detach} sleeping=${L.sleeping} rssMB=${L.rss}`)
    }
  }
}, 3600000)
