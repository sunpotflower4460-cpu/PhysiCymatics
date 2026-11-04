// M6-12 long-time sand motion study (skipped unless PHYSI_BENCH=1).
// Pure measurement: nothing here feeds back into the dynamics. The
// "activity map" (time-averaged |wtt| of the bare plate) is computed
// from a separate plate-only run and only compared with the emergent
// sand density afterwards.
import { it } from 'vitest'
// @ts-expect-error bench-only node import (runs under vitest node env)
import { writeFileSync, mkdirSync } from 'node:fs'
// @ts-expect-error bench-only node import (runs under vitest node env)
import { deflateSync } from 'node:zlib'
import { singleToneScore } from '../sequence'
import { TransientSandSim, createTransientSim, calibratedForceN } from '../transientSand'
import { TransientPlate } from '../plateTransient'

const proc = (globalThis as { process?: { env: Record<string, string | undefined> } }).process
const benchIt = proc?.env.PHYSI_BENCH ? it : it.skip

const GRID = 32
const MAX_STEPS = 80_000
const SNAP_TIMES = [0.1, 0.5, 1, 2, 3, 5]

interface Snap {
  t: number
  alive: number; sittingFrac: number; sleepingFrac: number
  meanSpeed: number; entropy: number
  radial: number[]; density: Float64Array
  pathMean: number
}

function snapshot(sim: TransientSandSim, path: Float64Array): Snap {
  const st = sim.state; const a = sim.plate.a
  const radial = new Array(10).fill(0)
  const density = new Float64Array(GRID * GRID)
  let alive = 0; let sitting = 0; let sleeping = 0; let speed = 0; let pathSum = 0
  for (let i = 0; i < st.id.length; i++) {
    if (!st.alive[i]) continue
    alive++
    if (st.sitting[i]) sitting++
    if (sim.sleeping[i]) sleeping++
    speed += Math.hypot(st.vx[i], st.vy[i], st.vz[i])
    pathSum += path[i]
    const r = Math.hypot(st.x[i], st.y[i]) / a
    radial[Math.min(9, (r * 10) | 0)]++
    const gx = Math.min(GRID - 1, Math.max(0, ((st.x[i] / a + 1) * GRID / 2) | 0))
    const gy = Math.min(GRID - 1, Math.max(0, ((st.y[i] / a + 1) * GRID / 2) | 0))
    density[gy * GRID + gx]++
  }
  let entropy = 0
  if (alive > 0) {
    for (let c = 0; c < density.length; c++) {
      const gx = c % GRID; const gy = (c / GRID) | 0
      const cx = (gx + 0.5) * 2 / GRID - 1; const cy = (gy + 0.5) * 2 / GRID - 1
      if (cx * cx + cy * cy > 1) continue
      const p = density[c] / alive
      if (p > 0) entropy -= p * Math.log(p)
    }
  }
  return {
    t: sim.t, alive,
    sittingFrac: alive ? sitting / alive : 0,
    sleepingFrac: alive ? sleeping / alive : 0,
    meanSpeed: alive ? speed / alive : 0,
    entropy, radial, density,
    pathMean: alive ? pathSum / alive : 0,
  }
}

function runLong(freq: number, n: number, seed: number, mode: 'hermite' | 'ode'): { snaps: Snap[]; sim: TransientSandSim; wallS: number } {
  const sim = createTransientSim(singleToneScore(freq, calibratedForceN(), 30), n, seed)
  if (mode === 'ode') sim.plate.accelMode = 'ode'
  const st = sim.state
  const path = new Float64Array(n)
  const px = Float64Array.from(st.x); const py = Float64Array.from(st.y); const pz = Float64Array.from(st.z)
  const dt = sim.dtNow()
  const snapSteps = new Map<number, number>()
  for (const ts of SNAP_TIMES) {
    const k = Math.round(ts / dt)
    if (k <= MAX_STEPS) snapSteps.set(k, ts)
  }
  const snaps: Snap[] = []
  const wall0 = Date.now()
  for (let s = 1; s <= MAX_STEPS; s++) {
    sim.step()
    for (let i = 0; i < n; i++) {
      if (st.alive[i]) {
        path[i] += Math.hypot(st.x[i] - px[i], st.y[i] - py[i], st.z[i] - pz[i])
      }
      px[i] = st.x[i]; py[i] = st.y[i]; pz[i] = st.z[i]
    }
    if (snapSteps.has(s)) snaps.push(snapshot(sim, path))
    if (s * dt >= 5) break
  }
  if (snaps.length === 0 || snaps[snaps.length - 1].t !== sim.t) snaps.push(snapshot(sim, path))
  return { snaps, sim, wallS: (Date.now() - wall0) / 1000 }
}

/** Time-averaged |wtt| over the bare plate after transients (measurement only). */
function activityMap(freq: number, mode: 'hermite' | 'ode'): Float64Array {
  const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), [0, 0])
  if (mode === 'ode') plate.accelMode = 'ode'
  const dt = 1 / (44 * freq)
  const warm = Math.min(4000, Math.round(0.2 / dt))
  for (let k = 0; k < warm; k++) plate.beginStep(k * dt, (k + 1) * dt)
  const act = new Float64Array(GRID * GRID)
  const out = new Float64Array(5)
  const a = plate.a
  const samples = 44 // one nominal drive period, 1 sample per step
  for (let k = 0; k < samples; k++) {
    const t0 = (warm + k) * dt
    plate.beginStep(t0, t0 + dt)
    const tm = t0 + dt / 2
    for (let gy = 0; gy < GRID; gy++) {
      for (let gx = 0; gx < GRID; gx++) {
        const x = ((gx + 0.5) * 2 / GRID - 1) * a
        const y = ((gy + 0.5) * 2 / GRID - 1) * a
        plate.surfAt(x, y, tm, out)
        act[gy * GRID + gx] += Math.abs(out[4]) / samples
      }
    }
  }
  return act
}

function pearson(d: Float64Array, act: Float64Array): number {
  // Correlation over in-disc cells between sand density and activity.
  // Negative = sand where the plate barely accelerates.
  const xs: number[] = []; const ys: number[] = []
  for (let gy = 0; gy < GRID; gy++) {
    for (let gx = 0; gx < GRID; gx++) {
      const cx = (gx + 0.5) * 2 / GRID - 1; const cy = (gy + 0.5) * 2 / GRID - 1
      if (cx * cx + cy * cy > 1) continue
      xs.push(d[gy * GRID + gx]); ys.push(act[gy * GRID + gx])
    }
  }
  const n = xs.length
  const mx = xs.reduce((s, v) => s + v, 0) / n; const my = ys.reduce((s, v) => s + v, 0) / n
  let sxy = 0; let sxx = 0; let syy = 0
  for (let i = 0; i < n; i++) {
    sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0
}

// Minimal grayscale PNG writer (no smoothing; nearest-cell upscale).
const CRC_TABLE = (() => {
  const t = new Int32Array(256)
  for (let n = 0; n < 256; n++) {
    let c = n
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1
    t[n] = c
  }
  return t
})()
function chunk(type: string, data: Uint8Array): Uint8Array {
  const out = new Uint8Array(data.length + 12)
  const dv = new DataView(out.buffer)
  dv.setUint32(0, data.length)
  for (let i = 0; i < 4; i++) out[4 + i] = type.charCodeAt(i)
  out.set(data, 8)
  let crc = 0xffffffff
  for (let i = 4; i < 8 + data.length; i++) crc = CRC_TABLE[(crc ^ out[i]) & 0xff] ^ (crc >>> 8)
  dv.setUint32(8 + data.length, (crc ^ 0xffffffff) >>> 0)
  return out
}
function writePng(path: string, grid: Float64Array, invert: boolean): void {
  const S = 256; const scale = S / GRID
  let mx = 0
  for (const v of grid) mx = Math.max(mx, v)
  const raw = new Uint8Array(S * (S + 1))
  for (let y = 0; y < S; y++) {
    raw[y * (S + 1)] = 0
    const gy = GRID - 1 - ((y / scale) | 0)
    for (let x = 0; x < S; x++) {
      const gx = (x / scale) | 0
      const cx = (gx + 0.5) * 2 / GRID - 1; const cy = (gy + 0.5) * 2 / GRID - 1
      let v = cx * cx + cy * cy > 1 ? 32 : 255 * (mx > 0 ? grid[gy * GRID + gx] / mx : 0)
      if (invert && cx * cx + cy * cy <= 1) v = 255 - v
      raw[y * (S + 1) + 1 + x] = v
    }
  }
  const ihdr = new Uint8Array(13)
  const dv = new DataView(ihdr.buffer)
  dv.setUint32(0, S); dv.setUint32(4, S)
  ihdr[8] = 8; ihdr[9] = 0
  const parts = [
    new Uint8Array([137, 80, 78, 71, 13, 10, 26, 10]),
    chunk('IHDR', ihdr),
    chunk('IDAT', deflateSync(raw)),
    chunk('IEND', new Uint8Array(0)),
  ]
  const total = parts.reduce((s, p) => s + p.length, 0)
  const out = new Uint8Array(total)
  let off = 0
  for (const p of parts) { out.set(p, off); off += p.length }
  writeFileSync(path, out)
}

function fmtSnap(s: Snap): string {
  return `t=${s.t.toFixed(3)}s alive=${s.alive} sit=${(s.sittingFrac * 100).toFixed(0)}% sleep=${(s.sleepingFrac * 100).toFixed(0)}% vMean=${s.meanSpeed.toExponential(2)} H=${s.entropy.toFixed(3)} path=${s.pathMean.toFixed(3)}m radial=[${s.radial.join(',')}]`
}

benchIt('M6-12 long-time sand motion: 300 grains, A vs B, 8 frequencies', () => {
  const freqs = [70, 152, 173.31, 271.3, 368.6, 500, 649, 1241.6]
  mkdirSync(new URL('../../../../docs/assets/m612', import.meta.url), { recursive: true })
  const assetDir = new URL('../../../../docs/assets/m612/', import.meta.url)
  for (const freq of freqs) {
    const act = activityMap(freq, 'hermite')
    writePng(new URL(`activity_${freq}hz.png`, assetDir).pathname, act, false)
    for (const mode of ['hermite', 'ode'] as const) {
      const { snaps, sim, wallS } = runLong(freq, 300, 42, mode)
      const last = snaps[snaps.length - 1]
      const prev = snaps.length > 1 ? snaps[snaps.length - 2] : null
      let l1 = 0
      if (prev && last.alive > 0) {
        for (let c = 0; c < last.density.length; c++) {
          l1 += Math.abs(last.density[c] / last.alive - prev.density[c] / Math.max(1, prev.alive))
        }
        l1 /= 2
      }
      const corr = pearson(last.density, act)
      console.log(`M612 f=${freq} mode=${mode} wall=${wallS.toFixed(1)}s bounce=${sim.nBounce} detach=${sim.nDetach} cap=${sim.nEventCap} corr(density,activity)=${corr.toFixed(3)} lastL1=${l1.toFixed(3)}`)
      for (const s of snaps) console.log(`  ${fmtSnap(s)}`)
      if (mode === 'hermite') {
        writePng(new URL(`density_${freq}hz.png`, assetDir).pathname, last.density, true)
      }
    }
  }
}, 3_600_000)

benchIt('M6-12 plate-only detach margin per frequency (explains frozen vs active)', () => {
  // Center drive excites only m=0 modes, so the wtt field is radial.
  // min over r and t of (g + wtt): if it never goes below 0, a seated
  // grain can never detach under this model (and with c_roll=0 it
  // then never moves horizontally again).
  const out = new Float64Array(5)
  for (const freq of [70, 152, 173.31, 271.3, 368.6, 500, 649, 1241.6]) {
    for (const mode of ['hermite', 'ode'] as const) {
      const plate = new TransientPlate(singleToneScore(freq, calibratedForceN(), 30), [0, 0])
      if (mode === 'ode') plate.accelMode = 'ode'
      const dt = 1 / (44 * freq)
      const nSeg = Math.min(Math.round(1 / dt), 60_000)
      let minMargin = Infinity; let minR = 0
      for (let k = 0; k < nSeg; k++) {
        plate.beginStep(k * dt, (k + 1) * dt)
        for (const tau of [k * dt, (k + 0.5) * dt, (k + 1) * dt]) {
          for (let ri = 0; ri <= 200; ri++) {
            const r = (ri / 200) * plate.a
            plate.surfAt(r, 0, tau, out)
            const m = 9.81 + out[4]
            if (m < minMargin) { minMargin = m; minR = r }
          }
        }
      }
      console.log(`M612-margin f=${freq} mode=${mode} min(g+wtt)=${minMargin.toFixed(3)} m/s2 at r/a=${(minR / plate.a).toFixed(2)} segs=${nSeg} -> ${minMargin < 0 ? 'detach possible' : 'detach IMPOSSIBLE everywhere'}`)
    }
  }
}, 1_800_000)

benchIt('M6-12 scale tendency: 1000 and 5000 grains, shorter time', () => {
  for (const [freq, n] of [[152, 1000], [152, 5000], [649, 1000]] as const) {
    const act = activityMap(freq, 'hermite')
    for (const mode of ['hermite', 'ode'] as const) {
      const sim = createTransientSim(singleToneScore(freq, calibratedForceN(), 30), n, 42)
      if (mode === 'ode') sim.plate.accelMode = 'ode'
      const dt = sim.dtNow()
      const steps = Math.min(Math.round(1 / dt), MAX_STEPS)
      const wall0 = Date.now()
      for (let s = 0; s < steps; s++) sim.step()
      const snap = snapshot(sim, new Float64Array(n))
      console.log(`M612-scale f=${freq} n=${n} mode=${mode} wall=${((Date.now() - wall0) / 1000).toFixed(1)}s corr=${pearson(snap.density, act).toFixed(3)} ${fmtSnap(snap)}`)
    }
  }
}, 3_600_000)
