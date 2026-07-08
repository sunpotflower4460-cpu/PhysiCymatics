// faraday.ts — TIER2: Faraday onset threshold + above-threshold pattern.
//
// PUT-IN     : pack threshold curves (a_c_ms2, band, λ, subharmonic) per depth,
//              measured by the c003 oracle (two-route 0.0013%). Log-f interp.
// EMERGED    : (i) the onset — below threshold the surface is SILENT; (ii) above
//              threshold, a subharmonic (f/2) standing-wave field on the unstable
//              band |k|=2π/λ. Isotropic random-phase seeding — no angle is put in.
// CLAIM-TIER : measured (threshold, λ, f/2). Amplitude is NOT set by linear
//              theory → the field is amplitude-normalized and labelled as such.
// FLOORS     : linear onset only; pattern TYPE (square/hex) is nonlinear (see
//              shape.ts / c005). Above ~5× threshold → out of linear scope.
import { thresholdCurve, type ThresholdPoint } from './pack'

export interface Threshold {
  a_c_ms2: number
  band_ms2: [number, number]
  lam_mm: number
  subharmonic: boolean
}

/** Log-frequency linear interpolation of the measured threshold curve. */
export function thresholdAt(freqHz: number, depth_m: number): Threshold {
  const c = thresholdCurve(depth_m)
  const f = Math.max(c[0].f_hz, Math.min(c[c.length - 1].f_hz, freqHz))
  let lo: ThresholdPoint = c[0]
  let hi: ThresholdPoint = c[c.length - 1]
  for (let i = 0; i < c.length - 1; i++) {
    if (f >= c[i].f_hz && f <= c[i + 1].f_hz) { lo = c[i]; hi = c[i + 1]; break }
  }
  const t = lo.f_hz === hi.f_hz ? 0
    : (Math.log(f) - Math.log(lo.f_hz)) / (Math.log(hi.f_hz) - Math.log(lo.f_hz))
  const lerp = (a: number, b: number) => a + (b - a) * t
  return {
    a_c_ms2: lerp(lo.a_c_ms2, hi.a_c_ms2),
    band_ms2: [lerp(lo.band_ms2[0], hi.band_ms2[0]), lerp(lo.band_ms2[1], hi.band_ms2[1])],
    lam_mm: lerp(lo.lam_mm, hi.lam_mm),
    subharmonic: lo.subharmonic && hi.subharmonic,
  }
}

/** Deterministic phases (mulberry32) so the pattern is reproducible per config. */
function mulberry32(seed: number) {
  let a = seed >>> 0
  return () => {
    a |= 0; a = (a + 0x6d2b79f5) | 0
    let t = Math.imul(a ^ (a >>> 15), 1 | a)
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

export interface FaradayView {
  active: boolean        // drive ≥ threshold ?
  ratio: number          // drive / threshold
  overLinear: boolean    // > 5× threshold (nonlinear, out of scope)
  lam_mm: number
  respHz: number         // f/2
}

export function faradayState(freqHz: number, depth_m: number, accel_ms2: number): FaradayView {
  const th = thresholdAt(freqHz, depth_m)
  const ratio = accel_ms2 / th.a_c_ms2
  return {
    active: accel_ms2 >= th.a_c_ms2,
    ratio,
    overLinear: ratio > 5,
    lam_mm: th.lam_mm,
    respHz: freqHz / 2,
  }
}

const NDIR = 12 // isotropic directions sampled on the unstable ring |k|=k_c

/**
 * Subharmonic standing-wave field on the unstable band, at time tSec.
 * windowMm = physical width the canvas spans. Field is amplitude-NORMALIZED
 * (linear theory does not fix the amplitude). Below threshold → all zeros.
 */
export function faradayField(
  freqHz: number, depth_m: number, accel_ms2: number, tSec: number,
  windowMm = 48, w = 160, h = 160, seed = 12345,
): { data: Float32Array; w: number; h: number; active: boolean } {
  const view = faradayState(freqHz, depth_m, accel_ms2)
  const data = new Float32Array(w * h)
  if (!view.active) return { data, w, h, active: false } // honest silence

  const kc = (2 * Math.PI) / view.lam_mm // 1/mm
  const rng = mulberry32(seed)
  const dirs = Array.from({ length: NDIR }, (_, j) => {
    const ang = (Math.PI * j) / NDIR // 0..π (k and −k identical for standing wave)
    return { kx: kc * Math.cos(ang), ky: kc * Math.sin(ang), phi: rng() * 2 * Math.PI }
  })
  // subharmonic temporal factor: surface oscillates at f/2
  const temporal = Math.cos(2 * Math.PI * (freqHz / 2) * tSec)

  let maxAbs = 0
  const raw = new Float32Array(w * h)
  for (let py = 0; py < h; py++) {
    for (let px = 0; px < w; px++) {
      const xmm = (px / (w - 1) - 0.5) * windowMm
      const ymm = (py / (h - 1) - 0.5) * windowMm
      let v = 0
      for (const d of dirs) v += Math.cos(d.kx * xmm + d.ky * ymm + d.phi)
      const val = v * temporal
      raw[py * w + px] = val
      const a = Math.abs(val)
      if (a > maxAbs) maxAbs = a
    }
  }
  for (let i = 0; i < raw.length; i++) data[i] = maxAbs > 0 ? raw[i] / maxAbs : 0
  return { data, w, h, active: true }
}
