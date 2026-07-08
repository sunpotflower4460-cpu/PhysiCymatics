// plate.ts — TIER1: circular free-edge plate forced modal response.
//
// PUT-IN     : pack modal table (m,n,f_hz,R(r/a)) + the pack forced-response
//              formula; drive point (r0,θ0), damping ζ (pack PLACEHOLDER 0.002).
// EMERGED    : the steady-state displacement field; its ZERO SET is where sand
//              collects (the Chladni node lines). The computation chooses the
//              shape — no pattern is drawn by hand, no mode is special-cased
//              (center drive kills m≥1 automatically because R_m≥1(0)=0).
// CLAIM-TIER : measured (c001 known-physics model; two-route eigenvalues 1.6e-10,
//              Leissa-anchored 0.34%). ζ is a PLACEHOLDER floor.
// Rendering-honesty: returns the physical field only; the canvas reads its zero
//              set. No easing, no exaggeration.
import { plateModes, plateZeta, eigenfrequenciesHz, type PlateMode } from './pack'

const modes = plateModes()

/** Detent frequencies (30), ascending. */
export const eigenfrequencies = (): number[] => eigenfrequenciesHz()

/** ζ placeholder value + status from the pack. */
export const zeta = () => plateZeta()

/** Linear interpolation of a 0..1 (r/a) radial table. */
function interpR(R: number[], rOverA: number): number {
  if (rOverA <= 0) return R[0]
  if (rOverA >= 1) return R[R.length - 1]
  const x = rOverA * (R.length - 1)
  const i = Math.min(Math.floor(x), R.length - 2)
  const frac = x - i
  return R[i] * (1 - frac) + R[i + 1] * frac
}

/** Per-mode forced participation amp_mn = R_mn(r0/a)/sqrt((f_mn²−f²)²+(2ζ f_mn f)²). */
function modeAmp(md: PlateMode, f: number, z: number, r0: number): number {
  const fmn = md.f_hz
  const denom = Math.sqrt(
    (fmn * fmn - f * f) ** 2 + (2 * z * fmn * f) ** 2,
  )
  return interpR(md.R, r0) / denom
}

export interface PlateField {
  data: Float32Array // normalized |field|, 0..1; NaN outside the disk
  w: number
  h: number
  maxAbs: number // pre-normalization max |field| (0 ⇒ nothing excited)
}

/**
 * Forced steady-state field over a square canvas containing the unit disk.
 * Drive point given in [0,1]² UI coords (0.5,0.5 = plate center).
 * Pixels outside the disk are marked NaN (transparent).
 */
export function forcedResponseField(
  freqHz: number,
  drivePosXY: [number, number],
  zetaValue: number,
  w = 160,
  h = 160,
): PlateField {
  // UI (0..1) → disk coords (−1..1); radius 1 = plate edge.
  const dx = (drivePosXY[0] - 0.5) * 2
  const dy = (drivePosXY[1] - 0.5) * 2
  const r0 = Math.min(1, Math.hypot(dx, dy))
  const theta0 = Math.atan2(dy, dx)

  const amps = modes.map((md) => modeAmp(md, freqHz, zetaValue, r0))

  const raw = new Float32Array(w * h)
  let maxAbs = 0
  for (let py = 0; py < h; py++) {
    for (let px = 0; px < w; px++) {
      const x = (px / (w - 1)) * 2 - 1
      const y = (py / (h - 1)) * 2 - 1
      const r = Math.hypot(x, y)
      const idx = py * w + px
      if (r > 1) {
        raw[idx] = NaN
        continue
      }
      const theta = Math.atan2(y, x)
      let val = 0
      for (let k = 0; k < modes.length; k++) {
        const md = modes[k]
        val += amps[k] * interpR(md.R, r) * Math.cos(md.m * (theta - theta0))
      }
      const a = Math.abs(val)
      raw[idx] = a
      if (a > maxAbs) maxAbs = a
    }
  }
  const data = new Float32Array(w * h)
  for (let i = 0; i < raw.length; i++) {
    data[i] = Number.isNaN(raw[i]) ? NaN : maxAbs > 0 ? raw[i] / maxAbs : 0
  }
  return { data, w, h, maxAbs }
}

/** Drive point [0,1]² → (r0/a, θ0) on the unit disk. */
function driveRadius(drivePosXY: [number, number]): number {
  const dx = (drivePosXY[0] - 0.5) * 2
  const dy = (drivePosXY[1] - 0.5) * 2
  return Math.min(1, Math.hypot(dx, dy))
}

/**
 * The strongest resonance this drive point can reach: the eigenmode maximizing
 * the on-resonance participation |R_mn(r0)| / (2ζ f_mn²). Returns that mode's
 * frequency and the field's absolute peak there — the reference the live field
 * is displayed against, so OFF-resonance drives read as (honestly) near-dark.
 */
export function referencePeak(
  drivePosXY: [number, number],
  zetaValue: number,
): { fRef: number; peak: number } {
  const r0 = driveRadius(drivePosXY)
  let best = -Infinity
  let fRef = modes[0].f_hz
  for (const md of modes) {
    const part = Math.abs(interpR(md.R, r0)) / (2 * zetaValue * md.f_hz * md.f_hz)
    if (part > best) { best = part; fRef = md.f_hz }
  }
  const peak = forcedResponseField(fRef, drivePosXY, zetaValue).maxAbs
  return { fRef, peak }
}
