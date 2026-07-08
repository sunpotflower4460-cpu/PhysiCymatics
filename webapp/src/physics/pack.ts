// pack.ts — typed access to the (trimmed) physics_pack. READ-ONLY.
//
// PUT-IN     : physics_pack/PhysiCymatics_physics_pack_v1.json (v1_4), subset.
// EMERGED    : nothing — faithful reader; no physics invented here.
// CLAIM-TIER : established (values are the pack's; provenance in the pack).
// Rendering-honesty (CLAUDE.md §7): this layer only reads pack numbers.
import raw from '../data/physics_pack.trimmed.json'

export interface PlateMode { m: number; n: number; f_hz: number; R: number[] }

interface Pack {
  schema: string
  plate_modal_shapes: {
    modes: PlateMode[]
    zeta: { value: number; status: string }
    formula: string
  }
  plate_detents: {
    circular_free_steel_24cm: {
      eigenfrequencies_hz: number[]
      modes: [number, number][]
      center_drive_reachable_hz: number[]
    }
  }
  fluid_water: {
    constants: { rho_kgm3: number; sigma_nm: number; nu_m2s: number }
    response: string
    threshold_curve_h20mm: ThresholdPoint[]
    threshold_curve_h5mm: ThresholdPoint[]
    app_default_440hz: { a_c_ms2: number; lam_mm: number; x_c_um: number }
  }
  pattern_shape: {
    capillary_regime: { f_example: number; result: string; tier: string; literature_match: string }
    mixed_regime: { f_example: number; result_snapshot: string; tier: string; literature_match: string }
    display_rule: string
  }
}

export interface ThresholdPoint {
  f_hz: number
  a_c_ms2: number
  band_ms2: [number, number]
  lam_mm: number
  subharmonic: boolean
  x_c_um: number
}

export const PACK = raw as unknown as Pack
export const SCHEMA = PACK.schema

// --- plate ---
export const plateModes = (): PlateMode[] => PACK.plate_modal_shapes.modes
export const plateZeta = () => PACK.plate_modal_shapes.zeta            // {value, status}
export const eigenfrequenciesHz = (): number[] =>
  PACK.plate_detents.circular_free_steel_24cm.eigenfrequencies_hz
export const centerDriveReachableHz = (): number[] =>
  PACK.plate_detents.circular_free_steel_24cm.center_drive_reachable_hz

// --- water ---
export const thresholdCurve = (depth_m: number): ThresholdPoint[] =>
  depth_m <= 0.01 ? PACK.fluid_water.threshold_curve_h5mm
                  : PACK.fluid_water.threshold_curve_h20mm
export const faradayResponse = () => PACK.fluid_water.response

// --- shape ---
export const patternShape = () => PACK.pattern_shape
