// shape.ts — Tier2.5: pattern-shape readout (what the computation SELECTED).
//
// PUT-IN     : pack.pattern_shape + the c005 direct-PDE snapshot.
// EMERGED    : the lattice symmetry c005 actually selected (capillary → square,
//              measured & Chen–Viñals cross-checked). Mixed/low-f is NOT yet
//              converged → shown as 未確定 (display_rule: show only what the
//              computation selected).
// CLAIM-TIER : capillary square = measured/established; mixed = frontier.
import snapshot from '../data/c005_snapshot.json'
import { patternShape } from './pack'

// Capillary vs gravity crossover for clean water: λ_c = 2π√(σ/(ρg)) ≈ 17 mm.
// Deep-capillary (square earned & Chen–Viñals-matched) only well below that.
const CAPILLARY_LAM_MM = 6.0

export interface ShapeView {
  regime: 'capillary' | 'mixed' | 'gravity'
  label: string
  tier: string
  note: string
  squarePairDeg?: number
  squareShare?: number
}

export function shapeAt(lamMm: number): ShapeView {
  const ps = patternShape()
  if (lamMm <= CAPILLARY_LAM_MM) {
    return {
      regime: 'capillary',
      label: '正方格子（SQUARE）',
      tier: 'measured / established (Chen-Viñals交差検証)',
      note: ps.capillary_regime.literature_match,
      squarePairDeg: (snapshot as { top2_pair_angle: number }).top2_pair_angle,
      squareShare: (snapshot as { top2_share: number }).top2_share,
    }
  }
  if (lamMm >= 25.0) {
    return {
      regime: 'gravity',
      label: '未確定（重力域・縞予想）',
      tier: 'frontier — CI収束待ち',
      note: 'Chen-Viñals予想は縞。c005での確定はCI夜間走行待ち。',
    }
  }
  return {
    regime: 'mixed',
    label: '未確定（混合域）',
    tier: 'frontier — CI収束待ち',
    note: ps.mixed_regime.literature_match,
  }
}
