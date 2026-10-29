// types.ts — sand simulation data contracts (no React / Canvas / Audio deps).
export interface ParticleState {
  x: Float64Array; y: Float64Array; z: Float64Array
  vx: Float64Array; vy: Float64Array; vz: Float64Array
  alive: Uint8Array; sitting: Uint8Array; id: Uint32Array
}

export interface SandParams {
  e: number; beta: number; cRoll: number; cRr: number
  spp: number; vStickFrac: number; bisectIters: number
}

export const DEFAULT_PARAMS: SandParams = {
  e: 0.75, beta: 0.8, cRoll: 0, cRr: 0.02, spp: 44, vStickFrac: 0.3, bisectIters: 14,
}

export interface SandConfig { gamma: number; count: number; seed: number }

export interface SandStats {
  t: number; medAbsW: number; trapped: number; ringEnrich: number
  sitting: number; alive: number; nBounce: number; stepsPerSec: number
}

/** Alive-particle positions for rendering (Float32 copy; physics stays Float64). */
export interface SandSnapshot {
  t: number; xy: Float32Array; aliveCount: number
  initialCount: number; stats: SandStats
}
