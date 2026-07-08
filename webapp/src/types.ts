// Config contract (ported from ui_shell/App.jsx, typed).
export interface Config {
  schema: string
  run_id: string
  media: 'plate' | 'fluid'
  drive: { type: string; pos_xy: [number, number]; freq_hz: number; accel_ms2: number }
  plate: {
    shape: string; size_m: number; thickness_m: number; material: string
    E_pa: number; rho_kgm3: number; nu: number; bc: string; damping_zeta: number
  }
  fluid: { depth_m: number; viscosity_pas: number; surface_tension_nm: number }
  particles: { kind: string; diameter_um: number; count: number }
  seed: number
  app_version: string
  physics_core: string
}
