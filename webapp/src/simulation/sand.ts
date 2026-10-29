// sand.ts — c002 bounce mechanics, TypeScript port of reference/sand_bounce.py.
//
// PUT-IN     : gravity; exact ballistic flight; impact on the moving tilted
//              surface (e, beta); sit / detach by support sign. Same update
//              order and bisection count as the Python reference.
// EMERGED    : (c002) transport toward low |W| — measured in stats(), never fed back.
// CLAIM-TIER : measured (within the c002 numeric model; port checked vs fixture).
// FLOORS     : no air / grain-grain / back-reaction; e/beta/vStick assumed;
//              single mode; JS RNG differs from NumPy (shared fixtures use arrays).
// Dynamics never references a target radius; measurement lives in stats().
import { G, SingleModeSurface } from './surface'
import { DEFAULT_PARAMS, type ParticleState, type SandParams, type SandSnapshot, type SandStats } from './types'

/** Deterministic RNG (mulberry32) for UI seeding. Not NumPy-compatible. */
export function mulberry32(seed: number): () => number {
  let a = seed >>> 0
  return () => {
    a = (a + 0x6d2b79f5) >>> 0
    let t = a
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

export function initialState(
  a: number, count: number, seed: number, initXY?: { x: ArrayLike<number>; y: ArrayLike<number> },
): ParticleState {
  const n = initXY ? initXY.x.length : count
  const x = new Float64Array(n); const y = new Float64Array(n)
  if (initXY) {
    x.set(initXY.x); y.set(initXY.y)
  } else {
    const rng = mulberry32(seed)
    for (let i = 0; i < n; i++) {
      const rr = a * Math.sqrt(rng())
      const th = 2 * Math.PI * rng()
      x[i] = rr * Math.cos(th); y[i] = rr * Math.sin(th)
    }
  }
  const z = new Float64Array(n).fill(5e-4)
  const id = new Uint32Array(n)
  for (let i = 0; i < n; i++) id[i] = i
  return {
    x, y, z, vx: new Float64Array(n), vy: new Float64Array(n), vz: new Float64Array(n),
    alive: new Uint8Array(n).fill(1), sitting: new Uint8Array(n), id,
  }
}

const S5 = new Float64Array(5)

export class SandSim {
  t = 0
  stepIndex = 0
  nBounce = 0
  constructor(
    readonly surface: SingleModeSurface,
    readonly state: ParticleState,
    readonly p: SandParams = DEFAULT_PARAMS,
  ) {}

  get dt(): number { return (2 * Math.PI) / this.surface.omega / this.p.spp }
  get vStick(): number { return this.p.vStickFrac * this.surface.A * this.surface.omega }
  get nAlive(): number { let s = 0; for (const v of this.state.alive) s += v; return s }

  step(): void {
    const { surface: s, state: st, p } = this
    const dt = this.dt; const t = this.t; const w = s.omega; const A = s.A
    const tNew = (this.stepIndex + 1) * dt
    const s_ = Math.sin(w * tNew); const c_ = Math.cos(w * tNew)
    const n = st.id.length
    for (let i = 0; i < n; i++) {
      if (!st.alive[i] || st.sitting[i]) continue
      const xf = st.x[i] + st.vx[i] * dt
      const yf = st.y[i] + st.vy[i] * dt
      const zf = st.z[i] + st.vz[i] * dt - 0.5 * G * dt * dt
      const vzf = st.vz[i] - G * dt
      s.surf(xf, yf, s_, c_, S5)
      if (zf >= S5[0]) {
        st.x[i] = xf; st.y[i] = yf; st.z[i] = zf; st.vz[i] = vzf
        continue
      }
      // Penetration: bisect contact time in [t, tNew] from step-start state.
      const x0 = st.x[i]; const y0 = st.y[i]; const z0 = st.z[i]
      const vx0 = st.vx[i]; const vy0 = st.vy[i]; const vz0 = st.vz[i]
      let t0 = t; let t1 = tNew
      for (let k = 0; k < p.bisectIters; k++) {
        const tm = 0.5 * (t0 + t1); const ddt = tm - t
        s.surf(x0 + vx0 * ddt, y0 + vy0 * ddt, Math.sin(w * tm), Math.cos(w * tm), S5)
        if (z0 + vz0 * ddt - 0.5 * G * ddt * ddt < S5[0]) t1 = tm; else t0 = tm
      }
      const tc = 0.5 * (t0 + t1); const ddt = tc - t
      const xc = x0 + vx0 * ddt; const yc = y0 + vy0 * ddt
      const vzc = vz0 - G * ddt
      s.surf(xc, yc, Math.sin(w * tc), Math.cos(w * tc), S5)
      const hc = S5[0]; const hx = S5[1]; const hy = S5[2]; const uz = S5[3]
      const nn = Math.sqrt(hx * hx + hy * hy + 1)
      const nx = -hx / nn; const ny = -hy / nn; const nz = 1 / nn
      const rvx = vx0; const rvy = vy0; const rvz = vzc - uz
      const vn = rvx * nx + rvy * ny + rvz * nz
      if (vn >= 0) continue // not a hit: keep step-start state (as c002b)
      const k1 = (1 + p.e) * vn
      const r2x = rvx - k1 * nx; const r2y = rvy - k1 * ny; const r2z = rvz - k1 * nz
      const vn2 = r2x * nx + r2y * ny + r2z * nz
      st.x[i] = xc; st.y[i] = yc; st.z[i] = hc + 1e-9
      st.vx[i] = vn2 * nx + p.beta * (r2x - vn2 * nx)
      st.vy[i] = vn2 * ny + p.beta * (r2y - vn2 * ny)
      st.vz[i] = vn2 * nz + p.beta * (r2z - vn2 * nz) + uz
      this.nBounce++
      if (Math.abs(vn2) < this.vStick) {
        st.sitting[i] = 1
        if (p.cRoll === 0) { st.vx[i] = 0; st.vy[i] = 0 }
      }
    }
    for (let i = 0; i < n; i++) {
      if (!st.alive[i] || !st.sitting[i]) continue
      s.surf(st.x[i], st.y[i], s_, c_, S5)
      st.z[i] = S5[0]
      const nAcc = Math.max(G - A * w * w * S5[4] * s_, 0)
      if (p.cRoll > 0) {
        st.vx[i] += -p.cRoll * nAcc * S5[1] * dt
        st.vy[i] += -p.cRoll * nAcc * S5[2] * dt
        const sp = Math.hypot(st.vx[i], st.vy[i]) + 1e-30
        const drag = Math.min(p.cRr * nAcc * dt, sp)
        st.vx[i] -= (drag * st.vx[i]) / sp
        st.vy[i] -= (drag * st.vy[i]) / sp
        st.x[i] += st.vx[i] * dt; st.y[i] += st.vy[i] * dt
      }
      if (A * w * w * S5[4] * s_ > G) {
        st.sitting[i] = 0
        st.vz[i] = S5[3]; st.z[i] = S5[0] + 1e-9
      }
    }
    for (let i = 0; i < n; i++) {
      if (st.alive[i] && Math.hypot(st.x[i], st.y[i]) > s.a) st.alive[i] = 0
    }
    this.t = tNew
    this.stepIndex++
  }

  advanceTo(tTarget: number): void {
    while (this.t < tTarget - 1e-15) this.step()
  }

  // ----- measurement only (reads state; never feeds back into step) -----
  stats(stepsPerSec = 0): SandStats {
    const { surface: s, state: st } = this
    const g: [number, number, number] = [0, 0, 0]
    const absW: number[] = []
    let ring = 0; let seated = 0; let alive = 0
    const band = 0.005
    for (let i = 0; i < st.id.length; i++) {
      if (!st.alive[i]) continue
      alive++
      s.wGrad(st.x[i], st.y[i], g)
      absW.push(Math.abs(g[0]))
      if (Number.isFinite(s.rNode) && Math.abs(Math.hypot(st.x[i], st.y[i]) - s.rNode) < band) ring++
      if (st.sitting[i]) seated++
    }
    absW.sort((a, b) => a - b)
    const med = absW.length ? absW[Math.floor(absW.length / 2)] : NaN
    const trapped = absW.length ? absW.filter((v) => s.gamma * v < 1).length / absW.length : NaN
    const expFrac = Number.isFinite(s.rNode) ? (4 * s.rNode * band) / (s.a * s.a) : NaN
    return {
      t: Math.round(this.t * 100) / 100, medAbsW: med, trapped,
      ringEnrich: alive && Number.isFinite(expFrac) ? ring / alive / expFrac : NaN,
      sitting: alive ? seated / alive : NaN, alive, nBounce: this.nBounce, stepsPerSec,
    }
  }

  snapshot(stepsPerSec = 0): SandSnapshot {
    const st = this.state
    const xy = new Float32Array(this.nAlive * 2)
    let k = 0
    for (let i = 0; i < st.id.length; i++) {
      if (!st.alive[i]) continue
      xy[k++] = st.x[i]; xy[k++] = st.y[i]
    }
    return { t: this.t, xy, aliveCount: k / 2, initialCount: st.id.length, stats: this.stats(stepsPerSec) }
  }
}

export function createSim(gamma: number, count: number, seed: number): SandSim {
  const surface = new SingleModeSurface(gamma)
  return new SandSim(surface, initialState(surface.a, count, seed))
}
