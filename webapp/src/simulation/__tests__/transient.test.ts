// M3+M6-1 TS gates: score/PCM timing, phase continuity, free-frequency path,
// steady-state response, center-drive symmetry, state continuity on switch,
// exact target time, residual time after contact, Python parity (M6 fixture).
import { describe, it, expect } from 'vitest'
import { ScoreSignal, renderScorePcm, singleToneScore, type PerformanceScore } from '../sequence'
import { TransientPlate } from '../plateTransient'
import { TransientSandSim, createTransientSim, calibratedForceN } from '../transientSand'
import { initialState } from '../sand'
import { G } from '../surface'
import parity from '../../data/truth/transient_m6_parity.json'

const SCORE: PerformanceScore = {
  schema: 'physicymatics.score.v1', duration_s: 0.8, initialPhase_rad: 0, audioGain: 0.2,
  events: [
    { start_s: 0, end_s: 0.4, frequency_hz: 152, forcePeak_N: 0.05, attack_s: 0.05, release_s: 0.05 },
    { start_s: 0.4, end_s: 0.8, frequency_hz: 649, forcePeak_N: 0.05, attack_s: 0.05, release_s: 0.05 },
  ],
}

describe('score + PCM', () => {
  it('phase is continuous across the event boundary', () => {
    const s = new ScoreSignal(SCORE)
    expect(Math.abs(s.u(0.4 - 1e-7) - s.u(0.4 + 1e-7))).toBeLessThan(1e-3)
    expect(s.freqAt(0.39)).toBe(152)
    expect(s.freqAt(0.41)).toBe(649)
  })

  it('PCM length and event onset align to samples', () => {
    const pcm = renderScorePcm(SCORE, 48000)
    expect(pcm.length).toBe(Math.round(0.8 * 48000))
    const s = new ScoreSignal(SCORE)
    const k = Math.round(0.4 * 48000)
    expect(Math.abs(pcm[k] - 0.2 * s.u(k / 48000))).toBeLessThan(1e-6)
    expect(Math.max(...pcm.map(Math.abs))).toBeLessThanOrEqual(0.2 + 1e-6)
  })

  it('arbitrary non-preset frequency passes through unchanged', () => {
    const s = new ScoreSignal(singleToneScore(173.31, 0.05, 10))
    expect(s.freqAt(1)).toBe(173.31)
  })
})

describe('transient plate', () => {
  it('center drive never excites m>0', () => {
    const pl = new TransientPlate(singleToneScore(152, 0.05, 10))
    const dt = 1 / (44 * 152)
    while (pl.t < 0.2) pl.beginStep(pl.t, pl.t + dt)
    pl.dofs.forEach((d, j) => {
      if (pl.modes[d.mode].m > 0) expect(pl.q[j]).toBe(0)
    })
  })

  it('off-resonance steady amplitude matches the analytic DOF response', () => {
    const F = calibratedForceN()
    const i152probe = new TransientPlate(singleToneScore(152, F, 30))
    const i152 = i152probe.modeIndex152()
    const f = i152probe.modes[i152].f_hz
    const fDrive = f * 1.05
    const pl = new TransientPlate(singleToneScore(fDrive, F, 30))
    const j = pl.dofs.findIndex((d) => d.mode === i152)
    const w = pl.omega[j]; const wd = 2 * Math.PI * fDrive
    const ampRef = (Math.abs(pl.phiDrive[j]) * F) / pl.mass[j]
      / Math.sqrt((w * w - wd * wd) ** 2 + (2 * pl.zeta * w * wd) ** 2)
    const dt = 1 / (44 * fDrive)
    const qs: number[] = []
    while (pl.t < 6) {
      pl.beginStep(pl.t, pl.t + dt)
      if (pl.t > 5) qs.push(pl.q[j])
    }
    const amp = (Math.max(...qs) - Math.min(...qs)) / 2
    expect(Math.abs(amp - ampRef) / ampRef).toBeLessThan(0.02)
  })

  it('frequency switch keeps plate and sand state continuous', () => {
    const sim = createTransientSim(SCORE, 30, 2)
    while (sim.t < 0.39) sim.step()
    const xBefore = Float64Array.from(sim.state.x)
    while (sim.t < 0.41) sim.step()
    let maxJump = 0
    for (let i = 0; i < 30; i++) {
      expect(Number.isFinite(sim.state.x[i])).toBe(true)
      maxJump = Math.max(maxJump, Math.abs(sim.state.x[i] - xBefore[i]))
    }
    expect(maxJump).toBeLessThan(0.01)
    expect(sim.nAlive).toBeLessThanOrEqual(30)
  })
})

describe('M6-1 time and contact', () => {
  it('advanceTo ends exactly at the target time', () => {
    const sim = createTransientSim(singleToneScore(152, 0.05, 10), 5, 1)
    sim.advanceTo(0.137)
    expect(Math.abs(sim.t - 0.137)).toBeLessThan(1e-12)
  })

  it('residual time after a flat-plate bounce is flown', () => {
    const plate = new TransientPlate(singleToneScore(152, 0, 10))
    const st = initialState(plate.a, 0, 0, { x: [0], y: [0] })
    st.z[0] = 0.02
    const sim = new TransientSandSim(plate, st)
    const dt = 0.12
    sim.step(dt)
    const tc = Math.sqrt((2 * 0.02) / G)
    const vUp = 0.75 * G * tc
    const rem = dt - tc
    expect(Math.abs(st.z[0] - (vUp * rem - 0.5 * G * rem * rem))).toBeLessThan(1e-6)
    expect(Math.abs(sim.t - dt)).toBeLessThan(1e-15)
  })

  it('setFrequency keeps phase and plate state', () => {
    const pl = new TransientPlate(singleToneScore(152, 0.05, 30))
    const dt = 1 / (44 * 152)
    while (pl.t < 0.2) pl.beginStep(pl.t, pl.t + dt)
    const qBefore = Float64Array.from(pl.q)
    const uBefore = pl.signal.u(pl.scoreTime(pl.t))
    pl.setFrequency(173.31, pl.t)
    expect(Math.max(...pl.q.map((v, i) => Math.abs(v - qBefore[i])))).toBe(0)
    expect(Math.abs(pl.signal.u(pl.scoreTime(pl.t)) - uBefore)).toBeLessThan(1e-9)
    expect(pl.signal.freqAt(pl.scoreTime(pl.t))).toBe(173.31)
  })

  it('matches the Python M6 parity fixture within tolerance', () => {
    const P = parity as unknown as {
      force_N: number; events: PerformanceScore['events']
      init_x: number[]; init_y: number[]; x: number[]; y: number[]
      alive: boolean[]; n_bounce: number; n_detach: number
    }
    const score: PerformanceScore = {
      schema: 'physicymatics.score.v1', duration_s: 0.5, initialPhase_rad: 0, audioGain: 0.15,
      events: P.events,
    }
    const plate = new TransientPlate(score)
    const st = initialState(plate.a, 0, 0, { x: P.init_x, y: P.init_y })
    // Cross-language physics parity is gated with sleep OFF on both sides
    // (fixture generated with sleep_enabled=False); sleep on/off parity is
    // gated separately per language in sleep.test.ts / test_m6_sleep.py.
    const sim = new TransientSandSim(plate, st, undefined, { enabled: false })
    sim.advanceTo(0.5)
    expect(Math.abs(sim.t - 0.5)).toBeLessThan(1e-12)
    const dx = Math.max(...P.x.map((v, i) => Math.abs(v - st.x[i])))
    const dy = Math.max(...P.y.map((v, i) => Math.abs(v - st.y[i])))
    // Pre-declared tolerance: identical algorithm in float64; only libm /
    // summation-order differences are expected over 0.5 s of collisions.
    expect(dx).toBeLessThan(1e-6)
    expect(dy).toBeLessThan(1e-6)
    expect(sim.nBounce).toBe(P.n_bounce)
    expect(sim.nDetach).toBe(P.n_detach)
    expect([...st.alive].map(Boolean)).toEqual(P.alive)
  })
})
