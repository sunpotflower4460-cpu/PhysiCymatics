// Golden: water threshold interpolation must match the c003 oracle sweep <0.1%.
import { describe, it, expect } from 'vitest'
import { thresholdAt, faradayState } from '../faraday'
import truth from '../../data/truth/c003_threshold_h20_bulk.json'

const pts = (truth as { points: { f_hz: number; a_c_ms2: number; lam_mm: number }[] }).points

describe('faraday threshold vs c003 oracle (h=20mm bulk)', () => {
  it('interpolated a_c matches the oracle at each curve node to < 0.1%', () => {
    for (const p of pts) {
      const th = thresholdAt(p.f_hz, 0.020)
      expect(Math.abs(th.a_c_ms2 - p.a_c_ms2) / p.a_c_ms2).toBeLessThan(1e-3)
      expect(Math.abs(th.lam_mm - p.lam_mm) / p.lam_mm).toBeLessThan(1e-3)
    }
  })

  it('response is subharmonic (f/2) and silent below threshold', () => {
    const th = thresholdAt(60, 0.020)
    expect(th.subharmonic).toBe(true)
    const below = faradayState(60, 0.020, th.a_c_ms2 * 0.5)
    const above = faradayState(60, 0.020, th.a_c_ms2 * 1.5)
    expect(below.active).toBe(false)
    expect(above.active).toBe(true)
    expect(above.respHz).toBeCloseTo(30, 6)
  })

  it('flags > 5× threshold as beyond linear scope', () => {
    const th = thresholdAt(60, 0.020)
    expect(faradayState(60, 0.020, th.a_c_ms2 * 6).overLinear).toBe(true)
  })
})
