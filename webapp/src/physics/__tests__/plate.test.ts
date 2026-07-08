// Golden: plate eigenfrequencies (from pack) must equal the c001 oracle rel<1e-6.
import { describe, it, expect } from 'vitest'
import { eigenfrequencies } from '../plate'
import truth from '../../data/truth/c001_eigenfrequencies.json'

describe('plate eigenfrequencies vs c001 oracle', () => {
  it('matches the oracle to rel < 1e-6', () => {
    const got = eigenfrequencies()
    const want = (truth as { eigenfrequencies_hz: number[] }).eigenfrequencies_hz
    expect(got.length).toBe(want.length)
    got.forEach((f, i) => {
      expect(Math.abs(f - want[i]) / want[i]).toBeLessThan(1e-6)
    })
  })

  it('center-drive first ring ≈ 152 Hz and 90.5 Hz present as a mode', () => {
    const eig = eigenfrequencies()
    expect(eig.some((f) => Math.abs(f - 152.0) < 1.0)).toBe(true)
    expect(eig.some((f) => Math.abs(f - 90.5) < 1.0)).toBe(true)
  })
})
