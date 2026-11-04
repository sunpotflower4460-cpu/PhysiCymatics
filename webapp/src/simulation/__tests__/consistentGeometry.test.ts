// M6-14: geometrically consistent plate surface (research path).
// The production geometry treats nodal eigenfunction samples as
// piecewise-constant bins for height but uses sampled derivatives
// for gradients — two different surfaces. The consistent path uses
// cubic Hermite interpolation on the same nodal data, so height,
// gradient, and Hessian all derive from ONE function.
import { describe, expect, it } from 'vitest'
import { singleToneScore } from '../sequence'
import { TransientPlate } from '../plateTransient'
import { calibratedForceN } from '../transientSand'

function excited(): TransientPlate {
  const plate = new TransientPlate(singleToneScore(152, calibratedForceN(), 30), [0.05, 0.03])
  let seed = 4242
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
  for (let j = 0; j < plate.ndof; j++) { plate.q[j] = rnd() * 2e-4; plate.qd[j] = rnd() * 0.3 }
  plate.beginStep(0, 1 / (44 * 152))
  return plate
}

const TAU = 1 / (44 * 152) / 2

describe('M6-14 consistent geometry', () => {
  it('radialInterp reproduces nodal values and slopes exactly', () => {
    const plate = excited()
    const out = new Float64Array(3)
    const R = (plate as unknown as { R: Float64Array[] }).R
    const dR = (plate as unknown as { dR: Float64Array[] }).dR
    for (let mi = 0; mi < R.length; mi++) {
      for (const k of [0, 1, 100, 255, 510, 511]) {
        plate.radialInterp(mi, k / 511, out)
        expect(out[0]).toBeCloseTo(R[mi][k], 15)
        expect(out[1]).toBeCloseTo(dR[mi][k], 12)
      }
    }
  })

  it('height and gradient are continuous across former bin edges', () => {
    const plate = excited()
    const o = new Float64Array(5)
    const th = 0.3; const cx = Math.cos(th); const cy = Math.sin(th)
    let maxJumpH = 0; let maxJumpG = 0
    let maxJumpHOld = 0
    for (let k = 1; k < 511; k++) {
      const re = (k / 511) * plate.a
      plate.surfAtConsistent((re - 1e-9) * cx, (re - 1e-9) * cy, TAU, o)
      const hl = o[0]; const gxl = o[1]; const gyl = o[2]
      plate.surfAtConsistent((re + 1e-9) * cx, (re + 1e-9) * cy, TAU, o)
      maxJumpH = Math.max(maxJumpH, Math.abs(o[0] - hl))
      maxJumpG = Math.max(maxJumpG, Math.hypot(o[1] - gxl, o[2] - gyl))
      plate.surfAt((re - 1e-9) * cx, (re - 1e-9) * cy, TAU, o); const hold = o[0]
      plate.surfAt((re + 1e-9) * cx, (re + 1e-9) * cy, TAU, o)
      maxJumpHOld = Math.max(maxJumpHOld, Math.abs(o[0] - hold))
    }
    expect(maxJumpHOld).toBeGreaterThan(5e-7) // old: ~um jumps exist
    expect(maxJumpH).toBeLessThan(1e-11)
    expect(maxJumpG).toBeLessThan(1e-6)
  })

  it('gradient and Hessian are the derivatives of its own height', () => {
    const plate = excited()
    const o = new Float64Array(5)
    const f = { h: 0, hx: 0, hy: 0, hxx: 0, hxy: 0, hyy: 0, ht: 0, htt: 0, htx: 0, hty: 0 }
    let seed = 31
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    const d = 2e-6
    const H = (x: number, y: number) => { plate.surfAtConsistent(x, y, TAU, o); return o[0] }
    const G = (x: number, y: number) => { plate.surfAtConsistent(x, y, TAU, o); return [o[1], o[2]] }
    let maxGradErr = 0; let maxHessErr = 0; let n = 0
    for (let k = 0; k < 400; k++) {
      const x = (rnd() * 2 - 1) * plate.a * 0.92
      const y = (rnd() * 2 - 1) * plate.a * 0.92
      if (Math.hypot(x, y) > plate.a * 0.92 || Math.hypot(x, y) < 1e-4) continue
      plate.surfAtConsistent(x, y, TAU, o)
      const fdX = (H(x + d, y) - H(x - d, y)) / (2 * d)
      const fdY = (H(x, y + d) - H(x, y - d)) / (2 * d)
      maxGradErr = Math.max(maxGradErr, Math.abs(fdX - o[1]), Math.abs(fdY - o[2]))
      plate.surfaceFieldConsistent(x, y, TAU, f)
      const [gxx1] = G(x + d, y); const [gxx0] = G(x - d, y)
      const [, gyy1] = G(x, y + d); const [, gyy0] = G(x, y - d)
      const [gxy1] = G(x, y + d); const [gxy0] = G(x, y - d)
      maxHessErr = Math.max(maxHessErr,
        Math.abs((gxx1 - gxx0) / (2 * d) - f.hxx),
        Math.abs((gyy1 - gyy0) / (2 * d) - f.hyy),
        Math.abs((gxy1 - gxy0) / (2 * d) - f.hxy))
      n++
    }
    expect(n).toBeGreaterThan(200)
    // FD truncation on a cubic-in-r interpolant: O(d^2 * 4th deriv).
    expect(maxGradErr).toBeLessThan(2e-5)
    expect(maxHessErr).toBeLessThan(0.5)
  })

  it('consistent surfAt differs from bin surfAt only at geometry scale', () => {
    const plate = excited()
    const o1 = new Float64Array(5); const o2 = new Float64Array(5)
    let seed = 7
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
    let maxDh = 0; let maxDg = 0
    for (let k = 0; k < 2000; k++) {
      const x = (rnd() * 2 - 1) * plate.a * 0.95
      const y = (rnd() * 2 - 1) * plate.a * 0.95
      if (Math.hypot(x, y) > plate.a * 0.95) continue
      plate.surfAt(x, y, TAU, o1)
      plate.surfAtConsistent(x, y, TAU, o2)
      maxDh = Math.max(maxDh, Math.abs(o1[0] - o2[0]))
      maxDg = Math.max(maxDg, Math.hypot(o1[1] - o2[1], o1[2] - o2[2]))
    }
    // Height differences are micrometre-scale (bin quantization),
    // gradient differences are the same order as gradients themselves.
    expect(maxDh).toBeGreaterThan(1e-7)
    expect(maxDh).toBeLessThan(2e-5)
    expect(maxDg).toBeGreaterThan(1e-4)
  })
})
