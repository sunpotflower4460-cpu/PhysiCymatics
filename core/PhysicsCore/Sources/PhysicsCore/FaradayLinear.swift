// FaradayLinear.swift — Faraday instability threshold for water (linear).
//
// PUT-IN     : linearized free-surface dynamics per wavenumber k under a
//              vertically oscillating gravity g(t)=g(1+Γcos ω_d t):
//                η_tt + 2γ_k η_t + [(g(t)k + σk³/ρ)tanh(kh)] η = 0.
//              water @20°C constants (established). damping is the stated
//              sub-model (bulk γ=2νk²; optional +bottom boundary layer).
// FORBIDDEN  : no "f/2", no threshold value, no wavelength is put in. We
//              integrate ONE drive period and ask if a cycle amplifies.
// EMERGED    : threshold Γ_c (bisection over drive amplitude), the subharmonic
//              response (from the multiplier sign), pattern scale = argmin_k.
// CLAIM-TIER : measured (within the linear/damping model). Cross-checks:
//              Liouville det(M)=exp(−2γT) and a weak-damping analytic threshold.
// FLOORS     : linearized (no saturation); damping is a sub-model; boundary-
//              layer coefficient is memory-based (YELLOW).
//
// Faithful line-by-line port of reference/faraday_linear.py. Golden-tested
// against the Python-generated JSON to <0.1%.
import Foundation

public enum DampingVariant: Sendable {
    case bulk
    case bulkBL   // + bottom boundary layer
}

public struct FaradayResult: Sendable {
    public let fDrive: Double
    public let h: Double
    public let kc: Double
    public let gc: Double
    public let aCms2: Double
    public let lamMm: Double
    public let lamDispMm: Double
    public let subharmonic: Bool
    public let gcAnalyticAtKc: Double
    public let liouvilleMaxDev: Double
}

public struct FaradayLinear: Sendable {
    public static let g = 9.81
    public var rho: Double
    public var sigma: Double
    public var nu: Double

    /// Defaults are the established water @20°C constants; a `PhysicsPack` can
    /// override them from `fluid_water.constants` for provenance parity.
    public init(rho: Double = 998.0, sigma: Double = 0.0728, nu: Double = 1.004e-6) {
        self.rho = rho
        self.sigma = sigma
        self.nu = nu
    }

    public init(pack: PhysicsPack) {
        let c = pack.fluidWater.constants
        self.init(rho: c.rhoKgm3, sigma: c.sigmaNm, nu: c.nuM2s)
    }

    private func tanhKH(_ k: Double, _ h: Double) -> Double {
        tanh(Swift.min(k * h, 20.0))
    }

    /// Undamped dispersion ω₀²(k) = (gk + σk³/ρ)·tanh(kh).
    public func om0Sq(_ k: Double, _ h: Double) -> Double {
        (Self.g * k + sigma * k * k * k / rho) * tanhKH(k, h)
    }

    public func gammaK(_ k: Double, _ h: Double, _ variant: DampingVariant) -> Double {
        let g0 = 2.0 * nu * k * k
        switch variant {
        case .bulk:
            return g0
        case .bulkBL:
            let kh = Swift.min(k * h, 20.0)
            let om0 = om0Sq(k, h).squareRoot()
            return g0 + k * (nu * om0 / 2.0).squareRoot() / sinh(2.0 * kh)
        }
    }

    /// Integrate the real equation over ONE drive period for two ICs (RK4),
    /// vectorized over the k-grid. Returns (s, trace, det, gamma, T).
    public func monodromy(_ k: [Double], _ h: Double, _ gam: [Double],
                          _ Gam: [Double], _ wD: Double,
                          nsteps: Int = 1200)
        -> (s: [Double], tr: [Double], det: [Double], T: Double) {
        let n = k.count
        let om0s = k.map { om0Sq($0, h) }
        var mod = [Double](repeating: 0.0, count: n)
        for i in 0..<n { mod[i] = Gam[i] * Self.g * k[i] * tanhKH(k[i], h) }
        let T = 2.0 * Double.pi / wD
        let dt = T / Double(nsteps)

        // state: E,V each carry two independent ICs → (E0,E1),(V0,V1)
        var e0 = [Double](repeating: 1.0, count: n)   // IC A: E=1,V=0
        var v0 = [Double](repeating: 0.0, count: n)
        var e1 = [Double](repeating: 0.0, count: n)   // IC B: E=0,V=1
        var v1 = [Double](repeating: 1.0, count: n)

        func accel(_ e: Double, _ v: Double, _ i: Double_i) -> Double {
            -(om0s[i.idx] + mod[i.idx] * cos(wD * i.t)) * e - 2.0 * gam[i.idx] * v
        }

        var t = 0.0
        for _ in 0..<nsteps {
            for i in 0..<n {
                // RK4 on the 2×2 linear system, per IC, per k.
                func stepIC(_ e: Double, _ v: Double) -> (Double, Double) {
                    let ctx0 = Double_i(idx: i, t: t)
                    let ctxH = Double_i(idx: i, t: t + 0.5 * dt)
                    let ctx1 = Double_i(idx: i, t: t + dt)
                    let k1e = v,                       k1v = accel(e, v, ctx0)
                    let k2e = v + 0.5 * dt * k1v,      k2v = accel(e + 0.5 * dt * k1e, v + 0.5 * dt * k1v, ctxH)
                    let k3e = v + 0.5 * dt * k2v,      k3v = accel(e + 0.5 * dt * k2e, v + 0.5 * dt * k2v, ctxH)
                    let k4e = v + dt * k3v,            k4v = accel(e + dt * k3e, v + dt * k3v, ctx1)
                    let ne = e + dt / 6.0 * (k1e + 2 * k2e + 2 * k3e + k4e)
                    let nv = v + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
                    return (ne, nv)
                }
                (e0[i], v0[i]) = stepIC(e0[i], v0[i])
                (e1[i], v1[i]) = stepIC(e1[i], v1[i])
            }
            t += dt
        }

        var s = [Double](repeating: 0.0, count: n)
        var tr = [Double](repeating: 0.0, count: n)
        var det = [Double](repeating: 0.0, count: n)
        for i in 0..<n {
            let trace = e0[i] + v1[i]
            let d = e0[i] * v1[i] - e1[i] * v0[i]
            let disc = trace * trace - 4.0 * d
            let lam: Double
            if disc >= 0 {
                lam = (abs(trace) + disc.squareRoot()) / 2.0
            } else {
                lam = Swift.max(d, 1e-300).squareRoot()
            }
            tr[i] = trace
            det[i] = d
            s[i] = log(Swift.max(lam, 1e-300)) / T
        }
        return (s, tr, det, T)
    }

    /// Route 2: weak-damping subharmonic threshold with detuning.
    public func gcAnalytic(_ k: Double, _ h: Double, _ wD: Double,
                           _ variant: DampingVariant) -> Double {
        let om0 = om0Sq(k, h).squareRoot()
        let gam = gammaK(k, h, variant)
        let delta = om0 - 0.5 * wD
        let mc = 2.0 * wD * (gam * gam + delta * delta).squareRoot()
        return mc / (Self.g * k * tanhKH(k, h))
    }

    /// k such that ω₀(k) = 2π·fWave, by bisection (mirrors scipy brentq range).
    public func kOfFreq(_ fWave: Double, _ h: Double) -> Double {
        let w = 2.0 * Double.pi * fWave
        let target = w * w
        var lo = 1e-2, hi = 1e6
        for _ in 0..<200 {
            let mid = 0.5 * (lo + hi)
            if om0Sq(mid, h) - target > 0 { hi = mid } else { lo = mid }
        }
        return 0.5 * (lo + hi)
    }

    private func bisectCurve(_ kg: [Double], _ h: Double, _ wD: Double,
                             _ variant: DampingVariant, iters: Int) -> [Double] {
        let n = kg.count
        let gam = kg.map { gammaK($0, h, variant) }
        var Ghi = kg.map { 1.6 * gcAnalytic($0, h, wD, variant) + 1e-4 }
        for _ in 0..<8 {
            let (s, _, _, _) = monodromy(kg, h, gam, Ghi, wD)
            var anyNeed = false
            for i in 0..<n where s[i] <= 0 { Ghi[i] *= 1.7; anyNeed = true }
            if !anyNeed { break }
        }
        var Glo = [Double](repeating: 0.0, count: n)
        for _ in 0..<iters {
            var Gm = [Double](repeating: 0.0, count: n)
            for i in 0..<n { Gm[i] = 0.5 * (Glo[i] + Ghi[i]) }
            let (s, _, _, _) = monodromy(kg, h, gam, Gm, wD)
            for i in 0..<n {
                if s[i] > 0 { Ghi[i] = Gm[i] } else { Glo[i] = Gm[i] }
            }
        }
        return (0..<n).map { 0.5 * (Glo[$0] + Ghi[$0]) }
    }

    private static func logspace(_ a: Double, _ b: Double, _ n: Int) -> [Double] {
        (0..<n).map { pow(10.0, a + (b - a) * Double($0) / Double(n - 1)) }
    }
    private static func linspace(_ a: Double, _ b: Double, _ n: Int) -> [Double] {
        (0..<n).map { a + (b - a) * Double($0) / Double(n - 1) }
    }
    private static func argmin(_ v: [Double]) -> Int {
        var j = 0
        for i in 1..<v.count where v[i] < v[j] { j = i }
        return j
    }

    /// Three-stage k refinement (the water tongue is razor-thin; a single log
    /// grid aliases it). Mirrors reference/faraday_linear.sweep_f.
    public func sweepF(_ fD: Double, _ h: Double,
                       _ variant: DampingVariant = .bulk) -> FaradayResult {
        let wD = 2.0 * Double.pi * fD
        let kHalf = kOfFreq(fD / 2.0, h)
        let kg = Self.logspace(-0.6, 0.6, 90).map { $0 * kHalf }
        let gc = bisectCurve(kg, h, wD, variant, iters: 14)
        let kctr = kg[Self.argmin(gc)]
        let kg2 = Self.linspace(0.955, 1.045, 120).map { $0 * kctr }
        let gc2 = bisectCurve(kg2, h, wD, variant, iters: 18)
        let kctr2 = kg2[Self.argmin(gc2)]
        let kg3 = Self.linspace(0.994, 1.006, 80).map { $0 * kctr2 }
        let gc3 = bisectCurve(kg3, h, wD, variant, iters: 22)
        let j = Self.argmin(gc3)

        let gam3 = kg3.map { gammaK($0, h, variant) }
        let Gtest = gc3.map { $0 * 1.0002 }
        let (_, tr, det, T) = monodromy(kg3, h, gam3, Gtest, wD)
        var liou = 0.0
        for i in 0..<kg3.count {
            let expected = exp(-2 * gam3[i] * T)
            liou = Swift.max(liou, abs(det[i] - expected) / expected)
        }
        return FaradayResult(
            fDrive: fD, h: h,
            kc: kg3[j], gc: gc3[j], aCms2: gc3[j] * Self.g,
            lamMm: 2e3 * Double.pi / kg3[j],
            lamDispMm: 2e3 * Double.pi / kHalf,
            subharmonic: tr[j] < 0,
            gcAnalyticAtKc: gcAnalytic(kg3[j], h, wD, variant),
            liouvilleMaxDev: liou
        )
    }
}

/// Tiny context carrier so the RK4 closures stay readable (index into the
/// k-grid + the current time).
private struct Double_i {
    let idx: Int
    let t: Double
}
