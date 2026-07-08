// FaradayLinearTests — golden comparison: Swift FaradayLinear ↔ Python JSON.
// Golden produced by reference/gen_golden.py from reference/faraday_linear.py.
// Acceptance (CLAUDE.md oracle map): FaradayLinear vs sweep/ext JSON within <0.1%.
import XCTest
@testable import PhysicsCore

final class FaradayLinearTests: XCTestCase {

    struct GoldenFaraday: Codable {
        let f_drive_hz: Double
        let h_m: Double
        let kc: Double
        let Gc: Double
        let a_c_ms2: Double
        let lam_mm: Double
        let lam_disp_mm: Double
        let subharmonic: Bool
        let Gc_analytic_at_kc: Double
        let liouville_maxdev: Double
        let tol_rel: Double
    }

    func golden() throws -> GoldenFaraday {
        let url = try XCTUnwrap(Bundle.module.url(forResource: "golden_faraday", withExtension: "json"))
        return try JSONDecoder().decode(GoldenFaraday.self, from: Data(contentsOf: url))
    }

    func testThresholdMatchesGoldenWithin0p1pct() throws {
        let g = try golden()
        let f = FaradayLinear()
        let r = f.sweepF(g.f_drive_hz, g.h_m, .bulk)

        func rel(_ a: Double, _ b: Double) -> Double { abs(a - b) / abs(b) }
        XCTAssertLessThan(rel(r.kc, g.kc), g.tol_rel, "kc")
        XCTAssertLessThan(rel(r.gc, g.Gc), g.tol_rel, "Gc")
        XCTAssertLessThan(rel(r.aCms2, g.a_c_ms2), g.tol_rel, "a_c")
        XCTAssertLessThan(rel(r.lamMm, g.lam_mm), g.tol_rel, "lambda")
        XCTAssertEqual(r.subharmonic, g.subharmonic)                 // P4: f/2
    }

    func testTwoRoutesAgree() throws {
        // Route 1 (monodromy bisection) vs route 2 (weak-damping analytic).
        let f = FaradayLinear()
        let r = f.sweepF(60.0, 0.020, .bulk)
        XCTAssertLessThan(abs(r.gc - r.gcAnalyticAtKc) / r.gc, 1e-3)
    }

    func testLiouvilleConservation() throws {
        // Integrator truth: det(M) = exp(-2γT).
        let f = FaradayLinear()
        let r = f.sweepF(60.0, 0.020, .bulk)
        XCTAssertLessThan(r.liouvilleMaxDev, 1e-12)
    }

    func testWavelengthTracksHalfFrequencyDispersion() throws {
        // P3: pattern scale = wavelength of the f/2 wave (~9 mm at 60 Hz).
        let f = FaradayLinear()
        let r = f.sweepF(60.0, 0.020, .bulk)
        XCTAssertLessThan(abs(r.lamMm - r.lamDispMm) / r.lamDispMm, 1e-3)
        XCTAssertGreaterThan(r.lamMm, 8.0)
        XCTAssertLessThan(r.lamMm, 9.5)
    }
}
