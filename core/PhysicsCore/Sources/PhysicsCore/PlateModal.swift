// PlateModal.swift — circular free-edge plate, forced linear response.
//
// PUT-IN     : the pack's modal table (m,n,f_hz,R(r/a)) and the pack's forced-
//              response formula (linear modal superposition with damping ζ).
// EMERGED    : which modes a given drive point excites (via R_mn at the drive
//              radius) and the resulting relative field — the sand pattern is
//              the zero set of this field, chosen by the computation, not drawn.
// CLAIM-TIER : measured (within the c001 known-physics model; eigenvalues are
//              two-route validated in the oracle at 1.58e-10).
// FLOORS     : linear small-deflection; ζ is a PLACEHOLDER (0.002, unmeasured);
//              field is display-normalized (absolute amplitude needs modal-mass
//              calibration, later); no air loading.
//
// Rendering-honesty note: `forcedResponseField` returns the physical field only.
// Callers render its zero crossings; they must not ease or exaggerate it.
import Foundation

public struct PlateModal: Sendable {
    public let modes: [PlateMode]
    /// Damping ratio ζ used in the resonance denominator (pack placeholder 0.002).
    public let zeta: Double

    public init(pack: PhysicsPack) {
        self.modes = pack.plateModalShapes.modes
        self.zeta = pack.plateModalShapes.zeta.value
    }

    public init(modes: [PlateMode], zeta: Double) {
        self.modes = modes
        self.zeta = zeta
    }

    /// Natural frequencies of all modes, ascending.
    public var eigenfrequenciesHz: [Double] {
        modes.map(\.fHz).sorted()
    }

    /// A center point drive excites only m==0 (ring) modes — the selection rule.
    /// Returned frequencies are rounded to 0.1 Hz to match the pack's detent list.
    public var centerDriveReachableHz: [Double] {
        modes.filter { $0.m == 0 }
            .map { ($0.fHz * 10).rounded() / 10 }
            .sorted()
    }

    /// Linear interpolation of a radial table on the 0..1 (r/a) grid.
    private func interpR(_ R: [Double], _ rOverA: Double) -> Double {
        let ngrid = R.count
        let x = max(0.0, min(1.0, rOverA)) * Double(ngrid - 1)
        let i = min(Int(x.rounded(.down)), ngrid - 2)
        let frac = x - Double(i)
        return R[i] * (1 - frac) + R[i + 1] * frac
    }

    /// Forced steady-state field along a ray at angle `theta`, for a point drive
    /// at (r0/a, theta0) oscillating at `fHz`. Pack formula:
    ///   amp_mn = R_mn(r0/a) / sqrt((f_mn²−f²)² + (2ζ f_mn f)²)
    ///   field(r,θ) = Σ amp_mn · R_mn(r/a) · cos(m(θ−θ0))
    public func forcedResponseField(fHz f: Double,
                                    driveROverA r0: Double = 0.0,
                                    driveTheta theta0: Double = 0.0,
                                    theta: Double = 0.0) -> [Double] {
        guard let ngrid = modes.first?.R.count else { return [] }
        var field = [Double](repeating: 0.0, count: ngrid)
        for md in modes {
            let fmn = md.fHz
            let Rr0 = interpR(md.R, r0)
            let denom = ((fmn * fmn - f * f) * (fmn * fmn - f * f)
                         + (2 * zeta * fmn * f) * (2 * zeta * fmn * f)).squareRoot()
            let amp = Rr0 / denom
            let cs = cos(Double(md.m) * (theta - theta0))
            let w = amp * cs
            for j in 0..<ngrid {
                field[j] += w * md.R[j]
            }
        }
        return field
    }
}
