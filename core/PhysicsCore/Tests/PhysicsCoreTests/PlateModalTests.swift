// PlateModalTests — golden comparison: Swift PlateModal ↔ Python/pack JSON.
// Golden produced by reference/gen_golden.py. Acceptance (CLAUDE.md oracle map):
// PlateModal eigenfrequencies vs c001 within rel <1e-6.
import XCTest
@testable import PhysicsCore

final class PlateModalTests: XCTestCase {

    struct GoldenPlate: Codable {
        let pack_schema: String
        let n_modes: Int
        let zeta: Double
        let eigenfrequencies_hz: [Double]
        let center_drive_reachable_hz: [Double]
        let forced_response: Forced
        struct Forced: Codable {
            let radius_index: [Int]
            let at_152hz: [Double]
            let at_300hz: [Double]
        }
    }

    func loadJSON<T: Decodable>(_ name: String, _ type: T.Type) throws -> T {
        let url = try XCTUnwrap(Bundle.module.url(forResource: name, withExtension: "json"))
        return try JSONDecoder().decode(T.self, from: Data(contentsOf: url))
    }

    func pack() throws -> PhysicsPack {
        let url = try XCTUnwrap(Bundle.module.url(forResource: "physics_pack", withExtension: "json"))
        return try PhysicsPack.load(from: url)
    }

    func testEigenfrequenciesMatchGolden() throws {
        let g = try loadJSON("golden_plate", GoldenPlate.self)
        let plate = PlateModal(pack: try pack())
        let eig = plate.eigenfrequenciesHz
        XCTAssertEqual(eig.count, g.n_modes)
        XCTAssertEqual(eig.count, g.eigenfrequencies_hz.count)
        for (got, want) in zip(eig, g.eigenfrequencies_hz) {
            XCTAssertEqual(got, want, accuracy: 1e-6 * max(1.0, abs(want)))
        }
    }

    func testCenterDriveReachableMatchesGolden() throws {
        let g = try loadJSON("golden_plate", GoldenPlate.self)
        let plate = PlateModal(pack: try pack())
        XCTAssertEqual(plate.centerDriveReachableHz, g.center_drive_reachable_hz)
        // P1: first center-reachable ring ≈ 152 Hz; 90.5 Hz (m=2) is NOT reachable.
        XCTAssertEqual(plate.centerDriveReachableHz.first ?? 0, 152.0, accuracy: 1.0)
        XCTAssertFalse(plate.centerDriveReachableHz.contains { abs($0 - 90.5) < 1.0 })
    }

    func testForcedResponseFieldMatchesGolden() throws {
        let g = try loadJSON("golden_plate", GoldenPlate.self)
        let plate = PlateModal(pack: try pack())
        for (freq, want) in [(151.99, g.forced_response.at_152hz),
                             (300.0, g.forced_response.at_300hz)] {
            let field = plate.forcedResponseField(fHz: freq, driveROverA: 0.0,
                                                  driveTheta: 0.0, theta: 0.0)
            for (k, idx) in g.forced_response.radius_index.enumerated() {
                XCTAssertEqual(field[idx], want[k],
                               accuracy: 1e-6 * max(1e-6, abs(want[k])),
                               "field mismatch at \(freq)Hz idx \(idx)")
            }
        }
    }
}
