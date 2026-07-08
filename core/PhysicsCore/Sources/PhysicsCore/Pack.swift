// Pack.swift — typed loader for the PhysiCymatics physics_pack.
//
// PUT-IN     : the committed physics_pack JSON (板のドレミ・水の定数・モード形).
// EMERGED    : nothing — this is a faithful reader, no physics invented here.
// CLAIM-TIER : established (values are the pack's; provenance lives in the pack).
// FLOORS     : decodes only the fields PhysicsCore consumes; the pack carries
//              more (rendering_oracle, sand, pattern_shape) added in later phases.
//
// Rendering-honesty note (CLAUDE.md §7): this layer only READS. It never
// rewrites, eases, or exaggerates a pack value.
import Foundation

/// One circular-plate modal shape: angular index m, radial index n, its natural
/// frequency, and the normalized radial table R(r/a) on a 0..1 grid (max|R|=1).
public struct PlateMode: Codable, Sendable {
    public let m: Int
    public let n: Int
    public let fHz: Double
    public let R: [Double]

    enum CodingKeys: String, CodingKey {
        case m, n, R
        case fHz = "f_hz"
    }
}

public struct PlateModalShapes: Codable, Sendable {
    public let modes: [PlateMode]
    public let zeta: Zeta

    public struct Zeta: Codable, Sendable {
        public let value: Double
        public let status: String
    }
}

public struct PlateDetents: Codable, Sendable {
    public let circle: CircularFreeSteel

    enum CodingKeys: String, CodingKey {
        case circle = "circular_free_steel_24cm"
    }

    public struct CircularFreeSteel: Codable, Sendable {
        public let eigenfrequenciesHz: [Double]
        public let centerDriveReachableHz: [Double]

        enum CodingKeys: String, CodingKey {
            case eigenfrequenciesHz = "eigenfrequencies_hz"
            case centerDriveReachableHz = "center_drive_reachable_hz"
        }
    }
}

public struct FluidWater: Codable, Sendable {
    public let constants: Constants

    public struct Constants: Codable, Sendable {
        public let rhoKgm3: Double
        public let sigmaNm: Double
        public let nuM2s: Double

        enum CodingKeys: String, CodingKey {
            case rhoKgm3 = "rho_kgm3"
            case sigmaNm = "sigma_nm"
            case nuM2s = "nu_m2s"
        }
    }
}

/// The subset of physics_pack that PhysicsCore consumes.
public struct PhysicsPack: Codable, Sendable {
    public let schema: String
    public let plateModalShapes: PlateModalShapes
    public let plateDetents: PlateDetents
    public let fluidWater: FluidWater

    enum CodingKeys: String, CodingKey {
        case schema
        case plateModalShapes = "plate_modal_shapes"
        case plateDetents = "plate_detents"
        case fluidWater = "fluid_water"
    }

    public static func load(from url: URL) throws -> PhysicsPack {
        let data = try Data(contentsOf: url)
        return try JSONDecoder().decode(PhysicsPack.self, from: data)
    }

    public static func load(from data: Data) throws -> PhysicsPack {
        try JSONDecoder().decode(PhysicsPack.self, from: data)
    }
}
