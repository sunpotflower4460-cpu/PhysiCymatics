// swift-tools-version:5.9
// PhysicsCore — PhysiCymatics Phase 2 physics kernel (Swift 6-ready).
// Reads the committed physics_pack and reproduces the Python oracles bit-for-bit
// within tolerance. Golden tests compare Swift output ↔ Python-generated JSON.
import PackageDescription

let package = Package(
    name: "PhysicsCore",
    platforms: [.macOS(.v12), .iOS(.v15)],
    products: [
        .library(name: "PhysicsCore", targets: ["PhysicsCore"]),
    ],
    targets: [
        .target(
            name: "PhysicsCore",
            swiftSettings: [.enableUpcomingFeature("StrictConcurrency")]
        ),
        .testTarget(
            name: "PhysicsCoreTests",
            dependencies: ["PhysicsCore"],
            resources: [
                .copy("Resources/physics_pack.json"),
                .copy("Resources/golden_plate.json"),
                .copy("Resources/golden_faraday.json"),
            ]
        ),
    ]
)
