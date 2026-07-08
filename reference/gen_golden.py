"""reference/gen_golden.py — generate golden JSON for the Swift PhysicsCore tests.

The Swift PhysicsCore package must reproduce these values from the SAME inputs
(physics_pack + the ported physics). This script is the single source of truth:
it reads the committed physics_pack and the pinned Python reference, and writes
the golden vectors the Swift XCTest suite compares against (CLAUDE.md Phase 2:
"golden tests: Swift出力↔Pythonオラクル JSON照合").

Run: python reference/gen_golden.py
Writes: core/PhysicsCore/Tests/PhysicsCoreTests/Resources/{golden_plate,golden_faraday}.json
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import faraday_linear as fl  # noqa: E402

PACK = json.load(open(os.path.join(
    ROOT, "physics_pack", "PhysiCymatics_physics_pack_v1.json")))
OUT_DIR = os.path.join(ROOT, "core", "PhysicsCore",
                       "Tests", "PhysicsCoreTests", "Resources")


def forced_response_field(modes, f_hz, zeta, r0=0.0, theta0=0.0, theta=0.0):
    """Pack formula: amp_mn = R_mn(r0)/sqrt((f_mn²−f²)²+(2ζ f_mn f)²);
    field(r,θ) = Σ amp_mn R_mn(r) cos(m(θ−θ0)). r0, r are r/a in [0,1]."""
    ngrid = len(modes[0]["R"])
    field = [0.0] * ngrid
    for md in modes:
        R = md["R"]
        fmn = md["f_hz"]
        # R at the drive radius r0 (linear interp on the 0..1 grid)
        x = r0 * (ngrid - 1)
        i = min(int(math.floor(x)), ngrid - 2)
        frac = x - i
        R_r0 = R[i] * (1 - frac) + R[i + 1] * frac
        denom = math.sqrt((fmn**2 - f_hz**2)**2 + (2 * zeta * fmn * f_hz)**2)
        amp = R_r0 / denom
        cs = math.cos(md["m"] * (theta - theta0))
        for j in range(ngrid):
            field[j] += amp * R[j] * cs
    return field


def gen_plate():
    pms = PACK["plate_modal_shapes"]
    modes = pms["modes"]
    zeta = pms["zeta"]["value"]                       # 0.002 placeholder
    det = PACK["plate_detents"]["circular_free_steel_24cm"]
    eig = sorted(m["f_hz"] for m in modes)
    reach = det["center_drive_reachable_hz"]
    idx = [0, 32, 64, 96, 128]                        # sampled radius indices

    def sample(f):
        fld = forced_response_field(modes, f, zeta, r0=0.0, theta0=0.0, theta=0.0)
        return [fld[i] for i in idx]

    golden = dict(
        pack_schema=PACK["schema"],
        n_modes=len(modes),
        zeta=zeta,
        eigenfrequencies_hz=eig,
        center_drive_reachable_hz=reach,
        forced_response=dict(
            note="center drive (r0=0, theta0=0), field sampled along theta=0",
            radius_index=idx,
            at_152hz=sample(151.99),                  # on the first ring resonance
            at_300hz=sample(300.0),                   # off resonance
        ),
    )
    return golden


def gen_faraday():
    r = fl.sweep_f(60.0, 0.020, "bulk")
    return dict(
        f_drive_hz=60.0, h_m=0.020, variant="bulk",
        kc=r["kc"], Gc=r["Gc"], a_c_ms2=r["a_c_ms2"],
        lam_mm=r["lam_mm"], lam_disp_mm=r["lam_disp_mm"],
        subharmonic=r["subharmonic"],
        Gc_analytic_at_kc=r["Gc_analytic_at_kc"],
        liouville_maxdev=r["liouville_maxdev"],
        tol_rel=1e-3,
    )


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    plate = gen_plate()
    faraday = gen_faraday()
    with open(os.path.join(OUT_DIR, "golden_plate.json"), "w") as f:
        json.dump(plate, f, indent=1)
    with open(os.path.join(OUT_DIR, "golden_faraday.json"), "w") as f:
        json.dump(faraday, f, indent=1)
    # copy the pack next to the golden so the Swift test bundle can load it
    with open(os.path.join(OUT_DIR, "physics_pack.json"), "w") as f:
        json.dump(PACK, f, ensure_ascii=False)
    print("wrote golden_plate.json:")
    print("  eig[:5] =", plate["eigenfrequencies_hz"][:5])
    print("  reach   =", plate["center_drive_reachable_hz"])
    print("  152Hz field sample =", [round(v, 4) for v in plate["forced_response"]["at_152hz"]])
    print("wrote golden_faraday.json:")
    print("  kc=%.4f Gc=%.6f a_c=%.5f lam=%.4f sub=%s" % (
        faraday["kc"], faraday["Gc"], faraday["a_c_ms2"],
        faraday["lam_mm"], faraday["subharmonic"]))


if __name__ == "__main__":
    main()
