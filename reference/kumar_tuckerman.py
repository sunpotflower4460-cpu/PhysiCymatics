"""reference/kumar_tuckerman.py — Phase-4 skeleton: full-viscous Faraday threshold.

PUT-IN     : the Kumar & Tuckerman (JFM 279:49, 1994) full-viscous Floquet
             formulation of the Faraday threshold — the exact linear-stability
             problem for a viscous fluid layer, WITHOUT the weak-damping
             (Benjamin–Ursell) approximation used in reference/faraday_linear.py.
WHY        : the weak-damping model is accurate for water (kh≫1, low ν) but
             over-estimates the threshold for high-viscosity fluids. The
             walking-droplet check (20 cSt oil, 80 Hz) reads γ_F 5.3 g vs the
             literature 4.2–4.3 g — an HONEST FLOOR recorded in
             oracles/c003_faraday/c003_walker_check.json. Closing that gap needs
             this full-viscous solver.
CLAIM-TIER : frontier (unimplemented). This module is the FRAME for the nightly
             CI job; the representative check below confirms the weak-damping
             baseline still reproduces the water oracle, so the frame is wired
             correctly before the heavy solver lands.
FLOORS     : `full_viscous_threshold` is not implemented — calling it raises,
             on purpose, so the nightly job fails loudly until the physics is in.

Oracle: oracles/c003_faraday/c003.py (weak-damping) + c003_walker_check.json (floor).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import faraday_linear as fl  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def weak_damping_threshold(f_drive, h, rho=998.0, sigma=0.0728, nu=1.004e-6):
    """Current (weak-damping) baseline for an arbitrary Newtonian fluid.
    This is the model the full-viscous KT version must eventually replace for
    high-viscosity fluids; for water it is already anchor-validated.

    reference/faraday_linear.py carries the fluid constants as module globals,
    so for a non-water fluid we override them for the call and restore after."""
    saved = (fl.RHO, fl.SIG, fl.NU)
    fl.RHO, fl.SIG, fl.NU = rho, sigma, nu
    try:
        return fl.sweep_f(f_drive, h, "bulk")
    finally:
        fl.RHO, fl.SIG, fl.NU = saved


def walker_oil_floor():
    """The documented high-viscosity floor (not a pass) — read straight from
    the committed oracle so the number cannot drift silently."""
    return json.load(open(os.path.join(
        ROOT, "oracles", "c003_faraday", "c003_walker_check.json")))


def full_viscous_threshold(f_drive, h, rho, sigma, nu):
    """Kumar–Tuckerman full-viscous Floquet threshold — NIGHTLY TODO (Phase 4).

    The plan: assemble the viscous Floquet operator (Chebyshev in depth,
    Fourier in time), find the marginal drive amplitude where the largest
    Floquet multiplier crosses |μ|=1, and confirm the 20 cSt / 80 Hz oil
    threshold lands at 4.2–4.3 g (vs the weak-damping 5.3 g floor)."""
    raise NotImplementedError(
        "Kumar-Tuckerman full-viscous solver is the Phase-4 nightly task; "
        "weak-damping baseline lives in weak_damping_threshold(). Floor: "
        "oracles/c003_faraday/c003_walker_check.json")


def representative_water_57hz():
    """Representative 1-shot: the frame reproduces the water oracle at 57 Hz."""
    r = weak_damping_threshold(57.0, 0.020)
    pack = json.load(open(os.path.join(
        ROOT, "physics_pack", "PhysiCymatics_physics_pack_v1.json")))
    row = next(x for x in pack["fluid_water"]["threshold_curve_h20mm"]
               if x["f_hz"] == 57.0)
    return dict(ours_a_c=r["a_c_ms2"], pack_a_c=row["a_c_ms2"],
                ours_lam=r["lam_mm"], pack_lam=row["lam_mm"])


if __name__ == "__main__":
    rep = representative_water_57hz()
    print("representative water 57Hz:", rep)
    floor = walker_oil_floor()
    print("walker floor:", floor["threshold"])
