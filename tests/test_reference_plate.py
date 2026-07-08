"""Phase-1 pinned tests for reference/plate_modal.py (c001 calibration).

These freeze the calibration numbers of the ruler: two independent numerical
routes must agree, the center-drive selection rule must hold, and Chladni's
exponent must emerge near 2. Acceptance line (CLAUDE.md oracle map): Swift
PlateModal eigenfrequencies vs c001 within rel <1e-6 — mirrored here for the
Python reference. Fast (~1 s)."""
import os, sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "reference"))
import plate_modal as pm  # noqa: E402


@pytest.fixture(scope="module")
def modes():
    return pm.mode_table(nu=0.30)


def test_two_routes_agree(modes):
    # Bessel-determinant route vs Rayleigh-Ritz energy route.
    max_err = max(d["rel_err"] for d in modes)
    assert max_err < 1e-6            # oracle recorded 1.58e-10


def test_center_drive_reaches_152hz(modes):
    reach = pm.center_drive_reachable_hz(modes)
    # P1: first center-reachable ring ≈ 152 Hz (single ring, radius 8.2 cm).
    assert any(abs(f - 152.0) < 1.0 for f in reach)
    assert abs(reach[0] - 152.0) < 1.0


def test_90hz_not_center_reachable(modes):
    # P1: 90.5 Hz is an (m=2) mode — a CENTER drive cannot excite it.
    reach = pm.center_drive_reachable_hz(modes)
    assert not any(abs(f - 90.5) < 1.0 for f in reach)


def test_center_drive_is_rings_only(modes):
    # Selection rule: a center point force only excites m=0 (ring) modes.
    E, rho, h, a, nu = (pm.STEEL[k] for k in ("E_pa", "rho", "h_m", "a_m", "nu"))
    import numpy as np
    D = E * h**3 / (12 * (1 - nu**2)); c0 = np.sqrt(D / (rho * h))
    reach = set(pm.center_drive_reachable_hz(modes))
    for d in modes:
        f = round(float(d["lam2"] / (2 * np.pi * a * a) * c0), 1)
        if f in reach and f <= 5000:
            assert d["m"] == 0


def test_chladni_exponent_emerges_near_two(modes):
    p, r2 = pm.chladni_exponent(modes)
    assert 2.0 <= p <= 2.3          # measured 2.155 (empirical law ≈ 2)
    assert r2 > 0.9


def test_topology_discovered_matches_labels(modes):
    # Measured nodal-diameter count must equal the angular index m for all modes.
    assert all(d["n_diameters_measured"] == d["m"] for d in modes)
