"""Phase-1 pinned tests for reference/faraday_linear.py (c003 threshold).

Freezes the Faraday calibration at 60 Hz / 20 mm water: the measured threshold,
the subharmonic response discovered from the multiplier sign, the Liouville
integrator check, and route-1↔route-2 agreement. Acceptance line (CLAUDE.md):
FaradayLinear vs sweep/ext JSON within <0.1%. Fast (~4 s)."""
import os, sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "reference"))
import faraday_linear as fl  # noqa: E402


@pytest.fixture(scope="module")
def sweep():
    return fl.sweep_f(60.0, 0.020, "bulk")


def test_threshold_kc(sweep):
    assert abs(sweep["kc"] - 729.8922293497279) / 729.8922293497279 < 1e-3


def test_threshold_acceleration_about_1ms2(sweep):
    # P2: clean water at 60 Hz waves at ~1 m/s². Measured a_c ≈ 1.105 m/s².
    assert abs(sweep["a_c_ms2"] - 1.1051504873992941) / 1.1051504873992941 < 1e-3
    assert 0.9 < sweep["a_c_ms2"] < 1.3


def test_response_is_subharmonic(sweep):
    # P4: response at exactly half the drive — discovered from trace < 0.
    assert sweep["subharmonic"] is True


def test_liouville_conservation(sweep):
    # Integrator truth: det(M) = exp(-2γT). Oracle acceptance < 1e-12.
    assert sweep["liouville_maxdev"] < 1e-12


def test_two_routes_agree(sweep):
    # Route 1 (monodromy bisection) vs route 2 (weak-damping analytic).
    dev = abs(sweep["Gc"] - sweep["Gc_analytic_at_kc"]) / sweep["Gc"]
    assert dev < 1e-3               # measured ~1.3e-5


def test_wavelength_tracks_half_frequency_dispersion(sweep):
    # P3: pattern scale = wavelength of the f/2 wave (≈9 mm at 60 Hz).
    dev = abs(sweep["lam_mm"] - sweep["lam_disp_mm"]) / sweep["lam_disp_mm"]
    assert dev < 1e-3
    assert 8.0 < sweep["lam_mm"] < 9.5
