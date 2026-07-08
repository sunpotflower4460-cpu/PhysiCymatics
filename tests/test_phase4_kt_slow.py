"""Phase-4 NIGHTLY skeleton: Kumar–Tuckerman full-viscous frame (slow).

Off the fast push gate; wired to nightly + manual. Confirms the frame is
correctly plumbed (weak-damping baseline reproduces the water oracle at 57 Hz),
the high-viscosity floor is recorded, and the full-viscous solver is a loud
NotImplemented TODO until the physics lands. Oracle: oracles/c003_faraday/c003.py."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "reference"))
import kumar_tuckerman as kt  # noqa: E402

pytestmark = pytest.mark.slow


def test_representative_water_57hz_matches_pack():
    rep = kt.representative_water_57hz()
    assert abs(rep["ours_a_c"] - rep["pack_a_c"]) / rep["pack_a_c"] < 1e-3
    assert abs(rep["ours_lam"] - rep["pack_lam"]) / rep["pack_lam"] < 1e-3


def test_walker_oil_floor_recorded():
    f = kt.walker_oil_floor()
    assert f["wavelength"]["verdict"] == "GREEN"
    assert f["threshold"]["verdict"] == "documented floor"
    lo, hi = f["threshold"]["lit_g"]
    assert lo <= 4.3 and hi >= 4.2                 # literature band
    assert f["threshold"]["ours_g"] > hi           # our weak-damping over-estimates


def test_full_viscous_is_explicit_todo():
    # The nightly frame must fail loudly until the KT solver is implemented.
    with pytest.raises(NotImplementedError):
        kt.full_viscous_threshold(80.0, 0.004, 960.0, 0.0206, 20e-6)
