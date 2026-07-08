"""Phase-4 NIGHTLY skeleton: c005 ensemble convergence frame (slow).

Off the fast push gate; wired to nightly + manual. Confirms the ensemble ENTRY
is defined (seed bundle × frequency sweep) and the representative capillary run
(440 Hz) reproduces a square lattice from the committed flagship output. Full
convergence (all seeds/freqs to ~1400 periods) is the nightly heavy job.
Oracle: oracles/c005_shape/c005.py."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "reference"))
import c005_ensemble as ens  # noqa: E402

pytestmark = pytest.mark.slow


def test_ensemble_entry_defined():
    configs = list(ens.iter_ensemble())
    assert len(configs) == len(ens.SEEDS) * len(ens.FREQ_SWEEP_HZ)
    assert 440.0 in ens.FREQ_SWEEP_HZ            # capillary square anchor present
    assert 40.0 in ens.FREQ_SWEEP_HZ             # mixed regime present
    assert all("seed" in c and "f_drive_hz" in c for c in configs)


def test_capillary_440hz_is_square():
    rep = ens.capillary_representative()
    assert ens.is_square(rep)
    assert abs(rep["top2_pair_angle"] - 90.0) < 5.0
    assert rep["temporal_peak_over_drive"] == 0.5   # f/2 preserved


def test_mixed_regime_left_square_is_frontier():
    # The mixed regime (40 Hz) leaves the square region as Chen–Viñals predicts;
    # exact hexagon/8-fold ID is the unconverged nightly target (frontier).
    import json
    mixed = json.load(open(os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "oracles", "c005_shape", "c005_mixed_result.json")))
    # leading angular pair sits well away from 90° (not a square lattice).
    lo, hi, _share = mixed["pair_bins_top"][0]
    assert not (lo <= 90.0 <= hi)
    assert abs((lo + hi) / 2.0 - 90.0) > 5.0
