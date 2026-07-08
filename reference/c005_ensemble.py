"""reference/c005_ensemble.py — Phase-4 skeleton: c005 ensemble convergence.

PUT-IN     : the Zhang–Viñals direct-PDE flagship (oracles/c005_shape/c005.py),
             run as an ENSEMBLE — a bundle of random seeds across a frequency
             sweep — so the selected pattern angle is a converged statistic, not
             a single-seed snapshot.
WHY        : the capillary-regime square (440 Hz) is confirmed single-seed and
             Chen–Viñals cross-validated; the mixed regime (40 Hz) is still
             converging (frontier). The ensemble is what turns "snapshot" into
             "measured with error bars".
CLAIM-TIER : frontier (entry only). This module defines the seed bundle and the
             frequency sweep — the ENTRY to the nightly job — and a representative
             check that the committed flagship output at 440 Hz is a square.
FLOORS     : running the full flagship to convergence (~1400+ drive periods per
             config) is the nightly heavy job; it is NOT run here. Only the entry
             (config list) and the pinned representative result are exercised.

Oracle: oracles/c005_shape/c005.py (flagship) + c005_result.json (440 Hz square).
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The ensemble ENTRY: a seed bundle over a frequency sweep spanning the
# capillary (high f), mixed (mid f) and toward gravity (low f) regimes.
SEEDS = (7, 11, 13, 17, 23)
FREQ_SWEEP_HZ = (40.0, 110.0, 220.0, 440.0)


def iter_ensemble():
    """Yield every (seed, f_drive) config the nightly job would run.
    Purely the entry — no simulation is launched here."""
    for f in FREQ_SWEEP_HZ:
        for s in SEEDS:
            yield dict(seed=s, f_drive_hz=f)


def capillary_representative():
    """Representative 1-shot: the committed flagship at 440 Hz selected a square
    lattice (top-2 pair angle ≈90°, f/2 preserved). Read from the oracle so the
    result is pinned, not re-derived."""
    return json.load(open(os.path.join(
        ROOT, "oracles", "c005_shape", "c005_result.json")))


def is_square(result, angle_tol_deg=5.0, share_min=0.5):
    """A capillary square: dominant angular pair near 90° with majority share."""
    return (abs(result["top2_pair_angle"] - 90.0) <= angle_tol_deg
            and result["top2_share"] >= share_min)


if __name__ == "__main__":
    configs = list(iter_ensemble())
    print(f"ensemble entry: {len(configs)} configs "
          f"({len(SEEDS)} seeds × {len(FREQ_SWEEP_HZ)} freqs)")
    rep = capillary_representative()
    print(f"440Hz representative: top2_pair={rep['top2_pair_angle']:.1f}° "
          f"share={rep['top2_share']:.3f} temporal={rep['temporal_peak_over_drive']} "
          f"-> square={is_square(rep)}")
