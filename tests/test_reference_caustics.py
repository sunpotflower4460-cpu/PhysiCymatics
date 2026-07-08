"""Phase-1 pinned tests for reference/caustics.py (c006 ray optics).

Freezes the optical calibration: the near-axial ray crossing depth must match
the paraxial anchor f = nR/(n−1) to ~0%, and the spherical aberration must be
NEGATIVE (marginal rays cross late — the soft deep glow). Acceptance line
(CLAUDE.md): Metal caustics vs profiles; paraxial crossing depth 0.000%. Fast."""
import os, sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "reference"))
import caustics as cs  # noqa: E402

A20 = 20e-6
A200 = 200e-6


def test_paraxial_focus_20um():
    # Oracle c006: A=20µm → f_analytic ≈ 41.066 cm.
    assert abs(cs.paraxial_focus_m(A20) * 100 - 41.066) < 0.01


def test_paraxial_focus_200um():
    # Oracle c006: A=200µm → f_analytic ≈ 4.107 cm.
    assert abs(cs.paraxial_focus_m(A200) * 100 - 4.107) < 0.01


def test_near_axial_crossing_matches_paraxial_20um():
    # THE calibration line: near-axial ray crossing vs nR/(n−1) = 0.000%.
    assert cs.paraxial_deviation_pct(A20) < 0.01


def test_near_axial_crossing_matches_paraxial_200um():
    assert cs.paraxial_deviation_pct(A200) < 0.01


@pytest.mark.parametrize("A", [A20, A200])
def test_spherical_aberration_is_negative(A):
    # Recorded fact: sinusoidal lens has NEGATIVE spherical aberration —
    # marginal rays cross LATER (deeper) than paraxial.
    assert cs.spherical_aberration_sign(A) == 1
    assert cs.marginal_focus_m(A, 0.4) > cs.paraxial_focus_m(A)
