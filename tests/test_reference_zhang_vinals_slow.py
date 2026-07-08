"""Phase-1 NIGHTLY skeleton for reference/zhang_vinals.py (c005 linear onset).

Marked `slow`: excluded from fast CI (`pytest -m "not slow"`), wired for the
nightly heavy-oracle run. This is the route-2 self-validation of the PDE solver
— the LINEAR onset of the Zhang–Viñals equations, whose growth rate must rise
monotonically with drive amplitude toward the Faraday threshold. The full
nonlinear pattern-selection flagship (which angle survives) stays in
oracles/c005_shape/c005.py and is added to nightly CI in Phase 4.

Growth values are deterministic (seed=7) and pinned to the oracle's `linear`
mode output. Runtime ~30 s."""
import os, sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "reference"))
import zhang_vinals as zv  # noqa: E402

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def growths():
    P = zv.setup(440.0)
    fc = zv.leading_threshold_scaled(P)
    return P, [zv.linear_growth_rate(P, fac * fc) for fac in (0.9, 1.0, 1.1)]


def test_setup_constants(growths):
    P, _ = growths
    assert abs(P["gam"] - 0.00634) < 1e-4
    assert abs(P["Gam0"] - 0.2462) < 1e-3
    assert abs(P["kc"] - 2954.8) < 1.0


def test_growth_monotone_in_drive(growths):
    _, g = growths
    # Stronger shaking → less decay / more growth, everywhere.
    assert g[0] < g[1] < g[2]


def test_growth_values_pinned(growths):
    _, g = growths
    for got, want in zip(g, (-0.00216, -0.00140, -0.00065)):
        assert abs(got - want) < 5e-5


def test_below_threshold_decays(growths):
    _, g = growths
    # Well below onset the ring band must decay.
    assert g[0] < 0
