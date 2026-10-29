"""M6-2 sleep parity (Python): sleep on/off must reproduce the same physics."""
import numpy as np
import pytest

from reference import plate_transient as pt
from reference import sand_bounce as sb


@pytest.fixture(scope="module")
def basis():
    return pt.build_basis()


def _run(basis, sleep, freq=152.0, t_end=0.3, n=30, force_scale=1.0):
    probe = pt.TransientPlate(basis, pt.ScoreSignal(
        [dict(start_s=0.0, end_s=1.0, frequency_hz=152.0, forcePeak_N=1.0,
              attack_s=0.0, release_s=0.0)]))
    i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
    F = float(probe.force_for_amplitude(i152, 860e-6)) * force_scale
    score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=freq,
                                  forcePeak_N=F, attack_s=0.05, release_s=0.0)])
    plate = pt.TransientPlate(basis, score)
    st, _ = sb.initial_state(basis["modes"][0], n, seed=7)
    sim = pt.TransientSandSim(plate, st, sleep_enabled=sleep)
    sim.advance_to(t_end)
    return sim, st


def test_sleep_parity_strong(basis):
    off, st_off = _run(basis, False)
    on, st_on = _run(basis, True)
    assert on.skipped_updates > 0
    assert abs(on.t - off.t) < 1e-12
    assert on.n_bounce == off.n_bounce
    assert on.n_detach == off.n_detach
    assert st_on.n_alive == st_off.n_alive
    assert np.max(np.abs(st_on.x - st_off.x)) < 1e-6
    assert np.max(np.abs(st_on.y - st_off.y)) < 1e-6
    assert np.array_equal(st_on.sitting, st_off.sitting)


def test_sleep_parity_weak_and_undriven(basis):
    for scale in (0.02, 0.0):
        off, st_off = _run(basis, False, force_scale=scale)
        on, st_on = _run(basis, True, force_scale=scale)
        assert np.max(np.abs(st_on.x - st_off.x)) < 1e-6
        assert st_on.n_alive == st_off.n_alive
        assert on.n_bounce == off.n_bounce


def test_wake_all_invalidates(basis):
    on, st = _run(basis, True, t_end=0.2, force_scale=0.02)
    on.wake_all()
    assert not on.sleeping.any()
    on.advance_to(0.3)
    assert np.all(np.isfinite(st.x))
