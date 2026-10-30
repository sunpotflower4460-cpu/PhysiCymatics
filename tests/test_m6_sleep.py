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


def test_residual_high_mode_no_missed_detach(basis):
    # Codex counterexample: 70 Hz score (undriven), residual 649 Hz mode with
    # q = 9.4/omega^2, seated particle at center. Sleep must not miss detach.
    def mk(sleep):
        score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=70.0,
                                      forcePeak_N=0.0, attack_s=0.0, release_s=0.0)])
        plate = pt.TransientPlate(basis, score)
        i649 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 1)
        j = next(i for i, d in enumerate(plate.dofs) if d["mode"] == i649)
        plate.q[j] = 9.4/plate.omega[j]**2
        st, _ = sb.initial_state(basis["modes"][0], 1,
                                   init_xy=(np.array([0.0]), np.array([0.0])))
        sim = pt.TransientSandSim(plate, st, sleep_enabled=sleep)
        plate.begin_step(0.0, 1e-9)
        h = float(plate.surf_at(np.array([0.0]), np.array([0.0]), 0.0)[0][0])
        st.z[0] = h; st.sitting[0] = True
        sim.t = plate.t
        return sim, st
    off, st_off = mk(False); on, st_on = mk(True)
    off.advance_to(0.05); on.advance_to(0.05)
    assert off.n_detach > 0
    assert on.n_detach == off.n_detach
    assert abs(float(st_on.z[0]) - float(st_off.z[0])) < 1e-6


def test_no_sitting_sleep_with_rolling(basis):
    # c_roll>0 seated particles move horizontally: sitting sleep is disabled.
    def run_roll(sleep):
        score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=152.0,
                                      forcePeak_N=0.276, attack_s=0.0, release_s=0.0)])
        plate = pt.TransientPlate(basis, score)
        st, _ = sb.initial_state(basis["modes"][0], 20, seed=5)
        params = sb.SandParams(c_roll=0.5, c_rr=0.02)
        sim = pt.TransientSandSim(plate, st, params, sleep_enabled=sleep)
        sim.advance_to(0.3)
        return sim, st
    off, st_off = run_roll(False)
    on, st_on = run_roll(True)
    sitting_sleeping = int((on.sleeping & st_on.sitting).sum())
    assert sitting_sleeping == 0
    assert np.max(np.abs(st_on.x - st_off.x)) < 1e-6
    assert np.max(np.abs(st_on.y - st_off.y)) < 1e-6


def test_event_cap_zero_in_normal_run(basis):
    on, _ = _run(basis, True)
    assert on.n_event_cap == 0


def test_wake_all_invalidates(basis):
    on, st = _run(basis, True, t_end=0.2, force_scale=0.02)
    on.wake_all()
    assert not on.sleeping.any()
    on.advance_to(0.3)
    assert np.all(np.isfinite(st.x))
