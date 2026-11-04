"""M6-11 research accel mode (Python parity at plate level + sim smoke).

'accel_mode == "ode"' evaluates qdd from the equation of motion applied
to the Hermite q/qd. Default "hermite" is the reference model and must
be untouched. Accuracy claims are anchored to the analytic free-decay
solution, not to TS/Python agreement.
"""
import numpy as np
import pytest

from reference import plate_transient as pt
from reference import sand_bounce as sb


@pytest.fixture(scope="module")
def basis():
    return pt.build_basis()


def _plate(basis, mode, freq=152.0, force=0.0):
    score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=freq,
                                  forcePeak_N=force, attack_s=0.0, release_s=0.0)])
    plate = pt.TransientPlate(basis, score, drive_xy=(0.05, 0.03))
    plate.accel_mode = mode
    rng = np.random.default_rng(222)
    plate.q[:] = rng.uniform(-1, 1, plate.ndof) * 2e-4
    plate.qd[:] = rng.uniform(-1, 1, plate.ndof) * 0.3
    return plate


def _analytic_wtt(plate, q0, qd0, x, y, t):
    w = plate.omega; z = plate.zeta
    wd = w * np.sqrt(1 - z * z)
    c = (qd0 + z * w * q0) / wd
    e = np.exp(-z * w * t)
    q = e * (q0 * np.cos(wd * t) + c * np.sin(wd * t))
    qd = e * (-z * w * (q0 * np.cos(wd * t) + c * np.sin(wd * t))
              + wd * (-q0 * np.sin(wd * t) + c * np.cos(wd * t)))
    qdd = -2 * z * w * qd - w * w * q
    # Surface weights via the same phi evaluation surf_at uses.
    phis = []
    for d in plate.dofs:
        phis.append(float(pt._phi_grad(plate.modes[d["mode"]], d["kind"],
                                       np.array([x]), np.array([y]))[0][0]))
    return float(qdd @ np.array(phis))


def test_ode_wtt_closer_to_analytic(basis):
    errs = {}
    for mode in ("hermite", "ode"):
        plate = _plate(basis, mode)
        q0, qd0 = plate.q.copy(), plate.qd.copy()
        dt = 1.0 / (44 * 152)
        max_e = 0.0
        seg_end = 0.0
        for k in range(1, 41):
            t = k * dt / 4
            while seg_end < t:
                plate.begin_step(seg_end, seg_end + dt)
                seg_end += dt
            for xy in ((0.03, -0.02), (-0.05, 0.04)):
                w = float(plate.surf_at(np.array([xy[0]]), np.array([xy[1]]), t)[4][0])
                max_e = max(max_e, abs(w - _analytic_wtt(plate, q0, qd0, xy[0], xy[1], t)))
        errs[mode] = max_e
    assert errs["hermite"] > 100.0
    assert errs["ode"] < errs["hermite"] / 5.0


def test_ode_lookahead_bounds_implemented_accel(basis):
    plate = _plate(basis, "ode", force=1.0)
    dt = 1.0 / (44 * 152)
    lc = plate.lookahead_ranges(0.0, 0.02, lambda t: dt)
    rng = np.random.default_rng(7)
    seg_end = 0.0
    for t in sorted(rng.uniform(0.0, 0.019, 60)):
        while seg_end < t:
            plate.begin_step(seg_end, seg_end + dt)
            seg_end += dt
        x, y = rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1)
        phis = np.array([float(pt._phi_grad(plate.modes[d["mode"]], d["kind"],
                                            np.array([x]), np.array([y]))[0][0])
                         for d in plate.dofs])
        lb = float(np.minimum(phis * lc["amin"], phis * lc["amax"]).sum())
        ub = float(np.maximum(phis * lc["amin"], phis * lc["amax"]).sum())
        w = float(plate.surf_at(np.array([x]), np.array([y]), float(t))[4][0])
        assert lb - 1e-6 <= w <= ub + 1e-6


def test_advance_to_grid_target_no_degenerate_step(basis):
    # M6-11R regression: advance_to(0.1) at 500 Hz used to create a
    # ~1e-15 s final segment whose Hermite acceleration is garbage,
    # spuriously detaching seated grains. Pure step() driving is the
    # reference.
    def run(driver):
        probe = pt.TransientPlate(basis, pt.ScoreSignal(
            [dict(start_s=0.0, end_s=1.0, frequency_hz=152.0, forcePeak_N=1.0,
                  attack_s=0.0, release_s=0.0)]))
        i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
        F = float(probe.force_for_amplitude(i152, 860e-6))
        score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=500.0,
                                      forcePeak_N=F, attack_s=0.05, release_s=0.0)])
        plate = pt.TransientPlate(basis, score)
        st, _ = sb.initial_state(basis["modes"][0], 30, seed=12345)
        sim = pt.TransientSandSim(plate, st, sleep_enabled=True)
        driver(sim)
        return sim

    def by_steps(sim):
        dt = sim.dt_now()
        for _ in range(round(0.1 / dt)):
            sim.step()

    stepped = run(by_steps)
    adv = run(lambda sim: sim.advance_to(0.1))
    assert adv.t == 0.1
    assert adv.n_detach == stepped.n_detach
    assert adv.n_bounce == stepped.n_bounce


def test_ode_wtt_lower_bound_and_finder(basis):
    plate = _plate(basis, "ode", force=1.0)
    dt = 1.0 / (44 * 152)
    plate.begin_step(0.0, dt)
    rng = np.random.default_rng(3)
    for _ in range(30):
        x, y = rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1)
        ta = float(rng.uniform(0.0, dt * 0.8))
        tb = ta + float(rng.uniform(0.2, 1.0)) * (dt - ta)
        lb = plate.wtt_lower_bound(x, y, ta, tb)
        for s in np.linspace(ta, tb, 21):
            w = float(plate.surf_at(np.array([x]), np.array([y]), float(s))[4][0])
            assert w >= lb - 1e-6
    # Finder agrees with a dense scan on segments where roots exist.
    seg_t = dt
    found = 0
    for _ in range(20):
        plate.begin_step(seg_t, seg_t + dt)
        for xy in ((0.02, 0.01), (-0.04, 0.03)):
            status, got = plate.first_detach_below(xy[0], xy[1], seg_t, seg_t + dt)
            ts = np.linspace(seg_t, seg_t + dt, 2001)
            ws = plate.surf_at(np.full(ts.size, xy[0]), np.full(ts.size, xy[1]), ts)[4]
            below = np.where(ws < -pt.G)[0]
            assert status != "unknown"
            if below.size == 0:
                assert status == "none"
            else:
                found += 1
                assert status == "found"
                assert abs(got - float(ts[below[0]])) < dt / 2000 + 2e-6
        seg_t += dt
        if found >= 3:
            break
    assert found > 0


def _rest_plate(basis, events, phase, mode="ode"):
    plate = pt.TransientPlate(basis, pt.ScoreSignal(events, initial_phase=phase),
                              drive_xy=(0.0, 0.0))
    plate.accel_mode = mode
    return plate


def test_m612r_first_detach_is_first_across_two_short_events(basis):
    # Mirrors the TS counterexample: two 10 N events (1-2 us, 18-19 us)
    # in one nominal step; the first root must win.
    plate = _rest_plate(basis, [
        dict(start_s=1e-6, end_s=2e-6, frequency_hz=152.0, forcePeak_N=10.0,
             attack_s=0.0, release_s=0.0),
        dict(start_s=18e-6, end_s=19e-6, frequency_hz=152.0, forcePeak_N=10.0,
             attack_s=0.0, release_s=0.0),
    ], -np.pi / 2)
    dt = 1.0 / (44 * 152)
    plate.begin_step(0.0, dt)
    w15 = float(plate.surf_at(np.array([0.02]), np.array([0.01]), 1.5e-6)[4][0])
    assert w15 < -pt.G
    status, t = plate.first_detach_below(0.02, 0.01, 0.0, dt)
    assert status == "found"
    assert 1e-6 <= t < 2e-6
    status0, _ = plate.first_detach_below(0.02, 0.01, 0.0, dt, budget_in=0)
    assert status0 == "unknown"


def test_m612r_tiny_segment_uses_ode_form_in_hermite_mode(basis):
    plate = _rest_plate(basis, [
        dict(start_s=10e-9, end_s=20e-9, frequency_hz=152.0, forcePeak_N=10.0,
             attack_s=0.0, release_s=0.0),
    ], np.pi / 2, mode="hermite")
    plate.begin_step(0.0, 9e-8)
    assert plate.ode_form_active()
    w = float(plate.surf_at(np.array([0.02]), np.array([0.01]), 15e-9)[4][0])
    assert abs(w) > 100


def test_default_mode_is_hermite(basis):
    plate = pt.TransientPlate(basis, pt.ScoreSignal(
        [dict(start_s=0.0, end_s=1.0, frequency_hz=152.0, forcePeak_N=0.0,
              attack_s=0.0, release_s=0.0)]))
    assert plate.accel_mode == "hermite"


def test_ode_sim_sleep_parity_and_detach(basis):
    def run(sleep):
        probe = pt.TransientPlate(basis, pt.ScoreSignal(
            [dict(start_s=0.0, end_s=1.0, frequency_hz=152.0, forcePeak_N=1.0,
                  attack_s=0.0, release_s=0.0)]))
        i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
        F = float(probe.force_for_amplitude(i152, 860e-6))
        score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=152.0,
                                      forcePeak_N=F, attack_s=0.05, release_s=0.0)])
        plate = pt.TransientPlate(basis, score)
        plate.accel_mode = "ode"
        st, _ = sb.initial_state(basis["modes"][0], 30, seed=7)
        sim = pt.TransientSandSim(plate, st, sleep_enabled=sleep)
        sim.advance_to(0.3)
        return sim, st

    off, st_off = run(False)
    on, st_on = run(True)
    assert on.n_bounce == off.n_bounce
    assert on.n_detach == off.n_detach
    assert off.n_detach > 0
    assert st_on.n_alive == st_off.n_alive
    assert np.max(np.abs(st_on.x - st_off.x)) < 1e-6
    assert np.all(np.isfinite(st_on.z))
