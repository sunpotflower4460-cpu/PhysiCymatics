"""M6-1A tests: exact target time, residual time after contact, hidden
mid-step contact, mid-step detach, no-penetration, dt convergence."""
import numpy as np
import pytest

from reference import plate_transient as pt
from reference import sand_bounce as sb


@pytest.fixture(scope="module")
def basis():
    return pt.build_basis()


def _sim(basis, freq=152.0, force=0.0, n=1, xy=(0.0, 0.0), spp=44):
    score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=freq,
                                  forcePeak_N=force, attack_s=0.0, release_s=0.0)])
    plate = pt.TransientPlate(basis, score)
    mode0 = basis["modes"][0]
    st, _ = sb.initial_state(mode0, n, init_xy=(np.full(n, xy[0]), np.full(n, xy[1])))
    params = sb.SandParams(spp=spp)
    return pt.TransientSandSim(plate, st, params), st


def test_advance_to_exact(basis):
    sim, _ = _sim(basis)
    sim.advance_to(0.137)
    assert abs(sim.t - 0.137) < 1e-12
    sim.advance_to(0.3)
    assert abs(sim.t - 0.3) < 1e-12
    # Single-mode sim too.
    mode = sb.make_mode(0, 0)
    st, _ = sb.initial_state(mode, 5, seed=1)
    s2 = sb.SandSim(sb.SingleModeSurface(mode, 80), st)
    s2.advance_to(0.137)
    assert abs(s2.t - 0.137) < 1e-12


def test_residual_time_flat_plate(basis):
    # Force=0 -> flat surface h=0. Drop from 2 cm in ONE long step; after the
    # bounce at tc the remaining time must be flown (old code froze at tc).
    sim, st = _sim(basis, force=0.0)
    st.z[0] = 0.02; st.vz[0] = 0.0
    dt = 0.12
    sim.step(dt)
    tc = np.sqrt(2*0.02/pt.G)
    v_hit = pt.G*tc
    v_up = 0.75*v_hit
    rem = dt - tc
    z_ref = v_up*rem - 0.5*pt.G*rem*rem
    vz_ref = v_up - pt.G*rem
    # A second bounce may occur if z_ref path recrosses; here rem keeps z>0 check.
    assert z_ref > 0
    assert abs(float(st.z[0]) - z_ref) < 1e-6
    assert abs(float(st.vz[0]) - vz_ref) < 1e-4
    assert sim.n_bounce >= 1 and abs(sim.t - dt) < 1e-15


def test_multiple_bounces_one_step(basis):
    sim, st = _sim(basis, force=0.0)
    st.z[0] = 0.001; st.vz[0] = 0.0
    sim.step(0.2)  # several bounces on a flat plate within one step
    assert sim.n_bounce >= 2
    assert float(st.z[0]) >= -1e-12


def test_fallback_remainder_never_below_surface(basis, monkeypatch):
    # Codex re-review repro: when the event cap is reached, fallback
    # substeps must keep checking contacts over their remainder instead
    # of ballistically accepting a below-surface endpoint. The cap is
    # lowered here so the fallback is really exercised (with the default
    # 512 rounds this particle sticks inside the main loop).
    monkeypatch.setattr(pt, "MAX_EVENT_ROUNDS", 2)
    sim, st = _sim(basis, force=0.0)
    st.z[0] = 0.001; st.vz[0] = -0.1
    sim.step(1.0)
    assert sim.n_event_cap > 0
    assert float(st.z[0]) >= -1e-12


def test_hidden_midstep_contact(basis):
    # Driven plate: in ONE long step (half period) the surface rises through
    # a hovering particle mid-step although it is above at both ends.
    plate_probe_score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=152.0,
                                             forcePeak_N=1.0, attack_s=0.0, release_s=0.0)])
    plate_probe = pt.TransientPlate(basis, plate_probe_score)
    i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
    F = plate_probe.force_for_amplitude(i152, 860e-6)
    sim, st = _sim(basis, force=float(F))
    # Warm up the PLATE ONLY to near steady state (cheap: no sand steps).
    pl_w = sim.plate
    _dt = 1.0/(44*152.0)
    while pl_w.t < 2.0:
        pl_w.begin_step(pl_w.t, pl_w.t + _dt)
    # Place particle just above the surface at center, zero velocity,
    # then take one step of ~half a period.
    f = 152.0; T = 1.0/f
    pl = sim.plate
    pl.begin_step(pl.t, pl.t + 1e-9)  # segment for surf_at probing
    h0 = float(pl.surf_at(np.array([0.0]), np.array([0.0]), pl.t)[0][0])
    st.x[0] = 0.0; st.y[0] = 0.0
    st.z[0] = h0 + 100e-6; st.vx[0] = st.vy[0] = st.vz[0] = 0.0
    st.sitting[0] = False; st.alive[0] = True
    sim.t = pl.t
    nb0 = sim.n_bounce
    sim.step(T/2)
    assert sim.n_bounce > nb0  # mid-step contact was detected and bounced


def test_midstep_detach(basis):
    plate_probe = pt.TransientPlate(basis, pt.ScoreSignal(
        [dict(start_s=0.0, end_s=60.0, frequency_hz=152.0, forcePeak_N=1.0,
              attack_s=0.0, release_s=0.0)]))
    i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
    F = float(plate_probe.force_for_amplitude(i152, 860e-6))
    sim, st = _sim(basis, force=F)
    pl_w = sim.plate
    _dt = 1.0/(44*152.0)
    while pl_w.t < 2.0:
        pl_w.begin_step(pl_w.t, pl_w.t + _dt)
    # Seat the particle at the center antinode; within one long step the
    # support must go negative and the particle must detach mid-step.
    pl = sim.plate
    pl.begin_step(pl.t, pl.t + 1e-9)
    h0 = float(pl.surf_at(np.array([0.0]), np.array([0.0]), pl.t)[0][0])
    st.x[0] = 0.0; st.y[0] = 0.0; st.z[0] = h0
    st.vx[0] = st.vy[0] = st.vz[0] = 0.0
    st.sitting[0] = True; st.alive[0] = True
    sim.t = pl.t
    nd0 = sim.n_detach
    sim.step(1.0/152.0)  # one full period
    assert sim.n_detach > nd0
    assert sim.last_detach_t is not None and pl.t - 1.0/152.0 - 1e-9 <= sim.last_detach_t <= pl.t + 1e-12


def test_no_penetration_and_convergence(basis):
    plate_probe = pt.TransientPlate(basis, pt.ScoreSignal(
        [dict(start_s=0.0, end_s=60.0, frequency_hz=152.0, forcePeak_N=1.0,
              attack_s=0.0, release_s=0.0)]))
    i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
    F = float(plate_probe.force_for_amplitude(i152, 860e-6))
    finals = {}
    for spp in (44, 88):
        score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=152.0,
                                      forcePeak_N=F, attack_s=0.0, release_s=0.0)])
        plate = pt.TransientPlate(basis, score)
        st, _ = sb.initial_state(basis["modes"][0], 30, seed=4)
        sim = pt.TransientSandSim(plate, st, sb.SandParams(spp=spp))
        worst_pen = 0.0
        while sim.t < 0.3 - 1e-15:
            sim.step()
            am = st.alive
            h,_,_,_,_ = plate.surf_at(st.x[am], st.y[am], sim.t)
            worst_pen = max(worst_pen, float(np.max(h - st.z[am])))
        assert worst_pen < 2e-6  # no unnatural penetration at step ends
        finals[spp] = (float(np.median(np.hypot(st.x[st.alive], st.y[st.alive]))),
                       int(st.n_alive), sim.n_bounce)
    # Halving dt keeps collective stats close (trajectories are chaotic;
    # compare population-level quantities, recorded as the convergence claim).
    assert abs(finals[44][0] - finals[88][0]) < 0.02
    assert abs(finals[44][1] - finals[88][1]) <= 5
