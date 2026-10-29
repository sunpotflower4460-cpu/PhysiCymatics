"""M3 transient plate tests: analytic steady state, decay, symmetry, continuity."""
import numpy as np
import pytest

from reference import plate_transient as pt
from reference import sand_bounce as sb


@pytest.fixture(scope="module")
def basis():
    return pt.build_basis()


def _plate(basis, freq, force, t_end=30.0, drive=(0.0, 0.0)):
    score = pt.ScoreSignal([dict(start_s=0.0, end_s=t_end, frequency_hz=freq,
                                  forcePeak_N=force, attack_s=0.0, release_s=0.0)])
    return pt.TransientPlate(basis, score, drive_xy=drive)


def test_basis_masses_and_dofs(basis):
    plate_mass = pt.RHO*pt.H_PLATE*np.pi*basis["a"]**2
    for md in basis["modes"]:
        assert 0 < md["M_kg"] < plate_mass
    n_modes = len(basis["modes"])
    n_mpos = sum(1 for md in basis["modes"] if md["m"] > 0)
    assert len(basis["dofs"]) == n_modes + n_mpos  # cos/sin both kept
    assert max(md["f_hz"] for md in basis["modes"]) <= 1300.0


def test_steady_state_amplitude_phase(basis):
    md0 = basis["modes"][0]  # 90.5 Hz is m=2; pick the m=0 152 Hz mode instead
    i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0)
    md = basis["modes"][i152]
    f, w = md["f_hz"], 2*np.pi*md["f_hz"]
    # Calibrate force for 100 um on resonance, then drive OFF resonance and
    # compare against the analytic single-DOF steady state for that DOF.
    probe = _plate(basis, f, 1.0)
    F = probe.force_for_amplitude(i152, 100e-6)
    f_drive = f*1.05
    pl = _plate(basis, f_drive, F, t_end=15.0)
    j = next(i for i, d in enumerate(pl.dofs) if d["mode"] == i152)
    wd = 2*np.pi*f_drive
    amp_ref = abs(pl.phi_drive[j])*F/pl.mass[j] / np.sqrt((w*w - wd*wd)**2 + (2*pl.zeta*w*wd)**2)
    # Integrate 6 s (>> decay time 1/(zeta*w) ~ 0.52 s), measure q amplitude.
    qs = []
    dt = 1.0/(44*f_drive)
    while pl.t < 6.0:
        pl.begin_step(pl.t, pl.t + dt)
        if pl.t > 5.0:
            qs.append(pl.q[j])
    amp = (max(qs) - min(qs))/2
    assert abs(amp - amp_ref)/amp_ref < 0.02


def test_free_decay(basis):
    pl = _plate(basis, 152.0, 0.0, t_end=5.0)
    pl.q[:] = 0.0; pl.q[0] = 1e-4
    e0 = float(np.sum(pl.qd**2 + (pl.omega*pl.q)**2))
    dt = 1.0/(44*152.0)
    while pl.t < 1.0:
        pl.begin_step(pl.t, pl.t + dt)
    e1 = float(np.sum(pl.qd**2 + (pl.omega*pl.q)**2))
    assert 0 < e1 < e0*np.exp(-2*pl.zeta*pl.omega[0]*1.0)*1.5


def test_center_drive_symmetry(basis):
    pl = _plate(basis, 152.0, 0.05, t_end=2.0, drive=(0.0, 0.0))
    dt = 1.0/(44*152.0)
    while pl.t < 0.3:
        pl.begin_step(pl.t, pl.t + dt)
    for j, d in enumerate(pl.dofs):
        if basis["modes"][d["mode"]]["m"] > 0:
            assert pl.q[j] == 0.0  # R_m>0(0)=0: never excited from center


def test_freq_switch_continuity(basis):
    ev = [dict(start_s=0.0, end_s=0.4, frequency_hz=152.0, forcePeak_N=0.05, attack_s=0.05, release_s=0.05),
          dict(start_s=0.4, end_s=0.8, frequency_hz=649.0, forcePeak_N=0.05, attack_s=0.05, release_s=0.05)]
    score = pt.ScoreSignal(ev)
    assert abs(score.u(0.4 - 1e-7) - score.u(0.4 + 1e-7)) < 1e-3  # phase integral
    pl = pt.TransientPlate(basis, score)
    mode0 = dict(m=0, a=basis["a"])
    st, _ = sb.initial_state({**basis["modes"][0]}, 30, seed=2)
    sim = pt.TransientSandSim(pl, st)
    sim.advance_to(0.39)
    x_before, q_before = st.x.copy(), pl.q.copy()
    sim.advance_to(0.41)
    assert np.all(np.isfinite(st.x)) and np.all(np.isfinite(pl.q))
    assert np.max(np.abs(st.x - x_before)) < 0.01  # no teleport across switch
    assert st.n_alive + st.n_lost == st.n_initial == 30


def test_batched_surface_parity(basis):
    # Stacked/einsum surf_at must match the per-DOF _phi_grad + Hermite loop.
    pl = _plate(basis, 152.0, 0.05, t_end=2.0, drive=(0.05, 0.01))
    dt = 1.0/(44*152.0)
    pl.begin_step(0.0, dt)
    pl.begin_step(dt, 2*dt)
    rng = np.random.default_rng(0)
    x = rng.uniform(-0.1, 0.1, 37); y = rng.uniform(-0.1, 0.1, 37)
    taus = dt + rng.random(37)*dt
    got = pl.surf_at(x, y, taus)
    ref = [np.zeros(37) for _ in range(5)]
    for k in range(37):
        q, qd, qdd = pl._hermite(float(taus[k]))
        for j, d in enumerate(pl.dofs):
            W, Wx, Wy = pt._phi_grad(pl.modes[d["mode"]], d["kind"],
                                     np.array([x[k]]), np.array([y[k]]))
            ref[0][k] += q[j]*W[0]; ref[1][k] += q[j]*Wx[0]; ref[2][k] += q[j]*Wy[0]
            ref[3][k] += qd[j]*W[0]; ref[4][k] += qdd[j]*W[0]
    for g_, r_ in zip(got, ref):
        assert np.max(np.abs(g_ - r_)) < 1e-12
    # Scalar-tau path too.
    got_s = pl.surf_at(x, y, 1.5*dt)
    q, qd, qdd = pl._hermite(1.5*dt)
    ref_s = [np.zeros(37) for _ in range(5)]
    for j, d in enumerate(pl.dofs):
        W, Wx, Wy = pt._phi_grad(pl.modes[d["mode"]], d["kind"], x, y)
        ref_s[0] += q[j]*W; ref_s[1] += q[j]*Wx; ref_s[2] += q[j]*Wy
        ref_s[3] += qd[j]*W; ref_s[4] += qdd[j]*W
    for g_, r_ in zip(got_s, ref_s):
        assert np.max(np.abs(g_ - r_)) < 1e-12


def test_surface_spatial_derivative(basis):
    pl = _plate(basis, 152.0, 0.05, t_end=2.0, drive=(0.05, 0.01))
    dt = 1.0/(44*152.0)
    pl.begin_step(0.0, dt)
    pl.begin_step(dt, 2*dt)
    x0, y0, tau, eps = 0.04, 0.02, 1.5*dt, 1e-5
    h, hx, hy, _, _ = pl.surf_at(np.array([x0]), np.array([y0]), tau)
    hx_num = (pl.surf_at(np.array([x0+eps]), np.array([y0]), tau)[0]
              - pl.surf_at(np.array([x0-eps]), np.array([y0]), tau)[0])/(2*eps)
    # Table lookup is piecewise-constant in r: allow O(table step) tolerance.
    assert abs(float(hx[0]) - float(hx_num[0])) < 0.05 + 0.2*abs(float(hx[0]))
