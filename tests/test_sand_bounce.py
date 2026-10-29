"""M1 sand tests: unit mechanics + c002 oracle parity + fixtures.

Fast tests only (short, small-N, deterministic). Long formation runs are
manual/slow and recorded in docs/sand_implementation_report.md.
"""
import hashlib, importlib.util, inspect, json, os

import numpy as np
import pytest

from reference import sand_bounce as sb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(ROOT, "tests", "fixtures", "sand")


@pytest.fixture(scope="module")
def mode():
    return sb.make_mode(0, 0)


def _oracle():
    spec = importlib.util.spec_from_file_location(
        "c002b", os.path.join(ROOT, "oracles", "c002_sand", "c002b.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


# ---------- single-body mechanics ----------
def test_free_flight_ballistic(mode):
    # High above the surface: no contact, z must follow z0 - g t^2 / 2.
    surf = sb.SingleModeSurface(mode, Gamma=80)
    st, _ = sb.initial_state(mode, 1, init_xy=(np.array([0.01]), np.array([0.0])))
    st.z[:] = 0.05
    sim = sb.SandSim(surf, st)
    for _ in range(10):
        sim.step()
    t = sim.t
    assert abs(float(st.z[0]) - (0.05 - 0.5*sb.G*t*t)) < 1e-12
    assert abs(float(st.vz[0]) - (-sb.G*t)) < 1e-12


def test_flat_surface_restitution_no_energy_gain(mode):
    # Gamma=0 -> flat still surface; drop from 2 cm, normal bounce must
    # return with speed <= e * impact speed (small-bisect tolerance).
    surf = sb.SingleModeSurface(mode, Gamma=0.0)
    assert surf.A == 0.0
    st, _ = sb.initial_state(mode, 1, init_xy=(np.array([0.02]), np.array([0.0])))
    st.z[:] = 0.02
    sim = sb.SandSim(surf, st)
    vz_before = None
    for _ in range(2000):
        prev_vz, prev_b = float(st.vz[0]), sim.n_bounce
        sim.step()
        if sim.n_bounce > prev_b:
            vz_before = prev_vz - sb.G*sim.dt  # approx pre-impact (within 1 step)
            break
    assert vz_before is not None and vz_before < 0
    assert 0 < float(st.vz[0]) <= 0.75*abs(vz_before) + 0.05


def test_surface_sample_consistency(mode):
    surf = sb.SingleModeSurface(mode, Gamma=80)
    assert abs(surf.f_hz - 152.0) < 0.5
    s = surf.sample(0.03, 0.01, 0.001)
    assert bool(s["inside"]) and not bool(surf.sample(0.2, 0, 0)["inside"])
    # wtt = -omega^2 * w for the single sinusoid.
    assert abs(float(s["wtt_ms2"]) + surf.omega**2*float(s["w_m"])) < 1e-12
    # signed shape: center and rim have opposite signs for mode (0,1).
    Wc, _, _ = surf.W_grad(np.array([0.0]), np.array([0.0]))
    Wr, _, _ = surf.W_grad(np.array([0.119]), np.array([0.0]))
    assert Wc[0]*Wr[0] < 0


def test_no_target_in_dynamics():
    src = inspect.getsource(sb.SandSim.step)
    assert "r_node" not in src and "r_nodes" not in src
    assert "node" not in src.lower()


# ---------- oracle parity (the M1 gate) ----------
def test_oracle_parity_short(mode):
    c002b = _oracle()
    rng = np.random.default_rng(7)
    Np, a = 50, 0.12
    rr = a*np.sqrt(rng.random(Np)); th = 2*np.pi*rng.random(Np)
    ix, iy = rr*np.cos(th), rr*np.sin(th)
    mode_o = c002b.make_mode(0, 0)
    ro = c002b.simulate(mode_o, Gamma=80, e=0.75, beta=0.8, Np=Np,
                        t_max=0.3, seed=999, init_xy=(ix, iy))
    rr_ = sb.simulate(mode, Gamma=80, e=0.75, beta=0.8, Np=Np,
                       t_max=0.3, seed=999, init_xy=(ix, iy))
    assert len(ro["x"]) == len(rr_["x"]) and ro["n_bounce"] == rr_["n_bounce"]
    assert np.max(np.abs(ro["x"] - rr_["x"])) < 1e-12
    assert np.max(np.abs(ro["y"] - rr_["y"])) < 1e-12
    assert np.array_equal(ro["sitting"], rr_["sitting"])


def test_fixture_matches_reference(mode):
    d = np.load(os.path.join(FIX, "m1_short.npz"))
    meta = json.load(open(os.path.join(FIX, "m1_short.json")))
    h = hashlib.sha256()
    for k in ("init_x", "init_y", "x", "y"):
        h.update(d[k].tobytes())
    assert h.hexdigest() == meta["arrays_sha256"]
    r = sb.simulate(mode, Gamma=80, e=0.75, beta=0.8, Np=int(len(d["init_x"])),
                    t_max=0.3, seed=999, init_xy=(d["init_x"], d["init_y"]))
    assert np.max(np.abs(r["x"] - d["x"])) < 1e-12
    assert r["n_bounce"] == int(d["n_bounce"])


def test_counts_and_gamma08_control(mode):
    st, _ = sb.initial_state(mode, 200, seed=3)
    sim = sb.SandSim(sb.SingleModeSurface(mode, Gamma=0.8), st)
    sim.advance_to(0.5)
    snap = sim.snapshot()
    assert snap["n_alive"] + snap["n_lost"] == snap["n_initial"] == 200
    assert np.all(np.isfinite(st.x)) and np.all(np.isfinite(st.z))
    # Under-threshold control: nearly everything trapped (|W|*Gamma<1).
    assert sim.stats()["trapped"] > 0.9
