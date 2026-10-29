"""scripts/gen_sand_fixtures.py — generate shared sand fixtures for port checks.

PUT-IN     : fixed initial arrays (seed 7, 50 particles), c002 oracle simulate.
EMERGED    : nothing — this only records oracle output for later comparison.
CLAIM-TIER : measured (numeric record of the c002 model run).
FLOORS     : fixture is a short (0.3 s) small-N record; not a formation claim.

Usage: python scripts/gen_sand_fixtures.py
Writes tests/fixtures/sand/m1_short.npz + .json (meta + sha256 of arrays).
"""
import hashlib, importlib.util, json, os, subprocess, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "tests", "fixtures", "sand")


def load_oracle():
    spec = importlib.util.spec_from_file_location(
        "c002b", os.path.join(ROOT, "oracles", "c002_sand", "c002b.py"))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def main():
    c002b = load_oracle()
    rng = np.random.default_rng(7)
    Np, a = 50, 0.12
    rr = a*np.sqrt(rng.random(Np)); th = 2*np.pi*rng.random(Np)
    ix, iy = rr*np.cos(th), rr*np.sin(th)
    mode = c002b.make_mode(0, 0)
    r = c002b.simulate(mode, Gamma=80, e=0.75, beta=0.8, Np=Np, t_max=0.3,
                       seed=999, init_xy=(ix, iy))
    os.makedirs(OUT, exist_ok=True)
    npz_path = os.path.join(OUT, "m1_short.npz")
    np.savez(npz_path, init_x=ix, init_y=iy, x=r["x"], y=r["y"],
             absW=r["absW"], sitting=r["sitting"],
             r_node=r["r_node"], A=r["A"], n_bounce=r["n_bounce"],
             hist_json=json.dumps([{k: (float(v) if isinstance(v, (float, np.floating)) else v)
                                     for k, v in h.items()} for h in r["hist"]]),
             params_json=json.dumps(r["params"]))
    h = hashlib.sha256()
    for k in ("init_x","init_y","x","y"):
        h.update(np.load(npz_path)[k].tobytes())
    try:
        commit = subprocess.check_output(["git","rev-parse","HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    meta = dict(schema="physicymatics.sand_fixture.v1", source_commit=commit,
                generator="scripts/gen_sand_fixtures.py", oracle="oracles/c002_sand/c002b.py",
                mode=dict(m=0, root_idx=0, f_hz=float(mode["f_hz"]), a_m=float(a),
                          ngrid=int(mode["ngrid"])),
                Gamma=80, e=0.75, beta=0.8, spp=44, v_stick_frac=0.3, bisect_iters=14,
                Np=Np, t_max_s=0.3, init_seed=7, arrays_sha256=h.hexdigest())
    with open(os.path.join(OUT, "m1_short.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(f"wrote {npz_path}\n  f={mode['f_hz']:.4f} Hz alive={len(r['x'])}/{Np} bounce={r['n_bounce']}")


if __name__ == "__main__":
    main()
