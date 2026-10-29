"""scripts/profile_m6.py — M6-1D cost breakdown (Python transient backend).

Measures separately: plate modal advance, surface evaluation, full sand step.
Usage: python scripts/profile_m6.py [--tmax 0.2] [--counts 500,2000,5000]
"""
import argparse, json, os, platform, sys, time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from reference import plate_transient as pt
from reference import sand_bounce as sb


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tmax", type=float, default=0.2)
    p.add_argument("--counts", default="500,2000,5000")
    p.add_argument("--output", default=None)
    a_ = p.parse_args()
    basis = pt.build_basis()
    score = pt.ScoreSignal([dict(start_s=0.0, end_s=60.0, frequency_hz=152.0,
                                  forcePeak_N=0.276, attack_s=0.0, release_s=0.0)])
    # 1) plate modal advance alone (particle-independent)
    plate = pt.TransientPlate(basis, score)
    dt = 1.0/(44*152.0)
    steps = int(round(a_.tmax/dt))
    t0 = time.time()
    for _ in range(steps):
        plate.begin_step(plate.t, plate.t + dt)
    plate_wall = time.time() - t0
    # 2) surface evaluation throughput at segment end state
    rng = np.random.default_rng(0)
    xs = rng.uniform(-0.11, 0.11, 5000); ys = rng.uniform(-0.11, 0.11, 5000)
    t0 = time.time()
    reps = 20
    for _ in range(reps):
        plate.surf_at(xs, ys, plate.t)
    surf_wall = (time.time() - t0)/reps
    rows = []
    for n in [int(v) for v in a_.counts.split(",")]:
        plate2 = pt.TransientPlate(basis, score)
        st, _ = sb.initial_state(basis["modes"][0], n, seed=1)
        sim = pt.TransientSandSim(plate2, st)
        t0 = time.time()
        sim.advance_to(a_.tmax)
        wall = time.time() - t0
        sit_frac = float(np.mean(st.sitting[st.alive])) if st.n_alive else 0.0
        snap_t0 = time.time()
        xy = np.stack([st.x[st.alive], st.y[st.alive]]).astype(np.float32)
        snap_wall = time.time() - snap_t0
        rows.append(dict(n=n, sim_time_s=a_.tmax, wall_s=round(wall, 3),
                         speed_x=round(a_.tmax/wall, 3), n_bounce=sim.n_bounce,
                         n_detach=sim.n_detach, n_event_cap=sim.n_event_cap,
                         sitting_frac=round(sit_frac, 3), alive=int(st.n_alive),
                         snapshot_copy_ms=round(snap_wall*1000, 3)))
        print(rows[-1], flush=True)
    out = dict(schema="physicymatics.m6_profile.v1", device=platform.platform(),
               plate_modal_advance=dict(steps=steps, wall_s=round(plate_wall, 3),
                                        us_per_step=round(plate_wall/steps*1e6, 1)),
               surface_eval=dict(particles=5000, wall_s_per_call=round(surf_wall, 4),
                                 ns_per_particle=round(surf_wall/5000*1e9, 1)),
               full_step_rows=rows)
    print(json.dumps({k: out[k] for k in ("plate_modal_advance", "surface_eval")}, indent=1))
    if a_.output:
        json.dump(out, open(a_.output, "w"), indent=1)


if __name__ == "__main__":
    main()
