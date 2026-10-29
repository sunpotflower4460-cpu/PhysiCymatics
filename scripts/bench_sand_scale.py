"""scripts/bench_sand_scale.py — particle-count scaling benchmark (M5 evidence).

Measures the single-mode vectorized backend at fixed physics dt for a short
fixed simulated time, then extrapolates the 600k cost from the MEASURED rate
(not a promise). Usage: python scripts/bench_sand_scale.py [--tmax 0.2]
"""
import argparse, json, os, platform, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from reference import sand_bounce as sb


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tmax", type=float, default=0.2)
    p.add_argument("--counts", default="500,2000,20000,100000")
    p.add_argument("--output", default=None)
    a_ = p.parse_args()
    mode = sb.make_mode(0, 0)
    rows = []
    for n in [int(v) for v in a_.counts.split(",")]:
        surf = sb.SingleModeSurface(mode, 80)
        st, _ = sb.initial_state(mode, n, seed=1)
        sim = sb.SandSim(surf, st)
        steps = int(round(a_.tmax/sim.dt))
        t0 = time.time()
        for _ in range(steps):
            sim.step()
        wall = time.time() - t0
        rate = steps*n/wall
        rows.append(dict(n=n, steps=steps, sim_time_s=a_.tmax, wall_s=round(wall, 3),
                         particle_steps_per_s=int(rate), alive=int(st.n_alive)))
        print(rows[-1])
    # Extrapolate 600k x 36 s at the largest measured rate (best case of the ladder).
    best = max(r["particle_steps_per_s"] for r in rows)
    steps_36 = int(round(36.0/rows[0]["sim_time_s"]*rows[0]["steps"]))
    est_s = 600000*steps_36/best
    out = dict(schema="physicymatics.sand_benchmark.v1",
               device=platform.platform(), python=sys.version.split()[0],
               backend="python-reference/sand_bounce.py float64, dt=T/44 at 152 Hz",
               rows=rows,
               estimate_600k_36s=dict(basis_rate_particle_steps_per_s=best,
                                      estimated_wall_hours=round(est_s/3600, 2),
                                      note="linear extrapolation from measured rate; memory/IO excluded"))
    print(json.dumps(out["estimate_600k_36s"], indent=1))
    if a_.output:
        json.dump(out, open(a_.output, "w"), indent=1)


if __name__ == "__main__":
    main()
