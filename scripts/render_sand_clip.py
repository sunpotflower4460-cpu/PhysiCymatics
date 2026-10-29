"""scripts/render_sand_clip.py — fixed-time sand clip + manifest (M1 CLI).

PUT-IN     : run.json config (physics + render), reference/sand_bounce.py.
EMERGED    : nothing visual is invented — each frame plots the simulated
             particle positions at that frame's exact physics time.
CLAIM-TIER : measured (within the c002 numeric model; see sand_bounce floors).
FLOORS     : M1 = single prescribed mode, silent (audio/score arrive in M3+);
             point markers are a position visualization, not real optics.

Usage:
  python scripts/render_sand_clip.py --config runs/example/run.json --output renders/example
  python scripts/render_sand_clip.py --gamma 80 --np 500 --tmax 24 --output renders/m1_gamma80
Outputs: clip.mp4 (if ffmpeg present), frames/*.png, history.json, manifest.json.
Frame k shows physics time k/fps (frame 0 = initial state). No easing.
"""
import argparse, hashlib, json, os, subprocess, sys, time, wave

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reference import sand_bounce as sb
from reference import plate_transient as pt


def draw_frame(path, st, a, r_node, size=(540, 960)):
    dpi = 100
    fig = plt.figure(figsize=(size[0]/dpi, size[1]/dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor("black"); fig.patch.set_facecolor("black")
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(a*np.cos(th), a*np.sin(th), color="#333333", lw=1)
    if np.isfinite(r_node):  # measurement overlay only (reported in manifest)
        ax.plot(r_node*np.cos(th), r_node*np.sin(th), color="#224422", lw=0.7, ls="--")
    m = st.alive
    ax.scatter(st.x[m], st.y[m], s=3, c="#f5e9c8", marker=".", linewidths=0)
    ax.set_xlim(-a*1.15, a*1.15); ax.set_ylim(-a*2.0, a*2.0)
    ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(path, dpi=dpi); plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--output", required=True)
    p.add_argument("--gamma", type=float, default=80.0)
    p.add_argument("--np", type=int, default=500, dest="npart")
    p.add_argument("--tmax", type=float, default=24.0)
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--e", type=float, default=0.75)
    p.add_argument("--beta", type=float, default=0.8)
    p.add_argument("--backend", choices=["single", "transient"], default="single")
    p.add_argument("--seg", type=float, default=12.0, help="transient: seconds per score event")
    p.add_argument("--drive-off", action="store_true", help="transient: drive at r=0.5a")
    p.add_argument("--width", type=int, default=540)
    p.add_argument("--height", type=int, default=960)
    p.add_argument("--ckpt-every", type=int, default=10, help="checkpoint every N frames (0=off)")
    p.add_argument("--resume", action="store_true", help="resume from output/checkpoint.npz")
    p.add_argument("--stop-after-frames", type=int, default=0, help="stop after N frames (checkpoint testing)")
    a_ = p.parse_args()
    cfg = dict(gamma=a_.gamma, npart=a_.npart, tmax_s=a_.tmax, fps=a_.fps,
               seed=a_.seed, e=a_.e, beta=a_.beta, width=a_.width, height=a_.height,
               backend=a_.backend, seg_s=a_.seg, drive_off=a_.drive_off)
    if a_.config:
        cfg.update(json.load(open(a_.config)))
        cfg = {**cfg}  # config file wins; recorded below as actual values
    os.makedirs(os.path.join(a_.output, "frames"), exist_ok=True)
    wall0 = time.time()
    mode = sb.make_mode(0, 0)
    score = None; audio_sr = None; n_audio = None
    if cfg["backend"] == "transient":
        basis = pt.build_basis()
        drive_xy = (0.06, 0.0) if cfg["drive_off"] else (0.0, 0.0)
        seg = cfg["seg_s"]
        freqs = [152.0, 649.0, 368.6]
        cfg["tmax_s"] = seg*len(freqs)
        probe = pt.TransientPlate(basis, pt.ScoreSignal(
            [dict(start_s=0.0, end_s=1.0, frequency_hz=152.0, forcePeak_N=1.0,
                  attack_s=0.0, release_s=0.0)]), drive_xy=drive_xy)
        i152 = next(i for i, md in enumerate(basis["modes"]) if md["m"] == 0 and md["root_idx"] == 0)
        force_N = probe.force_for_amplitude(i152, 860e-6)
        events = [dict(start_s=k*seg, end_s=(k+1)*seg, frequency_hz=f,
                       forcePeak_N=float(force_N), attack_s=0.1, release_s=0.1)
                  for k, f in enumerate(freqs)]
        score = pt.ScoreSignal(events, duration_s=cfg["tmax_s"])
        plate = pt.TransientPlate(basis, score, drive_xy=drive_xy)
        st, rng_state = sb.initial_state(mode, cfg["npart"], seed=cfg["seed"])
        sim = pt.TransientSandSim(plate, st, sb.SandParams(e=cfg["e"], beta=cfg["beta"]))
        r_node = float("nan")  # multi-mode: no single measurement circle drawn
        get_stats = lambda: dict(t=round(float(sim.t), 2), alive=int(st.n_alive),
                                 freq_hz=score.freq_at(sim.t))
    else:
        surf = sb.SingleModeSurface(mode, cfg["gamma"])
        st, rng_state = sb.initial_state(mode, cfg["npart"], seed=cfg["seed"])
        sim = sb.SandSim(surf, st, sb.SandParams(e=cfg["e"], beta=cfg["beta"]))
        r_node = mode["r_nodes"][0]*mode["a"] if mode["r_nodes"] else float("nan")
        get_stats = sim.stats
    init_hash = hashlib.sha256(st.x.tobytes() + st.y.tobytes()).hexdigest()
    # Frame k at t=k/fps for k=0..N-1, N=duration*fps (video length == duration).
    n_frames = int(round(cfg["tmax_s"]*cfg["fps"]))
    cfg_hash = hashlib.sha256(json.dumps(
        {k: cfg[k] for k in sorted(cfg) if k not in ()}, sort_keys=True, default=str).encode()).hexdigest()
    ckpt_path = os.path.join(a_.output, "checkpoint.npz")
    start_frame = 0
    hist = [get_stats()]

    def save_ckpt(next_frame):
        data = dict(x=st.x, y=st.y, z=st.z, vx=st.vx, vy=st.vy, vz=st.vz,
                    alive=st.alive, sitting=st.sitting,
                    sim_t=np.float64(sim.t), n_bounce=np.int64(sim.n_bounce),
                    next_frame=np.int64(next_frame), cfg_hash=cfg_hash,
                    hist_json=json.dumps(hist))
        if cfg["backend"] == "transient":
            data.update(q=sim.plate.q, qd=sim.plate.qd, plate_t=np.float64(sim.plate.t))
        else:
            data.update(step_index=np.int64(sim.step_index))
        np.savez(ckpt_path, **data)

    if a_.resume and os.path.exists(ckpt_path):
        ck = np.load(ckpt_path, allow_pickle=False)
        if str(ck["cfg_hash"]) != cfg_hash:
            raise SystemExit("checkpoint config mismatch; refusing to resume")
        for k in ("x", "y", "z", "vx", "vy", "vz", "alive", "sitting"):
            getattr(st, k)[:] = ck[k]
        sim.t = float(ck["sim_t"]); sim.n_bounce = int(ck["n_bounce"])
        if cfg["backend"] == "transient":
            sim.plate.q[:] = ck["q"]; sim.plate.qd[:] = ck["qd"]
            sim.plate.t = float(ck["plate_t"]); sim.plate._seg = None
        else:
            sim.step_index = int(ck["step_index"])
        start_frame = int(ck["next_frame"])
        hist = json.loads(str(ck["hist_json"]))
        print(f"resumed from frame {start_frame} at t={sim.t:.3f}s")
    if start_frame == 0:
        draw_frame(os.path.join(a_.output, "frames", "f00000.png"), st, mode["a"], r_node,
                   (cfg["width"], cfg["height"]))
        start_frame = 1
        if a_.ckpt_every:
            save_ckpt(1)
    for k in range(start_frame, n_frames):
        sim.advance_to(k/cfg["fps"])
        draw_frame(os.path.join(a_.output, "frames", f"f{k:05d}.png"), st, mode["a"], r_node,
                   (cfg["width"], cfg["height"]))
        if k % max(1, cfg["fps"]//2) == 0:
            hist.append(get_stats())
        if a_.ckpt_every and (k + 1) % a_.ckpt_every == 0:
            save_ckpt(k + 1)
        if a_.stop_after_frames and (k + 1) >= a_.stop_after_frames:
            save_ckpt(k + 1)
            print(f"stopped after {k+1} frames; checkpoint at {ckpt_path}")
            return
    hist.append(get_stats())
    if score is not None:
        audio_sr = 48000
        n_audio = int(round(cfg["tmax_s"]*audio_sr))
        pcm = np.array([0.15*score.u(k/audio_sr) for k in range(n_audio)])
        with wave.open(os.path.join(a_.output, "audio.wav"), "wb") as w:
            w.setnchannels(1); w.setsampwidth(2); w.setframerate(audio_sr)
            w.writeframes((np.clip(pcm, -1, 1)*32767).astype(np.int16).tobytes())
    clip = None
    if subprocess.run(["which", "ffmpeg"], capture_output=True).returncode == 0:
        clip = "clip.mp4"
        cmd = ["ffmpeg", "-y", "-framerate", str(cfg["fps"]),
               "-i", os.path.join(a_.output, "frames", "f%05d.png")]
        if score is not None:
            cmd += ["-i", os.path.join(a_.output, "audio.wav"), "-c:a", "aac"]
        cmd += ["-c:v", "libx264", "-pix_fmt", "yuv420p", os.path.join(a_.output, clip)]
        subprocess.run(cmd, check=True, capture_output=True)
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    is_tr = cfg["backend"] == "transient"
    manifest = dict(
        schema="physicymatics.sand_run.v1", source_commit=commit,
        implementation_backend=("python-reference/plate_transient.py (float64)" if is_tr
                                else "python-reference/sand_bounce.py (float64)"),
        model_name=("multi-mode transient plate + c002 contact" if is_tr
                    else "c002 single-mode bounce (prescribed sinusoid)"),
        claim_scope="measured within the numeric model; not a real-sand measurement",
        model_floors=["no air", "no grain-grain interaction", "no back-reaction on plate",
                      "e/beta/v_stick assumed"] + (["modes truncated at 1300 Hz",
                      "zeta=0.002 placeholder", "point-force drive idealization"] if is_tr
                      else ["single mode (0,1)", "silent clip (M1)"]),
        plate=dict(shape="circle", radius_m=float(mode["a"]), thickness_m=0.001, material="steel"),
        mode=dict(m=0, root_idx=0, f_hz=float(mode["f_hz"]), ngrid=int(mode["ngrid"]),
                  normalization="max|R|=1"),
        drive=(dict(score_events=[{k: e[k] for k in ("start_s", "end_s", "frequency_hz", "forcePeak_N")}
                                   for e in events],
                    score_hash=hashlib.sha256(json.dumps(events, sort_keys=True).encode()).hexdigest(),
                    drive_xy=list(drive_xy), audio="PCM from same score (gain 0.15)")
               if is_tr else
               dict(Gamma=cfg["gamma"], surfaceAccelerationPeak_ms2=cfg["gamma"]*sb.G,
                    amplitude_m=float(surf.A), audio="none (M1 silent)")),
        initial_state_hash=init_hash,
        final_state_hash=hashlib.sha256(st.x.tobytes() + st.y.tobytes() + st.z.tobytes()).hexdigest(),
        checkpoint_metadata=dict(every_frames=a_.ckpt_every, file="checkpoint.npz", resumed=bool(a_.resume)),
        rng_name="numpy PCG64 (default_rng)", seed=cfg["seed"],
        particle_count_initial=int(st.n_initial), particle_count_alive_final=int(st.n_alive),
        particle_count_rendered=int(st.n_alive), particle_interactions_enabled=False,
        physics_dt=("1/(44*f_drive) per event" if is_tr else float(sim.dt)),
        steps_per_period=44, surface_resolution=int(mode["ngrid"]),
        collision_solver=dict(e=cfg["e"], beta=cfg["beta"], v_stick_frac=0.3, bisect_iters=14),
        floating_point_precision="float64", simulated_duration_s=cfg["tmax_s"],
        wall_time_s=round(time.time()-wall0, 2),
        render_size=[cfg["width"], cfg["height"]], video_fps=cfg["fps"], frame_count=n_frames,
        frame_time_rule="frame k = physics time k/fps, k=0..duration*fps-1; frame 0 = initial state",
        measurement_overlay="dashed green circle = c001 node radius (measurement only, not dynamics)",
        codec_and_container=("H.264+AAC MP4 via ffmpeg" if is_tr else "H.264 MP4 via ffmpeg") if clip else "PNG sequence only",
        audio_sample_rate=audio_sr, audio_sample_count=n_audio,
        reproduction_command=" ".join(sys.argv),
        validation="see tests/ (parity + transient) and history.json")
    json.dump(manifest, open(os.path.join(a_.output, "manifest.json"), "w"), indent=1)
    json.dump(hist, open(os.path.join(a_.output, "history.json"), "w"), indent=1)
    print(f"done: frames={n_frames} alive={st.n_alive}/{st.n_initial} "
          f"bounce={sim.n_bounce} wall={manifest['wall_time_s']}s clip={clip}")
    print(f"final stats: {get_stats()}")


if __name__ == "__main__":
    main()
