"""reference/plate_transient.py — multi-mode transient plate + sand (M3 reference).

PUT-IN (linear modal plate, point force drive):
  w(x,y,t) = sum_j q_j(t) Phi_j(x,y)
  qdd_j + 2*zeta*omega_j*qd_j + omega_j^2*q_j = Phi_j(drive) * F(t) / M_j
  M_j = integral(rho*h*Phi_j^2 dA), Phi dimensionless (max|R|=1), q in metres.
  m>0 modes keep BOTH degenerate cos/sin directions as separate DOFs.
  F(t) = forcePeak_N * u(t), u from the shared score (phase = integral of f).
  Sand: same contact mechanics as reference/sand_bounce.py, with the surface
  (w, wt, wtt, wx, wy) evaluated from modal state; inside a sand step the
  modal state is the cubic-Hermite interpolant of the RK4 endpoint states
  (a state evaluation, defined here — not a visual easing).
EMERGED    : nothing is asserted about patterns; responses are computed and
  measured. Weak / unchanged / unformed outcomes are valid results.
CLAIM-TIER : measured (within this numeric model). Port/analytic checks in
  tests/test_plate_transient.py are not real-plate measurements.
FLOORS:
  retained modes truncated at f <= 1300 Hz (truncation error unquantified
  beyond the mode-count comparison in tests); rigid-body plate motion and
  real support mechanics excluded; zeta=0.002 is the pack PLACEHOLDER
  (unmeasured); point-force drive is an idealization (no shaker calibration);
  sand floors inherited from c002 (no air / grain-grain / back-reaction,
  e/beta assumed). Multi-mode stick threshold: v_stick = 0.3 * v_scale with
  v_scale = sqrt(sum_j qd_j^2 + (omega_j q_j)^2) (max|Phi|=1), which reduces
  exactly to c002's 0.3*A*omega for one steady sinusoid.
"""
import numpy as np

from reference import sand_bounce as sb

G = sb.G
RHO, H_PLATE = 7850.0, 1e-3
ZETA = 0.002  # pack PLACEHOLDER (unmeasured)


def build_basis(ngrid=512, f_max=1300.0, m_max=4):
    """Mode basis: radial tables + modal masses. DOFs expand cos/sin for m>0."""
    modes = []
    for m in range(0, m_max + 1):
        for k in range(4):
            md = sb.make_mode(m, k, ngrid=ngrid)
            if md["f_hz"] > f_max:
                break
            integ = np.trapezoid(md["R"]**2 * md["rg"], md["rg"])
            ang = 2*np.pi if m == 0 else np.pi
            M = RHO*H_PLATE * ang * md["a"]**2 * integ
            modes.append(dict(m=m, root_idx=k, f_hz=md["f_hz"], a=md["a"],
                              rg=md["rg"], R=md["R"], dR=md["dR"], M_kg=float(M),
                              r_nodes=md["r_nodes"]))
    modes.sort(key=lambda d: d["f_hz"])
    dofs = []
    for i, md in enumerate(modes):
        dofs.append(dict(mode=i, kind="cos"))
        if md["m"] > 0:
            dofs.append(dict(mode=i, kind="sin"))
    return dict(modes=modes, dofs=dofs, ngrid=ngrid, a=modes[0]["a"])


def _phi_grad(mode, kind, x, y):
    """Phi, dPhi/dx, dPhi/dy (per metre) for one DOF. Same table lookup as c002."""
    a, m, NG1 = mode["a"], mode["m"], len(mode["R"]) - 1
    r = np.hypot(x, y)
    idx = np.minimum((r*(NG1/a)).astype(np.int64), NG1)
    R, dR = mode["R"][idx], mode["dR"][idx]/a
    rs = np.where(r < 1e-12, 1e-12, r)
    c, s_ = x/rs, y/rs
    if m == 0:
        return R, dR*c, dR*s_
    th = np.arctan2(y, x)
    cm, sm = np.cos(m*th), np.sin(m*th)
    if kind == "cos":
        W, Wr, Wt = R*cm, dR*cm, -m*R*sm/rs
    else:
        W, Wr, Wt = R*sm, dR*sm, m*R*cm/rs
    return W, Wr*c - Wt*s_, Wr*s_ + Wt*c


class ScoreSignal:
    """u(t) = envelope(t)*sin(theta(t)), theta = 2*pi*int f dt (piecewise const f)."""

    def __init__(self, events, duration_s=None, initial_phase=0.0):
        self.events = [dict(e) for e in events]
        self.duration_s = duration_s if duration_s is not None else max(e["end_s"] for e in events)
        self._theta0 = []
        th = float(initial_phase)
        for e in self.events:
            self._theta0.append(th)
            th += 2*np.pi*e["frequency_hz"]*(e["end_s"] - e["start_s"])

    def event_at(self, t):
        for i, e in enumerate(self.events):
            if e["start_s"] <= t < e["end_s"]:
                return i, e
        return None, None

    def u(self, t):
        i, e = self.event_at(t)
        if e is None:
            return 0.0
        env = 1.0
        atk, rel = e.get("attack_s", 0.0), e.get("release_s", 0.0)
        if atk > 0 and t < e["start_s"] + atk:
            u_ = (t - e["start_s"])/atk; env = min(env, u_*u_*(3-2*u_))
        if rel > 0 and t > e["end_s"] - rel:
            u_ = (e["end_s"] - t)/rel; env = min(env, u_*u_*(3-2*u_))
        theta = self._theta0[i] + 2*np.pi*e["frequency_hz"]*(t - e["start_s"])
        return env*np.sin(theta)

    def force_peak_at(self, t):
        _, e = self.event_at(t)
        return e["forcePeak_N"] if e else 0.0

    def freq_at(self, t):
        _, e = self.event_at(t)
        return e["frequency_hz"] if e else None


class TransientPlate:
    """Modal state (q, qd) integrated by RK4 substeps; Hermite inside a step."""

    def __init__(self, basis, score, drive_xy=(0.0, 0.0), zeta=ZETA):
        self.basis = basis; self.score = score; self.zeta = zeta
        self.modes = basis["modes"]; self.dofs = basis["dofs"]
        self.ndof = len(self.dofs)
        self.omega = np.array([2*np.pi*self.modes[d["mode"]]["f_hz"] for d in self.dofs])
        self.mass = np.array([self.modes[d["mode"]]["M_kg"] for d in self.dofs])
        self.q = np.zeros(self.ndof); self.qd = np.zeros(self.ndof)
        self.t = 0.0
        # Stacked tables for batched surface evaluation (same values/lookup).
        self._R_stack = np.stack([md["R"] for md in self.modes])
        self._dR_stack = np.stack([md["dR"] for md in self.modes])
        self._dof_mode = np.array([d["mode"] for d in self.dofs])
        self._dof_sin = np.array([d["kind"] == "sin" for d in self.dofs])
        self._mode_m = np.array([md["m"] for md in self.modes])
        self.set_drive(*drive_xy)
        self._seg = None  # (t0, t1, q0, qd0, qdd0, q1, qd1, qdd1)

    def set_drive(self, x, y):
        self.drive_xy = (float(x), float(y))
        self.phi_drive = np.array([
            float(_phi_grad(self.modes[d["mode"]], d["kind"],
                            np.array([x]), np.array([y]))[0][0]) for d in self.dofs])

    def _accel(self, q, qd, t):
        F = self.score.force_peak_at(t) * self.score.u(t)
        return self.phi_drive*F/self.mass - 2*self.zeta*self.omega*qd - self.omega**2*q

    def _rk4(self, q, qd, t, h):
        def f(q_, qd_, t_): return qd_, self._accel(q_, qd_, t_)
        k1q, k1v = f(q, qd, t)
        k2q, k2v = f(q + 0.5*h*k1q, qd + 0.5*h*k1v, t + 0.5*h)
        k3q, k3v = f(q + 0.5*h*k2q, qd + 0.5*h*k2v, t + 0.5*h)
        k4q, k4v = f(q + h*k3q, qd + h*k3v, t + h)
        return (q + h*(k1q + 2*k2q + 2*k3q + k4q)/6,
                qd + h*(k1v + 2*k2v + 2*k3v + k4v)/6)

    def begin_step(self, t0, t1):
        """Advance modal state t0->t1 (substeps: omega_max*h <= 0.3); store segment."""
        assert abs(self.t - t0) < 1e-12
        q0, qd0 = self.q.copy(), self.qd.copy()
        qdd0 = self._accel(q0, qd0, t0)
        h_max = 0.3/self.omega.max()
        n = max(1, int(np.ceil((t1 - t0)/h_max)))
        h = (t1 - t0)/n
        q, qd, t = q0, qd0, t0
        for _ in range(n):
            q, qd = self._rk4(q, qd, t, h); t += h
        self.q, self.qd, self.t = q, qd, t1
        qdd1 = self._accel(q, qd, t1)
        self._seg = (t0, t1, q0, qd0, qdd0, q, qd, qdd1)
        # Conservative surface-height bound over this segment (any x, y, tau):
        # |Phi_j| <= 1 and Hermite |q_j(s)| <= max(|q0|,|q1|) + 0.105*h*(|qd0|+|qd1|).
        hseg = t1 - t0
        self.h_bound = float(np.sum(np.maximum(np.abs(q0), np.abs(q))
                                     + 0.105*hseg*(np.abs(qd0) + np.abs(qd))))

    def _hermite(self, tau):
        t0, t1, q0, qd0, qdd0, q1, qd1, qdd1 = self._seg
        h = t1 - t0; s = (tau - t0)/h
        s2, s3 = s*s, s*s*s
        q = ((2*s3 - 3*s2 + 1)*q0 + (s3 - 2*s2 + s)*h*qd0
             + (-2*s3 + 3*s2)*q1 + (s3 - s2)*h*qd1)
        qd = ((6*s2 - 6*s)*q0 + (3*s2 - 4*s + 1)*h*qd0
              + (-6*s2 + 6*s)*q1 + (3*s2 - 2*s)*h*qd1)/h
        qdd = ((12*s - 6)*q0 + (6*s - 4)*h*qd0
               + (-12*s + 6)*q1 + (6*s - 2)*h*qd1)/(h*h)
        return q, qd, qdd

    def surf_at(self, x, y, tau):
        """(h, hx, hy, uz, wtt) at time tau inside the current segment.

        tau may be a scalar or an array broadcastable with x/y (per-particle
        contact times during bisection). Radial index / angle terms are
        computed once per call and shared across modes.
        """
        x = np.asarray(x, float); y = np.asarray(y, float)
        t0, t1, q0, qd0, _a, q1, qd1, _b = self._seg
        hseg = t1 - t0; s = (np.asarray(tau, float) - t0)/hseg
        s2, s3 = s*s, s*s*s
        c_q0 = 2*s3 - 3*s2 + 1; c_qd0 = (s3 - 2*s2 + s)*hseg
        c_q1 = -2*s3 + 3*s2; c_qd1 = (s3 - s2)*hseg
        d_q0 = (6*s2 - 6*s)/hseg; d_qd0 = 3*s2 - 4*s + 1
        d_q1 = (-6*s2 + 6*s)/hseg; d_qd1 = 3*s2 - 2*s
        a_q0 = (12*s - 6)/(hseg*hseg); a_qd0 = (6*s - 4)/hseg
        a_q1 = (-12*s + 6)/(hseg*hseg); a_qd1 = (6*s - 2)/hseg
        # Batched evaluation in particle chunks: stacked tables + einsum.
        # Identical values/lookup to _phi_grad (parity-tested).
        a_plate = self.basis["a"]; NG1 = self.basis["ngrid"] - 1
        r_all = np.hypot(x, y)
        out = [np.zeros_like(r_all) for _ in range(5)]
        CH = 50000
        for lo in range(0, r_all.size, CH):
            sl = slice(lo, lo + CH)
            r = r_all[sl]
            idx = np.minimum((r*(NG1/a_plate)).astype(np.int64), NG1)
            rs = np.where(r < 1e-12, 1e-12, r)
            c, sn = x[sl]/rs, y[sl]/rs
            cm = [np.ones_like(r), c]; sm = [np.zeros_like(r), sn]
            for m in (2, 3, 4):
                cm.append(cm[m-1]*c - sm[m-1]*sn)
                sm.append(sm[m-1]*c + cm[m-1]*sn)
            Rsel = self._R_stack[:, idx]              # (nmode, n)
            dRsel = self._dR_stack[:, idx]/a_plate
            mm = self._mode_m
            cmM = np.stack([cm[m] for m in mm]); smM = np.stack([sm[m] for m in mm])
            # Per-mode cos/sin shapes and gradients, then gather DOF rows.
            Wc = np.where((mm == 0)[:, None], Rsel, Rsel*cmM)
            Ws = Rsel*smM
            Wrc = np.where((mm == 0)[:, None], dRsel, dRsel*cmM)
            Wrs = dRsel*smM
            Wtc = np.where((mm == 0)[:, None], 0.0, -mm[:, None]*Rsel*smM/rs)
            Wts = np.where((mm == 0)[:, None], 0.0, mm[:, None]*Rsel*cmM/rs)
            def gather(cos_arr, sin_arr):
                g = cos_arr[self._dof_mode]
                g = np.where(self._dof_sin[:, None], sin_arr[self._dof_mode], g)
                return g
            Phi = gather(Wc, Ws)
            Wr = gather(Wrc, Wrs); Wt = gather(Wtc, Wts)
            Phix = Wr*c - Wt*sn; Phiy = Wr*sn + Wt*c
            if np.ndim(s) == 0:
                q = c_q0*q0 + c_qd0*qd0 + c_q1*q1 + c_qd1*qd1
                qd = d_q0*q0 + d_qd0*qd0 + d_q1*q1 + d_qd1*qd1
                qdd = a_q0*q0 + a_qd0*qd0 + a_q1*q1 + a_qd1*qd1
                out[0][sl] = q @ Phi; out[1][sl] = q @ Phix; out[2][sl] = q @ Phiy
                out[3][sl] = qd @ Phi; out[4][sl] = qdd @ Phi
            else:
                ss = s[sl]; s2c, s3c = ss*ss, ss*ss*ss
                cq0 = 2*s3c - 3*s2c + 1; cqd0 = (s3c - 2*s2c + ss)*hseg
                cq1 = -2*s3c + 3*s2c; cqd1 = (s3c - s2c)*hseg
                dq0 = (6*s2c - 6*ss)/hseg; dqd0 = 3*s2c - 4*ss + 1
                dq1 = (-6*s2c + 6*ss)/hseg; dqd1 = 3*s2c - 2*ss
                aq0 = (12*ss - 6)/(hseg*hseg); aqd0 = (6*ss - 4)/hseg
                aq1 = (-12*ss + 6)/(hseg*hseg); aqd1 = (6*ss - 2)/hseg
                q = cq0*q0[:, None] + cqd0*qd0[:, None] + cq1*q1[:, None] + cqd1*qd1[:, None]
                qd = dq0*q0[:, None] + dqd0*qd0[:, None] + dq1*q1[:, None] + dqd1*qd1[:, None]
                qdd = aq0*q0[:, None] + aqd0*qd0[:, None] + aq1*q1[:, None] + aqd1*qd1[:, None]
                out[0][sl] = np.einsum("jn,jn->n", q, Phi)
                out[1][sl] = np.einsum("jn,jn->n", q, Phix)
                out[2][sl] = np.einsum("jn,jn->n", q, Phiy)
                out[3][sl] = np.einsum("jn,jn->n", qd, Phi)
                out[4][sl] = np.einsum("jn,jn->n", qdd, Phi)
        return tuple(out)

    def v_scale(self):
        q, qd = (self._seg[5], self._seg[6]) if self._seg else (self.q, self.qd)
        return float(np.sqrt(np.sum(qd**2 + (self.omega*q)**2)))

    def force_bounds(self):
        f_max = max((e["forcePeak_N"] for e in self.score.events), default=0.0)
        return np.abs(self.phi_drive) * f_max / self.mass

    def energy_per_dof(self):
        q, qd = (self._seg[5], self._seg[6]) if self._seg else (self.q, self.qd)
        return np.hypot(qd, self.omega * q)

    def force_for_amplitude(self, mode_index, target_A):
        """Steady on-resonance force giving surface amplitude target_A (max|Phi|=1)."""
        j = next(i for i, d in enumerate(self.dofs) if d["mode"] == mode_index)
        return target_A * self.mass[j] * 2*self.zeta * self.omega[j]**2 / abs(self.phi_drive[j])


CONTACT_SAMPLES = 8       # gap samples per remaining interval (first-contact search)
EVENT_BISECT_ITERS = 30   # contact-time bisection iterations
MAX_EVENT_ROUNDS = 8      # max event rounds per step before documented fallback
SLEEP_WINDOW_S = 0.5
SLEEP_ACCEL_MARGIN = 0.25   # m/s^2 (see docs/m6_2_event_driven_design.md)
SLEEP_HEIGHT_MARGIN = 2e-5  # m


class TransientSandSim:
    """Sand on the transient surface, event-accurate inside each step (M6-1).

    Per step, each particle carries a local time cursor. Flying particles:
    first contact in (cursor, t1] is searched by sampling the gap
    z_ballistic - h at CONTACT_SAMPLES fractions of the remaining interval,
    then bisecting the first sign-changing subinterval; after a bounce the
    remaining time is processed in further rounds (multiple contacts per
    step are handled). Sitting particles (c_roll=0): wtt at fixed (x,y) is
    LINEAR in time inside a Hermite segment, so the detach time
    (wtt crossing -g) is solved exactly from the two endpoint values and
    flight resumes from that instant.
    Detection floor (honest): a contact whose gap dips below zero and back
    entirely between two adjacent samples (width < remaining/8, with no
    sample <= 0) can be missed. This is NOT complete detection.
    If MAX_EVENT_ROUNDS is exhausted, remaining particles are advanced to
    t1 without further event checks and counted in n_event_cap.
    c_roll>0 sitting keeps the legacy whole-step rolling update (floor).
    """

    def __init__(self, plate, state, params=None, sleep_enabled=True):
        self.plate = plate; self.state = state
        self.p = params or sb.SandParams()
        self.t = 0.0; self.n_bounce = 0
        self.n_detach = 0; self.n_event_cap = 0
        self.last_detach_t = None
        # M6-2 sleep (computation management only; mirrors TS transientSand).
        self.sleep_enabled = sleep_enabled
        n = state.pid.size
        self.sleeping = np.zeros(n, bool)
        self.wake_t = np.zeros(n)
        self.sleep_t0 = np.zeros(n)
        self.eval_cooldown = np.zeros(n, dtype=np.int32)
        self.skipped_updates = 0
        self.sleep_evaluations = 0
        self.wake_count = 0
        self.force_buf = plate.force_bounds()

    def wake_all(self):
        st = self.state
        idx = np.where(self.sleeping & st.alive & ~st.sitting)[0]
        if idx.size:
            d = self.t - self.sleep_t0[idx]
            st.x[idx] += st.vx[idx]*d; st.y[idx] += st.vy[idx]*d
            st.z[idx] += st.vz[idx]*d - 0.5*G*d*d; st.vz[idx] -= G*d
        self.sleeping[:] = False

    def refresh_force_bounds(self):
        self.force_buf = self.plate.force_bounds()

    def materialize_all(self):
        st, pl = self.state, self.plate
        idx = np.where(self.sleeping & st.alive)[0]
        if not idx.size:
            return
        fly = idx[~st.sitting[idx]]; sit = idx[st.sitting[idx]]
        if fly.size:
            d = self.t - self.sleep_t0[fly]
            st.x[fly] += st.vx[fly]*d; st.y[fly] += st.vy[fly]*d
            st.z[fly] += st.vz[fly]*d - 0.5*G*d*d; st.vz[fly] -= G*d
            self.sleep_t0[fly] = self.t
        if sit.size and pl._seg is not None:
            h,_,_,_,_ = pl.surf_at(st.x[sit], st.y[sit], self.t)
            st.z[sit] = h
            self.sleep_t0[sit] = self.t

    def dt_now(self):
        f = self.plate.score.freq_at(self.t) or self.plate.score.freq_at(self.t + 1e-9) or 152.0
        return 1.0/(self.p.spp * max(f, 20.0))

    def step(self, dt=None):
        p, pl, st = self.p, self.plate, self.state
        e, beta = p.e, p.beta
        if dt is None:
            dt = self.dt_now()
        t = self.t; t_new = t + dt
        pl.begin_step(t, t_new)
        v_stick = p.v_stick_frac * pl.v_scale()
        x, y, z, vx, vy, vz = st.x, st.y, st.z, st.vx, st.vy, st.vz
        sitting, alive = st.sitting, st.alive
        cur = np.full(st.pid.size, t)
        if self.sleep_enabled:
            energy = pl.energy_per_dof()
        # Sleeping particles: skip unless due this step; due ones are restored.
        if self.sleeping.any():
            not_due = self.sleeping & alive & (self.wake_t > t_new)
            self.skipped_updates += int(not_due.sum())
            cur[not_due] = t_new
            due = np.where(self.sleeping & alive & (self.wake_t <= t_new))[0]
            if due.size:
                self.sleeping[due] = False
                self.wake_count += int(due.size)
                fly_d = due[~sitting[due]]; sit_d = due[sitting[due]]
                if fly_d.size:
                    d = self.wake_t[fly_d] - self.sleep_t0[fly_d]
                    x[fly_d] += vx[fly_d]*d; y[fly_d] += vy[fly_d]*d
                    z[fly_d] += vz[fly_d]*d - 0.5*G*d*d; vz[fly_d] -= G*d
                if sit_d.size:
                    h,_,_,_,_ = pl.surf_at(x[sit_d], y[sit_d], self.wake_t[sit_d])
                    z[sit_d] = h
                cur[due] = np.maximum(self.wake_t[due], t)

        def ballistic(idx, tau):
            d = tau - cur[idx]
            return (x[idx] + vx[idx]*d, y[idx] + vy[idx]*d,
                    z[idx] + vz[idx]*d - 0.5*G*d*d, vz[idx] - G*d)

        for _round in range(MAX_EVENT_ROUNDS):
            # --- flying ---
            fly = alive & ~sitting & (cur < t_new - 1e-15)
            if fly.any():
                fi_all = np.where(fly)[0]
                # Conservative skip: ballistic minimum above the segment's
                # surface-height bound => contact is impossible this round.
                rem_all = t_new - cur[fi_all]
                z_end = z[fi_all] + vz[fi_all]*rem_all - 0.5*G*rem_all*rem_all
                z_min = np.minimum(z[fi_all], z_end)
                tv = np.clip(vz[fi_all]/G, 0.0, rem_all)
                z_min = np.minimum(z_min, z[fi_all] + vz[fi_all]*tv - 0.5*G*tv*tv)
                safe = z_min > pl.h_bound
                if safe.any():
                    ok = fi_all[safe]
                    xf, yf, zf, vzf = ballistic(ok, np.full(ok.size, t_new))
                    x[ok] = xf; y[ok] = yf; z[ok] = zf; vz[ok] = vzf
                    cur[ok] = t_new
                fi = fi_all[~safe]
                rem = t_new - cur[fi] if fi.size else rem_all[:0]
                gaps = []
                for k in range(1, CONTACT_SAMPLES + 1):
                    tau = cur[fi] + rem*(k/CONTACT_SAMPLES)
                    xs, ys, zs, _ = ballistic(fi, tau)
                    hs,_,_,_,_ = pl.surf_at(xs, ys, tau)
                    gaps.append(zs - hs)
                gaps = np.stack(gaps)  # (S, n)
                hit_any = (gaps <= 0).any(axis=0)
                # no contact: ballistic to t_new
                ok = fi[~hit_any]
                if ok.size:
                    xf, yf, zf, vzf = ballistic(ok, np.full(ok.size, t_new))
                    x[ok] = xf; y[ok] = yf; z[ok] = zf; vz[ok] = vzf
                    cur[ok] = t_new
                if hit_any.any():
                    idx = fi[hit_any]
                    first = np.argmax(gaps[:, hit_any] <= 0, axis=0)  # 0-based sample idx
                    rem_i = t_new - cur[idx]
                    lo = cur[idx] + rem_i*(first/CONTACT_SAMPLES)
                    hi = cur[idx] + rem_i*((first + 1)/CONTACT_SAMPLES)
                    x0 = x[idx].copy(); y0 = y[idx].copy(); z0 = z[idx].copy()
                    vx0 = vx[idx].copy(); vy0 = vy[idx].copy(); vz0 = vz[idx].copy()
                    c0 = cur[idx].copy()
                    for _ in range(EVENT_BISECT_ITERS):
                        tm = 0.5*(lo + hi); d = tm - c0
                        xm = x0 + vx0*d; ym = y0 + vy0*d
                        zm = z0 + vz0*d - 0.5*G*d*d
                        hm,_,_,_,_ = pl.surf_at(xm, ym, tm)
                        below = zm <= hm
                        hi = np.where(below, tm, hi); lo = np.where(below, lo, tm)
                    tc = 0.5*(lo + hi); d = tc - c0
                    xc = x0 + vx0*d; yc = y0 + vy0*d
                    vzc = vz0 - G*d
                    hc, hx, hy, uz, _w = pl.surf_at(xc, yc, tc)
                    nn = np.sqrt(hx*hx + hy*hy + 1.0)
                    nx_, ny_, nz_ = -hx/nn, -hy/nn, 1.0/nn
                    rvx, rvy, rvz = vx0, vy0, vzc - uz
                    vn = rvx*nx_ + rvy*ny_ + rvz*nz_
                    hit = vn < 0
                    rvx2 = rvx - (1+e)*vn*nx_; rvy2 = rvy - (1+e)*vn*ny_; rvz2 = rvz - (1+e)*vn*nz_
                    vn2 = rvx2*nx_ + rvy2*ny_ + rvz2*nz_
                    tx, ty, tz = rvx2-vn2*nx_, rvy2-vn2*ny_, rvz2-vn2*nz_
                    wi = idx[hit]
                    x[wi] = xc[hit]; y[wi] = yc[hit]; z[wi] = hc[hit] + 1e-9
                    vx[wi] = vn2[hit]*nx_[hit] + beta*tx[hit]
                    vy[wi] = vn2[hit]*ny_[hit] + beta*ty[hit]
                    vz[wi] = vn2[hit]*nz_[hit] + beta*tz[hit] + uz[hit]
                    cur[wi] = tc[hit]
                    self.n_bounce += int(hit.sum())
                    stick = hit & (np.abs(vn2) < v_stick)
                    si = idx[stick]; sitting[si] = True
                    if p.c_roll == 0.0:
                        vx[si] = 0.0; vy[si] = 0.0
                    # Non-hit (vn>=0 at tc): separating contact; place at tc state.
                    nh = idx[~hit]
                    if nh.size:
                        x[nh] = xc[~hit]; y[nh] = yc[~hit]
                        z[nh] = hc[~hit] + 1e-9; vz[nh] = vzc[~hit]
                        cur[nh] = tc[~hit]
            # --- sitting ---
            sit = alive & sitting & (cur < t_new - 1e-15)
            if sit.any():
                si = np.where(sit)[0]
                if p.c_roll > 0.0:
                    # Legacy whole-step rolling (documented floor for c_roll>0).
                    hs, hx, hy, uz, wtt = pl.surf_at(x[si], y[si], t_new)
                    z[si] = hs
                    Nacc = np.maximum(G + wtt, 0.0)
                    vx[si] += -p.c_roll*Nacc*hx*dt; vy[si] += -p.c_roll*Nacc*hy*dt
                    sp = np.hypot(vx[si], vy[si]) + 1e-30
                    drag = np.minimum(p.c_rr*Nacc*dt, sp)
                    vx[si] -= drag*vx[si]/sp; vy[si] -= drag*vy[si]/sp
                    x[si] += vx[si]*dt; y[si] += vy[si]*dt
                    det = wtt < -G
                    di = si[det]; sitting[di] = False
                    vz[di] = uz[det]; z[di] = hs[det] + 1e-9
                    cur[si] = t_new
                else:
                    _,_,_,_, w0 = pl.surf_at(x[si], y[si], cur[si])
                    hs1,_,_,_, w1 = pl.surf_at(x[si], y[si], np.full(si.size, t_new))
                    det_now = w0 < -G
                    det_later = (~det_now) & (w1 < -G)
                    stay = ~(det_now | det_later)
                    if stay.any():
                        sj = si[stay]; z[sj] = hs1[stay]; cur[sj] = t_new
                    det = det_now | det_later
                    if det.any():
                        dj = si[det]
                        w0d, w1d = w0[det], w1[det]
                        frac = np.where(w1d != w0d, (-G - w0d)/np.where(w1d != w0d, w1d - w0d, 1.0), 0.0)
                        frac = np.clip(frac, 0.0, 1.0)
                        frac = np.where(det_now[det], 0.0, frac)
                        td = cur[dj] + frac*(t_new - cur[dj])
                        hd,_,_, uzd, _ = pl.surf_at(x[dj], y[dj], td)
                        sitting[dj] = False
                        z[dj] = hd + 1e-9; vz[dj] = uzd
                        cur[dj] = td
                        self.n_detach += int(dj.size)
                        if dj.size:
                            self.last_detach_t = float(td[0])
            if not ((alive & (cur < t_new - 1e-15)).any()):
                break
        # Fallback after the event cap: repeated inelastic micro-bounces
        # accumulate in finite time (inelastic collapse); the physical limit
        # is contact. Seat remaining flying particles on the surface for the
        # rest of the step (documented regularization, counted in n_event_cap)
        # instead of ballistically penetrating the plate.
        left = alive & (cur < t_new - 1e-15)
        if left.any():
            li = np.where(left)[0]
            self.n_event_cap += int(li.size)
            hs,_,_, uzl, _ = pl.surf_at(x[li], y[li], t_new)
            z[li] = hs; vz[li] = uzl
            sitting[li] = True
            if p.c_roll == 0.0:
                vx[li] = 0.0; vy[li] = 0.0
            cur[li] = t_new
        # Sleep evaluation for particles that finished this step awake.
        if self.sleep_enabled:
            self.eval_cooldown = np.maximum(self.eval_cooldown - 1, 0)
            done = np.where(alive & ~self.sleeping & (self.eval_cooldown == 0)
                            & (cur >= t_new - 1e-15))[0]
            if done.size:
                self.sleep_evaluations += int(done.size)
                k = 1 + 2*pl.zeta
                wE = pl.omega * k
                # |Phi_j(x,y)| per DOF for the finished subset.
                abs_phi = np.empty((pl.ndof, done.size))
                for j, d in enumerate(pl.dofs):
                    W, _, _ = _phi_grad(pl.modes[d["mode"]], d["kind"], x[done], y[done])
                    abs_phi[j] = np.abs(W)
                sit_d = done[sitting[done]]; fly_d = done[~sitting[done]]
                pos = {int(v): i for i, v in enumerate(done)}
                if sit_d.size:
                    cols = np.array([pos[int(v)] for v in sit_d])
                    s0 = ((self.force_buf + wE*energy)[:, None] * abs_phi[:, cols]).sum(axis=0)
                    s1 = ((wE*self.force_buf)[:, None] * abs_phi[:, cols]).sum(axis=0)
                    ok = (G - s0) > SLEEP_ACCEL_MARGIN
                    if ok.any():
                        idx = sit_d[ok]
                        delta = np.where(s1[ok] > 0,
                                          np.minimum((G - SLEEP_ACCEL_MARGIN - s0[ok])/np.where(s1[ok] > 0, s1[ok], 1.0), SLEEP_WINDOW_S),
                                          SLEEP_WINDOW_S)
                        good = delta > 8*dt
                        idx = idx[good]
                        self.sleeping[idx] = True
                        self.wake_t[idx] = t_new + delta[good]
                        self.sleep_t0[idx] = t_new
                if fly_d.size:
                    h_base = float((energy/pl.omega).sum())
                    h_rate = float((self.force_buf/pl.omega).sum())
                    z0 = z[fly_d]; vz0 = vz[fly_d]
                    def safe(d):
                        z_end = z0 + vz0*d - 0.5*G*d*d
                        tv = np.clip(vz0/G, 0.0, d)
                        z_min = np.minimum(np.minimum(z0, z_end), z0 + vz0*tv - 0.5*G*tv*tv)
                        return z_min > h_base + h_rate*d + SLEEP_HEIGHT_MARGIN
                    lo = np.zeros(fly_d.size); hi = np.full(fly_d.size, SLEEP_WINDOW_S)
                    ok0 = safe(lo)
                    full = safe(hi)
                    for _ in range(40):
                        mid = 0.5*(lo + hi)
                        sm = safe(mid)
                        lo = np.where(sm, mid, lo); hi = np.where(sm, hi, mid)
                    delta = np.where(full, SLEEP_WINDOW_S, lo)
                    # Edge-exit root (straight xy from inside the disk).
                    px, py, pvx, pvy = x[fly_d], y[fly_d], vx[fly_d], vy[fly_d]
                    A = pvx*pvx + pvy*pvy
                    B = 2*(px*pvx + py*pvy)
                    C = px*px + py*py - pl.basis["a"]**2
                    disc = B*B - 4*A*C
                    with np.errstate(invalid="ignore", divide="ignore"):
                        root = np.where((A > 0) & (disc >= 0),
                                        (-B + np.sqrt(np.maximum(disc, 0)))/np.where(A > 0, 2*A, 1.0),
                                        np.inf)
                    delta = np.where(root > 0, np.minimum(delta, root), delta)
                    good = ok0 & (delta > 8*dt)
                    idx = fly_d[good]
                    self.sleeping[idx] = True
                    self.wake_t[idx] = t_new + delta[good]
                    self.sleep_t0[idx] = t_new
                # Failed evaluations back off 16 steps (management only).
                failed = done[~self.sleeping[done]]
                self.eval_cooldown[failed] = 16
        # Edge loss: sleeping flyers still hold their sleep-start xy (inside
        # the disk), so they are not falsely culled; particles that finished
        # this step outside are culled now (same semantics as TS).
        off = alive & (np.hypot(x, y) > pl.basis["a"])
        alive[off] = False
        self.sleeping[off] = False
        self.t = t_new

    def advance_to(self, t_target):
        """Advance to EXACTLY t_target (final step shortened; never overshoots)."""
        while self.t < t_target - 1e-15:
            self.step(min(self.dt_now(), t_target - self.t))
        self.materialize_all()
