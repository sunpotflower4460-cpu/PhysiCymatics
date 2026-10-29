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

    def force_for_amplitude(self, mode_index, target_A):
        """Steady on-resonance force giving surface amplitude target_A (max|Phi|=1)."""
        j = next(i for i, d in enumerate(self.dofs) if d["mode"] == mode_index)
        return target_A * self.mass[j] * 2*self.zeta * self.omega[j]**2 / abs(self.phi_drive[j])


class TransientSandSim:
    """Sand on the transient surface. Same contact order as sand_bounce.SandSim;
    detach when wtt < -g (support N<0); stick when |vn2| < 0.3*v_scale."""

    def __init__(self, plate, state, params=None):
        self.plate = plate; self.state = state
        self.p = params or sb.SandParams()
        self.t = 0.0; self.n_bounce = 0

    def dt_now(self):
        f = self.plate.score.freq_at(self.t) or self.plate.score.freq_at(self.t + 1e-9) or 152.0
        return 1.0/(self.p.spp * max(f, 20.0))

    def step(self):
        p, pl, st = self.p, self.plate, self.state
        e, beta = p.e, p.beta
        dt = self.dt_now(); t = self.t; t_new = t + dt
        pl.begin_step(t, t_new)
        v_stick = p.v_stick_frac * pl.v_scale()
        x, y, z, vx, vy, vz = st.x, st.y, st.z, st.vx, st.vy, st.vz
        sitting, alive = st.sitting, st.alive
        fly = alive & ~sitting
        if fly.any():
            fi = np.where(fly)[0]
            xf = x[fi] + vx[fi]*dt; yf = y[fi] + vy[fi]*dt
            zf = z[fi] + vz[fi]*dt - 0.5*G*dt*dt; vzf = vz[fi] - G*dt
            hs,_,_,_,_ = pl.surf_at(xf, yf, t_new)
            pen = zf < hs; ok = fi[~pen]
            x[ok] = xf[~pen]; y[ok] = yf[~pen]; z[ok] = zf[~pen]; vz[ok] = vzf[~pen]
            if pen.any():
                idx = fi[pen]
                t0 = np.full(idx.size, t); t1 = np.full(idx.size, t_new)
                for _ in range(p.bisect_iters):
                    tm = 0.5*(t0+t1); ddt = tm - t
                    xm = x[idx] + vx[idx]*ddt; ym = y[idx] + vy[idx]*ddt
                    zm = z[idx] + vz[idx]*ddt - 0.5*G*ddt*ddt
                    hm,_,_,_,_ = pl.surf_at(xm, ym, tm)
                    below = zm < hm
                    t1 = np.where(below, tm, t1); t0 = np.where(below, t0, tm)
                tc = 0.5*(t0+t1); ddt = tc - t
                xc = x[idx] + vx[idx]*ddt; yc = y[idx] + vy[idx]*ddt
                vzc = vz[idx] - G*ddt
                hc, hx, hy, uz, _wtt = pl.surf_at(xc, yc, tc)
                nn = np.sqrt(hx*hx + hy*hy + 1.0)
                nx_, ny_, nz_ = -hx/nn, -hy/nn, 1.0/nn
                rvx, rvy, rvz = vx[idx], vy[idx], vzc - uz
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
                self.n_bounce += int(hit.sum())
                stick = hit & (np.abs(vn2) < v_stick)
                si = idx[stick]; sitting[si] = True
                vx[si] = 0.0; vy[si] = 0.0
        sit = alive & sitting
        if sit.any():
            si = np.where(sit)[0]
            hs, hx, hy, uz, wtt = pl.surf_at(x[si], y[si], t_new)
            z[si] = hs
            if p.c_roll > 0.0:
                Nacc = np.maximum(G + wtt, 0.0)
                vx[si] += -p.c_roll*Nacc*hx*dt; vy[si] += -p.c_roll*Nacc*hy*dt
                sp = np.hypot(vx[si], vy[si]) + 1e-30
                drag = np.minimum(p.c_rr*Nacc*dt, sp)
                vx[si] -= drag*vx[si]/sp; vy[si] -= drag*vy[si]/sp
                x[si] += vx[si]*dt; y[si] += vy[si]*dt
            det = wtt < -G
            di = si[det]; sitting[di] = False
            vz[di] = uz[det]; z[di] = hs[det] + 1e-9
        off = alive & (np.hypot(x, y) > pl.basis["a"])
        alive[off] = False
        self.t = t_new

    def advance_to(self, t_target):
        while self.t < t_target - 1e-15:
            self.step()
