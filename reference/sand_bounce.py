"""reference/sand_bounce.py — c002 bounce mechanics as a reusable state update.

PUT-IN (all textbook mechanics, same as oracles/c002_sand/c002b.py):
  gravity; exact ballistic flight between surface contacts; impact on the
  moving tilted surface (normal restitution e, tangential factor beta);
  sit when the rebound normal speed is tiny (|vn2| < v_stick_frac*A*omega);
  detach when the required support goes negative (N < 0). The mode shape
  enters ONLY through surface height / velocity / acceleration / slope.
EMERGED (in c002, measured — not assumed here):
  transport of grains toward low |W|, ring enrichment at the c001 node
  radius, Gamma window (0.8 frozen / 15 weak / 40 / 80 strong), edge losses.
CLAIM-TIER : measured (within the c002 numeric model). This module is a
  faithful extraction of the c002b state update; agreement with the oracle
  is a port check, NOT a measurement of real sand.
FLOORS (inherited from c002, do not soften):
  no air; no grain-grain interaction (dilute limit); no back-reaction of
  sand mass on the plate mode; e=0.75 / beta=0.8 / v_stick_frac=0.3 are
  ASSUMED (unmeasured) values; beta is a simplification of Coulomb friction;
  single prescribed sinusoidal mode (no transient multi-mode plate);
  angular grains only in this extraction (c_roll=0); formation curves are
  not saturated at 24 s in c002.

Dynamics code below never references a target pattern: the mode shape is
used only to evaluate the surface. Radii used for *measurement* live in
SandSim.stats() / measure(), which only read state.

Fidelity contract: make_mode() and SandSim.step() replicate the arithmetic
of oracles/c002_sand/c002b.py (same 4096-point signed table, same fast
table lookup, same update order, same bisection count), so short-run
parity with the oracle can be checked numerically. The oracle file itself
is not modified.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy import special, optimize

G = 9.81


# ---------- c001/c002 mode machinery (verbatim numerics from c002b) ----------
def bc_matrix(lam, m, nu):
    x = lam
    Jm, Jp, Jpp = special.jv(m, x), special.jvp(m, x, 1), special.jvp(m, x, 2)
    Im, Ip, Ipp = special.iv(m, x), special.ivp(m, x, 1), special.ivp(m, x, 2)
    s = Im if Im != 0 else 1.0
    M_J = x*x*Jpp + nu*(x*Jp - m*m*Jm)
    M_I = (x*x*Ipp + nu*(x*Ip - m*m*Im)) / s
    V_J = -x**3*Jp - (1-nu)*m*m*(x*Jp - Jm)
    V_I = ( x**3*Ip - (1-nu)*m*m*(x*Ip - Im)) / s
    return np.array([[M_J, M_I],[V_J, V_I]]), s


def det_free(lam, m, nu):
    A,_ = bc_matrix(lam, m, nu); return A[0,0]*A[1,1]-A[0,1]*A[1,0]


def det_roots(m, nu, lam_max=25.0):
    xs = np.arange(0.3, lam_max, 0.005)
    v = np.array([det_free(x,m,nu) for x in xs]); out=[]
    for i in range(len(xs)-1):
        if v[i]*v[i+1] < 0:
            out.append(optimize.brentq(det_free, xs[i], xs[i+1], args=(m,nu), xtol=1e-12))
    return out


def make_mode(m, root_idx, nu=0.30, E=2.0e11, rho=7850.0, h=1e-3, a=0.12, ngrid=4096):
    lam = det_roots(m, nu)[root_idx]
    A,s = bc_matrix(lam, m, nu)
    r0 = A[0] if np.hypot(*A[0])>=np.hypot(*A[1]) else A[1]
    ac, bs = r0[1], -r0[0]; bc = bs/s
    rg = np.linspace(0,1,ngrid)
    R  = ac*special.jv(m, lam*rg) + bc*special.iv(m, lam*rg)
    nrm = np.max(np.abs(R)); R/=nrm
    dR = (ac*lam*special.jvp(m, lam*rg,1) + bc*lam*special.ivp(m, lam*rg,1))/nrm
    D = E*h**3/(12*(1-nu**2)); f = lam**2/(2*np.pi*a*a)*np.sqrt(D/(rho*h))
    zr = [ (rg[i]+rg[i+1])/2 for i in range(int(0.05*ngrid), ngrid-1) if R[i]*R[i+1]<0 ]
    return dict(m=m, lam=lam, rg=rg, R=np.ascontiguousarray(R),
                dR=np.ascontiguousarray(dR), f_hz=f, a=a, r_nodes=zr, ngrid=ngrid)


# ---------- surface contract (signed, SI, time-dependent) ----------
class SingleModeSurface:
    """w(x,y,t) = A * W(x,y) * sin(omega*t + phase0), all terms from one formula.

    W / dW come from the same 4096-point signed radial table and the same
    fast lookup as c002b (index = min(int(r*(ngrid-1)/a), ngrid-1); no
    interpolation). Gradients wrt physical distance include the 1/a factor
    because the table coordinate is r/a and dR is stored wrt r/a.
    """

    def __init__(self, mode, Gamma, phase0=0.0):
        self.mode = mode
        self.a = mode["a"]; self.m = mode["m"]
        self.NG1 = mode["ngrid"] - 1
        self.Rt, self.dRt = mode["R"], mode["dR"]
        self.omega = 2*np.pi*mode["f_hz"]
        self.f_hz = mode["f_hz"]
        self.Gamma = float(Gamma)
        self.A = Gamma*G/self.omega**2
        self.phase0 = float(phase0)

    def W_grad(self, xp, yp):
        r = np.hypot(xp, yp)
        idx = np.minimum((r*(self.NG1/self.a)).astype(np.int64), self.NG1)
        R, dR = self.Rt[idx], self.dRt[idx]/self.a
        rs = np.where(r < 1e-12, 1e-12, r)
        c, s_ = xp/rs, yp/rs
        if self.m == 0:
            return R, dR*c, dR*s_
        cm = np.cos(self.m*np.arctan2(yp, xp)); sm = np.sin(self.m*np.arctan2(yp, xp))
        W = R*cm; Wr = dR*cm; Wt = -self.m*R*sm/rs
        return W, Wr*c - Wt*s_, Wr*s_ + Wt*c

    def _surf(self, xp, yp, s_, c_):
        # Same return order as c002b's surf(): (h, hx, hy, uz, W).
        Wp, Wxp, Wyp = self.W_grad(xp, yp)
        return self.A*Wp*s_, self.A*Wxp*s_, self.A*Wyp*s_, self.A*self.omega*Wp*c_, Wp

    def sample(self, x_m, y_m, t_s):
        """Surface sample at (x_m, y_m, t_s). Scalars or arrays broadcast.

        Returns dict with inside, w_m, wt_ms, wtt_ms2, wx, wy (SI units;
        wx/wy dimensionless slopes wrt metres).
        """
        xp = np.asarray(x_m, dtype=np.float64)
        yp = np.asarray(y_m, dtype=np.float64)
        th = self.omega*t_s + self.phase0
        s_, c_ = np.sin(th), np.cos(th)
        Wp, Wxp, Wyp = self.W_grad(xp, yp)
        return dict(
            inside=np.hypot(xp, yp) <= self.a,
            w_m=self.A*Wp*s_,
            wt_ms=self.A*self.omega*Wp*c_,
            wtt_ms2=-self.A*self.omega**2*Wp*s_,
            wx=self.A*Wxp*s_,
            wy=self.A*Wyp*s_,
        )


@dataclass
class ParticleState:
    x: np.ndarray; y: np.ndarray; z: np.ndarray
    vx: np.ndarray; vy: np.ndarray; vz: np.ndarray
    alive: np.ndarray; sitting: np.ndarray; pid: np.ndarray

    @property
    def n_initial(self): return int(self.pid.size)
    @property
    def n_alive(self): return int(self.alive.sum())
    @property
    def n_lost(self): return self.n_initial - self.n_alive

    def copy(self):
        return ParticleState(**{k: getattr(self, k).copy() for k in
                             ("x","y","z","vx","vy","vz","alive","sitting","pid")})


def initial_state(mode, Np, seed=1, ic="uniform", init_xy=None):
    """Initial particles. Uniform/midslope/antinode sampling copies c002b."""
    rng = np.random.default_rng(seed)
    a = mode["a"]; Rt = mode["R"]
    if init_xy is not None:
        x, y = np.asarray(init_xy[0], dtype=np.float64).copy(), np.asarray(init_xy[1], dtype=np.float64).copy()
        Np = len(x)
    else:
        if ic == "uniform":
            rr = a*np.sqrt(rng.random(Np)); th = 2*np.pi*rng.random(Np)
        elif ic == "midslope":
            rr = np.clip(0.35*a + 0.02*a*rng.standard_normal(Np), 0, 0.999*a)
            th = 2*np.pi*rng.random(Np)
        else:
            i_anti = int(np.argmax(np.abs(Rt))); r_anti = mode["rg"][i_anti]*a
            rr = np.clip(r_anti + 0.02*a*rng.standard_normal(Np), 0, 0.999*a)
            th = 2*np.pi*rng.random(Np)
        x, y = rr*np.cos(th), rr*np.sin(th)
    return ParticleState(
        x=x, y=y, z=np.full(Np, 5e-4),
        vx=np.zeros(Np), vy=np.zeros(Np), vz=np.zeros(Np),
        alive=np.ones(Np, bool), sitting=np.zeros(Np, bool),
        pid=np.arange(Np, dtype=np.uint32)), rng.bit_generator.state


# ---------- state update: one fixed dt, c002b order preserved ----------
@dataclass
class SandParams:
    e: float = 0.75
    beta: float = 0.8
    c_roll: float = 0.0
    c_rr: float = 0.02
    spp: int = 44
    v_stick_frac: float = 0.3
    bisect_iters: int = 14


class SandSim:
    """Fixed-step sim. step() advances exactly dt = T/spp of physics time."""

    def __init__(self, surface, state, params=None, t=0.0):
        self.surface = surface
        self.state = state
        self.p = params or SandParams()
        self.t = float(t)
        self.step_index = int(round(self.t / self.dt)) if self.dt else 0
        self.n_bounce = 0

    @property
    def dt(self):
        return (2*np.pi/self.surface.omega) / self.p.spp

    @property
    def v_stick(self):
        return self.p.v_stick_frac*self.surface.A*self.surface.omega

    def step(self, dt_override=None):
        p, s, st = self.p, self.surface, self.state
        e, beta, c_roll, c_rr = p.e, p.beta, p.c_roll, p.c_rr
        w, A, t = s.omega, s.A, self.t
        dt = self.dt if dt_override is None else float(dt_override)
        x, y, z = st.x, st.y, st.z
        vx, vy, vz = st.vx, st.vy, st.vz
        sitting, alive = st.sitting, st.alive
        v_stick = self.v_stick
        # Default fixed-step path keeps c002b's (step_index+1)*dt exactly;
        # shortened final steps (advance_to) use t+dt.
        t_new = (self.step_index + 1)*self.dt if dt_override is None else t + dt
        s_, c_ = np.sin(w*t_new), np.cos(w*t_new)
        fly = alive & ~sitting
        if fly.any():
            fi = np.where(fly)[0]
            xf = x[fi] + vx[fi]*dt
            yf = y[fi] + vy[fi]*dt
            zf = z[fi] + vz[fi]*dt - 0.5*G*dt*dt
            vzf = vz[fi] - G*dt
            hs,_,_,_,_ = s._surf(xf, yf, s_, c_)
            pen = zf < hs
            ok = fi[~pen]
            x[ok] = xf[~pen]; y[ok] = yf[~pen]; z[ok] = zf[~pen]; vz[ok] = vzf[~pen]
            if pen.any():
                idx = fi[pen]
                t0 = np.full(idx.size, t); t1 = np.full(idx.size, t_new)
                for _ in range(p.bisect_iters):
                    tm = 0.5*(t0+t1); ddt = tm - t
                    xm = x[idx] + vx[idx]*ddt
                    ym = y[idx] + vy[idx]*ddt
                    zm = z[idx] + vz[idx]*ddt - 0.5*G*ddt*ddt
                    hm,_,_,_,_ = s._surf(xm, ym, np.sin(w*tm), np.cos(w*tm))
                    below = zm < hm
                    t1 = np.where(below, tm, t1); t0 = np.where(below, t0, tm)
                tc = 0.5*(t0+t1); ddt = tc - t
                xc = x[idx] + vx[idx]*ddt
                yc = y[idx] + vy[idx]*ddt
                zc = z[idx] + vz[idx]*ddt - 0.5*G*ddt*ddt
                vzc = vz[idx] - G*ddt
                sc_, cc_ = np.sin(w*tc), np.cos(w*tc)
                hc, hx, hy, uz, Wc = s._surf(xc, yc, sc_, cc_)
                nn = np.sqrt(hx*hx + hy*hy + 1.0)
                nx_, ny_, nz_ = -hx/nn, -hy/nn, 1.0/nn
                rvx, rvy, rvz = vx[idx], vy[idx], vzc - uz
                vn = rvx*nx_ + rvy*ny_ + rvz*nz_
                hit = vn < 0
                rvx2 = rvx - (1+e)*vn*nx_
                rvy2 = rvy - (1+e)*vn*ny_
                rvz2 = rvz - (1+e)*vn*nz_
                vn2 = rvx2*nx_ + rvy2*ny_ + rvz2*nz_
                tx, ty, tz = rvx2-vn2*nx_, rvy2-vn2*ny_, rvz2-vn2*nz_
                rvx3 = vn2*nx_ + beta*tx
                rvy3 = vn2*ny_ + beta*ty
                rvz3 = vn2*nz_ + beta*tz
                wi = idx[hit]
                x[wi] = xc[hit]; y[wi] = yc[hit]; z[wi] = hc[hit] + 1e-9
                vx[wi] = rvx3[hit]; vy[wi] = rvy3[hit]; vz[wi] = rvz3[hit] + uz[hit]
                self.n_bounce += int(hit.sum())
                stick = hit & (np.abs(vn2) < v_stick)
                si = idx[stick]
                sitting[si] = True
                if c_roll == 0.0:
                    vx[si] = 0.0; vy[si] = 0.0
        sit = alive & sitting
        if sit.any():
            si = np.where(sit)[0]
            hs, hx, hy, uz, Ws = s._surf(x[si], y[si], s_, c_)
            z[si] = hs
            Nacc = np.maximum(G - A*w*w*Ws*s_, 0.0)
            if c_roll > 0.0:
                vx[si] += -c_roll*Nacc*hx*dt
                vy[si] += -c_roll*Nacc*hy*dt
                sp = np.hypot(vx[si], vy[si]) + 1e-30
                drag = np.minimum(c_rr*Nacc*dt, sp)
                vx[si] -= drag*vx[si]/sp
                vy[si] -= drag*vy[si]/sp
                x[si] += vx[si]*dt; y[si] += vy[si]*dt
            det = (A*w*w*Ws*s_) > G
            di = si[det]
            sitting[di] = False
            vz[di] = uz[det]; z[di] = hs[det] + 1e-9
        off = alive & (np.hypot(x, y) > s.a)
        alive[off] = False
        self.t = t_new
        self.step_index += 1

    def advance_to(self, t_target):
        """Advance to EXACTLY t_target (final step shortened; never overshoots)."""
        while self.t < t_target - 1e-15:
            rem = t_target - self.t
            if rem >= self.dt - 1e-18:
                self.step()
            else:
                self.step(rem)
                self.step_index = int(round(self.t/self.dt))

    def snapshot(self):
        st = self.state
        return dict(t_s=self.t, step_index=self.step_index,
                    n_initial=st.n_initial, n_alive=st.n_alive,
                    n_lost=st.n_lost, n_bounce=self.n_bounce,
                    state=st.copy())

    # ----- measurement only (reads state; never feeds back into step) -----
    def stats(self):
        s, st = self.surface, self.state
        am = st.alive
        r_node = s.mode["r_nodes"][0]*s.a if s.mode["r_nodes"] else np.nan
        band = 0.005
        exp_frac = 4*r_node*band/(s.a*s.a) if np.isfinite(r_node) else np.nan
        Wp,_,_ = s.W_grad(st.x[am], st.y[am])
        r_al = np.hypot(st.x[am], st.y[am])
        ring = float(np.mean(np.abs(r_al - r_node) < band)) if np.isfinite(r_node) else np.nan
        return dict(t=round(float(self.t), 2),
                    med_absW=float(np.median(np.abs(Wp))),
                    trapped=float(np.mean(s.Gamma*np.abs(Wp) < 1.0)),
                    ring_enrich=float(ring/exp_frac) if np.isfinite(r_node) else float("nan"),
                    sitting=float(np.mean(st.sitting[am])),
                    alive=int(am.sum()))


def simulate(mode, Gamma, e=0.75, beta=0.8, Np=500, t_max=6.0, seed=1,
             ic="uniform", init_xy=None, c_roll=0.0, c_rr=0.02, spp=44,
             v_stick_frac=0.3, log_every=0.5, bisect_iters=14):
    """Drop-in equivalent of c002b.simulate(); same return keys as the oracle."""
    surface = SingleModeSurface(mode, Gamma)
    state, _ = initial_state(mode, Np, seed=seed, ic=ic, init_xy=init_xy)
    sim = SandSim(surface, state, SandParams(e=e, beta=beta, c_roll=c_roll,
                    c_rr=c_rr, spp=spp, v_stick_frac=v_stick_frac,
                    bisect_iters=bisect_iters))
    hist = []; next_log = 0.0
    nsteps = int(round(t_max/sim.dt))
    for _ in range(nsteps):
        sim.step()
        if sim.t >= next_log:
            hist.append(sim.stats()); next_log += log_every
    st = sim.state; am = st.alive
    Wp,_,_ = surface.W_grad(st.x[am], st.y[am])
    r_node = mode["r_nodes"][0]*mode["a"] if mode["r_nodes"] else np.nan
    return dict(x=st.x[am].copy(), y=st.y[am].copy(), absW=np.abs(Wp),
                sitting=st.sitting[am].copy(), hist=hist,
                n_bounce=sim.n_bounce, A=surface.A, r_node=r_node,
                params=dict(Gamma=Gamma, e=e, beta=beta, Np=int(state.n_initial),
                            t_max=t_max, ic=ic, c_roll=c_roll, c_rr=c_rr,
                            f_hz=mode["f_hz"], m=mode["m"]))
