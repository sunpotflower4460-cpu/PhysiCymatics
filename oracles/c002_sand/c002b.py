# c002b: bounce + (optional) rolling — WHERE do grains gather? Let physics answer.
# ==============================================================================
# PUT-IN (all textbook mechanics, nothing about nodes):
#   gravity; exact ballistic flight; impact on moving tilted surface
#   (restitution e, tangential factor beta); sit when rebound tiny;
#   detach when support N<0;
#   NEW: seated grains may ROLL on the instantaneous incline:
#        a = -c_roll*(N/m)*grad(h) - c_rr*(N/m)*v_hat   (c_roll=5/7 solid sphere)
#   Two grain types: ANGULAR (c_roll=0, static friction holds) vs ROUND (c_roll=5/7).
# The mode shape W enters ONLY through surface height/velocity/accel/slope.
# MEASURED: direction of net transport, ring occupancy at c001 node radius,
#   trapped fraction, formation curve, edge losses. Direction is NOT assumed.
# ==============================================================================
import numpy as np, json, os, sys, time
from scipy import special, optimize

OUT = os.path.dirname(os.path.abspath(__file__))
G = 9.81

# ---------- c001 mode machinery (validated there) ----------
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

# ---------- simulator ----------
def simulate(mode, Gamma, e, beta, Np, t_max, seed, ic="uniform", init_xy=None,
             c_roll=0.0, c_rr=0.02, spp=44, v_stick_frac=0.3, log_every=0.5,
             bisect_iters=14):
    rng = np.random.default_rng(seed)
    a = mode["a"]; m = mode["m"]; NG1 = mode["ngrid"]-1
    Rt, dRt = mode["R"], mode["dR"]
    w = 2*np.pi*mode["f_hz"]
    A = Gamma*G/w**2
    T = 2*np.pi/w; dt = T/spp
    v_stick = v_stick_frac*A*w
    def W_grad(xp, yp):
        r = np.hypot(xp, yp)
        idx = np.minimum((r*(NG1/a)).astype(np.int64), NG1)   # fast table lookup
        R, dR = Rt[idx], dRt[idx]/a
        rs = np.where(r < 1e-12, 1e-12, r)
        c, s_ = xp/rs, yp/rs                                   # cos th, sin th
        if m == 0:
            return R, dR*c, dR*s_
        cm = np.cos(m*np.arctan2(yp, xp)); sm = np.sin(m*np.arctan2(yp, xp))
        W = R*cm; Wr = dR*cm; Wt = -m*R*sm/rs
        return W, Wr*c - Wt*s_, Wr*s_ + Wt*c
    def surf(xp, yp, s_, c_):
        Wp, Wxp, Wyp = W_grad(xp, yp)
        return A*Wp*s_, A*Wxp*s_, A*Wyp*s_, A*w*Wp*c_, Wp
    # initial condition
    if init_xy is not None:
        x, y = init_xy[0].copy(), init_xy[1].copy(); Np = len(x)
        rr = np.hypot(x, y); th = np.arctan2(y, x)
    elif ic == "uniform":
        rr = a*np.sqrt(rng.random(Np)); th = 2*np.pi*rng.random(Np)
    elif ic == "midslope":  # off node AND off antinode: fair audit-3 stress
        rr = np.clip(0.35*a + 0.02*a*rng.standard_normal(Np), 0, 0.999*a)
        th = 2*np.pi*rng.random(Np)
    else:  # all at strongest antinode ring (stress test: result must not be in IC)
        i_anti = int(np.argmax(np.abs(Rt))); r_anti = mode["rg"][i_anti]*a
        rr = np.clip(r_anti + 0.02*a*rng.standard_normal(Np), 0, 0.999*a)
        th = 2*np.pi*rng.random(Np)
    if init_xy is None:
        x, y = rr*np.cos(th), rr*np.sin(th)
    z = np.full(Np, 5e-4); vx = np.zeros(Np); vy = np.zeros(Np); vz = np.zeros(Np)
    sitting = np.zeros(Np, bool); alive = np.ones(Np, bool)
    r_node = mode["r_nodes"][0]*a if mode["r_nodes"] else np.nan
    band = 0.005
    exp_frac = 4*r_node*band/(a*a) if np.isfinite(r_node) else np.nan
    t = 0.0; hist = []; n_bounce = 0; next_log = 0.0
    nsteps = int(round(t_max/dt))
    for step in range(1, nsteps+1):
        t_new = step*dt
        s_, c_ = np.sin(w*t_new), np.cos(w*t_new)
        fly = alive & ~sitting
        if fly.any():
            fi = np.where(fly)[0]
            xf = x[fi] + vx[fi]*dt
            yf = y[fi] + vy[fi]*dt
            zf = z[fi] + vz[fi]*dt - 0.5*G*dt*dt
            vzf = vz[fi] - G*dt
            hs,_,_,_,_ = surf(xf, yf, s_, c_)
            pen = zf < hs
            ok = fi[~pen]
            x[ok] = xf[~pen]; y[ok] = yf[~pen]; z[ok] = zf[~pen]; vz[ok] = vzf[~pen]
            if pen.any():
                idx = fi[pen]
                t0 = np.full(idx.size, t); t1 = np.full(idx.size, t_new)
                for _ in range(bisect_iters):
                    tm = 0.5*(t0+t1); ddt = tm - t
                    xm = x[idx] + vx[idx]*ddt
                    ym = y[idx] + vy[idx]*ddt
                    zm = z[idx] + vz[idx]*ddt - 0.5*G*ddt*ddt
                    hm,_,_,_,_ = surf(xm, ym, np.sin(w*tm), np.cos(w*tm))
                    below = zm < hm
                    t1 = np.where(below, tm, t1); t0 = np.where(below, t0, tm)
                tc = 0.5*(t0+t1); ddt = tc - t
                xc = x[idx] + vx[idx]*ddt
                yc = y[idx] + vy[idx]*ddt
                zc = z[idx] + vz[idx]*ddt - 0.5*G*ddt*ddt
                vzc = vz[idx] - G*ddt
                sc_, cc_ = np.sin(w*tc), np.cos(w*tc)
                hc, hx, hy, uz, Wc = surf(xc, yc, sc_, cc_)
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
                n_bounce += int(hit.sum())
                stick = hit & (np.abs(vn2) < v_stick)
                si = idx[stick]
                sitting[si] = True
                if c_roll == 0.0:
                    vx[si] = 0.0; vy[si] = 0.0      # angular grain: static friction holds
        sit = alive & sitting
        if sit.any():
            si = np.where(sit)[0]
            hs, hx, hy, uz, Ws = surf(x[si], y[si], s_, c_)
            z[si] = hs
            Nacc = np.maximum(G - A*w*w*Ws*s_, 0.0)             # N/m = g + h_tt
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
        off = alive & (np.hypot(x, y) > a)
        alive[off] = False
        t = t_new
        if t >= next_log:
            am = alive
            Wp,_,_ = W_grad(x[am], y[am])
            r_al = np.hypot(x[am], y[am])
            ring = float(np.mean(np.abs(r_al - r_node) < band)) if np.isfinite(r_node) else np.nan
            hist.append(dict(t=round(t,2),
                             med_absW=float(np.median(np.abs(Wp))),
                             trapped=float(np.mean(Gamma*np.abs(Wp) < 1.0)),
                             ring_enrich=float(ring/exp_frac) if np.isfinite(r_node) else np.nan,
                             sitting=float(np.mean(sitting[am])),
                             alive=int(am.sum())))
            next_log += log_every
    am = alive
    Wp,_,_ = W_grad(x[am], y[am])
    return dict(x=x[am].copy(), y=y[am].copy(), absW=np.abs(Wp),
                sitting=sitting[am].copy(), hist=hist, n_bounce=n_bounce, A=A,
                r_node=r_node,
                params=dict(Gamma=Gamma, e=e, beta=beta, Np=Np, t_max=t_max, ic=ic,
                            c_roll=c_roll, c_rr=c_rr, f_hz=mode["f_hz"], m=mode["m"]))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--run", default="pilot")
    p.add_argument("--tmax", type=float, default=6.0)
    p.add_argument("--np", type=int, default=500, dest="npart")
    p.add_argument("--gamma", type=float, default=4.0)
    p.add_argument("--croll", type=float, default=0.0)
    p.add_argument("--ic", default="uniform")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--e", type=float, default=0.75)
    p.add_argument("--beta", type=float, default=0.8)
    a_ = p.parse_args()
    t0 = time.time()
    mode = make_mode(0, 0)
    print(f"mode(0,1): f={mode['f_hz']:.1f} Hz, r_node/a={mode['r_nodes'][0]:.4f}")
    r = simulate(mode, Gamma=a_.gamma, e=a_.e, beta=a_.beta, Np=a_.npart,
                 t_max=a_.tmax, seed=a_.seed, ic=a_.ic, c_roll=a_.croll)
    for h in r["hist"][::max(1,len(r["hist"])//8)]:
        print(h)
    print(f"[{a_.run}] bounces={r['n_bounce']}, A={r['A']*1e6:.1f}um, wall={time.time()-t0:.1f}s")
    np.savez(f"{OUT}/c002b_{a_.run}.npz", x=r["x"], y=r["y"], absW=r["absW"],
             sitting=r["sitting"], r_node=r["r_node"],
             hist_json=json.dumps(r["hist"]), params_json=json.dumps(r["params"]))
