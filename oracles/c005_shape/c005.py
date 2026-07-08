# c005: EARN THE SHAPE — direct simulation of the Zhang-Vinals equations.
# ==============================================================================
# PUT-IN (established, source captured: Topaz-Silber arXiv:nlin/0111039 quoting
#   Zhang & Vinals JFM 336:301 (1997), derived from Navier-Stokes,
#   quasi-potential approx for weakly viscous, semi-infinite fluid):
#     (dt - g*Lap) h   - D[Phi] = F(h,Phi)          [eq 19]
#     (dt - g*Lap) Phi - (Gam0*Lap - G(t)) h = G(h,Phi)   [eq 20]
#   F, G = the full quadratic+cubic terms, eqs (21)-(22), transcribed verbatim.
#   D = |k| in Fourier space. Scalings eq (24): time by omega_drive;
#   gamma=2 nu k^2/w, Gam0=sigma k^3/(rho w^2), G0=g k/w^2, f=a k/w^2,
#   G0+Gam0=1/4. Water @ 440 Hz: k_c=2971/m -> gamma=0.0064, Gam0=0.2462.
# FORBIDDEN: no pattern named, no angular seeding (isotropic noise IC),
#   no term that knows about 90 or 60 degrees.
# FAIRNESS: periodic square domain quantizes angles. Ring radius ~18 lattice
#   units -> unstable band covers shells R^2 in [318,333] giving pair angles
#   incl. ~8,27,45,56,63,70,90 deg. 60-ish IS representable (63.4 via (8,16)).
# MEASURED: which angular structure survives; temporal spectrum (f/2 check);
#   linear onset vs c003 threshold (route-2 for this new code).
# ==============================================================================
import numpy as np, json, os, sys, time
try:
    from scipy import fft as sfft
    F2 = lambda a: sfft.rfft2(a, workers=2)
    IF2 = lambda A, N: sfft.irfft2(A, s=(N, N), workers=2)
except Exception:
    F2 = np.fft.rfft2
    IF2 = lambda A, N: np.fft.irfft2(A, s=(N, N))

OUT = os.path.dirname(os.path.abspath(__file__))
G_SI, RHO, SIG, NU = 9.81, 998.0, 0.0728, 1.004e-6

def setup(f_drive=440.0, N=128, R=18.0):
    w = 2*np.pi*f_drive
    from scipy import optimize
    kc = optimize.brentq(lambda k: (G_SI*k + SIG*k**3/RHO) - (w/2)**2, 1e-2, 1e6)
    gam = 2*NU*kc*kc/w
    G0 = G_SI*kc/w**2
    Gam0 = 0.25 - G0
    L = 2*np.pi*R                      # k_c = 1 sits at lattice radius R
    ki = np.fft.fftfreq(N, d=1.0/N)    # integer lattice
    kj = np.fft.rfftfreq(N, d=1.0/N)
    KX = ki[:, None]/R; KY = kj[None, :]/R
    K2 = KX*KX + KY*KY
    KABS = np.sqrt(K2)
    kmax = (N//2)*(2.0/3.0)/R
    deal = (KABS <= kmax)
    return dict(f=f_drive, w=w, kc=kc, gam=gam, G0=G0, Gam0=Gam0, N=N, R=R,
                KX=KX, KY=KY, K2=K2, KABS=KABS, deal=deal)

def rhs(hh, ph, t, P, fdrive, nonlinear=True):
    """spectral in, spectral out. hh,ph = rFFT2 of h,Phi."""
    KX, KY, K2, KABS, deal = P["KX"], P["KY"], P["K2"], P["KABS"], P["deal"]
    G0, Gam0 = P["G0"], P["Gam0"]; N = P["N"]
    F = lambda A: IF2(A, N)
    h   = F(hh)
    DP  = F(KABS*ph)
    dh_lin = KABS*ph
    dp_lin = -(Gam0*K2 + G0 - fdrive*np.cos(t))*hh
    if not nonlinear:
        return deal*dh_lin, deal*dp_lin
    Px  = F(1j*KX*ph); Py = F(1j*KY*ph)
    hx  = F(1j*KX*hh); hy = F(1j*KY*hh)
    lapP= F(-K2*ph)
    q1h = F2(h*DP)
    Dq1 = F(KABS*q1h)
    Fh = (-(1j*KX*F2(h*Px) + 1j*KY*F2(h*Py))
          - 0.5*K2*F2(h*h*DP)
          - KABS*q1h
          + KABS*(F2(h*Dq1) + 0.5*F2(h*h*lapP)))
    grad2 = hx*hx + hy*hy
    Greal = 0.5*DP*DP - 0.5*(Px*Px + Py*Py) - DP*(h*lapP + Dq1)
    Gh = (F2(Greal)
          - 0.5*Gam0*(1j*KX*F2(hx*grad2) + 1j*KY*F2(hy*grad2)))
    return deal*(dh_lin + Fh), deal*(dp_lin + Gh)

def step_ifrk4(hh, ph, t, dt, P, fdrive, E, E2, nonlinear=True):
    k1h, k1p = rhs(hh, ph, t, P, fdrive, nonlinear)
    a_h, a_p = E2*(hh + 0.5*dt*k1h), E2*(ph + 0.5*dt*k1p)
    k2h, k2p = rhs(a_h, a_p, t+0.5*dt, P, fdrive, nonlinear)
    b_h, b_p = E2*hh + 0.5*dt*k2h, E2*ph + 0.5*dt*k2p
    k3h, k3p = rhs(b_h, b_p, t+0.5*dt, P, fdrive, nonlinear)
    c_h, c_p = E*hh + dt*E2*k3h, E*ph + dt*E2*k3p
    k4h, k4p = rhs(c_h, c_p, t+dt, P, fdrive, nonlinear)
    hh = E*hh + dt/6.0*(E*k1h + 2*E2*(k2h+k3h) + k4h)
    ph = E*ph + dt/6.0*(E*k1p + 2*E2*(k2p+k3p) + k4p)
    hh[0,0] = 0.0; ph[0,0] = 0.0                  # gauge / volume conservation
    return hh, ph

def run(P, fdrive, n_periods, seed=7, amp0=1e-3, spp=48, nonlinear=True,
        log_every=10, tag="run"):
    N = P["N"]; gam = P["gam"]
    dt = 2*np.pi/spp
    E  = np.exp(-gam*P["K2"]*dt)
    E2 = np.exp(-gam*P["K2"]*dt*0.5)
    rng = np.random.default_rng(seed)
    ki = np.fft.fftfreq(N, d=1.0/N); kj = np.fft.rfftfreq(N, d=1.0/N)
    I2 = ki[:,None]**2 + kj[None,:]**2
    band = (I2 >= 317) & (I2 <= 334)
    ck = f"{OUT}/c005_{tag}_ckpt.npz"
    if os.path.exists(ck):
        d = np.load(ck); hh = d["hh"]; ph = d["ph"]; t0 = float(d["t"])
        print(f"[{tag}] resumed at t={t0:.1f} ({t0/(2*np.pi):.0f} periods)", flush=True)
    else:
        h0 = amp0*rng.standard_normal((N, N))
        hh = F2(h0)*P["deal"]; hh[0,0]=0
        hh = hh*band                      # isotropic ring-band seed (radial filter only)
        ph = np.zeros_like(hh); t0 = 0.0
    ringlog = []
    t = t0
    for step in range(1, n_periods*spp+1):
        hh, ph = step_ifrk4(hh, ph, t, dt, P, fdrive, E, E2, nonlinear)
        t = t0 + step*dt
        if step % spp == 0:
            amp = np.sqrt(float(np.sum((np.abs(hh)**2)[band])))
            if (step//spp) % log_every == 0:
                hr = IF2(hh, N)
                print(f"[{tag}] per={(t/(2*np.pi)):4.0f} bandnorm={amp:.3e} rms_h={float(np.sqrt(np.mean(hr*hr))):.3e}", flush=True)
            ringlog.append(amp)
        if step % (spp*20) == 0:
            np.savez(ck, hh=hh, ph=ph, t=t)
    np.savez(ck, hh=hh, ph=ph, t=t)
    return hh, ph, np.array(ringlog)

def shell_masks(P):
    """lattice shells covering the unstable band around |k|=1 (R^2 317..334)."""
    N, R = P["N"], P["R"]
    ki = np.fft.fftfreq(N, d=1.0/N)
    I2 = (ki[:,None]**2 + ki[None,:]**2)
    masks = {}
    for s in range(317, 335):
        m = np.isclose(I2, s)
        if m.any(): masks[s] = m
    return masks

def angle_spectrum(hh, P, topn=12):
    N = P["N"]
    ki = np.fft.fftfreq(N, d=1.0/N); kj = np.fft.rfftfreq(N, d=1.0/N)
    A, B = np.meshgrid(ki, kj, indexing="ij")
    I2 = A*A + B*B
    band = (I2 >= 317) & (I2 <= 334) & (B >= 0) & ~((B==0)&(A<0))
    aa, bb, pw = A[band], B[band], (np.abs(hh)**2)[band]
    th = (np.degrees(np.arctan2(bb, aa))) % 180.0
    order = np.argsort(pw)[::-1]
    peaks = [(float(th[i]), float(pw[i]/ (pw[order[0]]+1e-300))) for i in order[:topn]]
    return peaks, float(pw.sum())

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "linear"
    P = setup(f_drive=float(sys.argv[4]) if len(sys.argv) > 4 else 440.0)
    fc0 = 2*P["gam"]
    print(f"setup: f=440Hz kc={P['kc']:.1f} gamma={P['gam']:.5f} "
          f"Gam0={P['Gam0']:.4f} G0={P['G0']:.5f}  f_c(leading)={fc0:.5f} "
          f"(c003: a_c=32.8 m/s^2 -> f_scaled={32.8*P['kc']/P['w']**2:.5f})")
    t0 = time.time()
    if mode == "linear":
        # route-2 validation of THIS code: growth sign must flip at c003 threshold
        for fac in (0.9, 1.0, 1.1):
            hh, ph, log = run(P, fac*fc0, 60, spp=48, nonlinear=False,
                              amp0=1e-6, log_every=1000, tag=f"lin{fac}")
            gmeas = np.log(log[-1]/log[19]) / (2*np.pi*(len(log)-20))
            print(f"  f/f_c={fac}: growth per unit time = {gmeas:+.5f} "
                  f"(expect sign {'+' if fac>1 else '-' if fac<1 else '~0'})", flush=True)
        print(f"wall={time.time()-t0:.0f}s")
    elif mode == "flagship":
        fac = float(sys.argv[2]) if len(sys.argv) > 2 else 1.30
        nper = int(sys.argv[3]) if len(sys.argv) > 3 else 800
        hh, ph, log = run(P, fac*fc0, nper, spp=36, nonlinear=True,
                          amp0=3e-3, seed=7, log_every=20, tag=sys.argv[5] if len(sys.argv) > 5 else "main")
        np.savez(f"{OUT}/c005_{sys.argv[5] if len(sys.argv) > 5 else chr(109)+chr(97)+chr(105)+chr(110)}.npz", hh=hh, ph=ph, log=log)
        peaks, tot = angle_spectrum(hh, P)
        print("angle peaks (deg, relpow):", [(round(a,1), round(p,3)) for a,p in peaks])
        print(f"wall={time.time()-t0:.0f}s")
