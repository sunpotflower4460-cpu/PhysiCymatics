# c003: Faraday instability — shake water, WATCH what grows. Analog-style, in silico.
# ==============================================================================
# PUT-IN (real water, real mechanics; nothing about f/2 or thresholds):
#   - linearized free-surface dynamics per wavenumber k under vertically
#     oscillating gravity  g(t) = g*(1 + Gamma*cos(w_d t)):
#         eta_tt + 2*gamma_k*eta_t + [ (g(t)*k + sigma*k^3/rho) * tanh(kh) ] * eta = 0
#     (Benjamin-Ursell 1954 structure; established derivation)
#   - water @20C: rho=998, sigma=0.0728, nu=1.004e-6  (established values)
#   - damping submodel (THE stated simplification, two variants -> honesty band):
#       bulk:    gamma = 2*nu*k^2                     (weak-damping free surface)
#       bulk+BL: + k*sqrt(nu*omega0/2)/sinh(2kh)      (bottom boundary layer;
#                 coefficient from memory -> YELLOW anchor, needs source check)
# FORBIDDEN: no "f/2", no threshold values, no selected wavelength anywhere in
#   the dynamics. We integrate one drive period and ask the physics: does a
#   cycle amplify? (monodromy of the REAL equation = the analog question)
# MEASURED: growth/decay -> threshold by bisection; response type from the
#   Floquet multiplier SIGN (negative = sign flip each period = subharmonic,
#   discovered not assumed); wavelength = argmin of measured threshold over k;
#   long noise-seeded runs + FFT for the experiment-style demonstration.
# CROSS-CHECK: Liouville det(M)=exp(-2*gamma*T) (integrator truth);
#   weak-damping analytic asymptote Gc = 4*gamma*omega0/(g*k*tanh(kh)).
# ==============================================================================
import numpy as np, json, os, sys, time
from scipy import optimize

OUT = os.path.dirname(os.path.abspath(__file__))
G, RHO, SIG, NU = 9.81, 998.0, 0.0728, 1.004e-6   # water @ ~20C (established)

def om0_sq(k, h):
    return (G*k + SIG*k**3/RHO)*np.tanh(np.minimum(k*h, 20.0))

def gamma_k(k, h, variant):
    g0 = 2.0*NU*k*k
    if variant == "bulk":
        return g0
    kh = np.minimum(k*h, 20.0)
    om0 = np.sqrt(om0_sq(k, h))
    return g0 + k*np.sqrt(NU*om0/2.0)/np.sinh(2.0*kh)   # + bottom BL (memory-based)

def monodromy(k, h, Gam, w_d, variant, nsteps=1200):
    """Integrate the real equation over ONE drive period for two independent
    initial conditions; return growth rate s, trace, det (all vectorized in k).
    This is the analog question in minimal form: does one shake-cycle amplify?"""
    om0s = om0_sq(k, h)
    gam  = gamma_k(k, h, variant)
    mod  = Gam*G*k*np.tanh(np.minimum(k*h, 20.0))   # only gravity is modulated
    T = 2*np.pi/w_d; dt = T/nsteps
    E = np.stack([np.ones_like(k), np.zeros_like(k)])   # (2 ICs, Nk)
    V = np.stack([np.zeros_like(k), np.ones_like(k)])
    def rhs(E, V, t):
        return V, -(om0s + mod*np.cos(w_d*t))*E - 2.0*gam*V
    t = 0.0
    for i in range(nsteps):
        k1e, k1v = rhs(E, V, t)
        k2e, k2v = rhs(E+0.5*dt*k1e, V+0.5*dt*k1v, t+0.5*dt)
        k3e, k3v = rhs(E+0.5*dt*k2e, V+0.5*dt*k2v, t+0.5*dt)
        k4e, k4v = rhs(E+dt*k3e, V+dt*k3v, t+dt)
        E = E + dt/6.0*(k1e+2*k2e+2*k3e+k4e)
        V = V + dt/6.0*(k1v+2*k2v+2*k3v+k4v)
        t += dt
    tr  = E[0] + V[1]
    det = E[0]*V[1] - E[1]*V[0]
    disc = tr*tr - 4.0*det
    lam = np.where(disc >= 0, (np.abs(tr)+np.sqrt(np.maximum(disc, 0)))/2.0,
                   np.sqrt(np.maximum(det, 1e-300)))
    s = np.log(np.maximum(lam, 1e-300))/T
    return s, tr, det, gam, T

def Gc_analytic(k, h, variant, w_d):
    """weak-damping subharmonic threshold WITH detuning (route 2 cross-check):
       M_c = 2*w_d*sqrt(gamma^2 + delta^2),  delta = omega0(k) - w_d/2"""
    om0 = np.sqrt(om0_sq(k, h)); gam = gamma_k(k, h, variant)
    delta = om0 - 0.5*w_d
    Mc = 2.0*w_d*np.sqrt(gam*gam + delta*delta)
    return Mc/(G*k*np.tanh(np.minimum(k*h, 20.0)))

def k_of_freq(f_wave, h):
    w = 2*np.pi*f_wave
    fn = lambda k: om0_sq(k, h) - w*w
    return optimize.brentq(fn, 1e-2, 1e6)

def bisect_curve(kg, h, w_d, variant, iters=20):
    """measured neutral curve: bisect drive amplitude per k (vectorized)."""
    Ghi = 1.6*Gc_analytic(kg, h, variant, w_d) + 1e-4
    for _ in range(8):                        # ensure bracket (expand if theory low)
        s,_,_,_,_ = monodromy(kg, h, Ghi, w_d, variant)
        need = s <= 0
        if not need.any(): break
        Ghi[need] *= 1.7
    Glo = np.zeros_like(kg)
    for _ in range(iters):
        Gm = 0.5*(Glo+Ghi)
        s,_,_,_,_ = monodromy(kg, h, Gm, w_d, variant)
        up = s > 0
        Ghi = np.where(up, Gm, Ghi); Glo = np.where(up, Glo, Gm)
    return 0.5*(Glo+Ghi)

def sweep_f(f_d, h, variant):
    """three-stage k refinement: the water tongue is razor-thin (dk/k ~ 0.5%),
       a single log grid ALIASES it (caught by route-2 disagreement — recorded)."""
    w_d = 2*np.pi*f_d
    k_half = k_of_freq(f_d/2.0, h)            # only sets the SCAN RANGE
    kg = k_half*np.logspace(-0.6, 0.6, 90)    # stage 1: wide honest scan
    Gc = bisect_curve(kg, h, w_d, variant, iters=14)
    kctr = kg[int(np.argmin(Gc))]
    kg2 = kctr*np.linspace(0.955, 1.045, 120) # stage 2: zoom on measured basin
    Gc2 = bisect_curve(kg2, h, w_d, variant, iters=18)
    kctr2 = kg2[int(np.argmin(Gc2))]
    kg3 = kctr2*np.linspace(0.994, 1.006, 80) # stage 3: resolve the tongue core
    Gc3 = bisect_curve(kg3, h, w_d, variant, iters=22)
    j = int(np.argmin(Gc3))
    s, tr, det, gam, T = monodromy(kg3, h, Gc3*1.0002, w_d, variant)
    liou = float(np.max(np.abs(det - np.exp(-2*gam*T))/np.exp(-2*gam*T)))
    return dict(f_d=f_d, h=h, variant=variant,
                kc=float(kg3[j]), Gc=float(Gc3[j]),
                lam_mm=float(2e3*np.pi/kg3[j]),
                lam_disp_mm=float(2e3*np.pi/k_half),
                subharmonic=bool(tr[j] < 0),
                Gc_analytic_at_kc=float(Gc_analytic(kg3[j:j+1], h, variant, w_d)[0]),
                liouville_maxdev=liou,
                kgrid=kg3.tolist(), Gc_curve=Gc3.tolist(),
                kgrid_wide=kg.tolist(), Gc_wide=Gc.tolist())

if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    t0 = time.time()
    if part == "smoke":
        r = sweep_f(60.0, 0.020, "bulk")
        print(json.dumps({k_: v for k_, v in r.items() if k_ not in ("kgrid","Gc_curve")}, indent=1))
        print(f"wall={time.time()-t0:.1f}s")
    elif part == "sweep":
        fs = [10, 14, 20, 28, 40, 57, 80, 113, 160, 200]
        res = []
        for h in (0.005, 0.020):
            for variant in ("bulk", "bulk_bl"):
                for f in fs:
                    r = sweep_f(float(f), h, variant)
                    res.append(r)
                    print(f"f={f:>3} h={h*1000:.0f}mm {variant:7s}: "
                          f"a_c={r['Gc']:.3f} g, lam={r['lam_mm']:.2f} mm "
                          f"(disp {r['lam_disp_mm']:.2f}), sub={r['subharmonic']}, "
                          f"liou={r['liouville_maxdev']:.1e}", flush=True)
        json.dump(res, open(f"{OUT}/c003_sweep.json", "w"))
        print(f"wall={time.time()-t0:.1f}s")
    elif part == "demo":
        # tongue map + noise-seeded analog-style traces + figures + result.json
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        res = json.load(open(f"{OUT}/c003_sweep.json"))
        pick = [r for r in res if r["f_d"]==57 and r["h"]==0.020 and r["variant"]=="bulk"][0]
        f_d, h = pick["f_d"], pick["h"]; w_d = 2*np.pi*f_d
        kc, Gc = pick["kc"], pick["Gc"]
        # --- tongue map at f_d (fine LINEAR k-grid: razor tongues, see AUDIT dead-end)
        k_half = k_of_freq(f_d/2, h)
        kg = k_half*np.linspace(0.90, 1.85, 340)
        Gg = np.linspace(0.0, 1.0, 64)
        S = np.zeros((Gg.size, kg.size)); TR = np.zeros_like(S)
        for i, Gm in enumerate(Gg):
            s, tr, _,_,_ = monodromy(kg, h, np.full_like(kg, Gm), w_d, "bulk")
            S[i] = s; TR[i] = tr
        har_mask = (S > 0) & (TR >= 0)
        harm_thresh = float(Gg[np.where(har_mask.any(axis=1))[0][0]]) if har_mask.any() else None
        # --- noise-seeded long runs at k=kc: below / at / above threshold
        rng = np.random.default_rng(3)
        Gams = np.array([0.7*Gc, 1.0*Gc, 1.3*Gc])
        om0s = om0_sq(np.full(3, kc), h); gam = gamma_k(np.full(3, kc), h, "bulk")
        mod = Gams*G*kc*np.tanh(min(kc*h, 20.0))
        nper, spp = 400, 800
        T = 2*np.pi/w_d; dt = T/spp
        e = 1e-9*rng.standard_normal(3); v = 1e-9*rng.standard_normal(3)*w_d
        trace = np.zeros((3, nper*spp)); t = 0.0
        def rhs2(e, v, t):
            return v, -(om0s + mod*np.cos(w_d*t))*e - 2*gam*v
        for i in range(nper*spp):
            k1e,k1v = rhs2(e,v,t)
            k2e,k2v = rhs2(e+0.5*dt*k1e, v+0.5*dt*k1v, t+0.5*dt)
            k3e,k3v = rhs2(e+0.5*dt*k2e, v+0.5*dt*k2v, t+0.5*dt)
            k4e,k4v = rhs2(e+dt*k3e, v+dt*k3v, t+dt)
            e = e + dt/6*(k1e+2*k2e+2*k3e+k4e)
            v = v + dt/6*(k1v+2*k2v+2*k3v+k4v)
            t += dt
            trace[:, i] = e
        # renormalization-free growth check via log-envelope on |analytic-ish|
        tt = np.arange(nper*spp)*dt
        # spectra of the last 300 periods
        seg = trace[:, spp*100:]
        segt = tt[spp*100:]
        win = np.hanning(seg.shape[1])
        fr = np.fft.rfftfreq(seg.shape[1], dt)
        peaks = []
        for j in range(3):
            sp = np.abs(np.fft.rfft(seg[j]*win))
            pj = fr[np.argmax(sp[1:])+1]
            peaks.append(float(pj))
        # figures --------------------------------------------------------
        fig = plt.figure(figsize=(13.5, 9))
        ax = fig.add_subplot(2, 2, 1)
        for h_, variant, cst in ((0.020,"bulk","o-"), (0.020,"bulk_bl","s--"),
                                 (0.005,"bulk","^-"), (0.005,"bulk_bl","v--")):
            rr = [r for r in res if r["h"]==h_ and r["variant"]==variant]
            ax.plot([r["f_d"] for r in rr], [r["Gc"] for r in rr], cst, ms=4,
                    label=f"h={h_*1000:.0f}mm {variant}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("drive frequency f [Hz]"); ax.set_ylabel("threshold a_c [g]")
        ax.set_title("measured Faraday threshold for water\n(band between damping variants = honest uncertainty)", fontsize=9)
        ax.legend(fontsize=7)
        ax = fig.add_subplot(2, 2, 2)
        rr = [r for r in res if r["h"]==0.020 and r["variant"]=="bulk"]
        ax.plot([r["f_d"] for r in rr], [r["lam_mm"] for r in rr], "o", label="measured argmin_k of threshold")
        ax.plot([r["f_d"] for r in rr], [r["lam_disp_mm"] for r in rr], "-",
                label="dispersion at f/2 (not put in)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("drive frequency f [Hz]"); ax.set_ylabel("wavelength [mm]")
        ax.set_title("emergent pattern scale = dispersion at HALF the drive", fontsize=9)
        ax.legend(fontsize=7)
        ax = fig.add_subplot(2, 2, 3)
        for j, lab in enumerate(("0.7 a_c (below)", "1.0 a_c (neutral)", "1.3 a_c (above)")):
            ax.plot(segt if False else tt, np.log10(np.abs(trace[j])+1e-300), lw=0.7, label=lab)
        ax.set_xlabel("time [s]"); ax.set_ylabel("log10 |surface amplitude|")
        ax.set_title("analog-style: shake and watch (noise-seeded)", fontsize=9)
        ax.legend(fontsize=7)
        ax = fig.add_subplot(2, 2, 4)
        sp = np.abs(np.fft.rfft(seg[2]*win)); sp/=sp.max()
        ax.plot(fr, sp, lw=0.8)
        ax.axvline(f_d, color="gray", ls=":", label=f"drive {f_d:.0f} Hz")
        ax.axvline(f_d/2, color="r", ls="--", label=f"f/2 = {f_d/2:.1f} Hz")
        ax.set_xlim(0, 2.2*f_d); ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("spectrum (norm.)")
        ax.set_title(f"response spectrum above threshold: peak at {peaks[2]:.2f} Hz", fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(f"{OUT}/c003_threshold.png", dpi=135); plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 5))
        sub = np.ma.masked_where(~((S>0)&(TR<0)), S)
        har = np.ma.masked_where(~((S>0)&(TR>=0)), S)
        ax.pcolormesh(kg, Gg, np.ma.log10(sub+1e-12), cmap="Blues", shading="auto")
        ax.pcolormesh(kg, Gg, np.ma.log10(har+1e-12), cmap="Reds", shading="auto")
        k_full = k_of_freq(f_d, h)
        ax.axvline(k_half, color="b", ls=":", lw=1, label="k(f/2) dispersion")
        ax.axvline(k_full, color="r", ls=":", lw=1, label="k(f) dispersion")
        Gan = Gc_analytic(kg, h, "bulk", w_d)
        ax.plot(kg, np.where(Gan <= Gg.max(), Gan, np.nan), "g--", lw=1.2,
                label="route-2 analytic (weak damping)")
        ax.plot(pick["kgrid"], pick["Gc_curve"], "k-", lw=1.6, label="measured tongue tip (stage-3)")
        ax.set_ylim(0, Gg.max())
        ax.set_xlabel("wavenumber k [1/m]"); ax.set_ylabel("drive amplitude [g]")
        ax.set_title(f"instability tongues @ f={f_d:.0f} Hz, water h=20mm\nblue = subharmonic (f/2), red = harmonic (f) — structure discovered, not drawn", fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        fig.tight_layout(); fig.savefig(f"{OUT}/c003_tongues.png", dpi=135); plt.close(fig)

        out = dict(
          experiment="c003_faraday_threshold (analog-style time-domain)",
          water=dict(rho=RHO, sigma=SIG, nu=NU, note="established 20C values"),
          method="integrate real linearized surface dynamics per k over one drive period; ask if a cycle amplifies; bisect drive amplitude; minimize over k; response type from multiplier sign; noise-seeded long runs for the experiment-style demonstration",
          demo=dict(f_drive_hz=f_d, h_m=h, a_c_g=Gc, k_c=kc, lam_mm=pick["lam_mm"],
                    spectrum_peak_hz=peaks[2], f_half=f_d/2,
                    peak_over_drive=peaks[2]/f_d,
                    below_at_neutral_above=[0.7, 1.0, 1.3]),
          validation=dict(liouville_maxdev=pick["liouville_maxdev"],
                          analytic_weak_damping_at_kc=pick["Gc_analytic_at_kc"],
                          measured_Gc=pick["Gc"],
                          deviation_pct=100*abs(pick["Gc"]-pick["Gc_analytic_at_kc"])/pick["Gc"]),
          unasked=dict(harmonic_tongue_threshold_g=harm_thresh,
                       harmonic_over_subharmonic=(harm_thresh/Gc if harm_thresh else None),
                       note="first instability is subharmonic everywhere scanned; harmonic tongue opens only at much higher drive"),
          dead_end_recorded="single log k-grid (4.8% spacing) ALIASED the razor tongue (dk/k~0.5% for water) -> apparent a_c 5.7x too high; caught by route-2 disagreement; fixed by 3-stage zoom; same family as e014 denominator-2 / FFT power-of-2 traps",
          table=[{k_: r[k_] for k_ in ("f_d","h","variant","Gc","lam_mm","lam_disp_mm","subharmonic")} for r in res],
        )
        json.dump(out, open(f"{OUT}/c003_result.json", "w"), indent=1)
        print(json.dumps(out["demo"], indent=1))
        print(json.dumps(out["validation"], indent=1))
        print(f"wall={time.time()-t0:.1f}s — saved figures + c003_result.json")
