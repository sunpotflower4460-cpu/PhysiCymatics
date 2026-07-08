"""reference/faraday_linear.py — c003 Faraday threshold, tidy importable core.

PUT-IN     : linearized free-surface dynamics per wavenumber k under a
             vertically oscillating gravity g(t)=g(1+Γ cos ω_d t):
                 η_tt + 2γ_k η_t + [(g(t)k + σk³/ρ) tanh(kh)] η = 0
             (Benjamin–Ursell 1954 structure; established derivation).
             water @20°C: ρ=998, σ=0.0728, ν=1.004e-6 (established).
             damping submodel is the STATED simplification (two variants →
             honesty band): bulk γ=2νk², optional +bottom boundary layer.
FORBIDDEN  : no "f/2", no threshold value, no selected wavelength is put in.
             We integrate ONE drive period and ask the real equation whether a
             cycle amplifies (monodromy = the analog question in minimal form).
EMERGED    : threshold Γ_c by bisection; response type from the Floquet
             multiplier SIGN (negative trace ⇒ subharmonic, discovered);
             pattern scale = argmin_k of the measured threshold.
CLAIM-TIER : measured (within the stated linear/damping model). Cross-checks:
             Liouville det(M)=exp(−2γT) (integrator truth) and a weak-damping
             analytic threshold (route 2).
FLOORS     : linearized (no saturation/pattern selection); damping is a
             sub-model, not first-principles; boundary-layer coefficient is
             memory-based (YELLOW until source-checked).

Faithful refactor of oracles/c003_faraday/c003.py (linear/threshold part only;
the nonlinear demo + figures stay in the oracle).
"""
import numpy as np
from scipy import optimize

G, RHO, SIG, NU = 9.81, 998.0, 0.0728, 1.004e-6   # water @ ~20°C (established)


def om0_sq(k, h):
    """Undamped dispersion ω₀²(k) = (gk + σk³/ρ) tanh(kh)."""
    return (G * k + SIG * k**3 / RHO) * np.tanh(np.minimum(k * h, 20.0))


def gamma_k(k, h, variant="bulk"):
    g0 = 2.0 * NU * k * k
    if variant == "bulk":
        return g0
    kh = np.minimum(k * h, 20.0)
    om0 = np.sqrt(om0_sq(k, h))
    return g0 + k * np.sqrt(NU * om0 / 2.0) / np.sinh(2.0 * kh)   # +bottom BL


def monodromy(k, h, Gam, w_d, variant="bulk", nsteps=1200):
    """Integrate the real equation over ONE drive period for two ICs (RK4,
    vectorized in k). Returns growth rate s, trace, det, γ, T."""
    om0s = om0_sq(k, h)
    gam = gamma_k(k, h, variant)
    mod = Gam * G * k * np.tanh(np.minimum(k * h, 20.0))   # only gravity modulated
    T = 2 * np.pi / w_d
    dt = T / nsteps
    E = np.stack([np.ones_like(k), np.zeros_like(k)])
    V = np.stack([np.zeros_like(k), np.ones_like(k)])

    def rhs(E, V, t):
        return V, -(om0s + mod * np.cos(w_d * t)) * E - 2.0 * gam * V

    t = 0.0
    for _ in range(nsteps):
        k1e, k1v = rhs(E, V, t)
        k2e, k2v = rhs(E + 0.5 * dt * k1e, V + 0.5 * dt * k1v, t + 0.5 * dt)
        k3e, k3v = rhs(E + 0.5 * dt * k2e, V + 0.5 * dt * k2v, t + 0.5 * dt)
        k4e, k4v = rhs(E + dt * k3e, V + dt * k3v, t + dt)
        E = E + dt / 6.0 * (k1e + 2 * k2e + 2 * k3e + k4e)
        V = V + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        t += dt
    tr = E[0] + V[1]
    det = E[0] * V[1] - E[1] * V[0]
    disc = tr * tr - 4.0 * det
    lam = np.where(disc >= 0, (np.abs(tr) + np.sqrt(np.maximum(disc, 0))) / 2.0,
                   np.sqrt(np.maximum(det, 1e-300)))
    s = np.log(np.maximum(lam, 1e-300)) / T
    return s, tr, det, gam, T


def Gc_analytic(k, h, w_d, variant="bulk"):
    """Route 2: weak-damping subharmonic threshold with detuning,
       M_c = 2 ω_d sqrt(γ² + δ²), δ = ω₀(k) − ω_d/2."""
    om0 = np.sqrt(om0_sq(k, h))
    gam = gamma_k(k, h, variant)
    delta = om0 - 0.5 * w_d
    Mc = 2.0 * w_d * np.sqrt(gam * gam + delta * delta)
    return Mc / (G * k * np.tanh(np.minimum(k * h, 20.0)))


def k_of_freq(f_wave, h):
    w = 2 * np.pi * f_wave
    return optimize.brentq(lambda k: om0_sq(k, h) - w * w, 1e-2, 1e6)


def _bisect_curve(kg, h, w_d, variant, iters=20):
    Ghi = 1.6 * Gc_analytic(kg, h, w_d, variant) + 1e-4
    for _ in range(8):
        s, *_ = monodromy(kg, h, Ghi, w_d, variant)
        need = s <= 0
        if not need.any():
            break
        Ghi[need] *= 1.7
    Glo = np.zeros_like(kg)
    for _ in range(iters):
        Gm = 0.5 * (Glo + Ghi)
        s, *_ = monodromy(kg, h, Gm, w_d, variant)
        up = s > 0
        Ghi = np.where(up, Gm, Ghi)
        Glo = np.where(up, Glo, Gm)
    return 0.5 * (Glo + Ghi)


def sweep_f(f_d, h, variant="bulk"):
    """Three-stage k refinement. The water tongue is razor-thin (Δk/k~0.5%);
    a single log grid ALIASES it (recorded trap) — so we zoom in three stages."""
    w_d = 2 * np.pi * f_d
    k_half = k_of_freq(f_d / 2.0, h)              # sets the SCAN RANGE only
    kg = k_half * np.logspace(-0.6, 0.6, 90)
    Gc = _bisect_curve(kg, h, w_d, variant, iters=14)
    kctr = kg[int(np.argmin(Gc))]
    kg2 = kctr * np.linspace(0.955, 1.045, 120)
    Gc2 = _bisect_curve(kg2, h, w_d, variant, iters=18)
    kctr2 = kg2[int(np.argmin(Gc2))]
    kg3 = kctr2 * np.linspace(0.994, 1.006, 80)
    Gc3 = _bisect_curve(kg3, h, w_d, variant, iters=22)
    j = int(np.argmin(Gc3))
    s, tr, det, gam, T = monodromy(kg3, h, Gc3 * 1.0002, w_d, variant)
    liou = float(np.max(np.abs(det - np.exp(-2 * gam * T)) / np.exp(-2 * gam * T)))
    return dict(
        f_d=f_d, h=h, variant=variant,
        kc=float(kg3[j]), Gc=float(Gc3[j]),
        a_c_ms2=float(Gc3[j] * G),
        lam_mm=float(2e3 * np.pi / kg3[j]),
        lam_disp_mm=float(2e3 * np.pi / k_half),
        subharmonic=bool(tr[j] < 0),
        Gc_analytic_at_kc=float(Gc_analytic(kg3[j:j + 1], h, w_d, variant)[0]),
        liouville_maxdev=liou,
    )


if __name__ == "__main__":
    r = sweep_f(60.0, 0.020, "bulk")
    for kk in ("kc", "Gc", "a_c_ms2", "lam_mm", "lam_disp_mm", "subharmonic",
               "Gc_analytic_at_kc", "liouville_maxdev"):
        print(f"  {kk:18s} = {r[kk]}")
