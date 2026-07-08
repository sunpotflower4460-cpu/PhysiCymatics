"""reference/zhang_vinals.py — c005 Zhang–Viñals PDE, linear-onset skeleton.

PUT-IN     : the Zhang & Viñals (JFM 336:301, 1997) quasi-potential equations
             for a weakly-viscous semi-infinite fluid, driven at f_drive:
                 (∂t − γ∇²) h   = D[Φ] + F(h,Φ)
                 (∂t − γ∇²) Φ   = (Γ₀∇² − G(t)) h + G(h,Φ)
             scalings eq (24); water @440 Hz → k_c=2971/m, γ=0.0064, Γ₀=0.2462.
FORBIDDEN  : no pattern named, no angular seeding (isotropic ring-band noise IC),
             no term that knows about 60° or 90°.
EMERGED    : (linear part, exercised here) growth/decay of the ring band; its
             SIGN must flip exactly at the c003 Faraday threshold — this is the
             route-2 self-validation of the PDE solver.
CLAIM-TIER : measured (within the quasi-potential model). The full nonlinear
             pattern selection (which angle survives) is the FLAGSHIP run and is
             heavy — kept in oracles/c005_shape/c005.py and exercised only by the
             `slow` nightly test, not the fast CI.
FLOORS     : weakly-viscous quasi-potential approx; periodic square domain
             quantizes representable angles; linear onset here does not resolve
             saturated pattern competition.

This module refactors the LINEAR onset path of the oracle into importable
functions. The nonlinear flagship stays in the oracle by design (see CLAUDE.md
Phase-1/Phase-4 split).
"""
import numpy as np

try:
    from scipy import fft as _sfft
    _F2 = lambda a: _sfft.rfft2(a, workers=2)
    _IF2 = lambda A, n: _sfft.irfft2(A, s=(n, n), workers=2)
except Exception:                                    # pragma: no cover
    _F2 = np.fft.rfft2
    _IF2 = lambda A, n: np.fft.irfft2(A, s=(n, n))

G_SI, RHO, SIG, NU = 9.81, 998.0, 0.0728, 1.004e-6


def setup(f_drive=440.0, N=128, R=18.0):
    """Build the spectral operators. k_c=1 sits at lattice radius R."""
    from scipy import optimize
    w = 2 * np.pi * f_drive
    kc = optimize.brentq(lambda k: (G_SI * k + SIG * k**3 / RHO) - (w / 2)**2,
                         1e-2, 1e6)
    gam = 2 * NU * kc * kc / w
    G0 = G_SI * kc / w**2
    Gam0 = 0.25 - G0
    ki = np.fft.fftfreq(N, d=1.0 / N)
    kj = np.fft.rfftfreq(N, d=1.0 / N)
    KX = ki[:, None] / R
    KY = kj[None, :] / R
    K2 = KX * KX + KY * KY
    KABS = np.sqrt(K2)
    deal = (KABS <= (N // 2) * (2.0 / 3.0) / R)
    return dict(f=f_drive, w=w, kc=kc, gam=gam, G0=G0, Gam0=Gam0, N=N, R=R,
                KX=KX, KY=KY, K2=K2, KABS=KABS, deal=deal)


def _rhs_linear(hh, ph, t, P, fdrive):
    """Linear (Mathieu-per-k) part only: no term knows any angle."""
    dh = P["KABS"] * ph
    dp = -(P["Gam0"] * P["K2"] + P["G0"] - fdrive * np.cos(t)) * hh
    return P["deal"] * dh, P["deal"] * dp


def _step(hh, ph, t, dt, P, fdrive, E, E2):
    k1h, k1p = _rhs_linear(hh, ph, t, P, fdrive)
    k2h, k2p = _rhs_linear(E2 * (hh + 0.5 * dt * k1h), E2 * (ph + 0.5 * dt * k1p),
                           t + 0.5 * dt, P, fdrive)
    k3h, k3p = _rhs_linear(E2 * hh + 0.5 * dt * k2h, E2 * ph + 0.5 * dt * k2p,
                           t + 0.5 * dt, P, fdrive)
    k4h, k4p = _rhs_linear(E * hh + dt * E2 * k3h, E * ph + dt * E2 * k3p,
                           t + dt, P, fdrive)
    hh = E * hh + dt / 6.0 * (E * k1h + 2 * E2 * (k2h + k3h) + k4h)
    ph = E * ph + dt / 6.0 * (E * k1p + 2 * E2 * (k2p + k3p) + k4p)
    hh[0, 0] = 0.0
    ph[0, 0] = 0.0
    return hh, ph


def linear_growth_rate(P, fdrive, n_periods=60, spp=48, seed=7, amp0=1e-6):
    """Integrate the LINEAR PDE from an isotropic ring-band seed; return the
    measured growth rate (per unit scaled time) of the k_c=1 band."""
    N = P["N"]
    dt = 2 * np.pi / spp
    E = np.exp(-P["gam"] * P["K2"] * dt)
    E2 = np.exp(-P["gam"] * P["K2"] * dt * 0.5)
    rng = np.random.default_rng(seed)
    ki = np.fft.fftfreq(N, d=1.0 / N)
    kj = np.fft.rfftfreq(N, d=1.0 / N)
    I2 = ki[:, None]**2 + kj[None, :]**2
    band = (I2 >= 317) & (I2 <= 334)              # ring shell around |k|=1
    hh = _F2(amp0 * rng.standard_normal((N, N))) * P["deal"]
    hh[0, 0] = 0
    hh = hh * band
    ph = np.zeros_like(hh)
    t = 0.0
    log = []
    for step in range(1, n_periods * spp + 1):
        hh, ph = _step(hh, ph, t, dt, P, fdrive, E, E2)
        t = step * dt
        if step % spp == 0:
            log.append(np.sqrt(float(np.sum((np.abs(hh)**2)[band]))))
    log = np.array(log)
    return np.log(log[-1] / log[19]) / (2 * np.pi * (len(log) - 20))


def leading_threshold_scaled(P):
    """Scaled linear onset amplitude f_c = 2γ (leading Mathieu tongue)."""
    return 2 * P["gam"]


if __name__ == "__main__":
    P = setup(440.0)
    fc = leading_threshold_scaled(P)
    print(f"setup: kc={P['kc']:.1f}  gamma={P['gam']:.5f}  Gam0={P['Gam0']:.4f}  "
          f"f_c(leading)={fc:.5f}")
    for fac in (0.9, 1.0, 1.1):
        g = linear_growth_rate(P, fac * fc)
        print(f"  f/f_c={fac}: growth={g:+.5f} "
              f"(expect {'+' if fac > 1 else '-' if fac < 1 else '~0'})")
