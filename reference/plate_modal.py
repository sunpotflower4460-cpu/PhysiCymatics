"""reference/plate_modal.py — c001 free-edge circular plate, tidy importable core.

PUT-IN     : Kirchhoff-Love plate equation D∇⁴w = ρh ω² w, free-edge boundary
             conditions, isotropic material (steel default), linear point-drive
             modal-participation model.
EMERGED    : eigenvalue ladder λ_mn (TWO independent numerical routes agree),
             nodal topology counts (discovered, not asserted), Chladni exponent
             p in f ∝ (m+2n)^p (fit), center-drive selection rule (only m=0 rings).
CLAIM-TIER : measured (within a known-physics model; two independent routes).
             This CALIBRATES the ruler — it is not an emergence experiment.
FLOORS     : linear small-deflection (von Kármán excluded); no air loading /
             radiation damping; point-force drive (real shaker has finite
             contact area); damping ζ not measured here.

This is a faithful refactor of oracles/c001_plate/c001.py: same physics, same
two routes, exposed as importable functions with a vectorized determinant scan
so the pinned tests stay fast. Nothing about node targets or pattern names lives
here — the topology is only ever *measured* from the computed mode shape.
"""
import numpy as np
from scipy import special, optimize, linalg
from numpy.polynomial import polynomial as P
from numpy.polynomial import legendre as L

# ---- default steel plate (24 cm diameter, 1 mm thick) ----------------------
STEEL = dict(E_pa=2.0e11, rho=7850.0, h_m=1.0e-3, a_m=0.12, nu=0.30)


# ---------------- Route 1: Bessel boundary determinant ----------------------
# w = [A J_m(kr) + B I_m(kr)] cos(mθ), λ = k a. Free edge at r=a:
#   M_r  ∝ R'' + ν(R'/r − m²R/r²) = 0
#   V_r  ∝ Kirchhoff effective shear = 0
# The I_m column is scaled by 1/I_m(λ) to avoid overflow (zeros unchanged).
def _bc_matrix(lam, m, nu):
    x = lam
    Jm, Jp, Jpp = special.jv(m, x), special.jvp(m, x, 1), special.jvp(m, x, 2)
    Im, Ip, Ipp = special.iv(m, x), special.ivp(m, x, 1), special.ivp(m, x, 2)
    s = Im if Im != 0 else 1.0
    M_J = x * x * Jpp + nu * (x * Jp - m * m * Jm)
    M_I = (x * x * Ipp + nu * (x * Ip - m * m * Im)) / s
    V_J = -x**3 * Jp - (1 - nu) * m * m * (x * Jp - Jm)
    V_I = (x**3 * Ip - (1 - nu) * m * m * (x * Ip - Im)) / s
    return np.array([[M_J, M_I], [V_J, V_I]]), s


def _det_free_vec(x, m, nu):
    """Vectorized free-edge determinant over an array of λ (fast grid scan)."""
    Jm, Jp, Jpp = special.jv(m, x), special.jvp(m, x, 1), special.jvp(m, x, 2)
    Im, Ip, Ipp = special.iv(m, x), special.ivp(m, x, 1), special.ivp(m, x, 2)
    s = np.where(Im != 0, Im, 1.0)
    M_J = x * x * Jpp + nu * (x * Jp - m * m * Jm)
    M_I = (x * x * Ipp + nu * (x * Ip - m * m * Im)) / s
    V_J = -x**3 * Jp - (1 - nu) * m * m * (x * Jp - Jm)
    V_I = (x**3 * Ip - (1 - nu) * m * m * (x * Ip - Im)) / s
    return M_J * V_I - M_I * V_J


def _det_free(lam, m, nu):
    A, _ = _bc_matrix(lam, m, nu)
    return A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]


def det_roots(m, nu, lam_max=42.0, step=0.005):
    """λ ladder for angular index m via sign-change bracketing + brentq."""
    xs = np.arange(0.3, lam_max, step)
    vals = _det_free_vec(xs, m, nu)
    roots = []
    fin = np.isfinite(vals)
    for i in range(len(xs) - 1):
        if fin[i] and fin[i + 1] and vals[i] * vals[i + 1] < 0:
            roots.append(optimize.brentq(_det_free, xs[i], xs[i + 1],
                                         args=(m, nu), xtol=1e-12))
    return roots


def mode_radial(lam, m, nu, n=2001):
    """Null vector → radial mode shape R(r) on [0,1], normalized max|R|=1."""
    A, s = _bc_matrix(lam, m, nu)
    r0 = A[0] if np.hypot(*A[0]) >= np.hypot(*A[1]) else A[1]
    a_coef, b_coef = r0[1], -r0[0] / s
    rg = np.linspace(0, 1, n)
    R = a_coef * special.jv(m, lam * rg) + b_coef * special.iv(m, lam * rg)
    return rg, R / np.max(np.abs(R))


# ---------------- Route 2: Rayleigh–Ritz energy method ----------------------
# Independent of the BC algebra above (free BCs are natural in the energy
# functional), so agreement validates Route 1. Basis R_j = r^m P_j(2r²−1).
def ritz_eigs(m, nu, N=14, nq=120):
    xg, wg = np.polynomial.legendre.leggauss(nq)
    r = 0.5 * (xg + 1.0)
    w = 0.5 * wg
    u = P.Polynomial([-1.0, 0.0, 2.0])
    rm = P.Polynomial([0.0] * m + [1.0]) if m > 0 else P.Polynomial([1.0])
    basis = []
    for j in range(N):
        c = np.zeros(j + 1); c[j] = 1.0
        basis.append(rm * P.Polynomial(L.leg2poly(c))(u))
    Rv = np.array([b(r) for b in basis])
    R1 = np.array([b.deriv(1)(r) for b in basis])
    R2 = np.array([b.deriv(2)(r) for b in basis])
    invr = 1.0 / r
    Lop = R2 + R1 * invr - (m * m) * Rv * invr**2
    Aop = R1 * invr - (m * m) * Rv * invr**2
    Bop = (R1 - Rv * invr) * invr
    K = np.zeros((N, N)); M = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            integ = (Lop[i] * Lop[j]
                     - (1 - nu) * (R2[i] * Aop[j] + R2[j] * Aop[i]
                                   - 2 * m * m * Bop[i] * Bop[j]))
            K[i, j] = K[j, i] = np.sum(w * integ * r)
            M[i, j] = M[j, i] = np.sum(w * Rv[i] * Rv[j] * r)
    evals = linalg.eigh(K, M, eigvals_only=True)
    evals = evals[evals > 1e-6]
    return np.sqrt(np.sqrt(evals))


# ---------------- node topology: code DISCOVERS (m, n) ----------------------
def count_nodes(rg, R, m):
    core = (rg > 0.03) & (rg < 0.995)
    s = np.sign(R[core]); s = s[s != 0]
    n_circles = int(np.sum(s[1:] * s[:-1] < 0))
    if m == 0:
        n_diam = 0
    else:
        th = np.linspace(0, 2 * np.pi, 720, endpoint=False)
        ss = np.sign(np.cos(m * th)); ss = ss[ss != 0]
        crossings = int(np.sum(ss[1:] * ss[:-1] < 0)) + (1 if ss[0] * ss[-1] < 0 else 0)
        n_diam = crossings // 2
    return n_diam, n_circles


# ---------------- assemble the mode table (both routes) ---------------------
def mode_table(nu=0.30, m_max=8, roots_per_m=4):
    modes = []
    for m in range(0, m_max):
        roots = det_roots(m, nu)
        ritz = ritz_eigs(m, nu)
        for k, lam in enumerate(roots[:roots_per_m]):
            j = int(np.argmin(np.abs(ritz - lam)))
            md, nc = count_nodes(*mode_radial(lam, m, nu), m)
            modes.append(dict(m=m, root_index=k, lam_det=float(lam),
                              lam_ritz=float(ritz[j]),
                              rel_err=abs(lam - ritz[j]) / lam,
                              lam2=float(lam * lam),
                              n_diameters_measured=md, n_circles_measured=nc))
    modes.sort(key=lambda d: d["lam_det"])
    return modes


def eigenfrequencies_hz(modes, plate=STEEL, f_max=5000.0):
    """Dimensional detents f = λ²/(2π a²) · sqrt(D/ρh) for a given plate."""
    E, rho, h, a, nu = (plate["E_pa"], plate["rho"], plate["h_m"],
                        plate["a_m"], plate["nu"])
    D = E * h**3 / (12 * (1 - nu**2))
    c0 = np.sqrt(D / (rho * h))
    out = []
    for d in modes:
        f = d["lam2"] / (2 * np.pi * a * a) * c0
        if f <= f_max:
            out.append(dict(m=d["m"], n=d["n_circles_measured"], f_hz=round(float(f), 1)))
    out.sort(key=lambda z: z["f_hz"])
    return out


def chladni_exponent(modes):
    """Fit f ∝ (m+2n)^p through λ² vs (m+2n); returns (p, R²)."""
    pts = [(d["n_diameters_measured"] + 2 * d["n_circles_measured"], d["lam2"])
           for d in modes
           if (d["n_diameters_measured"] + 2 * d["n_circles_measured"]) >= 2]
    x = np.log([p[0] for p in pts]); y = np.log([p[1] for p in pts])
    p_fit, c_fit = np.polyfit(x, y, 1)
    yhat = p_fit * x + c_fit
    r2 = 1 - np.sum((y - yhat)**2) / np.sum((y - np.mean(y))**2)
    return float(p_fit), float(r2)


def center_drive_reachable_hz(modes, plate=STEEL, f_max=5000.0):
    """Modes a CENTER point drive can excite = m==0 rings only (selection rule)."""
    m0 = [d for d in modes if d["m"] == 0]
    return [z["f_hz"] for z in eigenfrequencies_hz(m0, plate, f_max)]


if __name__ == "__main__":
    modes = mode_table(0.30)
    max_err = max(d["rel_err"] for d in modes)
    p, r2 = chladni_exponent(modes)
    reach = center_drive_reachable_hz(modes)
    print(f"modes={len(modes)}  cross-route max rel err={max_err:.2e}")
    print(f"Chladni p={p:.3f} (R²={r2:.4f})")
    print(f"center-drive reachable Hz (first 4): {reach[:4]}")
