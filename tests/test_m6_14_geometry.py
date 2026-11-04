"""M6-14A: candidate evaluation for a geometrically consistent plate
surface. The implemented height uses piecewise-constant radial bins
while contact gradients use sampled eigenfunction derivatives — two
different surfaces. Candidates are compared against the analytic
Bessel eigenfunction (rebuilt here from sand_bounce primitives, the
independent definition the tables were sampled from).

Candidate B (chosen for the TS research path): cubic Hermite
interpolation on the nodal grid using the tabulated values AND
slopes; its gradient is by construction the derivative of its own
height. This test quantifies B vs the analytic truth and vs the
current bin scheme, and prints the selection evidence.
"""
import numpy as np
import pytest
from scipy import special

from reference import plate_transient as pt
from reference import sand_bounce as sb


@pytest.fixture(scope="module")
def basis():
    return pt.build_basis()


def _analytic_R(md, rho):
    """Analytic normalized eigenfunction and d/drho at arbitrary rho,
    rebuilt exactly as sand_bounce.make_mode does (same lam/coeffs/
    normalization), but evaluated off-grid."""
    m = md["m"]
    # Recover root index by matching f_hz is fragile; instead rebuild
    # from the mode's own table metadata via make_mode signature by
    # scanning root indices for a frequency match.
    nu = 0.30
    roots = sb.det_roots(m, nu)
    best = None
    for ri, lam in enumerate(roots):
        probe = sb.make_mode(m, ri, ngrid=8)
        if abs(probe["f_hz"] - md["f_hz"]) < 1e-6 * md["f_hz"]:
            best = (lam, probe)
            break
    assert best is not None
    lam, _ = best
    A, s = sb.bc_matrix(lam, m, nu)
    r0 = A[0] if np.hypot(*A[0]) >= np.hypot(*A[1]) else A[1]
    ac, bs = r0[1], -r0[0]
    bc = bs / s
    Rg = ac * special.jv(m, lam * np.linspace(0, 1, 512)) + bc * special.iv(m, lam * np.linspace(0, 1, 512))
    nrm = np.max(np.abs(Rg))
    R = (ac * special.jv(m, lam * rho) + bc * special.iv(m, lam * rho)) / nrm
    dR = (ac * lam * special.jvp(m, lam * rho, 1) + bc * lam * special.ivp(m, lam * rho, 1)) / nrm
    return R, dR


def _hermite_interp(md, rho):
    """Candidate B: cubic Hermite on nodal grid (values+slopes)."""
    R = np.asarray(md["R"]); dR = np.asarray(md["dR"])
    n = len(R) - 1
    x = np.clip(rho, 0.0, 1.0) * n
    k = np.minimum(np.floor(x).astype(int), n - 1)
    u = x - k
    drho = 1.0 / n
    h00 = 2*u**3 - 3*u**2 + 1; h10 = u**3 - 2*u**2 + u
    h01 = -2*u**3 + 3*u**2;    h11 = u**3 - u**2
    val = h00*R[k] + h10*drho*dR[k] + h01*R[k+1] + h11*drho*dR[k+1]
    d00 = 6*u**2 - 6*u; d10 = 3*u**2 - 4*u + 1
    d01 = -6*u**2 + 6*u; d11 = 3*u**2 - 2*u
    der = (d00*R[k] + d10*drho*dR[k] + d01*R[k+1] + d11*drho*dR[k+1]) / drho
    dd00 = 12*u - 6; dd10 = 6*u - 4; dd01 = -12*u + 6; dd11 = 6*u - 2
    dder = (dd00*R[k] + dd10*drho*dR[k] + dd01*R[k+1] + dd11*drho*dR[k+1]) / drho**2
    return val, der, dder


def _bin_scheme(md, rho):
    """Candidate A (current): piecewise-constant bin lookup."""
    R = np.asarray(md["R"]); dR = np.asarray(md["dR"])
    n = len(R) - 1
    idx = np.minimum((np.clip(rho, 0, 1) * n).astype(int), n)
    return R[idx], dR[idx]


def test_geometry_candidates_vs_analytic(basis):
    rho = np.linspace(0, 1, 20001)
    rows = []
    for md in basis["modes"]:
        Rt, dRt = _analytic_R(md, rho)
        # Sanity: analytic at nodes reproduces the table.
        nodes = np.linspace(0, 1, len(md["R"]))
        Rn, _ = _analytic_R(md, nodes)
        assert np.max(np.abs(Rn - np.asarray(md["R"]))) < 1e-9
        Ra, dRa = _bin_scheme(md, rho)
        Rb, dRb, _ = _hermite_interp(md, rho)
        rows.append((
            md["m"], md["f_hz"],
            float(np.max(np.abs(Ra - Rt))), float(np.max(np.abs(Rb - Rt))),
            float(np.max(np.abs(dRa - dRt))), float(np.max(np.abs(dRb - dRt))),
        ))
    for r in rows:
        print(f"m={r[0]} f={r[1]:7.1f}  maxErr R: bin={r[2]:.3e} hermite={r[3]:.3e}   "
              f"maxErr dR: bin={r[4]:.3e} hermite={r[5]:.3e}")
    # Selection criteria: Hermite must be dramatically closer to the
    # analytic eigenfunction in both value and derivative.
    for r in rows:
        assert r[3] < r[2] / 20, r
        assert r[5] < r[4] / 5, r
        assert r[3] < 1e-5, r


def test_hermite_is_c1_continuous_and_self_consistent(basis):
    md = basis["modes"][3]
    n = len(md["R"]) - 1
    eps = 1e-9
    max_jump = 0.0; max_djump = 0.0
    for k in range(1, n):
        rho_e = k / n
        vl = _hermite_interp(md, np.array([rho_e - eps]))
        vr = _hermite_interp(md, np.array([rho_e + eps]))
        max_jump = max(max_jump, abs(float(vr[0][0] - vl[0][0])))
        max_djump = max(max_djump, abs(float(vr[1][0] - vl[1][0])))
    assert max_jump < 1e-7
    assert max_djump < 1e-4
    # Derivative of the interpolation == FD of the interpolation.
    rho = np.linspace(0.001, 0.999, 500)
    _, der, _ = _hermite_interp(md, rho)
    d = 1e-6
    fd = (_hermite_interp(md, rho + d)[0] - _hermite_interp(md, rho - d)[0]) / (2 * d)
    assert np.max(np.abs(fd - der)) < 5e-3
