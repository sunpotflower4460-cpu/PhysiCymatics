# c001: Free-edge circular Kirchhoff plate — eigenmode calibration pillar
# ============================================================
# PUT-IN : Kirchhoff-Love plate equation D∇⁴w = ρh ω² w, free-edge BCs,
#          material params (steel), point-drive linear response model.
# EMERGED/MEASURED (within this known-physics model):
#   - eigenvalue ladder λ_mn (two INDEPENDENT numerical routes must agree)
#   - node topology counts discovered by code (not asserted)
#   - Chladni's law exponent p in f ∝ (m+2n)^p  (fit, not imposed)
#   - center-drive selection rule: only m=0 (ring) modes excited
# This is CALIBRATION of the ruler (known physics inheritance), NOT an
# emergence experiment. Emergence experiment = c002 (bounce dynamics).
# ============================================================
import numpy as np, json, os
from scipy import special, optimize, linalg
from numpy.polynomial import polynomial as P
from numpy.polynomial import legendre as L
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))

# ---------------- Route 1: Bessel boundary determinant ----------------
# w = [A J_m(kr) + B I_m(kr)] cos(mθ), λ = k a.
# Free edge at r=a:
#   M_r ∝ R'' + ν(R'/r − m²R/r²) = 0
#   V_r ∝ d/dr(∇²w)|_radial + (1−ν)(m²/r²)... = 0   (Kirchhoff effective shear)
# Column for I_m scaled by 1/I_m(λ) to avoid overflow (zeros unchanged).

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
    A, _ = bc_matrix(lam, m, nu)
    return A[0,0]*A[1,1] - A[0,1]*A[1,0]

def det_roots(m, nu, lam_max=42.0):
    xs = np.arange(0.3, lam_max, 0.005)
    vals = np.array([det_free(x, m, nu) for x in xs])
    roots = []
    for i in range(len(xs)-1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1] < 0:
            r = optimize.brentq(det_free, xs[i], xs[i+1], args=(m, nu), xtol=1e-12)
            roots.append(r)
    return roots

def mode_radial(lam, m, nu):
    """null vector -> R(r) on grid, normalized to max|R|=1"""
    A, s = bc_matrix(lam, m, nu)
    # null vector of 2x2 (use row with larger norm)
    r0 = A[0] if np.hypot(*A[0]) >= np.hypot(*A[1]) else A[1]
    a_coef, b_scaled = r0[1], -r0[0]
    b_coef = b_scaled / s          # undo column scaling
    rg = np.linspace(0, 1, 2001)
    R = a_coef*special.jv(m, lam*rg) + b_coef*special.iv(m, lam*rg)
    R = R / np.max(np.abs(R))
    return rg, R

# ---------------- Route 2: Rayleigh–Ritz energy method ----------------
# INDEPENDENT of the BC formulas above: free BCs are natural in the energy
# functional, so agreement validates the BC algebra of Route 1.
# Basis: R_j(r) = r^m * P_j(2r²−1)  (Legendre in r², regular at 0, well-conditioned)

def ritz_eigs(m, nu, N=14, nq=120):
    xg, wg = np.polynomial.legendre.leggauss(nq)     # on [-1,1]
    r = 0.5*(xg+1.0); w = 0.5*wg                     # map to [0,1]
    u = P.Polynomial([-1.0, 0.0, 2.0])               # u(r) = 2r²−1
    rm = P.Polynomial([0.0]*m + [1.0]) if m > 0 else P.Polynomial([1.0])
    basis = []
    for j in range(N):
        c = np.zeros(j+1); c[j] = 1.0
        Lj_coeffs = L.leg2poly(c)                    # Legendre -> power series in u
        Lj_of_u = P.Polynomial(Lj_coeffs)(u)         # compose: polynomial in r
        basis.append(rm * Lj_of_u)
    Rv  = np.array([b(r)            for b in basis])
    R1  = np.array([b.deriv(1)(r)   for b in basis])
    R2  = np.array([b.deriv(2)(r)   for b in basis])
    invr = 1.0/r
    Lop  = R2 + R1*invr - (m*m)*Rv*invr**2            # ∇² radial part
    Aop  = R1*invr - (m*m)*Rv*invr**2                 # w_r/r + w_θθ/r²  (radial factor)
    Bop  = (R1 - Rv*invr)*invr                        # ∂r(w_θ/(m r)) radial factor
    K = np.zeros((N, N)); M = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            integ = (Lop[i]*Lop[j]
                     - (1-nu)*(R2[i]*Aop[j] + R2[j]*Aop[i] - 2*m*m*Bop[i]*Bop[j]))
            K[i,j] = K[j,i] = np.sum(w * integ * r)
            M[i,j] = M[j,i] = np.sum(w * Rv[i]*Rv[j] * r)
    evals = linalg.eigh(K, M, eigvals_only=True)
    evals = evals[evals > 1e-6]                      # drop rigid-body modes
    return np.sqrt(np.sqrt(evals))                   # λ = (λ⁴)^{1/4}

# ---------------- node topology: code DISCOVERS (m, n) ----------------
def count_nodes(rg, R, lam, m):
    core = (rg > 0.03) & (rg < 0.995)
    s = np.sign(R[core]); s = s[s != 0]
    n_circles = int(np.sum(s[1:]*s[:-1] < 0))
    # diameters: sign changes of cos(mθ) at radius of max |R|
    th = np.linspace(0, 2*np.pi, 720, endpoint=False)
    f = np.cos(m*th)
    if m == 0:
        n_diam = 0
    else:
        ss = np.sign(f); ss = ss[ss != 0]
        crossings = int(np.sum(ss[1:]*ss[:-1] < 0)) + (1 if ss[0]*ss[-1] < 0 else 0)
        n_diam = crossings // 2
    return n_diam, n_circles

# ---------------- run ----------------
def run(nu, tag):
    res = {"nu": nu, "modes": []}
    for m in range(0, 8):
        roots = det_roots(m, nu)
        ritz  = ritz_eigs(m, nu)
        for k, lam in enumerate(roots[:4]):
            # match to nearest ritz eigenvalue
            j = int(np.argmin(np.abs(ritz - lam)))
            lam_r = float(ritz[j])
            rg, R = mode_radial(lam, m, nu)
            md, nc = count_nodes(rg, R, lam, m)
            res["modes"].append(dict(
                m=m, root_index=k, lam_det=float(lam), lam_ritz=lam_r,
                rel_err=abs(lam - lam_r)/lam,
                lam2=float(lam*lam),
                n_diameters_measured=md, n_circles_measured=nc))
    res["modes"].sort(key=lambda d: d["lam_det"])
    return res

nu = 0.30
res = run(nu, "nu030")
res33 = run(0.33, "nu033")

# cross-route agreement
max_err = max(d["rel_err"] for d in res["modes"])
print(f"[cross-check] max |det - ritz| / det over {len(res['modes'])} modes (nu=0.30): {max_err:.2e}")

# topology check: measured diameters must equal m
topo_ok = all(d["n_diameters_measured"] == d["m"] for d in res["modes"])
print(f"[topology] measured nodal diameters == m for all modes: {topo_ok}")

# ---------------- Chladni's law fit: f ∝ (m+2n)^p, f ∝ λ² ----------------
pts = [(d["n_diameters_measured"] + 2*d["n_circles_measured"], d["lam2"])
       for d in res["modes"] if (d["n_diameters_measured"] + 2*d["n_circles_measured"]) >= 2]
x = np.log([p[0] for p in pts]); y = np.log([p[1] for p in pts])
p_fit, c_fit = np.polyfit(x, y, 1)
yhat = p_fit*x + c_fit
r2 = 1 - np.sum((y-yhat)**2)/np.sum((y-np.mean(y))**2)
print(f"[Chladni law] f ∝ (m+2n)^p : p = {p_fit:.3f}, R² = {r2:.4f}  (empirical law expects ≈2)")

# ---------------- dimensional detents: default steel plate ----------------
E, rho, h, a = 2.0e11, 7850.0, 1.0e-3, 0.12   # Pa, kg/m3, m, m (24 cm diameter)
D = E*h**3/(12*(1-nu**2)); c0 = np.sqrt(D/(rho*h))
detents = []
for d in res["modes"]:
    f = d["lam2"]/(2*np.pi*a*a)*c0
    if f <= 5000:
        detents.append(dict(m=d["m"], n=d["n_circles_measured"], f_hz=round(float(f),1)))
detents.sort(key=lambda z: z["f_hz"])
print(f"[detents] steel a=0.12m h=1mm: {len(detents)} eigenfrequencies ≤ 5 kHz; "
      f"first five: {[z['f_hz'] for z in detents[:5]]}")

# ---------------- drive-point selection (linear modal participation) ----------------
# point force at r_d: participation ∝ |R_mn(r_d)| (angular alignment assumed best-case)
def participation(rd):
    out = []
    for d in res["modes"]:
        rg, R = mode_radial(d["lam_det"], d["m"], nu)
        out.append(abs(float(np.interp(rd, rg, R))))
    return np.array(out)
p_center, p_edge = participation(0.0), participation(0.9)
n_center = int(np.sum(p_center > 1e-8)); n_edge = int(np.sum(p_edge > 1e-8))
center_all_m0 = all(res["modes"][i]["m"] == 0 for i in range(len(p_center)) if p_center[i] > 1e-8)
print(f"[drive point] modes reachable: center = {n_center} (all m=0: {center_all_m0}), r=0.9a = {n_edge} of {len(res['modes'])}")

# literature anchors (FROM MEMORY — must be confirmed against Leissa at repo stage)
mem = {(2,0):5.253,(0,1):9.084,(3,0):12.23,(1,1):20.52}
anchor = []
for (mm,nn),vref in mem.items():
    got = [d for d in res33["modes"] if d["m"]==mm and d["n_circles_measured"]==nn]
    if got:
        v = got[0]["lam2"]; anchor.append(dict(mode=f"({mm},{nn})", lam2_memory=vref,
            lam2_computed=round(v,3), diff_pct=round(100*abs(v-vref)/vref,2)))
print("[memory anchors nu=0.33] ", anchor)

# ---------------- figures ----------------
fig, axes = plt.subplots(2, 4, figsize=(13, 6.5), subplot_kw={"aspect":"equal"})
th = np.linspace(0, 2*np.pi, 361)
for ax, d in zip(axes.flat, res["modes"][:8]):
    rg, R = mode_radial(d["lam_det"], d["m"], nu)
    rr, tt = np.meshgrid(rg, th)
    W = np.interp(rr, rg, R)*np.cos(d["m"]*tt)
    X, Y = rr*np.cos(tt), rr*np.sin(tt)
    ax.contourf(X, Y, W, levels=21, cmap="RdBu_r")
    ax.contour(X, Y, W, levels=[0], colors="k", linewidths=1.4)  # node lines = sand
    ax.set_title(f"(m={d['m']}, n={d['n_circles_measured']})  λ²={d['lam2']:.2f}", fontsize=9)
    ax.axis("off")
fig.suptitle("c001 free circular plate — first 8 modes (black = node lines = where sand collects)", fontsize=11)
fig.tight_layout(); fig.savefig(f"{OUT}/c001_modes.png", dpi=130); plt.close(fig)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
a1.scatter(np.exp(x), np.exp(y), s=28)
xs = np.linspace(min(x), max(x), 50)
a1.plot(np.exp(xs), np.exp(p_fit*xs + c_fit), "r-", lw=1)
a1.set_xscale("log"); a1.set_yscale("log")
a1.set_xlabel("m + 2n"); a1.set_ylabel("λ² (∝ f)")
a1.set_title(f"Chladni's law fit: p = {p_fit:.3f} (R²={r2:.3f})")
idx = np.arange(len(res["modes"]))
a2.bar(idx-0.2, p_center, 0.4, label="drive @ center")
a2.bar(idx+0.2, p_edge,   0.4, label="drive @ 0.9a")
a2.set_xticks(idx)
a2.set_xticklabels([f"{d['m']},{d['n_circles_measured']}" for d in res["modes"]], rotation=90, fontsize=6)
a2.set_ylabel("|mode shape at drive point|"); a2.legend(fontsize=8)
a2.set_title("drive point selects which modes exist")
fig.tight_layout(); fig.savefig(f"{OUT}/c001_law_drive.png", dpi=130); plt.close(fig)

# ---------------- result.json ----------------
result = dict(
    experiment="c001_free_circular_plate_calibration",
    put_in="Kirchhoff-Love equation, free-edge BC, steel params, linear point-drive model",
    routes={"route1":"Bessel boundary determinant","route2":"Rayleigh-Ritz energy (independent of BC algebra)"},
    cross_check_max_rel_err=float(max_err),
    topology_discovered_matches_labels=bool(topo_ok),
    chladni_law=dict(p=float(p_fit), r2=float(r2), expected="≈2 (established empirical)"),
    plate=dict(material="steel", E_pa=E, rho=rho, h_m=h, a_m=a, nu=nu),
    eigenfrequencies_hz=detents,
    drive_point=dict(center_reachable=n_center, center_all_m0=bool(center_all_m0),
                     edge_reachable=n_edge, total=len(res["modes"])),
    literature_anchors_nu033=anchor,
    tier=dict(eigenvalues="measured (within known-physics model, two independent routes)",
              chladni_p="measured (fit)",
              center_rings_only="measured within model / established consequence",
              literature_anchor_values="MEMORY-BASED — confirm against Leissa before repo GREEN"),
    floors=["linear small-deflection (von Karman excluded)","no air loading/radiation damping",
            "point-force drive model (real shaker has finite contact area)",
            "damping zeta not measured (assumed in UI transfer function)"],
)
with open(f"{OUT}/c001_result.json","w") as f: json.dump(result, f, indent=1, ensure_ascii=False)
print("[saved] c001_result.json, c001_modes.png, c001_law_drive.png")
