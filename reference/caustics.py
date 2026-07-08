"""reference/caustics.py — c006 water-surface caustics, tidy importable core.

PUT-IN     : geometric ray optics. A collimated beam travels downward from air
             (n=1) into water (n=1.333) through a sinusoidal surface
             h(x)=A cos(kx), λ=9 mm. Vector Snell refraction at the interface.
EMERGED    : near-axial focal depth (rays crossing the optical axis), the
             paraxial anchor f = nR/(n−1) with crest curvature R = 1/(A k²),
             and the SIGN of spherical aberration (measured from the ray fan).
CLAIM-TIER : measured (geometric-optics model). The near-axial crossing depth
             agrees with nR/(n−1) to ~0% — this is the calibration line.
FLOORS     : geometric optics only (no diffraction / wave caustic structure);
             single interface; monochromatic; 2-D cross-section.

Faithful refactor of the c006 oracle's ray-tracing core. Two recorded traps are
respected here: (1) the deflection sign — light must go air→water (downward),
crest → converging; the water→air / crest→trough mirror focuses at the wrong
depth (~+25%). (2) near-axial focus is measured as depth BELOW THE CREST VERTEX
(z=A), extrapolated slope→0 — not argmax-of-intensity, which is biased by the
fold caustics that out-shine the cusp.
"""
import numpy as np

N_WATER = 1.333
LAM_M = 9.0e-3
K = 2 * np.pi / LAM_M


def crest_radius(A):
    """Radius of curvature of h(x)=A cos(kx) at the crest x=0: R = 1/(A k²)."""
    return 1.0 / (A * K * K)


def paraxial_focus_m(A, n=N_WATER):
    """Single-surface paraxial focus below the crest: f = nR/(n−1)."""
    return n * crest_radius(A) / (n - 1.0)


def _snell(d, Nrm, r):
    """Vector Snell refraction. d,Nrm unit; r = n_incident / n_transmit."""
    Nrm = np.asarray(Nrm, float)
    cosi = -np.dot(d, Nrm)
    if cosi < 0:
        Nrm = -Nrm
        cosi = -np.dot(d, Nrm)
    sin2t = r * r * (1.0 - cosi * cosi)
    if sin2t > 1.0:
        return None                      # total internal reflection
    t = r * d + (r * cosi - np.sqrt(1.0 - sin2t)) * Nrm
    return t / np.hypot(*t)


def axis_crossing_z(x, A, n=N_WATER):
    """z-coordinate where the downward ray entering at horizontal x crosses the
    optical axis (X=0). Light goes air→water (r = 1/n)."""
    h = A * np.cos(K * x)
    hp = -A * K * np.sin(K * x)
    Nrm = np.array([-hp, 1.0])
    Nrm = Nrm / np.hypot(*Nrm)                    # upward surface normal
    t = _snell(np.array([0.0, -1.0]), Nrm, 1.0 / n)
    if t is None or t[0] == 0:
        return np.nan
    s = -x / t[0]
    return h + t[1] * s


def near_axial_focus_m(A, n=N_WATER, slope_max=1e-3):
    """Measured near-axial focal depth below the crest, by extrapolating the
    ray crossing z→ x=0 (quadratic in x). Should match paraxial_focus_m ~0%."""
    x_at = slope_max / (A * K * K)                # x where |slope| ≈ slope_max
    xs = np.linspace(0.2, 1.0, 5) * x_at
    z = np.array([axis_crossing_z(x, A, n) for x in xs])
    z0 = np.polyfit(xs**2, z, 1)[-1]              # crossing z as x→0
    return A - z0                                 # depth below the crest (z=A)


def paraxial_deviation_pct(A, n=N_WATER, slope_max=1e-3):
    fp = paraxial_focus_m(A, n)
    fr = near_axial_focus_m(A, n, slope_max)
    return 100.0 * abs(fr - fp) / fp


def marginal_focus_m(A, u, n=N_WATER):
    """Axis-crossing depth below crest for a ray of normalized aperture u∈(0,1],
    where u=1 maps to a quarter wavelength (x = u·λ/4). Reveals aberration."""
    x = u * LAM_M / 4.0
    return A - axis_crossing_z(x, A, n)


def spherical_aberration_sign(A, n=N_WATER):
    """+1 if marginal rays cross LATER (deeper) than paraxial = NEGATIVE
    spherical aberration (soft deep glow); −1 otherwise."""
    fp = paraxial_focus_m(A, n)
    fm = marginal_focus_m(A, 0.4, n)
    return 1 if fm > fp else -1


if __name__ == "__main__":
    for A_um in (20, 200):
        A = A_um * 1e-6
        print(f"A={A_um}µm  R={crest_radius(A)*1e3:7.3f}mm  "
              f"f_paraxial={paraxial_focus_m(A)*100:8.4f}cm  "
              f"near-axial dev={paraxial_deviation_pct(A):.5f}%  "
              f"aberration sign={spherical_aberration_sign(A):+d} "
              f"(marginal f={marginal_focus_m(A,0.4)*100:.3f}cm)")
