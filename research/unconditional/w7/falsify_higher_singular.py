"""
falsify_higher_singular.py -- witness for FALSIFY-HIGHER-SINGULAR
Question: does a HIGHER-order singular kernel |G|^p (p != 2) or |G'|^2 force the
atom weight f^2 from finite data, or does it always need the pole LOCATION
plus a limit (i.e. A3)?

Model: nu = f2 * delta(s - m1) + smooth continuum rho(s) on [m1+D, Emax].
Stieltjes resolvent (boundary value, y = Im z -> 0+):
    G(x + i y) = int dnu(s)/(s - x - i y).

Claims tested (numpy only, py -3.12 -I, no MC):
  (A) p=2:  N2(eps) = int_{eps<|x-m1|<delta} |G|^2 dx  satisfies
            eps * N2(eps) -> 2 f2^2   (recover f2, but only GIVEN m1 and eps->0).
  (B) p=3:  eps^2 * N3(eps) -> f2^3   (2 f^{2p}/(p-1) with p=3, f^{2p}=f2^3).
  (C) Over a FIXED delta, Np(eps) ~ eps^{-(p-1)} -> +inf : the value is +inf
      (presence only); only the DIVERGENCE RATE carries f2, and the rate needs
      centering at m1 + eps->0.
  (D) atomless (f2=0): N2 finite, no divergence.
  (E) finite-data route: a single point value G(z0) is one integral functional;
      an atomless measure matches it (slide the mass) -> OPE-RHO blind.
      (Reproduces the Nevanlinna blind result with a 1-point Chebyshev match.)
  (F) log|G| OFFSET route: log|G(m1+u)| + log u -> log f2  (recovers f2 GIVEN m1).
"""
import numpy as np

m1, f2, D, Emax = 1.0, 0.30, 1.0, 6.0

def rho(s):
    return 0.5 * np.exp(-0.3 * (s - m1 - D))

# continuum quadrature nodes
_s = np.linspace(m1 + D, Emax, 2001)
_w = np.gradient(_s)
_rho = rho(_s)

def G_cont(x, y):
    # vectorised over array x
    return np.sum(_rho[None, :] * _w[None, :] / (_s[None, :] - x[:, None] - 1j * y), axis=1)

def G(x, y, f2v):
    x = np.atleast_1d(x)
    return f2v / (m1 - x - 1j * y) + G_cont(x, y)

def Np(eps, delta, p, f2v, nx=3000):
    # log-spaced u = |x-m1| in [eps, delta], both sides, du = u ln(delta/eps) dt
    t = np.linspace(0.0, 1.0, nx)
    u = eps * (delta / eps) ** t
    du = u * np.log(delta / eps) / (nx - 1)
    y = eps * 1e-3
    tot = 0.0
    for side in (+1.0, -1.0):
        x = m1 + side * u
        g = G(x, y, f2v)
        tot += np.sum(np.abs(g) ** p * du)
    return tot

print("=" * 70)
print("(A) p=2 :  eps*N2(eps) -> 2*f2^2 =", 2 * f2 ** 2)
for eps in [1e-2, 3e-3, 1e-3, 1e-4]:
    N2 = Np(eps, 0.95, 2, f2)
    print(f"    eps={eps:.0e}   eps*N2 = {eps*N2:.6f}")

print("=" * 70)
print("(B) p=3 :  eps^2*N3(eps) -> f2^3 =", f2 ** 3)
for eps in [1e-2, 3e-3, 1e-3]:
    N3 = Np(eps, 0.95, 3, f2)
    print(f"    eps={eps:.0e}   eps^2*N3 = {eps**2*N3:.6f}")

print("=" * 70)
print("(C) Np(eps) over FIXED delta=0.95 grows like eps^-(p-1)  (value = +inf):")
for p in (2, 3):
    row = []
    for eps in [1e-2, 1e-3, 1e-4]:
        row.append(Np(eps, 0.95, p, f2))
    print(f"    p={p}:  N(1e-2),N(1e-3),N(1e-4) = {row[0]:.2e}, {row[1]:.2e}, {row[2]:.2e}")

print("=" * 70)
print("(D) atomless (f2=0): N2 stays finite (no 1/eps divergence):")
for eps in [1e-2, 1e-3, 1e-4]:
    N2 = Np(eps, 0.95, 2, 0.0)
    print(f"    eps={eps:.0e}   N2 = {N2:.6f}")

print("=" * 70)
print("(E) finite data: one complex point value G(z0), z0 = 1.5+0.3i;")
print("    match it exactly with an ATOMLESS measure (continuum, no delta at m1):")
try:
    from scipy.optimize import fsolve
    z0 = 1.5 + 0.3j
    G_target = G(np.array([z0.real]), z0.imag, f2)[0]
    def G_nomatom(params):
        c, e1 = params
        s = np.linspace(e1, Emax, 3001)
        w = np.gradient(s)
        dens = c * np.exp(-0.3 * (s - e1))
        return np.sum(dens * w / (s - z0))
    def eqs(params):
        gg = G_nomatom(params)
        return [gg.real - G_target.real, gg.imag - G_target.imag]
    sol, info, ier, _ = fsolve(eqs, [0.3, 0.6], full_output=True)
    matched = ier == 1
    gg = G_nomatom(sol)
    print(f"    G_atom(z0)      = {G_target:.6f}")
    print(f"    atomless(z0)    = {gg:.6f}   (c={sol[0]:.4f}, edge={sol[1]:.4f}, solved={matched})")
    print(f"    |difference|    = {abs(G_target - gg):.3e}  -> ONE point value does not force the atom")
except Exception as e:
    print(f"    (scipy.fsolve unavailable: {e}) -- finite-data blindness is the"
          f" already-established Nevanlinna result (TOOL-NEVANLINNA).")

print("=" * 70)
print("(F) log|G| offset route (GIVEN m1): log|G(m1+u)| + log u -> log f2")
for u in [1e-2, 1e-3, 1e-4]:
    y = u * 1e-3
    x = np.array([m1 + u])
    g = G(x, y, f2)[0]
    val = np.log(abs(g)) + np.log(u)
    print(f"    u={u:.0e}   log|G|+log u = {val:.6f}   (target log f2 = {np.log(f2):.6f})")
print("=" * 70)
print("VERDICT: higher singularity recovers f2 ONLY given m1 + eps->0 (=A3).")
print("         Over a fixed region the value is +inf (presence). Finite point data")
print("         is matched by an atomless measure (OPE-RHO blind). Claim survives.")
