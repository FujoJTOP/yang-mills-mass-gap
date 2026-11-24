"""
wi_su3.py -- witness for WI-SU3.md
Weight-invariance (WI) atom test:  A = pi / Im(tau_zero)
of  Z(tau) = w0 + sum_i c_i e^{-E_i tau}   (full 0++ correlator incl. the VEV^2 term w0).

Demonstrates the finite-L atom signature (CONSTRUCT-IR-ANALYTIC-NUM section 4):
  ATOM (one isolated E1)      ->  A = m1 EXACTLY, INVARIANT under overlap reweighting
  CUT  (comparable levels)    ->  A MOVES with the weight distribution
"Reweighting" = change the overlaps c_i (which operator / how smeared), NOT the energies E_i.

Synthetic spectral model, NOT SU(3) data. numpy only.
Run:  py -3.12 -I wi_su3.py
"""
import numpy as np

w0    = 0.5          # VEV^2 term |<tr F^2>|^2  (>0 : needed for a zero)
m1    = 1.0          # isolated atom mass (the lightest 0++)
theta = 2.0          # two-glueball threshold 2*m1
K, Delta = 60, 0.12
sk  = theta + Delta*np.arange(K)
rho = np.exp(-(sk - theta)/0.5)

def A_of(z):
    im = abs(z.imag)
    return np.pi/im if im > 1e-9 else np.nan

def spread(xs):
    xs = [x for x in xs if np.isfinite(x)]
    return float(max(xs) - min(xs)) if xs else float('nan')

# ---- atom + two-glueball continuum ----
def Z_and_dZ(zeta, a, cs):
    z, dz = w0 + 0j, 0j
    ea = np.exp(-m1*zeta); z += a*ea; dz += a*(-m1)*ea
    for s, r in zip(sk, rho):
        c = cs*r*Delta; e = np.exp(-s*zeta); z += c*e; dz += c*(-s)*e
    return z, dz

def atom_guess(a):
    return (np.log(max(a, 1e-12)/w0)/m1) + 1j*np.pi/m1   # zero of w0 + a e^{-m1 z}

def newton(a, cs, z):
    for _ in range(300):
        f, fp = Z_and_dZ(z, a, cs)
        if abs(fp) < 1e-14: break
        step = f/fp; z -= step
        if abs(step) < 1e-13: break
    return z

# ---- discretised cut: two comparable levels (no isolated lightest) ----
def Z2(zeta, b, d, E1=1.0, E2=1.35):
    e1, e2 = np.exp(-E1*zeta), np.exp(-E2*zeta)
    return w0 + b*e1 + d*e2, b*(-E1)*e1 + d*(-E2)*e2

def newton2(b, d, z):
    for _ in range(300):
        f, fp = Z2(z, b, d)
        if abs(fp) < 1e-14: break
        step = f/fp; z -= step
        if abs(step) < 1e-13: break
    return z

def nearest_upper(zs):
    """nearest-to-real-axis zero with Im>0."""
    cand = [z for z in zs if np.isfinite(z.imag) and z.imag > 1e-6]
    return min(cand, key=lambda z: z.imag) if cand else None

print("="*70)
print("WI-SU3 witness: weight-invariance of A distinguishes atom from cut")
print("="*70)

# [1] ATOM only: invariant under overlap reweighting
As = []
print("\n[1] ATOM only (cs=0): reweight atom overlap a -> A INVARIANT")
for a in [0.3, 1.0, 5.0]:
    z = newton(a, 0.0, atom_guess(a)); A = A_of(z); As.append(A)
    print(f"    a={a:<4}  Im(tau)={z.imag:.6f}  A=pi/Im={A:.6f}")
print(f"    spread(A) = {spread(As):.2e}      (atom: ~0)")

# [2] ATOM + CUT: moves with continuum-overlap weight
As = []
print("\n[2] ATOM + two-glueball CUT: vary continuum overlap cs -> A MOVES")
for cs in [0.0, 0.3, 0.8, 2.0]:
    z = newton(1.0, cs, atom_guess(1.0)); A = A_of(z); As.append(A)
    print(f"    cs={cs:<4}  Im(tau)={z.imag:.6f}  A=pi/Im={A:.6f}")
print(f"    spread(A) = {spread(As):.2e}      (cut: >0)")

# [3] CUT (two comparable levels): moves with weights, no invariant A
As = []
print("\n[3] CUT (two comparable levels 1.00, 1.35): reweight (b,d) -> A MOVES")
for (b, d) in [(1.0, 0.0), (1.0, 0.3), (0.3, 1.0), (1.0, 1.0)]:
    z0 = complex(np.log(max(b,1e-9)/w0)/1.0, np.pi/1.0)
    z = newton2(b, d, z0)
    # guard: if Newton left the physical window, scan a bounded grid
    if not (0 < z.imag < 6 and -6 < z.real < 6):
        grid = [newton2(b, d, complex(re, im))
                for re in np.linspace(0.1, 3.0, 40) for im in np.linspace(0.2, 5.0, 40)]
        z = nearest_upper(grid)
    A = A_of(z); As.append(A)
    print(f"    (b,d)=({b},{d})  Im(tau)={z.imag:.6f}  A={A:.6f}")
print(f"    spread(A) = {spread(As):.2e}      (cut: large)")

print("\n" + "="*70)
print("VERDICT  atom -> A invariant (spread ~0) | cut -> A moves (spread >0).")
print("Finite-L atom signature; PROTOCOL = vary the 0++ operator at fixed beta,")
print("fit Z(tau), track A.  Clay 0%; no QED. Synthetic model, not SU(3) data.")
print("="*70)
