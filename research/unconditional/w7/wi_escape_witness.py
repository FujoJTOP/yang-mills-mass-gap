"""
wi_escape_witness.py -- does weight-invariance (WI) escape the OPE-RHO no-go?

OPE-RHO theorem: the support floor is usc, every finite measurable INTEGRAL datum
(moments, tail, threshold, operator bounds = weak-continuous lsc integrals) is blind
to the atom (Delta_min = 0 in every integral class).

WI signature (CONSTRUCT-IR-ANALYTIC-NUM): A = pi/Im(tau_zero), tau_zero = nearest
complex-tau zero of  Z(tau) = w0 + sum_i f_i^2 exp(-E_i tau);  A is INVARIANT under
reweighting for an ISOLATED atom, but MOVES for a continuum.

Structural claim tested here: WI reads the ANALYTIC SINGULARITY of the Laplace
transform -- a simple POLE (isolated atom) vs a BRANCH CUT (continuum).  The POLE
LOCATION Im(tau)=pi/E is weight-independent; only the residue f^2 (via Re(tau)) shifts.
That is exactly A3 (f1^2 = -Res of the Stieltjes transform).  So WI is OUT of the
OPE-RHO *integral* class (it is not an integral) but INTO the A3 class (exact reading
that costs the spectrum = H1).

numpy only.  py -3.12 -I wi_escape_witness.py
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)

def Z(tau, w0, levels):
    s = complex(w0)
    for E, a in levels:
        s += a * np.exp(-E * tau)
    return s

def dZ(tau, w0, levels):
    s = 0j
    for E, a in levels:
        s += a * (-E) * np.exp(-E * tau)
    return s

def newton(seed, w0, levels, it=200):
    tau = complex(seed)
    for _ in range(it):
        fp = dZ(tau, w0, levels)
        if abs(fp) < 1e-300:
            break
        tau = tau - Z(tau, w0, levels) / fp
    return tau

def nearest_zero(w0, levels):
    """Zero with smallest |Im| > 0, Re in a physical window (search + Newton)."""
    best = None
    Es = [E for E, _ in levels]
    Emax = max(Es)
    for x in np.linspace(-4.0, 8.0, 31):
        for y in np.linspace(0.15 * np.pi / Emax, 4.0 * np.pi / min(Es), 60):
            for sgn in (+1.0, -1.0):
                z = newton(x + 1j * sgn * y, w0, levels)
                if abs(Z(z, w0, levels)) < 1e-7 and abs(z.imag) > 1e-6 \
                        and -1.0 < z.real < 12.0:
                    if best is None or abs(z.imag) < abs(best.imag):
                        best = z
    return best

def A_of(w0, levels):
    z = nearest_zero(w0, levels)
    if z is None:
        return float('nan'), None
    return np.pi / abs(z.imag), z

print("=" * 78)
print("1) ISOLATED ATOM  ->  pole at s=-E : Im(tau) fixed, A = E weight-invariant")
print("=" * 78)
E, w0 = 1.0, 1.0
print("exact zero: tau = [ln(f^2/w0) - i*pi]/E  =>  Im = -pi/E independent of f^2")
print(f"{'f^2':>8} {'Re(tau)':>12} {'Im(tau)':>12} {'A=pi/|Im|':>14}")
for f2 in (0.05, 1.0, 20.0):
    Re = np.log(f2 / w0) / E
    Im = -np.pi / E
    print(f"{f2:8.3g} {Re:12.6f} {Im:12.6f} {np.pi/abs(Im):14.6f}")
print("=> weight only shifts Re(tau); A = pi/|Im| = E EXACTLY, independent of f^2.\n")

print("=" * 78)
print("2) CONTINUUM (two close levels) -> nearest zero MOVES under reweighting")
print("=" * 78)
w0 = 1.0
print(f"{'(a,b) at E=(1.0,1.3)':>24} {'nearest Re':>12} {'nearest Im':>12} {'A=pi/|Im|':>14}")
for a, b in ((1.0, 1.0), (5.0, 1.0), (1.0, 5.0), (1.0, 0.3)):
    A, z = A_of(w0, [(1.0, a), (1.3, b)])
    if z is None:
        print(f"{str((a,b)):>24} {'--- no zero found ---':>40}")
    else:
        print(f"{str((a,b)):>24} {z.real:12.6f} {z.imag:12.6f} {A:14.6f}")
print("=> A MOVES with the weight split; no weight-invariant (cut, not pole).\n")

print("=" * 78)
print("3) DENSITY: a discretised continuum keeps MOVING as N grows (finite-L cut)")
print("=" * 78)
w0 = 1.0
print(f"{'N levels on [1,1.3]':>22} {'nearest Im':>12} {'A=pi/|Im|':>14}")
for N in (2, 3, 5, 8):
    Es = np.linspace(1.0, 1.3, N)
    levels = [(E, 1.0) for E in Es]
    A, z = A_of(w0, levels)
    if z is None:
        print(f"{N:22d} {'--- none ---':>26}")
    else:
        print(f"{N:22d} {z.imag:12.6f} {A:14.6f}")
print("=> as a 'continuum' is resolved (N up), the nearest-zero A is NOT weight/density")
print("   invariant: it tracks the discrete level ladder = the gap's isolation.\n")

print("=" * 78)
print("VERDICT (structural)")
print("=" * 78)
print("""
WI reads the singularity TYPE of the Laplace transform:
  isolated atom  = simple POLE at s=-E  -> Im(tau)=pi/E weight-invariant  (A=E exact)
  continuum      = BRANCH CUT           -> zeros dense/moving            (A not invariant)
Reading 'pole vs cut' = reading the Stieltjes pole = A3 (f1^2 = -Res Pi).
  => WI is OUT of the OPE-RHO *integral* class (it is not an lsc integral), so the
     OPE-RHO theorem does not *forbid* it.
  => but WI does NOT escape the wall: it is a re-encoding of A3, and at finite L it
     measures ISOLATION of the first level = the gap = H1/I2.
""")
