"""
savvidy_dterm.py  --  SAVVIDY-DTERM witness  [V-here]
numpy only, py -3.12 -I.

Exact objects (SLOT-SD-FIXPOINT / INVENT-VARIANCE-INEQ):
    G''(0) = S - D,   S = E[h_n - <K^{-1}c,c>|0],   D = Var(f|0) - E[<K^{-1}c,c>|0].
Savvidy cell:  Hess = H_+ - |lam_-| P_-  (P_- rank one),  D = D_+ + E[|<phi_-|c>|^2]/(g|B|)
with D_+ <= 0 by Brascamp-Lieb on the positive fiber.  Hence D <= E[|<phi_-|c>|^2]/(g|B|),
and  D < S  <=>  G''>0  is implied by  E[|<phi_-|c>|^2]/(g|B|) < S_bare.

Toy RG closure: S=kappa (bare), D=g^2 Var(phi^2), G''=kappa - g^2 Var(phi^2), D/S=g^2Var/kappa.

Checks:
 (1) D_+ <= 0 on a positive (log-concave) fiber  [BL, exact];
 (2) Savvidy term > 0 (negative mode);
 (3) D/S table (O(1)) and the sign flip at kappa=soft;
 (4) IR scale separation  |<phi_-|c>|^2 ~ (k/sqrt(gB))^2 -> 0.
"""
import numpy as np

print("=" * 70)
print("SAVVIDY-DTERM:  G''=S-D ;  D <= E|<phi_-|c>|^2/(gB)  =>  D<S <=> G''>0")
print("=" * 70)

# ---------- (1) D_+ <= 0 on the positive fiber (Brascamp-Lieb) ----------
def dplus(B, k, h=0.3, n=400001):
    y = np.linspace(-20.0, 20.0, n)
    w = np.exp(-(0.5 * B * y**2 + 0.25 * k * y**4)); w /= np.trapezoid(w, y)
    f = h * y**2                      # f = d_a V at a=0
    K = B + 3.0 * k * y**2            # fiber Hessian (>=0)
    c = 2.0 * h * y                   # cross-Hessian d2_{ay} V
    Var_f = np.trapezoid(w * f**2, y) - np.trapezoid(w * f, y)**2
    Ec2K = np.trapezoid(w * c**2 / K, y)
    return Var_f - Ec2K

print("\n(1) positive log-concave fiber: D_+ = Var(f|0) - E[c^2/K]   (BL theorem => <= 0)")
for B in [0.5, 1.0, 2.0]:
    for k in [0.0, 0.5, 1.0]:
        print(f"    B={B:>4} k={k:>4}:  D_+ = {dplus(B, k):+.4f}")

# ---------- (2)+(3) Savvidy term and toy D/S ----------
def var_phi2(c, lam, n=400001):
    x = np.linspace(-30.0, 30.0, n)
    w = np.exp(0.5 * c * x**2 - 0.25 * lam * x**4); w /= np.trapezoid(w, x)
    m2 = np.trapezoid(w * x**2, x)
    m4 = np.trapezoid(w * x**4, x)
    return m4 - m2**2

print("\n(2)+(3) Savvidy term D=g^2 Var(phi^2)>0 ;  D/S = g^2Var/kappa  (g=1).  V=-0.5c phi^2+0.25 lam phi^4")
hdr = f"    {'c':>4} {'lam':>5} {'Var(phi^2)':>11} {'g^2Var':>8}   D/S@kap=1.5   2.0   3.0"
print(hdr)
for c, lam in [(1.0, 1.0), (1.0, 0.5), (1.0, 2.0), (2.0, 1.0), (4.0, 1.0), (1.0, 0.2)]:
    Var = var_phi2(c, lam)
    row = f"    {c:>4} {lam:>5} {Var:>11.3f} {Var:>8.3f}"
    for kap in (1.5, 2.0, 3.0):
        row += f"   {Var/kap:>5.2f}"
    print(row)

print("\n    sign flip: D/S = soft/kappa  (soft = g^2 Var(phi^2));  D<S (G''>0) iff kappa > soft:")
soft = 1.0
for ratio in [3.0, 1.5, 1.05, 1.0, 0.9]:
    kap = ratio * soft
    gpp = kap - soft
    var_a = 1.0 / gpp if gpp > 1e-9 else float("inf")
    print(f"      kappa/soft={ratio:>4}: D/S={1.0/ratio:>5.3f}  Var(a)={var_a:>8.3f}  G''={gpp:+.3f}")

# ---------- (4) IR scale separation ----------
print("\n(4) IR scale separation: direct overlap |<phi_-|c>|^2 ~ (k/sqrt(gB))^2 -> 0 (deep IR); O(1) at crossover")
for r in [0.1, 0.3, 0.5, 1.0, 2.0]:
    ov = r**2 / (1.0 + r**2)          # smooth model: 0 at k=0, ->1
    print(f"    k/sqrt(gB)={r:>4}: overlap^2 ~ {ov:>6.3f}")

print("\nVerdict: deep IR  D/S -> 0 (floor survives, uniform);  crossover k~sqrt(gB): D/S = O(1) competition;")
print("         uniform  D<S  at beta=2.5  is OPEN = H1 (softened).   Clay 0%; no QED.")
