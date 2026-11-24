"""
construct_dual_stabilized.py  --  CONSTRUCT-DUAL-STABILIZED witness  [V-num]

Decide the atom question for the TRAPPED vortex ring in the dual abelian-Higgs.

Model (units xi = 1, A = 1):
    E(R) = A R ( ln(8R/xi) - 7/4 ) + B xi^2 / R
    E'(R) = A ( ln(8R/xi) - 3/4 ) - B xi^2 / R^2
    E''(R) = A / R + 2 B xi^2 / R^3

Questions:
 (1) Does B > 0 give a GENUINE KL atom (isolated stable ring at finite R* => delta in rho_trF2)?
 (2) Does the atom survive as B -> 0 (a small-B limit of a physical atom), or only for B >= B_c ?
 (3) Does the MINIMAL dual GL / YM have any B > 0 ?
     (the only power-1/R term of the minimal theory is the (D-2) Luescher term, NEGATIVE)

Labels [V-num]/[T-model]/[C]/[O]/[NG].  Nothing here proves the gap.  No QED.
Run:  py -3.12 -I construct_dual_stabilized.py
"""
import numpy as np

xi = 1.0
A  = 1.0          # line-tension coefficient (~ sigma)

def E(R, B):   return A*R*(np.log(8.0*R/xi) - 7.0/4.0) + B*xi**2/R
def dE(R, B):  return A*(np.log(8.0*R/xi) - 3.0/4.0) - B*xi**2/R**2
def d2E(R, B): return A/R + 2.0*B*xi**2/R**3

def find_min(B, lo=0.05, hi=200.0):
    """Root of E'(R)=0 on (lo,hi) by bisection; E' is increasing (E''>0 for B>0)."""
    f = lambda R: dE(R, B)
    a, b = lo, hi
    if f(a) * f(b) > 0:
        return None
    for _ in range(200):
        m = 0.5*(a+b)
        if f(a)*f(m) <= 0: b = m
        else: a = m
    return 0.5*(a+b)

print("="*74)
print("CONSTRUCT-DUAL-STABILIZED : trapped vortex ring, atom or artifact?")
print("="*74)

# ---- (0) minimal ring B=0 : monotone => continuum (baseline, reproduces prior) ----
Rs = np.linspace(xi, 40*xi, 8000)
print("\n[0] BASELINE  B = 0  (minimal ring, no trapping)")
print(f"    min E'(R) on [xi,40xi] = {dE(Rs,0.0).min():+.4f}  (>0 => monotone, no min => continuum)")

# ---- (1) B > 0 : unique interior minimum, strictly convex => isolated level = ATOM ----
print("\n[1] B > 0 : E'' = A/R + 2B xi^2/R^3 > 0 (strictly convex) => at most one min;")
print("    E'(R->0+) -> -inf, E'(R->inf) -> +inf => exactly one root R* (stable => delta = atom).")
print(f"    {'B':>6} {'R*':>9} {'R*>xi?':>7} {'m*=E(R*)':>10} {'m*/sqrt(sigma)':>14} {'atom?':>7}")
rows=[]
for B in [0.5, 1.0, 1.3294, 2.0, 3.0, 5.0, 10.0]:
    Rst = find_min(B)
    mst = E(Rst, B)
    phys = Rst > xi
    rows.append((B,Rst,mst,phys))
    print(f"    {B:6.3f} {Rst:9.4f} {str(phys):>7} {mst:10.4f} {mst/np.sqrt(A):14.4f} "
          f"{('YES' if phys else 'no(<-xi)'):>7}")

# threshold B_c where R* = xi (atom first becomes physical / thin-core valid)
B_c = dE(xi, 0.0) * xi**2   # solve A(ln8-3/4) = B_c/xi^2
print(f"\n    B_c (R*=xi):  A(ln8 - 3/4) = {A*(np.log(8.0)-0.75):.4f}  =>  B_c = {B_c:.4f}")
print("    => a PHYSICAL atom (R* > xi, thin core valid) needs B > B_c ~ 1.33 A (order-one trapping).")

# ---- (2) B -> 0 : the minimum retreats BELOW xi (outside thin-core validity) ----
print("\n[2] B -> 0+ : where does R* go?")
for B in [0.5, 0.2, 0.05, 0.01]:
    Rst = find_min(B, lo=1e-3, hi=10.0)
    print(f"    B={B:5.3f}:  R*={Rst:.4f} xi  ({'physical' if Rst>xi else 'BELOW xi -> thin-core invalid'})")
R0 = (xi/8.0)*np.exp(0.75)
print(f"    limit B->0+: R* -> (xi/8)e^(3/4) = {R0:.4f} xi  <  xi  (unphysical).")
print("    => the atom is NOT a small-B limit of a physical state; B must be an O(A) input.")

# ---- (3) does the MINIMAL dual GL / YM have B > 0 ? ----
print("\n[3] The minimal 1/R coefficient.  Ring energy = 2 pi R sigma + (curvature log) + c/R + ...")
print("    The only power 1/R term of the minimal (massless, London-limit) dual string is the")
print("    free transverse-fluctuation (Luescher) term:")
D = 4
B_luescher = -(D-2)*np.pi/6.0
print(f"      B_Luescher = -(D-2) pi/6 = -(4-2) pi/6 = {B_luescher:+.4f}   (NEGATIVE)")
print("    A NEGATIVE B makes E'(R) even MORE positive on R>=xi => still monotone => NO min => NO atom.")
print("    Physical reading: a neutral closed flux tube SHRINKS (no stable neutral ring in minimal YM).")
print("    A positive B (>B_c) must come from NON-MINIMAL input: a repulsive core /")
print("    a positive massive-mode Casimir / a non-minimal tr F^2|phi|^2 coupling.")

# ---- (4) even in model B, the atom WEIGHT f1^2 is not fixed by the minimal theory ----
print("\n[4] In model B the atom is a genuine delta at m*=E(R*) with weight w(R*)^2 = f1^2 > 0,")
print("    BUT the tr F^2 overlap w is set by the SAME non-minimal tr F^2|phi|^2 coupling that")
print("    supplies B => f1^2 magnitude is model-dependent (= H1), exactly as for the dual Higgs.")

print("\nVERDICT: with +B/R and B>B_c ~1.33A the trapped ring IS a genuine KL atom (isolated stable")
print("  ring, delta in rho, f1^2>0).  But (i) B->0 puts R* below xi (not a physical limit);")
print("  (ii) the MINIMAL dual GL has B = B_Luescher < 0 (no atom; neutral tube shrinks);")
print("  (iii) the atom weight is fixed only by the non-minimal coupling.  So the atom is REAL")
print("  in the modified model but requires non-minimal input = H1.  Clay 0%; no QED.")
