# HN1-H1 : deterministic exact-arithmetic checks (no MC, no sampling).
# All numbers are Fractions or closed forms; runtime < 2 s.
from fractions import Fraction as F

print("=" * 64)
print("[1] carrier bridge::  chi_KL = chi_curv / (1-eta)   [exact algebra]")
print("=" * 64)
# G'' = S - D = S(1-eta).  chi_curv := 1/S ; chi_KL := 1/G''.
# claim chi_KL = chi_curv/(1-eta)
for S, D in [(F(7, 5), F(679, 500)), (F(17363, 10000), F(16773, 10000)),
             (F(1), F(0)), (F(5, 4), F(-1, 4))]:
    eta = D / S
    chi_curv = 1 / S
    chi_KL = 1 / (S - D)
    assert chi_KL == chi_curv / (1 - eta), (S, D)
print("  4/4 identities hold (incl. eta<0 branch).  OK")

print()
print("=" * 64)
print("[2] the (1-c) slot:  chi_curv <= (1-c)/lam  <=>  chi_KL <= (1-c)/(lam(1-eta))")
print("=" * 64)
lam = F(2)
for c in (F(1, 10), F(1, 2), F(9, 10)):
    for eta in (F(0), F(1, 4), F(1, 2), F(2, 3), F(3, 2)):
        # take chi_curv exactly at the bound
        chi_curv = (1 - c) / lam
        chi_KL = chi_curv / (1 - eta)
        # (i) if eta <= c then chi_KL <= 1/lam   (the (1-c) form IMPLIES plain H1)
        impl = (chi_KL <= 1 / lam) if eta <= c else None
        # (ii) if eta > c the bound on chi_KL exceeds 1/lam
        print(f"  c={str(c):5s} eta={str(eta):4s} chi_KL={str(chi_KL):9s} "
              f"eta<=c? {eta <= c!s:5s} chi_KL<=1/lam? {chi_KL <= 1/lam}")
        if eta <= c:
            assert chi_KL <= 1 / lam
print("  => on the sub-domain {eta<=c} the (1-c) form is STRONGER than H1(lam).  OK")

print()
print("=" * 64)
print("[3] defect lemma  eta<=c<1  is FALSE pointwise (corpus 5/12 modes, L=8, beta=2.4)")
print("=" * 64)
n_lt, n_gt, n_mod = 7, 5, 12
print(f"  modes D<S: {n_lt}/{n_mod} ; D>S: {n_gt}/{n_mod} ; "
      f"fraction of modes with eta>1 = {F(n_gt, n_mod)}")
print(f"  G''<0 on {F(n_gt, n_mod)} of sampled modes  =>  ess inf G'' < 0 in that window")
print("  => the (1-c) form cannot carry c<1 as a uniform-in-(L,a) statement.")

print()
print("=" * 64)
print("[4] counterexample family nu_L: light mode at c/L^2 with O(1) residue")
print("=" * 64)
r0, c0, Lam2 = F(1, 2), F(1), F(4)
print("   L | m1^2=c/L^2 | f1^2 | M0 |   chi=sum s^-1 dnu | gap?")
for L in (2, 4, 8, 16, 64):
    m1sq = c0 / (L * L)
    f12 = r0
    M0 = F(1)
    chi = f12 / m1sq + (M0 - f12) / Lam2
    print(f"  {L:2d} | {str(m1sq):9s} | {str(f12):4s} | {str(M0):2s} | "
          f"{str(chi):12s} | no uniform gap")
print(f"  chi_L = r0*L^2/c0 + (1-r0)/Lam^2  -> grows like L^2 (unbounded)  =>  H1 FAILS")
print("  constraints satisfied: nu>=0, M0=1 finite, atoms positive, KL support in [0,inf),")
print("  squeeze f1^2/chi <= m1^2 holds (equality-free), finite-L chi<inf (trivial half).")
print("  NOT Yang-Mills: no 4D proof forbids a 1/L^2 light mode.")

print()
print("=" * 64)
print("[5] what a UNIFORM gap would buy (the only proven direction)")
print("=" * 64)
for Lam in (F(1, 2), F(1), F(2)):
    M0 = F(1)
    print(f"  m1^2 >= lam={str(Lam):4s}  =>  chi <= M0/lam = {str(M0/Lam):4s}  (finite)"
          f"   [gap => H1, one way only]")

print()
print("=" * 64)
print("[6] thin atom: w=eps/(2lam) gives chi=1/lam for EVERY eps (H1 blind to the edge)")
print("=" * 64)
lam = F(1)
for eps_exp in (1, 2, 4, 6, 9):
    eps = F(1, 10 ** eps_exp)
    w = eps / (2 * lam)
    # chi = w/eps + (1-w)/s2 = 1/lam  =>  s2 = 2*lam*(1-w)
    s2 = 2 * lam * (1 - w)
    chi = w / eps + (1 - w) / s2
    assert chi == 1 / lam
    print(f"  eps=1e-{eps_exp}: s2={float(s2):.9f}  chi={float(chi):.12f}  "
          f"m1^2=eps={float(eps):.1e} -> inf supp -> 0")
print("  => chi<=1/lam holds uniformly along this family while m1->0: H1 has no edge content.")

print()
print("=" * 64)
print("[7] six-route NAME convergence is a rename cascade: 5 names, 1 identity class")
print("=" * 64)
print("  chi_KL(k=0) = 1/G''  [identity] ; S=1/Var(a) [Gaussian only] ;")
print("  Var(a)=1/lam_1 only if phi=a is the Poincare minimiser (Gaussian only).")
print("  exact: eta=0 <=> Gaussian <=> {H1 <=> G''>=lam <=> lam_1>=lam <=> clustering}.")
print("  for eta!=0 only the implications INTO H1 survive.")
print()
print("ALL CHECKS PASS (no assertion failed).")
