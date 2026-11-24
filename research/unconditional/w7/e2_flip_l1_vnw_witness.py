"""
e2_flip_l1_vnw_witness.py  —  [V-here] for E2-FLIP-L1-VNW.md
Pure-Python (no deps).  Witnesses the EXACT distinction that decides the L-1 attack:

  vNW isolation is RELATIVE:  E2 > E1  (pointwise, homogeneous in the mass scale).
  Temple needs an ABSOLUTE certified floor:  E2 >= C  with C a NUMBER independent of m1.
  By single-scale covariance (SSC)  E2 = rho_true * m1,  so "E2 >= C uniformly" IS "m1 >= C/rho_true",
  i.e. the gap.  And a certified RELATIVE ratio E2 >= (1+c) E1 fed to Temple is a tautology
  (output = m1), carrying NO absolute floor.

Model (two-level, exact):  |psi> = sqrt(1-s)|psi1> + sqrt(s)|psi2>,  true E1=m, E2=rt*m.
  R   = <H>   = (1-s)E1 + s E2
  Var = <H^2>-R^2 = s(1-s)(E2-E1)^2
  Temple(R,Var,rho) = R - Var/(rho-R),  valid iff R < rho <= E2.
"""

def trial(s, m, rt):
    E1, E2 = m, rt * m
    R = (1.0 - s) * E1 + s * E2
    Var = s * (1.0 - s) * (E2 - E1) ** 2
    return R, Var, E1, E2

def temple(R, Var, rho):
    return R - Var / (rho - R)

line = "=" * 74

print(line); print("(1) RELATIVE isolation: certified E2 >= (1+c)E1  ->  rho_slot = (1+c)*m"); print(line)
s, rt, c = 0.02, 1.5, 0.4
print(f"true ratio rt={rt}, certified (1+c)={1+c}  (valid since {1+c} <= {rt})")
print(" m      R        Var        rho_slot   Temple_out   out/m")
for m in [0.1, 1.0, 3.0, 10.0, 100.0]:
    R, Var, E1, E2 = trial(s, m, rt)
    rho = (1 + c) * m
    out = temple(R, Var, rho)
    print(f"{m:6.2f} {R:9.4f} {Var:10.6f} {rho:9.4f} {out:11.4f} {out/m:8.5f}")
print("  => out/m is CONSTANT (independent of m): Temple output is a fixed multiple of m.")
print("     A homogeneous output is NO absolute floor (m may still -> 0).")

print(); print(line); print("  (1b) certified ratio EQUALS true ratio (c = rt-1): Temple reads 'm >= m'"); print(line)
c2 = rt - 1.0
for m in [1.0, 7.0]:
    R, Var, E1, E2 = trial(s, m, rt)
    rho = (1 + c2) * m
    print(f" m={m:5.2f}  rho_slot={rho:8.4f}  Temple_out={temple(R,Var,rho):8.4f}  ( == m exactly )")
print("  => RELATIVE isolation -> Temple is an exact tautology m1 >= m1: zero absolute information.")

print(); print(line); print("(2) ABSOLUTE isolation: certified E2 >= C (fixed number, indep of m)"); print(line)
C = 1.5; rt = 1.5; s = 0.02
print("Hypothesis rho=C <= E2=rt*m  is itself the gap lower bound m >= C/rt.")
print(" m      E2        C<=E2?   R<C?   Temple_out")
for m in [0.5, 1.0, 1.01, 2.0, 5.0]:
    R, Var, E1, E2 = trial(s, m, rt)
    ok = (C <= E2 + 1e-15)
    out = temple(R, Var, C) if (R < C) else float('nan')
    print(f"{m:6.2f} {E2:9.4f}  {str(ok):6s}  {str(R<C):6s} {out:11.4f}")
print(f"  => absolute floor appears only for m >= C/rt = {C/rt:.3f}, which was ASSUMED to state C<=E2.")
print("     'E2 >= C uniformly' == 'm1 >= C/rt' (the gap) by SSC.  CIRCULAR.")

print(); print(line); print("(3) FREE MASSLESS acid test: no isolated 0++, so E1=E2=0"); print(line)
m, rt, s = 0.0, 1.5, 0.02
R, Var, E1, E2 = trial(s, m, rt)
print(f"  free: E1={E1}, E2={E2}, R={R}, Var={Var}")
print("  ABSOLUTE tool needs C <= E2 = 0: any C>0 FAILS -> returns nothing  (acid-clean: fails on free field).")
print("  RELATIVE tool reads E2 >= (1+c)E1 as 0 >= 0 (trivial); NO tower exists -> vacuous, no gap certified.")
print("  => free field gives NO isolation; the tool cannot manufacture a gap there.  Acid test PASSES.")
print("  => BUT the acid test does NOT separate pointwise (real) from uniform (would-be): SSC does.")

print(); print(line); print("(4) SCALE COVARIANCE: relative statement invariant under m -> lambda*m"); print(line)
s, rt, c = 0.02, 1.5, 0.4
for lam in [0.01, 1.0, 100.0]:
    m = lam
    R, Var, E1, E2 = trial(s, m, rt)
    rho = (1 + c) * m
    out = temple(R, Var, rho)
    print(f"  lambda={lam:8.3f}: Temple_out/lambda = {out/lam:.5f}  (identical for all lambda)")
print("  => under m->lambda*m the whole relative statement is unchanged: no scale-invariant floor.")

print(); print(line)
print("CONCLUSION  [T given SSC]:")
print("  vNW gives RELATIVE isolation E2>E1 (pointwise, scale-covariant, no uniform margin).")
print("  Temple needs E2 >= C with C an ABSOLUTE constant; by E2=rt*m1 that is m1 >= C/rt = the gap.")
print("  A certified RELATIVE ratio alone makes Temple output = m1 (tautology) -> no floor.")
print("  Therefore uniformizing vNW = (the ratio rho) AND (the absolute scale) = H*.  IRREDUCIBLE.")
print(line)
