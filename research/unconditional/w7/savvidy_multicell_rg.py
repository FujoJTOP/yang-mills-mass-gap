"""
savvidy_multicell_rg.py  --  witness for research/unconditional/SAVVIDY-MULTICELL-RG.md

Does the RG-accumulated Savvidy softening
    S = int_{k_IR}^{k_UV} (dk/k) g^2(k) Var(phi^2)
stay sub-critical (kappa_bare > S) across all cells, using asymptotic freedom?

One-loop AF (pure SU(N)):  g^2(k) = C_RG / ln(k/Lam),  C_RG = 24 pi^2/11 (SU(2)).
Then  S = C_RG * Var * int du/u = C_RG * Var * ln( ln(k_UV/Lam)/ln(k_IR/Lam) ),  u=ln(k/Lam).

Cheap, single process, numpy only.  Labels [T]/[V-num].  No QED.
"""
import numpy as np

C_RG = 24 * np.pi**2 / 11.0   # SU(2), matches lean C_RG = 24 pi^2 / 11
print(f"C_RG (SU(2)) = 24*pi^2/11 = {C_RG:.4f}")

def g2(k, Lam, C=C_RG):
    return C / np.log(k / Lam)

def S_closed(k_ir, k_uv, Lam, Var, C=C_RG):
    return C * Var * np.log(np.log(k_uv / Lam) / np.log(k_ir / Lam))

def S_numeric(k_ir, k_uv, Lam, Var, C=C_RG, n=4001):
    k = np.geomspace(k_ir, k_uv, n)
    integrand = (1.0 / k) * g2(k, Lam, C) * Var
    return np.trapezoid(integrand, k)

Lam = 1.0
k_uv = 100.0 * Lam
Var = 1.0

print("\n(A) closed form vs numeric  [T]")
for k_ir in [10.0, 5.0, 2.0]:
    sc = S_closed(k_ir, k_uv, Lam, Var)
    sn = S_numeric(k_ir, k_uv, Lam, Var)
    print(f"  k_IR={k_ir:5.1f}Lam  closed={sc:10.4f}  numeric={sn:10.4f}  rel.diff={abs(sc-sn)/sc:.2e}")

print("\n(B) S(k_IR) grows LOG-LOG as k_IR -> Lam (diverges at the Landau pole)  [V-num]")
print(f"  fixed: k_UV=100 Lam, Var=1, C_RG={C_RG:.3f}")
print("   k_IR/Lam   ln(k_IR/Lam)   S(k_IR)")
for r in [10.0, 5.0, 3.0, 2.0, 1.5, 1.2, 1.1, 1.05, 1.01]:
    u = np.log(r)
    print(f"   {r:6.2f}   {u:11.4f}   {S_closed(r, k_uv, Lam, Var):10.3f}")
print("  => S finite for any k_IR>Lam, but grows without bound (ln(1/ln k_IR)) as k_IR->Lam.")

print("\n(C) the softening is IR-DOMINATED  [T]")
print("  per-decade contribution dS/d ln k = g^2(k)*Var ~ 1/ln(k/Lam), LARGEST at small k:")
for r in [50.0, 10.0, 3.0, 1.5]:
    print(f"   k={r:5.1f}Lam : g^2*Var = {g2(r,Lam)*Var:7.3f}")
print("  => ~all of S accumulates in the last few IR cells near Lam.")

print("\n(D) sub-criticality margin kappa_bare - S  [V-num]")
print("  Var in the SAVVIDY-CELL range (0.41..8.51); natural kappa_bare ~ O(1) in IR units:")
for Var_v in [0.41, 1.0, 8.51]:
    S20 = S_closed(2.0, k_uv, Lam, Var_v)
    S11 = S_closed(1.1, k_uv, Lam, Var_v)
    print(f"   Var={Var_v:5.2f}: S(k_IR=2Lam)={S20:8.2f}, S(k_IR=1.1Lam)={S11:8.2f}"
          f"   (kappa_bare=1 margin: {1.0-S20:8.2f} / {1.0-S11:8.2f})")
print("  => S is O(10-100)*Var, FAR above a natural O(1) kappa_bare; margin NEGATIVE for natural kappa_bare.")
print("     Sub-criticality is NOT generic; it needs kappa_bare large enough, decided at k_IR~Lam.")

print("\nVERDICT: accumulation is IR-marginal (log-log), IR-dominated, DIVERGES at the Landau pole;")
print("sub-criticality kappa_bare>S is not perturbative -- decided at k_IR~Lam where g^2->inf")
print("(perturbation theory breaks down) = H1.  No QED.")
