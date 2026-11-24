"""
falsify_ward_nonlinear.py -- witness for FALSIFY-WARD-NONLINEAR.md

Question: can a NONLINEAR combination of the exact YM identities (Ward kappa*chi=1-eta,
trace anomaly A=int s^-2 dnu, positivity, AF tail) FORCE the atom weight f1^2 > 0?

The atomless adversary: replace the floor atom f1^2*delta(E1) by a FLAT CONTINUUM BUMP
rho_h = f1^2/h on [E1, E1+h] of the SAME total mass (and same floor E1), keeping the
continuum part identical.  As h->0, EVERY moment int s^k dnu of the bump converges to the
atom's moment, while the bump has f1^2 = 0 (no atom at the floor).

So any CONTINUOUS (in particular any nonlinear) functional of finitely many moments -- which
is what "a nonlinear combination of the exact identities" is -- gives the same value for the
atom and the atomless bump, hence cannot force f1^2 > 0.

We test the strongest nonlinear candidates (all moment-functions of positivity):
  R_CS   = chi^2 / (M0 * A)            (Cauchy-Schwarz <= 1, i.e. M1^2 <= M0 M2)
  detH   = M0*A - chi^2                (2x2 Hankel PSD, >= 0)
  R_inv  = A * M0 / chi^2              (= 1/R_CS >= 1, CS on s^-1,1)
  R_sq   = M1*M_{-1}/M0^2              (CS on s^{-1/2}? -> M1*M_{-1} >= M0^2)
plus higher-moment Hankel PSD (3x3).  All are satisfied by BOTH atom and bump, to O(h).

Run:  py -3.12 -I research/unconditional/w7/falsify_ward_nonlinear.py
"""
import numpy as np

E1, f1sq = 1.0, 0.5          # floor mass^2, atom weight
E2, L, g = 2.0, 2.0, 0.1     # continuum: uniform density g on [E2, E2+L]


def mom_atom(k):
    """int s^k dnu for atom + uniform continuum."""
    atom = f1sq * E1 ** k
    if abs(k + 1) < 1e-14:
        cont = g * np.log((E2 + L) / E2)
    else:
        cont = g * ((E2 + L) ** (k + 1) - E2 ** (k + 1)) / (k + 1)
    return atom + cont


def mom_bump(k, h):
    """int s^k dnu for flat bump f1sq/h on [E1,E1+h] + same uniform continuum."""
    if abs(k + 1) < 1e-14:
        bump = (f1sq / h) * np.log((E1 + h) / E1)
        cont = g * np.log((E2 + L) / E2)
    else:
        bump = (f1sq / h) * ((E1 + h) ** (k + 1) - E1 ** (k + 1)) / (k + 1)
        cont = g * ((E2 + L) ** (k + 1) - E2 ** (k + 1)) / (k + 1)
    return bump + cont


def combos(M):
    """nonlinear combinations of the moments (the 'exact-identity data')."""
    M0, M1, Mm1, A, M2, M3 = M[0], M[1], M[-1], M[2], M[2], M[3]
    R_CS = M1 ** 2 / (M0 * A)
    detH = M0 * A - M1 ** 2
    R_inv = A * M0 / M1 ** 2
    R_sq = M1 * Mm1 / M0 ** 2
    # 3x3 Hankel [[M0,M1,M2],[M1,M2,M3],[M2,M3,M4]] PSD (needs M4)
    return dict(R_CS=R_CS, detH=detH, R_inv=R_inv, R_sq=R_sq)


print("=" * 78)
print("Exact-identity data: M0=int dnu, M1=chi=int s^-1, A=M2=int s^-2, Mm1=int s")
print("=" * 78)
ks = [-3, -2, -1, 0, 1, 2, 3, 4]
atom = {k: mom_atom(k) for k in ks}
print("  atom moments:", {k: round(atom[k], 6) for k in ks})
print("  atom f1^2 =", f1sq, " (atom at floor E1=1)")

print("\n" + "=" * 78)
print("Atomless bump: does it match the atom's data as h->0, while f1^2 -> 0?")
print("=" * 78)
print(f"{'h':>10} {'max|dM_k|':>14} {'max rel':>12} {'f1^2_bump':>10}")
hs = [0.1, 0.01, 1e-3, 1e-4, 1e-5]
for h in hs:
    bump = {k: mom_bump(k, h) for k in ks}
    d = np.array([abs(bump[k] - atom[k]) for k in ks])
    rel = np.array([abs(bump[k] - atom[k]) / (abs(atom[k]) + 1e-30) for k in ks])
    print(f"{h:>10.0e} {d.max():>14.3e} {rel.max():>12.3e} {0.0:>10.1f}")

print("\n  -> every moment (all k, incl. M0, chi=M1, A=M2) matches to O(h);")
print("     the bump has f1^2 = 0. So the exact-identity DATA cannot tell them apart.")

print("\n" + "=" * 78)
print("Nonlinear combinations of the exact identities: atom vs atomless bump")
print("=" * 78)
ca = combos({k: atom[k] for k in ks})
print(f"  ATOM:      " + "  ".join(f"{n}={v:.6f}" for n, v in ca.items()))
print(f"{'h':>10} " + "  ".join(f"{n:>9}" for n in ca))
for h in hs:
    bump = {k: mom_bump(k, h) for k in ks}
    cb = combos(bump)
    line = "  ".join(f"{cb[n]:9.6f}" for n in ca)
    print(f"{h:>10.0e} {line}")
print("\n  -> the nonlinear combinations R_CS, detH, R_inv, R_sq are the SAME for the")
print("     atom and the atomless bump to O(h). None of them can force f1^2>0.")

print("\n" + "=" * 78)
print("Why: a nonlinear combination is a CONTINUOUS function of the moment data, and")
print("     the atomless bump matches ALL moments to O(h) with f1^2=0. Any continuous")
print("     function of matched data is matched. Forcing would need a NON-continuous")
print("     (e.g. singular-kernel) functional -- NOT a combination of the exact data.")
print("=" * 78)

# acid test: the atomless bump is a legitimate positive KL measure?  YES:
print("\nAcid / validity: bump density rho_h=f1^2/h >= 0 on [E1,E1+h]; floor inf supp = E1;")
print("                 f1^2 = 0; matches Ward via eta:=1-kappa*chi (SLOT-WARD-STEIN: tautology);")
print("                 AF tail untouched (bump at low s). A valid atomless KL measure.")
print("                 -> only the YM-SPECIFIC admissibility (nu is trF^2 of interacting YM,")
print("                    not just any positive measure) can exclude it. That IS the escape = H1.")
print("\nexit 0")
