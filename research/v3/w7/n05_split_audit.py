"""N05: is G'' = Delta_dyn + W_full an IDENTITY or an INEQUALITY, and is W_full >= N/12 a theorem?

1D marginal fixture (legitimate for the split: eq:split is per-mode algebra, ch05:230-240 / ch18:30-37).
Measure mu ~ exp(-V), V = S + Phi, Phi = a^2/12 (Phi'' = 1/6 = N/12 at N=2, Lemma B floor).
G = d_v S, H = d_v^2 S, Phi' = d_v Phi, Phi'' = d_v^2 Phi.
  G''      = E[H+Phi''] - Var(G+Phi')
  Delta_dyn= E[H] - Var(G)
  W_geo    = E[Phi''] - Var(Phi')
  X        = -2 Cov(G,Phi')
  W_full   = W_geo + X
Py 3.12, numpy only.  Output: n05_split_audit.log
"""
import numpy as np

N = 400001
x = np.linspace(-3.0, 3.0, N)
dx = x[1] - x[0]
PHI2 = 1.0 / 6.0          # N/12, N=2, ||v||=1


def moments(V, Phipp=PHI2):
    """V: callable a -> full free energy; Phi = a^2/12."""
    Phi = 0.5 * x * x / 6.0
    S = V(x) - Phi
    p = np.exp(-(V(x) - V(x).min()))
    p /= p.sum()
    # derivatives by finite differences on the grid
    Sp = np.gradient(S, dx)
    Spp = np.gradient(Sp, dx)
    Phip = x / 6.0
    G = Sp
    F = G + Phip
    EH, EHp2 = np.sum(p * Spp), Phipp
    EG, EF = np.sum(p * G), np.sum(p * F)
    vG = np.sum(p * (G - EG) ** 2)
    vF = np.sum(p * (F - EF) ** 2)
    vPhi = np.sum(p * (Phip - np.sum(p * Phip)) ** 2)
    cov = np.sum(p * (G - EG) * (Phip - np.sum(p * Phip)))
    d = dict(EH=EH, EG=EG, vG=vG, vF=vF, vPhi=vPhi, cov=cov,
             EPhi2=EHp2,
             Gpp=(EH + EHp2) - vF,
             Ddyn=EH - vG,
             Wgeo=EHp2 - vPhi,
             X=-2.0 * cov)
    d["Wfull"] = d["Wgeo"] + d["X"]
    d["resid"] = d["Gpp"] - (d["Ddyn"] + d["Wfull"])
    return d


print("=" * 78)
print("[1] IDENTITY TEST: residual  G'' - (Delta_dyn + W_geo + X)  over random smooth V")
print("=" * 78)
rng = np.random.default_rng(20261125)
worst = 0.0
for k in range(12):
    c = rng.uniform(-0.5, 1.5, 5)
    V = lambda a, c=c: sum(ci * a ** (i + 2) / (i + 2) for i, ci in enumerate(c)) + 0.5 * 0.3 * a * a
    m = moments(V)
    worst = max(worst, abs(m["resid"]))
    print(f"  V#{k:02d}  G''={m['Gpp']:+.9f}  Dy={m['Ddyn']:+.9f}  Wg={m['Wgeo']:+.9f} "
          f" X={m['X']:+.9f}  Wf={m['Wfull']:+.9f}  resid={m['resid']:+.3e}")
print(f"  => max |residual| = {worst:.3e}   (identity, no bounded-potential hypothesis used)")

print()
print("=" * 78)
print("[2] WITNESS: G'' > 0 with Delta_dyn < 0   (needs NO positivity of W_full)")
print("=" * 78)
best = None
win = []
for c0 in np.linspace(-1.5, 1.2, 28):          # negative c0 = double well (non-log-concave marginal)
    for g in np.linspace(0.03, 1.2, 24):
        V = lambda a, c0=c0, g=g: 0.5 * c0 * a * a + 0.25 * g * a ** 4
        m = moments(V)
        if m["Gpp"] > 0 and m["Ddyn"] < 0:
            win.append((c0, g, m))
            if best is None or m["Ddyn"] < best[2]["Ddyn"]:
                best = (c0, g, m)
print(f"  witnesses found: {len(win)} / 144")
if best:
    c0, g, m = best
    print(f"  V = {c0:.3f}/2 a^2 + {g:.3f}/4 a^4  (double well; marginal NOT log-concave)")
    print(f"    G''      = {m['Gpp']:+.6f}")
    print(f"    Delta_dyn= {m['Ddyn']:+.6f}   <0")
    print(f"    W_geo    = {m['Wgeo']:+.6f}")
    print(f"    X        = {m['X']:+.6f}")
    print(f"    W_full   = {m['Wfull']:+.6f}  (sign = {np.sign(m['Wfull']):+.0f})")
    print(f"    resid    = {m['resid']:+.3e}")
    print("  => G''>0 while Delta_dyn<0: the implication fails from the IDENTITY alone.")

print()
print("=" * 78)
print("[3] IS  W_full >= N/12  A THEOREM?  (scan the same family; floor=1/6=0.166667)")
print("=" * 78)
below, tot, mn = 0, 0, 1e9
xneg = 0
for c0 in np.linspace(0.02, 1.2, 12):
    for g in np.linspace(0.05, 2.0, 12):
        V = lambda a, c0=c0, g=g: 0.5 * c0 * a * a + 0.25 * g * a ** 4
        m = moments(V)
        tot += 1
        if m["Wfull"] < PHI2:
            below += 1
        if m["X"] < 0:
            xneg += 1
        mn = min(mn, m["Wfull"])
print(f"  points: {tot}   W_full < N/12: {below}   X<0: {xneg}   min W_full = {mn:+.6f}")
print("  E[Phi''] >= N/12 holds identically here (Phi''=1/6), yet W_full drops below N/12")
print("  => 'W_full >= N/12' is NOT implied by Lemma B (ch05:203-219 bounds W_geo only).")
print("=" * 78)
