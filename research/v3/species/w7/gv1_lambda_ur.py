"""
GV1 — lambda_Ur at beta=2.5 (SU(2)): extrapolation with error bars, verdict arithmetic,
      u-convention audit, free-field acid test.

Inputs : gv1_enum_n6.json  (exact m_n = A_n, n<=6, from gv1_enum.py)
Outputs: console (tee'd to gv1_lambda_ur.log)

Conventions (GU4-NECESSARY.md §4.5):
    c_Gamma  = phi^U(Gamma)/u^{|Gamma|} ,  lambda_Ur = limsup |c_Gamma|^{1/|Gamma|}
    m_n      = sum_{|Gamma|=n, p0 in Gamma} |c_Gamma|          (per-site shell mass)
    criterion: lambda_Ur * u < 1  <=>  lambda_Ur < 1/u
Usage: py -3.12 gv1_lambda_ur.py
"""
import json
import numpy as np
import mpmath as mp

mp.mp.dps = 30
BETA = mp.mpf("2.5")
N = 2
EDELTA = 20 * mp.e

# ------------------------------------------------------------------ [A] u audit
print("=" * 78)
print("GV1 [A] u-convention audit at beta_W = 2.5  (which u is 0.625?)")
print("=" * 78)
u_nom = BETA / 4
cands = [
    ("beta/4                     (corpus: P4:143, TT1:60)", u_nom),
    ("I_2(beta)/I_1(beta)        (standard double-expansion ratio)",
     mp.besseli(2, BETA) / mp.besseli(1, BETA)),
    ("2 I_2(beta/2)/I_1(beta/2)  (task's guess)",
     2 * mp.besseli(2, BETA / 2) / mp.besseli(1, BETA / 2)),
    ("2 I_2(beta)/beta           (fundamental char coeff, norm. by beta)",
     2 * mp.besseli(2, BETA) / BETA),
    ("I_1(beta)/I_0(beta)        (adjoint-ish ratio)",
     mp.besseli(1, BETA) / mp.besseli(0, BETA)),
]
for name, v in cands:
    print(f"  {name:56s} u = {mp.nstr(v, 8):>12s}   1/u = {mp.nstr(1/v, 8):>10s}"
          f"   ubar = eD*u = {mp.nstr(EDELTA*v, 8):>10s}")
print("  => 0.625 is beta/4 EXACTLY; no Bessel ratio gives 0.625 "
      "(0.507 / 0.588 / 1.021 / 0.771).")
print("  => u = 0.625 is a DEFINITION (nominal strong-coupling value), not a derived")
print("     Bessel coefficient. The verdict is re-checked below under all candidates.")

# ------------------------------------------------------------- [B] free-field acid test
print()
print("=" * 78)
print("GV1 [B] free-field acid test (GU4 CE-1: c_f->1 => u->1/N=0.5, ubar->27.18)")
print("=" * 78)
u_free = mp.mpf(1) / N
print(f"  u_free      = 1/N              = {mp.nstr(u_free, 10)}")
print(f"  ubar(free)  = eDelta*u_free    = {mp.nstr(EDELTA*u_free, 10)}   "
      f"(GU4 CE-1 quotes 27.18)  {'OK' if abs(EDELTA*u_free-27.18) < 0.01 else 'FAIL'}")
print(f"  ubar(2.5)   = eDelta*0.625     = {mp.nstr(EDELTA*u_nom, 10)}   "
      f"(corpus: 33.98)")
print("  In free field eps_T == 0 identically (no T-dependence) while ubar > 1")
print("  and (see [E]) lambda_Ur*u > 1 as well => the criterion lambda_Ur*u<1 is")
print("  SUFFICIENT-ONLY, exactly like ubar<1. Failing it is NOT evidence eps_T != 0.")

# ------------------------------------------------------------------ [C] exact table
print()
print("=" * 78)
print("GV1 [C] exact shell masses (gv1_enum_n6.json)")
print("=" * 78)
d = json.load(open("gv1_enum_n6.json"))
rows = d["rows"]
u = float(u_nom)
Lt = float(20 ** 20 / 19 ** 19)          # = 53.0006865, tree-counting constant (UR6:U2)
print("   n  #clusters      A_n=sum|mu|      avg|mu|   lam^sum=A_n^(1/n)  lam^sup")
lam_sum = []
for r in rows:
    print(f"  {r['n']:2d}  {r['count']:>10d}  {r['A']:>14d}  {r['avg_abs_mu']:9.4f}"
          f"   {r['lam_sum']:14.6f}   {r['lam_sup']:8.6f}")
    lam_sum.append(r["lam_sum"])
print(f"\n  animals growth ratio a_n/a_(n-1): " +
      ", ".join(f"{rows[i]['count']/rows[i-1]['count']:.2f}" for i in range(1, len(rows))))
print(f"  A_n growth ratio A_n/A_(n-1):     " +
      ", ".join(f"{rows[i]['A']/rows[i-1]['A']:.2f}" for i in range(1, len(rows))))

# --------------------------------------------------------------- [D] extrapolation
print()
print("=" * 78)
print("GV1 [D] extrapolation of lambda^sum = limsup m_n^(1/n)  (only n<=6 is exact)")
print("=" * 78)
ns = np.array([r["n"] for r in rows], dtype=float)
ls = np.array(lam_sum, dtype=float)
As = np.array([r["A"] for r in rows], dtype=float)
rho = As[1:] / As[:-1]
n_rho = ns[1:]


def lsq(Amat, y):
    x, *_ = np.linalg.lstsq(Amat, y, rcond=None)
    return x


est = []
# (i) constant (rigorous lower bound: limsup >= lam_6)
est.append(("constant = lam_6 (rigorous LOWER bound)", ls[-1]))
# (ii) lambda_n = a + b/n   (all n>=2, and last 3)
m = ns >= 2
a2 = lsq(np.column_stack([np.ones(m.sum()), 1 / ns[m]]), ls[m])
est.append(("lam_n = a + b/n      (n=2..6)", a2[0]))
m3 = ns >= 4
a3 = lsq(np.column_stack([np.ones(m3.sum()), 1 / ns[m3]]), ls[m3])
est.append(("lam_n = a + b/n      (n=4..6)", a3[0]))
# (iii) lambda_n = a + b/n + c/n^2  (exact 3-point solve n=4,5,6)
Aq = np.column_stack([np.ones(3), 1 / ns[3:], 1 / ns[3:] ** 2])
est.append(("lam_n = a+b/n+c/n^2  (n=4,5,6)", lsq(Aq, ls[3:])[0]))
# (iv) theta scan: lam_n = a - b/n^theta   (theta bounded away from 0: collinear design)
best = None
for th in np.arange(0.50, 3.0, 0.005):
    X = np.column_stack([np.ones(len(ns[1:])), -1 / ns[1:] ** th])
    co = lsq(X, ls[1:])
    pred = X @ co
    sse = float(((pred - ls[1:]) ** 2).sum())
    if best is None or sse < best[0]:
        best = (sse, th, co[0], co[1])
est.append((f"lam_n = a-b/n^theta (theta={best[1]:.3f} fit)", best[2]))
# (v) ratio method: rho_n = a + b/n
ar = lsq(np.column_stack([np.ones(len(n_rho)), 1 / n_rho]), rho)
est.append(("ratio rho_n=A_n/A_(n-1) = a+b/n", ar[0]))
# (vi) ratio theta scan
bestr = None
for th in np.arange(0.05, 3.0, 0.005):
    X = np.column_stack([np.ones(len(n_rho[1:])), -1 / n_rho[1:] ** th])
    co = lsq(X, rho[1:])
    sse = float(((X @ co - rho[1:]) ** 2).sum())
    if bestr is None or sse < bestr[0]:
        bestr = (sse, th, co[0], co[1])
est.append((f"ratio rho_n = a-b/n^theta (theta={bestr[1]:.3f} fit)", bestr[2]))

for name, v in est:
    print(f"  {name:44s} lambda^sum_inf = {v:8.3f}   *u = {v*u:8.3f}")
vals = np.array([v for _, v in est])
print(f"\n  spread of the 6 estimators: [{vals.min():.2f}, {vals.max():.2f}], "
      f"median {np.median(vals):.2f}")
print(f"  RIGOROUS: lambda^sum >= lam_6 = {ls[-1]:.6f}  (limsup >= any element)")
print(f"  RIGOROUS: lambda^sum <= L_t   = {Lt:.6f}  (tree identity, C01/C03; = 33.125*u)")
print(f"  point estimate (median of fits, n<=6 data): lambda^sum = {np.median(vals):.1f} "
      f"+- {0.5*(vals.max()-vals.min()):.1f}")
print(f"  consistent cross-check: C06's rate-ratio reading (see [E]) gives 33.42")

print("\n  [D2] sup-slot probe from gv1_enum [D] (exact mu, clusters containing a K_6):")
for sp in d.get("sup_probe", []):
    nn, mm = sp["n"], sp["abs_mu"]
    print(f"    n={nn:2d}: |mu|={mm:<8d} |mu|^(1/n)={mm**(1.0/nn):.6f}   "
          f"*u = {mm**(1.0/nn)*u:.6f}   {sp['tag']}")

# ----------------------------------------------- [E] GU4 §3.2 rate- vs value-ratio fix
print()
print("=" * 78)
print("GV1 [E] GU4 §3.2: 15.909 is a VALUE ratio, the RATE ratio is 15.909^(1/6)")
print("=" * 78)
val_ratio = 384 / 1.7 ** 6
rate_ratio = 384 ** (1 / 6) / 1.7
Lt = float(20 ** 20 / 19 ** 19)
print(f"  value ratio  384/1.7^6        = {val_ratio:.6f}   <- what GU4:167 used")
print(f"  rate  ratio  384^(1/6)/1.7    = {rate_ratio:.6f}   <- the one that enters lambda")
print(f"  L_t/15.90881 = {Lt/val_ratio:.4f}   (GU4:168 'lambda_Ur ~ 3.3315')")
print(f"  L_t/1.589633 = {Lt/rate_ratio:.4f}   (corrected)   *u = {Lt/rate_ratio*u:.4f}")
print(f"  => GU4:168 underestimates lambda_Ur by the factor "
      f"{val_ratio/rate_ratio:.3f} = 15.909^(5/6);")
print(f"     corrected reading lambda_Ur*u ~ {Lt/rate_ratio*u:.2f} sits on top of the")
print(f"     direct enumeration (lam_6 = {ls[-1]:.2f}, extrapolated ~{np.median(vals):.0f}).")

# ------------------------------------------------------------------ [F] verdict
print()
print("=" * 78)
print("GV1 [F] verdict arithmetic")
print("=" * 78)
lam6 = ls[-1]
sup6 = rows[-1]["lam_sup"]
sup8 = max([sup6] + [sp["abs_mu"] ** (1.0 / sp["n"]) for sp in d.get("sup_probe", [])])
for label, lam in (("lambda^sup @n=6 (full enum)", sup6),
                   ("lambda^sup floor @n=8 (probe)", sup8),
                   ("lambda^sum (per site) @n=6", lam6),
                   ("lambda^sum extrapolated", float(np.median(vals)))):
    print(f"  {label:26s} = {lam:8.4f}   *u = {lam*u:8.4f}   "
          f"vs 1: {'CLOSED' if lam*u > 1 else 'open':6s} by {lam*u:6.3f}x   "
          f"vs 1.6: {lam/1.6:6.3f}x   vs 0.5595 (ring 1): {lam/0.5595:7.3f}x")
print()
print("  under the other u-candidates (lambda values unchanged, only 1/u moves):")
for name, v in cands:
    print(f"    {name:56s} 1/u = {float(1/v):7.4f}   "
          f"lambda^sup={float(sup6):6.3f} -> {'CLOSED' if sup6 > float(1/v) else 'open'}"
          f"   lambda^sum={lam6:6.3f} -> {'CLOSED' if lam6 > float(1/v) else 'open'}")

# ------------------------------------------------------------------ [G] certificate
print()
print("=" * 78)
print("GV1 [G] why more n cannot REVERSE the verdict (named certificate)")
print("=" * 78)
print("  A finite-n enumeration gives LOWER bounds on a limsup, never an upper bound.")
print("  Ring 3 open needs lambda_Ur < 1.6 = an UPPER bound valid for ALL clusters.")
print("  => no n can open it; only a structural bound |c_Gamma| <= C^n with C < 1.6 can,")
print("     which is exactly what C06:187 hoped for and what n<=6 refutes.")
print("  Cost of the full enumeration (measured: n=6 -> 7.218768e6 sets, 177 s DP):")
a6 = rows[-1]["count"]
r = a6 / rows[-2]["count"]
for k in (7, 8, 9, 10):
    est_cnt = a6 * r ** (k - 6)
    print(f"    n={k:2d}: ~{est_cnt:9.3e} clusters  ->  level-set memory ~"
          f"{est_cnt*140/1e9:8.2f} GB   DP time ~{est_cnt*24.5e-6/60:9.1f} min")
print("  => n=7 is memory-marginal here (16.7 GB free, one 0.88 GB job running),")
print("     n>=8 needs an orderly-generation (Redelmeier) enumerator in C. Not done.")
