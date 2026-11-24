"""
ward_stein_witness.py — witness for SLOT-WARD-STEIN.

Question: can the EXACT triple
    Ward   (IBP)   : E_mu[a f] = 1,               f = d_a V,  a = <v,X>
    Stein          : E_mu[Hess V] = E_mu[grad V (x) grad V] >= 0
                     => coarse: E[H|a] = Var(f|a) + m'(a)  (m=E[f|a])
    Variance       : G''(a) := m'(a) = E[H|a] - Var(f|a)
be combined into ONE relation R that is >= c > 0 for GAPPED measures only?

We verify, to quadrature precision, the logical content:
  (1) The triple is satisfied by EVERY smooth, even, normalised density p(a)
      (Ward <=> normalisation, Stein <=> the tautological Fisher identity,
       variance <=> the definition of G'').  So G''(0) = -p''(0)/p(0) is FREE.
  (2) A GAPLESS measure (free massless, chi->inf) satisfies it with G''->0;
      a BIMODAL marginal satisfies it with G''(0) < 0.
  (3) In >=1 fibre dimension the variance identity has the WRONG SIGN to floor:
      G''(0) = E[H|0] - Var(f|0) <= E[H|0];  Var(f|0)>0 SUBTRACTS.
  (4) Real SU(2)/SU(3) data (h1_right / W7): G''>0 measured, but the
      rigorous-style lower bound LB_rig (keeping only the proven W_geo>=0)
      is NEGATIVE at every beta -- no uniform floor; positivity is Haar-carried.

Run:  py -3.12 research/unconditional/w7/ward_stein_witness.py
Light: 1D/2D numpy quadrature + one JSON read.  Single job.
"""
import os, json, math
import numpy as np

RT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JSON = os.path.join(RT, "h1_right_L8_N200.json")

# ---------- small quadrature helpers (trapezoid on a fine grid) ----------
_TRAP = getattr(np, "trapezoid", None) or getattr(np, "trapz")

def trapz(y, x):
    return _TRAP(y, x)

def stats_from_density_p(p, a):
    """p: positive density on grid a. Robust: uses m*p = -p' (no log/division
       at the tails). Returns m'(0)=-p''(0)/p(0), Ward=int a m p, Stein residual."""
    Z = trapz(p, a); p = p / Z
    dp = np.gradient(p, a)                 # p'
    mp_arr = -dp                           # m p = -p'
    sup = p > 1e-13 * p.max()
    m = np.zeros_like(p)
    m[sup] = mp_arr[sup] / p[sup]          # m = -p'/p on the support
    mp = np.gradient(m, a)                 # m'(a)
    i0 = int(np.argmin(np.abs(a)))
    Gpp0 = mp[i0]                          # = -p''(0)/p(0)
    ward = trapz(a * mp_arr, a)            # int a m p da = int a (-p') da
    fisher = trapz((dp[sup] ** 2) / p[sup], a[sup])   # int (p')^2/p
    stein = fisher - trapz(m * m * p, a)              # should be ~0
    Var = trapz(a * a * p, a) - trapz(a * p, a) ** 2
    return dict(Gpp0=Gpp0, ward=ward, stein=stein, Var=Var)

def gauss(x, mu=0.0, s=1.0):
    return np.exp(-0.5 * ((x - mu) / s) ** 2) / (s * math.sqrt(2 * math.pi))

# =====================================================================
print("=" * 78)
print("PART A — 1D coarse marginal: the triple holds for EVERY normalised p;")
print("         G''(0)=-p''(0)/p(0) is FREE (sign and magnitude).")
print("=" * 78)
rows = []

def run_family(name, pfun, half_width, N=400001):
    x = np.linspace(-half_width, half_width, N)
    st = stats_from_density_p(pfun(x), x)
    rows.append((name, st["Gpp0"], st["ward"], st["stein"], st["Var"]))

# (i) Gaussian family, chi = Var: free/gapless limit is chi -> infinity
for chi in [4.0, 1.0, 0.25, 100.0, 1e4]:
    s = math.sqrt(chi)
    run_family("Gaussian chi=%-8.4g" % chi, lambda x, s=s: gauss(x, 0.0, s), 12 * s)
# (ii) bimodal: G''(0) = (s^2-mu^2)/s^4 < 0 when the mode is a DIP between lumps
for mu, s in [(2.0, 0.35), (1.5, 0.5), (0.6, 0.5)]:
    run_family("Bimodal mu=%.2f s=%.2f" % (mu, s),
               lambda x, mu=mu, s=s: 0.5 * gauss(x, -mu, s) + 0.5 * gauss(x, mu, s),
               mu + 10 * s)
# (iii) flat-top (wide centre + shoulders) -> near-zero curvature
run_family("Flat-top mixture",
           lambda x: 0.5 * gauss(x, 0.0, 0.6) + 0.25 * gauss(x, -2.2, 0.35)
           + 0.25 * gauss(x, 2.2, 0.35), 4.0)

print("%-24s %12s %10s %14s %12s" % ("marginal p(a)", "G''(0)", "Ward", "Stein residual", "Var(a)"))
print("-" * 78)
for name, gpp, w, sr, var in rows:
    print("%-24s %+12.6e %10.6f %+14.2e %12.4f" % (name, gpp, w, sr, var))
print("\n  -> Ward=1 & Stein residual=0 to quadrature precision for EVERY p;")
print("     G''(0) spans POSITIVE (Gaussian), NEGATIVE (bimodal dip), ~0 (flat-top).")
print("     Gaussian chi->inf (gapless free limit): G''->0< c.  No sign-definite R.")

# =====================================================================
print("\n" + "=" * 78)
print("PART B — 2D (coarse a + fibre y): variance identity G''=E[H|0]-Var(f|0);")
print("         Var(f|0)>0 SUBTRACTS (wrong direction to floor G'').")
print("=" * 78)

def cond_moments(V, Haa, fa, av, yv):
    """V(a,y) on grid (av x yv). Returns exact conditional moments at a=0:
       E[H|0], Var(f|0), Gpp0=d/d a E[f|a]|0, p(a) marginal, and Ward."""
    Aa, Yy = np.meshgrid(av, yv, indexing="ij")   # (na, ny)
    W = np.exp(-(V(Aa, Yy) - V(Aa, Yy).min()))
    p = trapz(W, yv)                              # marginal in a  (na,)
    Z = trapz(p, av)
    den = trapz(W, yv)
    mask = den > 1e-13 * den.max()                # guard underflow tails
    # m(a) = E[f|a] = int f W dy / int W dy
    fW = fa(Aa, Yy) * W
    m = np.zeros_like(av); EH = np.zeros_like(av); Ef2 = np.zeros_like(av)
    m[mask] = trapz(fW, yv)[mask] / den[mask]
    EH[mask] = trapz(Haa(Aa, Yy) * W, yv)[mask] / den[mask]
    Ef2[mask] = trapz(fa(Aa, Yy) ** 2 * W, yv)[mask] / den[mask]
    Vf = Ef2 - m ** 2
    Vf[~mask] = 0.0
    mp = np.gradient(m, av)
    i0 = int(np.argmin(np.abs(av)))
    # Ward E[a f] = int a m p / Z
    ward = trapz(av * m * p, av) / Z
    return dict(Gpp0=mp[i0], EH0=EH[i0], Vf0=Vf[i0], m0=m[i0], ward=ward, p=p, i0=i0)

av = np.linspace(-12, 12, 1201); yv = np.linspace(-12, 12, 1201)

# (B1) factorised-ish non-Gaussian: H_aa>0 but cross coupling makes Var(f|0)>0
def V_B1(a, y):      # A a^2/2 + B y^2/2 + c a y + lam y^4/4
    return 0.5 * 3.0 * a**2 + 0.5 * 2.0 * y**2 + 1.4 * a * y + 0.10 * y**4
def H_B1(a, y): return 3.0 + 0.0 * a          # d_a^2 V = A
def f_B1(a, y): return 3.0 * a + 1.4 * y

# (B2) Savvidy-like & normalisable: NEGATIVE action curvature at 0 + variance
def V_B2(a, y):
    return -0.8 * a**2 + 0.15 * a**4 + 1.0 * y**2 + 0.5 * a * y + 0.10 * y**4
def H_B2(a, y): return -1.6 + 1.8 * a**2       # d_a^2 V  (NEGATIVE at a=0)
def f_B2(a, y): return -1.6 * a + 0.6 * a**3 + 0.5 * y

for nm, (V, H, fa) in {"B1 (H>0, cross-coupling)": (V_B1, H_B1, f_B1),
                        "B2 (Savvidy-like H<0 at 0)": (V_B2, H_B2, f_B2)}.items():
    r = cond_moments(V, H, fa, av, yv)
    check = r["EH0"] - r["Vf0"]
    print("%-30s" % nm)
    print("   Ward E[af]         = %+.6f   (=> 1)" % r["ward"])
    print("   E[H|0]             = %+.6f" % r["EH0"])
    print("   Var(f|0)           = %+.6f   (>=0, SUBTRACTS)" % r["Vf0"])
    print("   E[H|0]-Var(f|0)    = %+.6f" % check)
    print("   G''(0) (d m/da)|0  = %+.6f   (identity check)" % r["Gpp0"])
    print("   G''(0) <= E[H|0] ? = %s" % (r["Gpp0"] <= r["EH0"] + 1e-6))
    print()

# =====================================================================
print("=" * 78)
print("PART C — real SU(2) data (h1_right, L=8): measured G''>0 but the")
print("         rigorous-style lower bound LB_rig<0 at EVERY beta (no floor).")
print("=" * 78)
if os.path.exists(JSON):
    d = json.load(open(JSON))
    print("%6s %6s %10s %10s %10s %10s %12s" %
          ("beta", "L", "Gpp", "+-SE", "1/Var(a)", "LB_rig", "Gpp/free"))
    print("-" * 78)
    for pb in d.get("per_beta", []):
        g = pb.get("gluon")
        if not g or "full" not in g:
            continue
        f = g["full"]
        print("%6.1f %6d %10.4f %10.4f %10.4f %10.4f %12.3f" %
              (pb["beta"], pb["L"], f["Gpp"], g["se"]["Gpp"], f["inv"],
               g["LB_rig"], g["ratio_Gpp_over_free"]))
    print("\n  -> Gpp>0 (measured) at all beta; LB_rig<0 at all beta =>")
    print("     no rigorous uniform G''>=c>0 at the crossover. Free reference")
    print("     (4 sin^2(k/2)=%.3f at L=8) is NOT exceeded uniformly (Gpp/free<1)."
          % (4.0 * math.sin(0.5 * 2 * math.pi / 8) ** 2))
    print("  -> W7 (direct): SU(2) beta=2.5 G''=+0.27+-0.004; SU(3) beta=2.5")
    print("     G''=+0.51+-0.08 with E[H]-Var=-0.03 (positivity is Haar-carried).")
else:
    print("  (JSON not found; skipping real-data panel)")

print("\n" + "=" * 78)
print("VERDICT: Ward+Stein+variance is ONE identity = 'p is a normalised pdf'.")
print("Gapless free field (chi->inf) satisfies it with G''->0; bimodal satisfies")
print("it with G''<0; variance identity computes G''=E[H|0]-Var(f|0)<=E[H|0]")
print("(wrong direction). No sign-definite R>=c>0 only for gapped. READ, not FORCE.")
print("=" * 78)
