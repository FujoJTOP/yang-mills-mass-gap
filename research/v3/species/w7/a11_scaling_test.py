# -*- coding: utf-8 -*-
r"""
a11_scaling_test.py — A11: is "two points agree at xi/a ~ 1" weak evidence?

Question (single): at xi/a = 1/(m a) ~ 1 the A8 two-point agreement
   m/Lambda = 612.5 +- 50.2  (beta=6.00)   vs   709.4 +- 127.3  (beta=6.30)
is consistent at 0.71 sigma.  How much of that is accident?

Inputs are A8-MASS-LAMBDA.md's own four numbers (E).
Everything below is arithmetic + two explicit approximation formulas.  No new dynamics.

Formulas used (written out, leading order only):
  (F1) strong-coupling (SC) leading 0++ mass, 4D SU(N) Wilson, four-wall tube:
         m a = 4 ln(1/u) + O(u),  u = fundamental char. coeff ~ beta/(2N)  [or beta/(2N^2)]
       => m a = -4 ln beta + const  =>  d(m a)/d ln beta = -4  EXACTLY, independent of
       the constant inside the log (6/beta or 18/beta) -- so the slope is convention-free. [T]
       Absolute normalisation: m a = 4 ln(2N/beta) zero at beta = 2N = 6 (SU(3)); at
       beta >= 6 the leading formula is already outside its domain (sign flip). [M]
  (F2) two-loop running a(beta):  a*Lambda = exp(-1/(2 b0 g^2)) (2 b0 g^2)^(-b1/(2 b0^2)),
       g^2 = 2N/beta = 6/beta;  b0 = 11N/(48 pi^2), b1 = 34N^2/(768 pi^4).   [C]
       Scaling hypothesis m = const  =>  m a(beta) = const * a(beta).

Tags: [T] theorem/standard, [C] computed-with-formula, [E] corpus measurement,
      [M] model/estimate, [NG] not obtainable from the corpus.

Run:  py -3.12 research\v3\species\w7\a11_scaling_test.py
"""
import sys, math
import numpy as np
from scipy.stats import chi2 as chi2d

sys.stdout.reconfigure(encoding="utf-8")
out = []
def P(s=""):
    out.append(s)
    print(s)

def tag(t): return "[" + t + "]"

# ---------------------------------------------------------------- inputs [E]
N = 3
d1 = dict(beta=6.00, ma=1.0730, sma=0.0880, mLam=612.5, sLam=50.2)
d2 = dict(beta=6.30, ma=0.8860, sma=0.1590, mLam=709.4, sLam=127.3)
for d in (d1, d2):
    d["xi"] = 1.0 / d["ma"]

P("=" * 78)
P("A11 · two-point consistency at xi/a ~ 1: is it evidence for scaling?")
P("=" * 78)
P()
P("[A] inputs, straight from A8-MASS-LAMBDA.md:51-52  " + tag("E"))
P(f"    beta={d1['beta']:.2f}  m a={d1['ma']:.4f}+-{d1['sma']:.4f}  xi/a={d1['xi']:.4f}"
  f"   m/Lambda={d1['mLam']:.1f}+-{d1['sLam']:.1f}")
P(f"    beta={d2['beta']:.2f}  m a={d2['ma']:.4f}+-{d2['sma']:.4f}  xi/a={d2['xi']:.4f}"
  f"   m/Lambda={d2['mLam']:.1f}+-{d2['sLam']:.1f}")
P(f"    => both points have xi/a in [{d1['xi']:.2f},{d2['xi']:.2f}]: O(1), no scale separation.")
P()

# ------------------------------------------------------- [B] joint + chi^2
P("[B] joint value and chi^2 (1 dof) on m/Lambda  " + tag("C"))
w1, w2 = 1.0 / d1["sLam"] ** 2, 1.0 / d2["sLam"] ** 2
joint = (d1["mLam"] * w1 + d2["mLam"] * w2) / (w1 + w2)
sjoint = 1.0 / math.sqrt(w1 + w2)
diff = d1["mLam"] - d2["mLam"]
sdiff = math.hypot(d1["sLam"], d2["sLam"])
z = diff / sdiff
chi2 = z * z
p = chi2d.sf(chi2, 1)
P(f"    weights w_i=1/sigma_i^2 -> joint  m/Lambda = {joint:.1f} +- {sjoint:.1f}")
P(f"    difference = {diff:+.1f} +- {sdiff:.1f}  =  {z:+.3f} sigma")
P(f"    chi^2(1 dof) = {chi2:.4f}   p-value = {p:.4f}   (p = erfc(sqrt(chi^2/2)) = "
  f"{math.erfc(math.sqrt(chi2 / 2)):.4f})")
P(f"    verdict on TENSION only (no tension): p = {p:.2f} > 0.30 -> branch (a) is alive.")
P()

# ---------------------------------------------- [C] the SAME test as a ratio
P("[C] the identical test written as a ratio (this is what 'scaling' means)  " + tag("C"))
r = d2["ma"] / d1["ma"]
sr = r * math.hypot(d1["sma"] / d1["ma"], d2["sma"] / d2["ma"])
P(f"    measured   r_meas = m a(6.30)/m a(6.00) = {r:.4f} +- {sr:.4f}  ({100*sr/r:.1f}%)")
P(f"    m/Lambda consistency  <=>  r_meas = a(6.30)/a(6.00) (two-loop)  -- same test, same number.")
P()

# --------------------------------------------------- [D] three model ratios
P("[D] three predictions for r, all computed, none fitted to r  " + tag("C") + tag("M"))

# (D1) strong-coupling leading ABSOLUTE form, F1 at the two betas
sc_abs1 = 4.0 * math.log(2.0 * N / d1["beta"])
sc_abs2 = 4.0 * math.log(2.0 * N / d2["beta"])
P(f"    D1 SC absolute  m a = 4 ln(2N/beta):  beta=6.00 -> {sc_abs1:+.5f} ; "
  f"beta=6.30 -> {sc_abs2:+.5f}")
P(f"       ratio = undefined / sign flip: the leading SC formula has its zero exactly at")
P(f"       beta = 2N = {2*N}.  At these two betas the SC expansion is at/past its edge.  {tag('M')}")

# (D2) SC SLOPE (the convention-free content), anchored with ONE data number m a(6.00)
#      two equivalent readings of the same leading law m a = -4 ln beta + const
d_ma = 4.0 * math.log(d2["beta"] / d1["beta"])
pred_lin = d1["ma"] - d_ma            # additive in m a
pred_log = d1["ma"] * math.exp(-d_ma / d1["ma"])   # multiplicative in m a
r_sc_lin = pred_lin / d1["ma"]
r_sc_log = pred_log / d1["ma"]
P(f"    D2 SC slope  d(m a)/d ln beta = -4 (exact at leading order, F1):")
P(f"       Delta m a = 4 ln(6.30/6.00) = {d_ma:.5f}")
P(f"       additive   -> m a(6.30) = {pred_lin:.4f}  => r_SC = {r_sc_lin:.4f}")
P(f"       log-slope  -> m a(6.30) = {pred_log:.4f}  => r_SC = {r_sc_log:.4f}")
P(f"       (one anchor number used: m a(6.00); the slope -4 is parameter-free.)")

# (D3) two-loop scaling, F2
b0 = 11.0 * N / (48.0 * math.pi ** 2)
b1 = 34.0 * N ** 2 / (768.0 * math.pi ** 4)
def La(beta):
    g2 = 2.0 * N / beta
    return math.exp(-1.0 / (2 * b0 * g2)) * (2 * b0 * g2) ** (-b1 / (2 * b0 ** 2))
r_scal = La(6.30) / La(6.00)
P(f"    D3 two-loop scaling  a(6.30)/a(6.00) = {r_scal:.4f}"
  f"   (A8 quotes 0.713 from its own Lambda a values 1.7518e-3/1.2489e-3 = {1.2489/1.7518:.4f})")

# (D4) null: no beta-dependence at all (pure cutoff mass, xi/a frozen at O(1))
r_null = 1.0
P(f"    D4 null / cutoff  r = {r_null:.4f}  (m a set by the lattice only, no running)")
P()

# -------------------------------------------------------- [E] discrimination
P("[E] how far is each prediction from the datum, and can the datum separate them?  "
  + tag("M"))
rows = [
    ("SC slope (additive)", r_sc_lin),
    ("SC slope (log)", r_sc_log),
    ("two-loop scaling", r_scal),
    ("null (no running)", r_null),
]
P(f"    {'model':>22} {'r_pred':>8} {'r_pred-r_meas':>14} {'sigmas':>8} {'chi^2':>8} {'p':>7}")
for name, rp in rows:
    dd = rp - r
    zz = dd / sr
    c2 = zz * zz
    P(f"    {name:>22} {rp:8.4f} {dd:+14.4f} {zz:+8.3f} {c2:8.4f} {chi2d.sf(c2,1):7.4f}")
P()
sep_ss = abs(r_sc_log - r_scal)
P(f"    model separation |r_SC - r_scaling| = {sep_ss:.4f}  =  {sep_ss/sr:.3f} sigma of the")
P(f"    present measurement  =>  the present data cannot tell them apart (< 1 sigma).")
P(f"    sigma_r needed for 3 sigma separation: {sep_ss/3:.4f} ({100*(sep_ss/3)/r:.1f}%),")
P(f"    i.e. a factor {sr/(sep_ss/3):.1f} smaller than today's {sr:.3f}.")
P()
P("    Log-slope form (scheme-free, no anchor at all):")
s_meas = math.log(r) / math.log(d2["beta"] / d1["beta"])
s_sig = (sr / r) / math.log(d2["beta"] / d1["beta"])
s_scal = math.log(r_scal) / math.log(d2["beta"] / d1["beta"])
s_sc = -4.0 / d1["ma"]
P(f"      s := d ln(m a)/d ln beta = {s_meas:+.2f} +- {s_sig:.2f}   (error > 100%: slope unmeasured)")
P(f"      scaling predicts {s_scal:+.2f} ; SC predicts -4/(m a) = {s_sc:+.2f} ; null predicts 0.00")
P(f"      all three sit within |z| < 1.1 sigma of the datum.")
P()
P("    Power check -- run the SAME test on a world with zero scaling (r = 1 exactly):")
c2n = (r_null - r) ** 2 / sr ** 2
P(f"      it returns chi^2 = {c2n:.3f}, p = {chi2d.sf(c2n,1):.3f}  -- nearly the same p as the")
P(f"      scaling hypothesis ({chi2d.sf(((r_scal-r)/sr)**2,1):.3f}).")
P(f"    Likelihood ratios on the same two numbers:")
lr_scal_null = math.exp((c2n - ((r_scal - r) / sr) ** 2) / 2)
lr_sc_null = math.exp((c2n - ((r_sc_log - r) / sr) ** 2) / 2)
P(f"      L(scaling)/L(null)  = {lr_scal_null:.2f}   (needs ~10 for 'moderate')")
P(f"      L(SC slope)/L(null)= {lr_sc_null:.2f}   (also 'not worth mentioning')")
P(f"    Fluke probability, stated as asked:")
P(f"      if the lattice mass had NO beta-running at all (r = 1), the chance of landing")
P(f"      at least {abs(r-1)/sr:.2f} sigma away from r=1 is p = {chi2d.sf(c2n,1):.3f}  (~{100*chi2d.sf(c2n,1):.0f}%).")
P(f"      So 'the two points agree' is reproduced by a zero-scaling world in ~1 case out of 3.5;")
P(f"      the observed agreement is not rare information about scaling.")
P()

# ------------------------------------- [F] what would a real lever arm look like
P("[F] lever arm: where do the two hypotheses actually separate?  " + tag("M"))
P(f"    {'beta':>6} {'xi/a(SC)':>9} {'m a SC':>9} {'m a scal':>9} {'diff':>7} {'sigmas@8%':>10}")
ma0 = d1["ma"]
for bt in (6.00, 6.30, 6.60, 6.90, 7.20):
    msc = ma0 - 4.0 * math.log(bt / 6.00)
    msl = ma0 * La(bt) / La(6.00)
    xi = 1.0 / max(msc, 1e-9)
    dif = msc - msl
    # assume an 8% measurement error on m a (optimistic for a glueball at xi/a >~ 1.5)
    P(f"    {bt:6.2f} {xi:7.2f} {msc:9.4f} {msl:9.4f} {dif:+7.4f} {abs(dif)/(0.08*msc):10.2f}")
P("    8% is the *statistical* scale A8 already has at beta=6.00; at beta=6.30 the corpus has 18%.")
P("    Read: only past xi/a ~ 2 do the two laws pull apart by >= 3 sigma at fixed precision.")
P()

# solve for beta where each model reaches m a = 0.5  (xi/a = 2)
def solve(f, lo=6.0, hi=8.0):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.5:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
b_sc = solve(lambda b: ma0 - 4.0 * math.log(b / 6.0))
b_sl = solve(lambda b: ma0 * La(b) / La(6.00))
P(f"    beta at which xi/a = 2 (m a = 0.5):  SC says {b_sc:.2f} ; scaling says {b_sl:.2f}")
P(f"    => a discriminating pair is beta ~ {min(b_sc,b_sl):.1f} and {max(b_sc,b_sl):.1f}, both with")
P(f"       an INDEPENDENT scale (r0/a or sqrt(sigma) a) so xi/a is not circular.  {tag('NG')} in corpus")
P()

# ------------------------------------------------ [G] corpus sanity / caveats
P("[G] why the 8%/18% errors should not be read as gaussian  " + tag("E"))
P("    research/probe-m/REPORT.md:137-146 (beta=6.2, L=10 pilot): the SAME operator with")
P("    ladder level n=0,8,16,32 gives xi/a(t<=3) = 0.41, 0.99, 1.18, 1.42 -- a factor 3.5")
P("    spread at ONE beta.  At xi/a ~ 1 the 'mass' is not yet an operator-independent")
P("    observable, so 8% and 18% are lower bounds on the systematic error.")
P("    A8-MASS-LAMBDA.md:40 itself notes the truncation spread (one-loop vs two-loop) is 30%.")
P()
P("=" * 78)
P("VERDICT: branch (a).  p = %.2f (no tension), and the strong-coupling slope" % p)
P("reproduces r = 0.82-0.83 vs measured 0.826 (0.05 sigma), i.e. equally well as")
P("two-loop scaling (0.713, 0.76 sigma).  Two-point agreement at xi/a ~ 1 is WEAK")
P("evidence, not a scaling test.  Strong evidence needs two points with measured")
P("xi/a >= 2 (beta ~ 6.9-7.0) and an independent scale at each.")
P("=" * 78)

with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
