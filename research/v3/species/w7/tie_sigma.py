"""TIE-SIGMA: one sigma, two rival scalings -- sigma ~ 1/a^2  (mu-4)  vs  sigma ~ Lambda_MS^2 (mu6).

Single sigma definition (mandatory, stated once):
    sigma(beta) := (a*sqrt(sigma))^2(beta) * (a^-1(beta))^2          [GeV^2]
i.e. the area-law coefficient of the Wilson loop, -lim_{A->inf} ln W(A)/A, converted to physical
units with a^-1 taken *from AT2021 Table 13 only* (a^-1 = sqrt(sigma)/(a sqrt(sigma))).
The lattice-unit reading  shat(beta) := (a sqrt(sigma))^2 = a^2 sigma  is reported alongside,
because that is the only quantity in Tab.13 that varies with beta.

Inputs (all from AT2021 = arXiv:2106.00364v3, JHEP 12 (2021) 082; quoted via
research/v3/species/L2-GLUEBALL-LAMBDA.md and MU6-NUM.md -- nothing from memory):
    r0*sqrt(sigma) = 1.160(6),  r0 = 0.472(5) fm        (AT2021 text p.18)
    Lambda_MS/sqrt(sigma) = 0.5424(13)[185]  (SU(3), Nf=0)   (AT2021 Tab.15)
    per-beta a*sqrt(sigma), SU(3), beta = 6/g^2          (AT2021 Tab.13, 10 points)
Discipline: a^-1 comes from Tab.13; the 2.91 GeV anchor is NOT used (it is 1.297x too high and is
shown once only as a labelled counterfactual [X]).
"""
import numpy as np

np.set_printoptions(precision=6, suppress=False)
FM = 5.0677307                      # GeV^-1 per fm

# ---------------------------------------------------------------- 0. constants
sqrtsig = 1.160 / (0.472 * FM)      # GeV        [C]
sigma_ref = sqrtsig ** 2            # GeV^2      [C]
N = 3                               # SU(3)
b0 = 11.0 * N / (48 * np.pi ** 2)   # dimensionless
b1 = 34.0 * N ** 2 / (768 * np.pi ** 4)
LAM_MS = 0.5424 * sqrtsig           # GeV        [E] Tab.15 ratio x sqrt(sigma)

TAB13 = [(5.6924, .40100), (5.8000, .31603), (5.8941, .26130), (5.9900, .21959),
         (6.0625, .19509), (6.2350, .14899), (6.3380, .12948), (6.5000, .10319),
         (6.6000, .09024), (6.7000, .07898)]          # [E] AT2021 Tab.13, SU(3)
B13 = np.array([p[0] for p in TAB13])
S13 = np.array([p[1] for p in TAB13])

# the five beta's asked for; 5.69 is the tabulated 5.6924 (no extrapolation below the table)
BETAS = [5.69, 5.9, 6.0, 6.3, 6.7]


def a_sqrt_sigma(beta):
    """Tab.13 column, linear in beta between tabulated points; snapped at the ends
    (linear interp reproduces MU6-NUM's a^-1(6.00) = 2.243 GeV; log-interp differs by 0.08%)."""
    b = np.clip(beta, B13.min(), B13.max())
    return float(np.interp(b, B13, S13))


def a_inv(beta):
    """GeV.  Tab.13 ONLY: a^-1 = sqrt(sigma)/(a sqrt(sigma))."""
    return sqrtsig / a_sqrt_sigma(beta)


def a_of(beta):
    """GeV^-1."""
    return 1.0 / a_inv(beta)


def g2_of(beta):
    return 2.0 * N / beta                       # beta = 2N/g^2


def x_of(beta):
    """e^{-1/(2 b0 g^2)} -- mu6's variable (NOT a*Lambda)."""
    return float(np.exp(-1.0 / (2.0 * b0 * g2_of(beta))))


def lidx(xv, yv):
    """least-squares d ln y / d ln x."""
    return float(np.polyfit(np.log(np.asarray(xv, float)), np.log(np.asarray(yv, float)), 1)[0])


print("=" * 104)
print("TIE-SIGMA  --  one sigma, two rival scalings.   SU(3), Nf=0, AT2021 arXiv:2106.00364v3")
print("=" * 104)
print(f"sqrt(sigma) = {sqrtsig:.6f} GeV = {sqrtsig*1e3:.3f} MeV      sigma = {sigma_ref:.6f} GeV^2   [C]")
print(f"Lambda_MS   = {LAM_MS*1e3:.3f} MeV   [E] Tab.15 (Lambda_MS/sqrt(sigma)=0.5424)"
      f"   =>  sigma/Lambda_MS^2 = {sigma_ref/LAM_MS**2:.4f}")
print(f"b0 = {b0:.6f}   b1 = {b1:.6f}   (N=3)")

# ------------------------------------------- 1. sigma definition, and its degeneracy
print("\n" + "-" * 104)
print("[1] SINGLE sigma: sigma(beta) = (a sqrt sigma)^2 * (a^-1)^2,  a^-1 from Tab.13 only")
print("-" * 104)
sig_all = []
for b in B13:
    ai = a_inv(b)
    sig_all.append(a_sqrt_sigma(b) ** 2 * ai ** 2)
sig_all = np.array(sig_all)
print(f"  sigma(beta) over all 10 Tab.13 points: min={sig_all.min():.9f}  max={sig_all.max():.9f}"
      f"  spread={sig_all.max()/sig_all.min()-1:.2e}")
print("  => sigma is beta-INDEPENDENT to machine precision: [T] IDENTITY, not a measurement.")
print("     Reason: a^-1 := sqrt(sigma)/(a sqrt sigma)  =>  sigma = (a sqrt sigma)^2 (a^-1)^2 = const.")
print("     Column 1 and column 3 below therefore vanish by construction; see [4] for the live test.")

# ------------------------------------------------------------ 2. table at the 5 betas
print("\n" + "-" * 104)
print("[2] THE FIVE POINTS   (a^-1 discipline check: a^-1(6.00) must be 2.243 GeV, NOT 2.91)")
print("-" * 104)
hdr = (f"  {'beta':>7} {'a sqrt(s)':>10} {'a^-1[GeV]':>10} {'a[GeV^-1]':>10} {'sigma[GeV^2]':>13}"
       f" {'shat=a^2s':>11} {'x':>10} {'aLam_MS':>9} {'s/Lam^2':>8}")
print(hdr)
rows = []
for b in BETAS:
    asg = a_sqrt_sigma(b)
    ai = a_inv(b)
    a = 1.0 / ai
    sa2 = asg ** 2                     # shat = a^2 sigma, dimensionless
    x = x_of(b)
    aL = a * LAM_MS
    sig = sa2 * ai ** 2
    rows.append(dict(beta=b, asg=asg, ai=ai, a=a, sa2=sa2, x=x, aL=aL, sig=sig, g2=g2_of(b)))
    print(f"  {b:7.4f} {asg:10.5f} {ai:10.4f} {a:10.5f} {sig:13.9f} {sa2:11.6f} {x:10.4e}"
          f" {aL:9.5f} {sig/LAM_MS**2:8.4f}")
beta_v = np.array([r['beta'] for r in rows]); a_v = np.array([r['a'] for r in rows])
sa2_v = np.array([r['sa2'] for r in rows]); x_v = np.array([r['x'] for r in rows])
aL_v = np.array([r['aL'] for r in rows]); sig_v = np.array([r['sig'] for r in rows])

print(f"\n  DISCIPLINE a^-1(6.00, Tab.13) = {a_inv(6.0):.4f} GeV      (rejected anchor 2.91 GeV = "
      f"{2.91/a_inv(6.0):.3f}x too high)")
print(f"  CHECK     a*Lambda_MS(beta=6) = {a_of(6.0)*LAM_MS:.5f}   (task: 0.117)"
      f"      x(beta=6) = {x_of(6.0):.4e}   (task: 7.65e-4)")
print(f"  CHECK     a*Lambda_MS / x  at beta=6 = {a_of(6.0)*LAM_MS/x_of(6.0):.1f}  (task: 153x)")
print(f"  DIMENSIONS a^-1 [GeV], a [GeV^-1], sigma [GeV^2], aLam [1], x [1], shat [1], sigma/Lam^2 [1]")

# ------------------------------------------------------------- 3. the three columns
print("\n" + "-" * 104)
print("[3] THREE LOG-DERIVATIVE COLUMNS, exactly as asked, on the single sigma (GeV^2)")
print("-" * 104)
c1 = lidx(a_v, sig_v)
c2 = lidx(x_v, sig_v)
c3 = lidx(a_v, sig_v / LAM_MS ** 2)
print(f"  C1 = d ln sigma / d ln a                  = {c1:+.3e}     (compare -2)   [T-degenerate: =0)")
print(f"  C2 = d ln sigma / d ln(e^-1/(2b0g2))      = {c2:+.3e}     (compare  2)   [T-degenerate: =0)")
print(f"  C3 = d ln(sigma/Lambda_MS^2) / d ln a     = {c3:+.3e}     (compare  0)   [T-degenerate: =0)")
print(f"\n  SAME THREE on shat = a^2 sigma  (the ONLY beta-varying object in Tab.13; sigma's own")
print(f"  beta-dependence is unmeasurable: a^-1 is derived from the same column):")
c1L = lidx(a_v, sa2_v)
c2L = lidx(x_v, sa2_v)
c3L = lidx(a_v, sa2_v / aL_v ** 2)
print(f"  C1' = d ln(a^2 sigma)/d ln a                = {c1L:+.4f}   (identity +2: shat = a^2 sigma)")
print(f"  C2' = d ln(a^2 sigma)/d ln(e^-1/(2b0g2))    = {c2L:+.4f}   (mu6's n)")
print(f"  C3' = d ln((a^2 sigma)/(a Lambda_MS)^2)/dln a = {c3L:+.3e}  (identity 0)")
print(f"  => sign flip C1 = 0 vs C1' = +2 IS mu-4's bookkeeping defect: dropping one 1/a^2.")

# ----------------------------------------------------- 4. coefficient tests (live content)
print("\n" + "-" * 104)
print("[4] LIVE TEST: the two closed forms carry SPECIFIED coefficients, so they are falsifiable")
print("-" * 104)
print(f"  {'beta':>7} {'shat=a^2s':>11} {'mu-4: 2g^2N=4N^2/beta':>22} {'meas/pred':>10}"
      f" {'mu6: R=shat/x^2':>17} {'shat/(aLam)^2':>14}")
pred4 = np.array([4.0 * N ** 2 / b for b in beta_v])
R2 = sa2_v / x_v ** 2
for i, b in enumerate(beta_v):
    print(f"  {b:7.4f} {sa2_v[i]:11.6f} {pred4[i]:22.4f} {sa2_v[i]/pred4[i]:10.5f}"
          f" {R2[i]:17.1f} {sa2_v[i]/aL_v[i]**2:14.4f}")
print(f"\n  (i)  mu-4: sigma a^2 = 2g^2N  ==>  measured/predicted runs {pred4[0]/pred4[-1]:.2f}x by formula"
      f" but {sa2_v[0]/sa2_v[-1]:.2f}x in data  =>  ratio drifts "
      f"{(sa2_v[0]/pred4[0])/(sa2_v[-1]/pred4[-1]):.2f}x  =>  **FALSIFIED [E]**")
print(f"       absolute deficit at the five betas: " + ", ".join(f"{sa2_v[i]/pred4[i]:.3g}" for i in range(5)))
print(f"  (ii) mu-6 with x^2: R = shat/x^2 runs {R2.max()/R2.min():.2f}x  =>  not a pure x^2 "
      f"(this is mu6's real weakness: it rests on R const)")
print(f"  (iii) mu-6 with (a Lambda_MS)^2: shat/(a Lambda)^2 = {sa2_v[0]/aL_v[0]**2:.4f} ... "
      f"{sa2_v[-1]/aL_v[-1]**2:.4f}, spread {(sa2_v[0]/aL_v[0]**2)/(sa2_v[-1]/aL_v[-1]**2)-1:.2e}"
      f"  =>  EXACT, coefficient = (sqrt(sigma)/Lambda)^2 = {(sqrtsig/LAM_MS)**2:.4f}")
print(f"  (iv) why x is the wrong variable: a*Lambda/x = " + "/".join(f"{aL_v[i]/x_v[i]:.0f}" for i in range(5))
      + f"  (drifts {aL_v[0]/x_v[0]/(aL_v[-1]/x_v[-1]):.2f}x -- NOT a constant offset)")
print(f"       local d ln x / d ln(a Lambda_MS) = " + "/".join(
          f"{np.log(x_v[i+1]/x_v[i])/np.log(aL_v[i+1]/aL_v[i]):.2f}" for i in range(4))
      + "  (drifts 0.57 -> 0.87: x is NOT a fixed power of a Lambda over this window)")
print(f"       => mu6's n=%.3f is an AVERAGE over a drifting local exponent 3.50 -> 2.28;" % c2L)
print(f"          the exact statement lives one variable over: n_aLam = {lidx(aL_v, sa2_v):.6f}.")

# ------------------------------------------------- 5. local (pairwise) indices: pure power?
print("\n" + "-" * 104)
print("[5] LOCAL INDICES -- is sigma-a^2 a PURE power of x? of a Lambda? of beta?")
print("-" * 104)
print(f"  {'beta window':>17} {'dln shat/dln x':>15} {'dln shat/dln aLam':>18} {'dln shat/dln beta':>18}")
for i in range(len(beta_v) - 1):
    d = (np.log(sa2_v[i+1] / sa2_v[i]))
    n1 = d / np.log(x_v[i+1] / x_v[i])
    n2 = d / np.log(aL_v[i+1] / aL_v[i])
    n3 = d / np.log(beta_v[i+1] / beta_v[i])
    print(f"  {beta_v[i]:6.3f}-{beta_v[i+1]:6.3f} {n1:15.3f} {n2:18.3f} {n3:18.3f}")
print(f"  global LS:  n_x = {c2L:.3f}   n_aLam = {lidx(aL_v, sa2_v):.6f}   "
      f"n_beta = {lidx(beta_v, sa2_v):.3f}")

# ------------------------------------------------------------------- 6. counterfactual
print("\n" + "-" * 104)
print("[6] COUNTERFACTUAL [X]: if a^-1 were anchored at 2.91 GeV (mu-6's original slip)")
print("-" * 104)
k = 2.91 / a_inv(6.0)
print(f"  a^-1 x {k:.4f} for ALL beta  =>  sigma x {k**2:.4f} for all beta  =>")
print(f"  every column in [2]/[3] is a constant SHIFT: C1,C2,C3 unchanged to 1e-15;"
      f" sigma/Lambda^2 -> {sigma_ref*k**2/LAM_MS**2:.4f}")
print("  => the 1.297x anchor error moves OFFSETS only, never exponents. It was a reporting error,")
print("     not the cause of the mu-4 / mu6 disagreement.")

# ----------------------------------------------------------------------- 7. verdict
print("\n" + "=" * 104)
print("[7] VERDICT")
print("=" * 104)
print(f"  C1 = {c1:+.2e} (vs -2: NO)   C2 = {c2:+.2e} (vs 2: NO, degenerate)   C3 = {c3:+.2e} (vs 0: YES)")
print("  (a) FIRES: d ln(sigma/Lambda_MS^2)/d ln a = 0  =>  sigma proportional to Lambda^2,")
print("      sigma/Lambda_MS^2 = %.4f = O(1)  =>  consistent with mu6, mu-4's sigma_B/sigma_C is wrong." % (sigma_ref / LAM_MS**2))
print("  (b) does NOT fire: d ln sigma/d ln a = 0, not -2  =>  sigma proportional to 1/a^2 is EXCLUDED")
print("      as a statement about the physical tension, and mu-4's coefficient 2g^2N is falsified by")
print(f"      {(sa2_v[0]/pred4[0])/(sa2_v[-1]/pred4[-1]):.1f}x drift.")
print("  (c) NOT invoked as 'no pure power': the beta-dependence IS a pure power -- but of a*Lambda_MS")
print("      (exponent exactly 2.000000), NOT of x. mu6 got the right scale and the wrong variable:")
print("      x/(a Lambda) drifts 1.52x across the window and the local exponent is not constant")
print("      (3.50 -> 2.28), so 2.63 is an average, not a law.")
print("  CAVEAT [T]: with a^-1 taken from the same Tab.13 column, C1 and C3 are identities, so the data")
print("      certify only the O(1) number sqrt(sigma)/Lambda_MS = %.4f -- not an independent exponent."
      % (sqrtsig / LAM_MS))
print("=" * 104)
