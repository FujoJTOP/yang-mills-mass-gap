"""at_scalar.py -- T1 scalarization + T2 sharp feasibility for the BD criterion (star).

(star)   Cdot_t Hess V_t Cdot_t - (1/2) Cddot_t  >=  lamdot_t Cdot_t     (all fields, all t>0)

Setup fixed by the task:  C_t = t g^{-1}  =>  Cdot_t = g^{-1} (constant),  Cddot_t = 0.
In a g-orthonormal frame g^{-1} = identity, so for every g-unit v

(star)  <=>  <v, Hess_g V_t(phi) v>_g  >=  lamdot_t      for all phi in X, all t>0
        <=>  lamdot_t <= m(t) := inf_phi lambda_min( Hess_g V_t(phi) ).          [T1]

So (star) is a CEILING on lamdot_t, nothing else.  Define

        M(t) := int_0^t m(s) ds,        lam_t = int_0^t lamdot_s ds.

Since lamdot <= m pointwise and M(0)=lam(0)=0, we get lam_t <= M(t) POINTWISE.
Hence e^{-2 lam_t} >= e^{-2 M(t)} and:

    (T2 THEOREM)  There exists lamdot_t <= m(t) with lam_t -> +oo and
    int_0^inf e^{-2 lam_t} dt < inf   <=>   I := int_0^inf e^{-2 M(t)} dt < inf,
    and then the OPTIMAL choice is lam_t = M(t), gamma_max = 1/I.

Necessity: lam <= M => int e^{-2 lam} >= int e^{-2M} = I.  If I = infinity no
admissible lam exists at all (lambda can only be LOWERED below M, never raised).
Sufficiency: take lam = M.

Consequence: the three conditions of T2 collapse into ONE scalar condition on the
ceiling m.  In particular lam_t is NOT "free to grow": (star) bounds it from
ABOVE, so a low ceiling forbids growth.  Growth needs a NON-INTEGRABLE ceiling.

This script quantifies six candidate ceilings, all with the same integrable
negative part -a e^{-c t} (a = the measured negative part of Hess S_W, c = the
spectral gap of the heat semigroup).  Only the four last columns matter.
"""
import json
import numpy as np

OUT = {}
N, BETA = 3, 2.5
LAM1 = (N * N - 1) / (2.0 * N)      # spectral gap of -Delta_g on SU(N), g=-2tr(XY)
HALF_RIC = N / 8.0                  # (1/2) Ric_g = (N/8) g   for bi-invariant g
A_NEG = 4.8 * BETA / N              # negative part scale from w7_pl_hess (Hess S_W >= -(4.8/N) beta g)
C_RATE = LAM1

print("=" * 92)
print("T1  SCALARIZATION            N=%d  beta=%.2f   lambda_1=(N^2-1)/(2N)=%.4f   (1/2)Ric=N/8=%.4f"
      % (N, BETA, LAM1, HALF_RIC))
print("=" * 92)
print("  C_t = t g^{-1}  =>  Cdot = g^{-1},  Cddot = 0   =>   Cdot Hess V Cdot = Hess V,  -(1/2)Cddot = 0")
print("  (star)  <v,Hess_g V_t v>_g >= lamdot_t  for every g-unit v, every phi, every t>0")
print()
print("  DIMENSION of lamdot_t:  <v,Hess V_t v>_g is dimensionless (v is g-unit, V_t is an action),")
print("  so lamdot_t is DIMENSIONLESS = a curvature FLUX (curvature per unit flow time t).")
print("  t itself is heat time measured in units of g^{-1}, i.e. t = (length)^2 in the physical metric")
print("  g_phys = a^2 g.  Then lam_t is dimensionless, int e^{-2 lam} dt has dimension a^2, so")
print("  gamma = 1/int e^{-2 lam} scales like 1/a^2 = a physical mass^2.  Free field check below.")
print()

# ---------------------------------------------------------------- ceilings
def M_curve(mfun, t):
    """M(t) = int_0^t m, with M(0)=0 (cumulative trapezoid on a log grid)."""
    vals = np.array([mfun(s) for s in t])
    dt = np.diff(t)
    return np.concatenate([[0.0], np.cumsum(0.5 * (vals[1:] + vals[:-1]) * dt)])

t = np.concatenate([np.logspace(-6, 0, 4000), np.linspace(1.0, 60.0, 4001)[1:],
                    np.logspace(np.log10(60.0), 6, 4000)[1:]])

# NOTE on regularisation: BD require lam_t = int_0^t lamdot_s ds to exist, so lamdot must be
# integrable near 0.  A pure delta/t ceiling is therefore inadmissible (it diverges at t=0 for
# delta>=1/2 and at t=inf for delta<=1/2 -- in fact for NO delta).  The physically occurring
# "1/t" budgets are regularised: 1/(1+t) with a 1/(1+t)^2 regulator.  Row 6 models the flat
# phi^4 ceiling lamdot = 1/t - chi/t^2 in exactly that regularised form (PL-POLCHINSKI 2.3).
cases = [
    ("free field  S_W=0            (literal form, no source)",
     lambda s: 0.0 * s, "lamdot<=0"),
    ("literal (star):  m=-a e^{-ct}            ",
     lambda s: -A_NEG * np.exp(-C_RATE * s), "lamdot<=-a e^{-ct}"),
    ("literal + delta/(1+t), delta=0.40        ",
     lambda s: 0.40 / (1.0 + s) - A_NEG * np.exp(-C_RATE * s), "delta<1/2"),
    ("literal + delta/(1+t), delta=0.60        ",
     lambda s: 0.60 / (1.0 + s) - A_NEG * np.exp(-C_RATE * s), "delta>1/2"),
    ("MANIFOLD (3.100):  N/8 - a e^{-ct}       ",
     lambda s: HALF_RIC - A_NEG * np.exp(-C_RATE * s), "half-Ric source"),
    ("flat phi^4: 1/(1+t) - chi/(1+t)^2        ",
     lambda s: 1.0 / (1.0 + s) - (A_NEG / C_RATE) / (1.0 + s) ** 2, "1/t source, delta=1"),
    ("Ric>0 + delta/(1+t), delta=0.40          ",
     lambda s: HALF_RIC + 0.40 / (1.0 + s) - A_NEG * np.exp(-C_RATE * s), "both"),
]

print("=" * 92)
print("T2  FEASIBILITY OF (star):  lam_t <= M(t),  feasible <=> I = int_0^inf e^{-2M} < inf")
print("=" * 92)
print("%-42s %12s %12s %12s   %s" % ("ceiling m(t)", "M(inf)", "I=int e^-2M", "gamma_max", "verdict"))
print("-" * 92)
rows = []
for name, f, note in cases:
    M = M_curve(f, t)
    Minf = M[-1]
    integ = np.exp(-2.0 * M)
    I = np.trapezoid(integ, t)
    # CONVERGENCE TEST: I(T) must plateau (I(1e6)/I(1e2) -> 1); a power-law tail fails it.
    def tail(T0):
        mm = t >= T0
        return float(np.trapezoid(integ[mm], t[mm]))
    t2, t4 = tail(1e2), tail(1e4)
    ratio = t4 / t2 if t2 > 0 else 0.0
    conv = (ratio < 0.7)          # power-law tails with exponent > -1 fail this test
    gm = 1.0 / I if conv else float("inf")
    ok = conv and (Minf > 0)
    verdict = "FEASIBLE  gamma_max=%.4g" % gm if ok else ("DEAD (I(T) not convergent)" if not conv else "DEAD (M(inf)<=0)")
    print("%-42s %12.5f %12.4g %12s   %s" % (name, Minf, I, ("%.4g" % gm) if ok else "-", verdict))
    rows.append(dict(name=name, Minf=float(Minf), I=float(I), tailratio=float(ratio), conv=bool(conv),
                     gamma_max=float(gm) if np.isfinite(gm) else None, feasible=bool(ok)))
print()

# ---------------------------------------------------------------- free field exact
print("=" * 92)
print("FREE-FIELD ACID TEST (S_W = 0  =>  V_t = 0  =>  Hess V_t = 0, i.e. m(t) = 0 for all t)")
print("=" * 92)
print("  literal form:   lamdot <= 0  =>  lam_t <= 0  =>  int e^{-2 lam} >= int 1 = inf  =>  NO admissible lam.")
print("  Yet the TRUE answer for Haar on SU(N)^E is a LSI with gamma = lambda_1(Delta_g).  Table:")
print("     N     lambda_1 = (N^2-1)/(2N)     N/4 (yield of the half-Ric source)     ratio")
for n in (1, 2, 3, 4, 5, 10):
    l1 = (n * n - 1) / (2.0 * n) if n > 1 else 0.0
    print("   %3d %14.4f %30.4f %20s" % (n, l1, n / 4.0, ("%.2f" % (n / 4 / l1)) if l1 else "inf (Ric=0)"))
print("  => the half-Ric source N/4 reproduces the correct ORDER of gamma for every N, and gives")
print("     exactly 0 for N=1 (U(1), flat torus: Ric=0), where the criterion MUST fail because")
print("     the truth is the massless L^-2 behaviour.  The criterion is acid-test-clean in the")
print("     MANIFOLD form and acid-test-DEAD in the literal flat form.")
print()

# ---------------------------------------------------------------- sharp statement
print("=" * 92)
print("SHARP FORM OF T2 (answers 'is lam_t free to grow?')")
print("=" * 92)
print("  Write m(t) = mu - q(t) with q >= 0 integrable (q = the negative part, decaying).")
print("  Then M(t) = mu t - Q(t),  Q(t) = int_0^t q <= Q(inf) < inf, and")
print("      I = int_0^inf e^{-2 mu t + 2 Q(t)} dt.")
print("    * mu > 0 : M(t) -> +oo linearly, I = e^{2Q(inf)}/(2 mu) < inf  =>  FEASIBLE,")
print("               gamma_max = 2 mu e^{-2 Q(inf)}, optimal lam_t = mu t - Q(t).")
print("    * mu = 0 : M(t) = -Q(t) <= 0 bounded above  =>  lam_t bounded above  =>  lam_t -/-> +oo.")
print("    * mu < 0 : M(t) -> -oo  =>  lam_t -> -oo.")
print("  ALSO: if q is NOT integrable but mu=0, a 1/t-type ceiling still works: m = delta/t - q,")
print("  M ~ delta log t, I ~ int t^{-2 delta} converges iff delta > 1/2  (table row 4 vs 3).")
print("  => ANSWER: lam_t is NOT free to grow.  The ceiling is what it is; the criterion is")
print("     satisfiable iff the ceiling has a NON-INTEGRABLE PART: either a positive constant")
print("     floor mu > 0, or a delta/t with delta > 1/2.  In our problem the ONLY candidate for")
print("     mu is (1/2)Ric_g = N/8 (manifold form) -- the literal form has mu = 0 identically.")
print()

# ---------------------------------------------------------------- (lam, lamdot) plane
print("=" * 92)
print("FEASIBLE REGION IN THE (t, lam) PLANE  (upper envelope = M; the region is the hypograph)")
print("=" * 92)
print("  rows: literal ceiling (mu=0) vs manifold ceiling (mu=N/8).  columns: t.")
ts = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
Mlit = M_curve(cases[1][1], t)
Mman = M_curve(cases[4][1], t)
print("   t      m_lit(t)   M_lit(t)   lam_lit_max    m_man(t)   M_man(t)   lam_man_max")
for tt in ts:
    i = int(np.argmin(np.abs(t - tt)))
    print("  %5.2f %10.4f %10.4f %12.4f %12.4f %10.4f %12.4f"
          % (t[i], cases[1][1](t[i]), Mlit[i], Mlit[i], cases[4][1](t[i]), Mman[i], Mman[i]))
print("  -> literal: lam_lit_max converges to -a/c = %.4f  (a bounded ceiling, NO divergence)."
      % (-A_NEG / C_RATE))
print("  -> manifold: lam_man_max ~ (N/8) t  -> +oo.  Same integrable negative part, different source.")
print()
print("  gamma_max(manifold, exact grid integral) = %.6f" % rows[4]["gamma_max"])
print("  certified lower bound  gamma >= 2 mu exp(-2 Q(inf)),  Q(inf)=int_0^inf ||(Hess V_t)_-|| dt:")
print("     Q(inf) = a/c = %.4f  =>  gamma >= (N/4)exp(-2a/c) = %.6g" %
      (A_NEG / C_RATE, (N / 4.0) * np.exp(-2 * A_NEG / C_RATE)))
print("  NOTE the exponential smallness: gamma ~ exp(-2 int||(Hess V_t)_-||).  For SU(3), beta=2.5")
print("  the bound is gamma ~ %.3g -- POSITIVE but useless as a uniform bound.  The whole" %
      ((N / 4.0) * np.exp(-2 * A_NEG / C_RATE)))
print("  residual question of the programme is therefore NOT 'does (star) hold' but 'is")
print("  Q(inf) = int_0^inf ||(Hess V_t)_-|| dt bounded UNIFORMLY IN L'.")

OUT["params"] = dict(N=N, beta=BETA, lam1=LAM1, half_Ric=HALF_RIC, a_neg=A_NEG, c_rate=C_RATE)
OUT["cases"] = rows
OUT["t"] = t[::50].tolist()
OUT["M_literal"] = Mlit[::50].tolist()
OUT["M_manifold"] = Mman[::50].tolist()
with open("research/v3/species/w7/at_scalar.json", "w") as fh:
    json.dump(OUT, fh, indent=1)
print()
print("[saved] research/v3/species/w7/at_scalar.json")
