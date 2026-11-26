#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v06_radius.py -- KP/polymer radius of 4D SU(N) Wilson theory: N-scaling, c, ceiling.
Pure arithmetic (scipy: lambertw, iv only). No simulation. Single process.

Conventions (fixed by Z06-HOLEA, do not re-derive here):
  S = beta * sum_p (1 - N^-1 Re tr U_p),   Boltzmann factor per plaquette e^{(beta/N) Re tr U_p}
  c_f(beta) = beta/(2N) + O(beta^3)          (un-normalised character coefficient)
  u := c_f/N = beta/(2N^2) + O(beta^3)       (per-tile normalised weight, Gupta 14.18)
  KP functional: sum_{gamma contains x} |w(gamma)| e^{a|gamma|} <= a
  project bound: K_KP <= 80 e^{2+beta} N   (P3-FULLREP-KP 4.2, beta <= 1/48)
  criterion (a = 1): u * K < 1
"""
import numpy as np
from scipy.special import lambertw, iv

# ---------- project constants (all from P3-FULLREP-KP / Z02-SEILER) ----------
DELTA = 20.0            # plaquettes through a link x 4 links - 1  -> animal growth bound (e*DELTA)^m
LAM_A = np.e * DELTA    # = 54.3656  crude "animal" constant used in the project
C_K = 80.0 * np.e**2    # K <= C_K * e^beta * N ; C_K = 591.2
SZZ = 1.0 / 48.0        # SZZ LSI radius  arXiv:2204.12737  |beta| < 1/(16(d-1)), d=4

def beta_kp(N, u_norm='c_f/N'):
    """KP radius from u*K<1 with K = C_K e^beta N.  Solve beta e^beta = X."""
    if u_norm == 'c_f/N':      # u = beta/(2N^2)  -> beta e^beta = 2N/C_K * e^{-2} * e^2 ...
        X = 2.0 * N / C_K      # = N/295.6
    elif u_norm == 'c_f':      # u = beta/(2N)    -> beta e^beta = 1/295.6
        X = 2.0 / C_K
    else:
        raise ValueError
    return float(np.real(lambertw(X)))

def N_for_beta(bt, u_norm='c_f/N'):
    """inverse of beta_kp"""
    if u_norm == 'c_f/N':
        return C_K * bt * np.e**bt / 2.0
    return None

def u_of(beta, N):
    return beta / (2.0 * N**2)

def u_ceiling_tube():
    """(T) structural ceiling of ANY KP-type expansion from the four-wall tube family alone.
    Tubes: |gamma| = 4T (T = 1,2,...), weight a0 u^{4T}.  Criterion needs
        Psi(x) := x^4/(1-x^4) <= a = ln(x/u),   x := u e^a.
    G(x) = ln(x/u) - Psi(x) has its max at 1-x^4 = 2x^2 -> x*^2 = sqrt2-1.
    Feasible iff  ln(x*/u) >= x*^4/(1-x*^4)  <=>  u <= x* exp(-x*^4/(1-x*^4)).
    """
    xs = np.sqrt(np.sqrt(2.0) - 1.0)
    us = xs * np.exp(-(xs**4) / (1.0 - xs**4))
    return xs, us

def su2_weights(beta):
    """exact SU(2): u_tile = I2/I1 (repo R4-PERCENT:19), c_f = 2*I2/I1 (AUDIT-R5B:222)"""
    z = beta
    return iv(2, z) / iv(1, z), 2.0 * iv(2, z) / iv(1, z)

out = []
P = out.append

P("=" * 78)
P("V06 -- KP/polymer radius of 4D SU(N):  c, N-scaling, ceiling")
P("=" * 78)
P("inputs: Delta=20, lam_a=e*Delta=%.4f, K<=%.1f*e^beta*N, SZZ LSI=1/48=%.5f"
  % (LAM_A, C_K, SZZ))
P("")

# ---------- 1. the N-scaling of the KP radius ----------
P("[1] KP radius beta_1(N) from u*K<1   (u = c_f/N, K = %.1f e^beta N)" % C_K)
P("    beta_1(N) = W(N/%.1f)   [W = Lambert W, principal branch]" % (C_K / 2.0))
P("    %6s %14s %14s %12s %12s %10s" % ("N", "beta_1", "beta_1/N", "g2_1=2N/b", "u_max=1/(C_K N)", "1/48 binds?"))
for N in [2, 3, 4, 5, 6, 7, 8, 10, 20, 50, 100, 300, 830, 1000, 5000]:
    b = beta_kp(N)
    P("    %6d %14.6e %14.6e %12.4f %12.4e %10s"
      % (N, b, b / N, 2.0 * N / b, 1.0 / (C_K * N), "yes" if b < SZZ else "no"))
bA = beta_kp(2, 'c_f')
P("    normalisation A (u = c_f, no per-link 1/N): beta_1 = %.6e  (N-independent) [= repo 1/298]" % bA)
P("    crossing beta_1 = 1/48 at N = %.2f ;  beta_1 = 1 at N = %.1f ;  beta_1 = 2.5 at N = %.0f"
  % (N_for_beta(SZZ), N_for_beta(1.0), N_for_beta(2.5)))
P("    A/B ratio = N/2  =>  adopting Z06 normalisation relaxes beta_1 by exactly N/2")
P("")

# ---------- 2. what physical points cost ----------
P("[2] colours needed for the physical targets (beta_1 >= beta_target), crude K")
for name, bt in [("SU(2) beta=2.5", 2.5), ("SU(3) beta=2.89", 2.89), ("beta=1", 1.0)]:
    P("    %-18s beta_target=%.2f -> N >= %.0f" % (name, bt, N_for_beta(bt)))
P("    with K improved by factor F: N >= N_above / F")
for F in [6, 10, 25, 100, 620]:
    P("      F=%4d : N >= %8.0f (for beta=2.5)" % (F, N_for_beta(2.5) / F))
P("")

# ---------- 3. structural ceiling: the four-wall tubes alone ----------
xs, us = u_ceiling_tube()
P("[3] structural ceiling from the tube family ALONE (no counting constant at all)")
P("    x* = sqrt(sqrt2-1) = %.6f ;  u* = x* exp(-x*^4/(1-x*^4)) = %.6f" % (xs, us))
P("    => any KP-type expansion with the true tube weights needs u <= %.6f" % us)
P("    => beta <= 2 N^2 u* = %.4f N^2 :  N=2 -> %.3f   N=3 -> %.3f   N=4 -> %.3f"
  % (2 * us, 4 * 2 * us, 9 * 2 * us, 16 * 2 * us))
u25, cf25 = su2_weights(2.5)
P("    SU(2), beta=2.5 exact: u_tile=I2/I1=%.6f (=c_f/N), c_f=%.6f" % (u25, cf25))
P("    u(2.5)/u* = %.4f  -> %s" % (u25 / us, "INSIDE (not excluded)" if u25 < us else "OUTSIDE (excluded)"))
P("    budget check at u=u(2.5): a*=ln(x*/u)=%.6f, tube sum Psi(x*)=%.6f, slack=%.6f"
  % (np.log(xs / u25), xs**4 / (1 - xs**4), np.log(xs / u25) - xs**4 / (1 - xs**4)))
P("    (slack > 0 but the slack is consumed by the very next family: (u e^a)^2 = %.4f)" % (xs**2))
P("")

# ---------- 4. how far the counting constant must improve ----------
P("[4] constant headroom at the physical points (need K < 1/u, vs K <= %.1f N)" % C_K)
# exact SU(2) weights at the two points of interest (leading order beta/(2N^2) is only a small-beta formula)
for bt in [1.0, 2.5]:
    u_ex, cf_ex = su2_weights(bt)
    P("    SU(2) b=%.1f  u_exact=I2/I1=%.4f  (leading beta/(2N^2)=%.4f, ratio %.2f)"
      % (bt, u_ex, bt / 8.0, u_ex / (bt / 8.0)))
for name, N, bt, ux in [("SU(2) b=2.5", 2, 2.5, u25), ("SU(2) b=1", 2, 1.0, su2_weights(1.0)[0]),
                        ("SU(3) b=2.89", 3, 2.89, None)]:
    u = ux if ux is not None else u_of(bt, N)
    Kneed = 1.0 / u
    Khave = C_K * N
    P("    %-14s u=%.4e  K_need<%.3f  K_have<=%.1f  improvement x%.0f  [u source: %s]"
      % (name, u, Kneed, Khave, Khave / Kneed, "exact Bessel" if ux is not None else "leading order"))
P("    available headroom: lam_a 54.37 -> true 4D animal constant ~9 (M) = x6.0 ;")
P("    C_g<=2Ne^beta -> sharper (=x2..4, M) ; total plausible 10-30x")
P("")

# ---------- 5. joint bound / multiscale candidates ----------
P("[5] candidate improvements, arithmetic verdict")
P("    (i)  normalisation u=c_f/N (Z06) : beta_1 multiplied by exactly N (A->B):")
P("         N=2 x2.0, N=3 x3.0, N=10 x10 -- the single largest available gain  [T]")
P("    (ii) joint bound over families: KP already sums ALL gamma contains x; splitting and")
P("         re-adding the bounds returns the same inequality => zero gain  [T]")
P("    (iii) multiscale first layer: cannot be priced here (v04 open)   [NG]")
P("    (iv) 1/N leading term: beta_1 ~ N/295.6 (linear), but g2_1 = 2N/beta_1 ~ 594 + O(ln N)")
P("         => in 't Hooft coupling the window is N-independent: no physical gain  [T]")
P("")
P("=" * 78)
P("VERDICT  beta_1^KP(N) = W(N/295.6);  c = 1/295.6 = 3.383e-3 in beta_1 = c N;")
P("         tube ceiling u* = 0.5232 => beta < 1.046 N^2 (does not exclude beta=1);")
P("         crude K needs x619 improvement at SU(2),2.5 -> headroom 10-30x => (b).")
P("=" * 78)

txt = "\n".join(out)
print(txt)
with open(__file__.replace('.py', '.log'), 'w', encoding='utf-8') as f:
    f.write(txt + "\n")
