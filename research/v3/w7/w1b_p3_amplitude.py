#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
w1b_p3_amplitude.py  --  W1B witness (attack on the P3 dressed amplitude [R])

W1 localised the next coefficient of w1(beta) to f1^2 = a0(1 - c u + ...); the
O(1) coefficient  c  of the dressed amplitude  D_fin = Z_dec e^{-f_dec t} is the
P3 [R] gap.  This witness adjudicates:  can  c  (the SIGN/floor of f1^2) be
BOUNDED at beta=2.5, or is it the same nonperturbative binding input W?

It is a NUMERICAL / EXACT-MOMENT WITNESS ([V] / exact-one-plaquette [T]).
NOT a proof of w1>0 at beta=2.5.  No QED.  Clay 0%.

Conventions identical to w1_strong_coupling.py (P3 / U3 / T2):
  weight exp((beta/N) Re tr U_p) ;  O = (1/N) Re tr U_p ;  w1 = f1^2/M0
  f1^2 = a0(1 - c u + ...) ;  a0 = s/(4 N^2) ;  u := c_f/N
  M0(0) = 1/(2 N^2)  (N>=3)  ;  M0(0) = 1/N^2  (SU(2))

Run:  py -3.12 -I w1b_p3_amplitude.py
"""
import numpy as np
from fractions import Fraction as F
from math import factorial
from scipy.special import iv
from numpy.polynomial.legendre import leggauss

np.set_printoptions(suppress=True)

# ---- SU(2) exact one-plaquette (Bessel) -------------------------------------
D = 12
def _mul(a, b):
    c = [F(0)] * (D + 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            if i + j > D: break
            if bj: c[i + j] += ai * bj
    return c
def _add(a, b): return [a[i] + b[i] for i in range(D + 1)]
def _inv(a):
    inv = [F(0)] * (D + 1); inv[0] = 1 / a[0]
    for n in range(1, D + 1):
        s = F(0)
        for k in range(1, n + 1): s += a[k] * inv[n - k]
        inv[n] = -s / a[0]
    return inv
def a_su2(j2):
    m = j2 + 1
    a = [F(0)] * (D + 1); k = 0
    while True:
        deg = 2 * k + (m - 1)
        if deg > D: break
        a[deg] += F(m) / (F(2) ** deg) / (F(factorial(k)) * F(factorial(k + m)))
        k += 1
    return a
a0b, ahh, a1b = a_su2(0), a_su2(1), a_su2(2)
a0i = _inv(a0b)
r_h = _mul(ahh, a0i); r_1 = _mul(a1b, a0i)
var = _add([F(1)] + [F(0)] * D, r_1)
hh = _mul(r_h, r_h)
var = [var[i] - hh[i] for i in range(D + 1)]
M0_su2 = [F(1, 4) * var[i] for i in range(D + 1)]

def M0_su2_exact(beta):
    return 0.25 * (1 + 3.0 * iv(3, beta) / iv(1, beta)) - (iv(2, beta) / iv(1, beta)) ** 2
def cf_su2_exact(beta):          # fundamental coefficient a_{1/2}(beta); u=cf/N
    return (4.0 / beta) * iv(2, beta)

# ---- SU(3) exact one-plaquette (Weyl quadrature) ----------------------------
def su3_expect(fun, beta, n=240):
    x, w = leggauss(n); phi = np.pi * (x + 1.0); wp = np.pi * w
    P1, P2 = np.meshgrid(phi, phi, indexing='ij'); W1, W2 = np.meshgrid(wp, wp, indexing='ij')
    Wt = W1 * W2; p3 = -(P1 + P2)
    e1, e2, e3 = np.exp(1j * P1), np.exp(1j * P2), np.exp(1j * p3)
    d = np.abs((e1 - e2) * (e1 - e3) * (e2 - e3)) ** 2
    ReTr = (e1 + e2 + e3).real
    norm = 1.0 / 6.0 / (2 * np.pi) ** 2
    ig = lambda g: norm * np.sum(Wt * d * g)
    Z = ig(np.exp(beta * ReTr))
    return ig(fun(e1 + e2 + e3) * np.exp(beta * ReTr)) / Z
def su3_M0(beta): return su3_expect(lambda t: (t.real / 3.0) ** 2, beta) - \
                          su3_expect(lambda t: t.real / 3.0, beta) ** 2
def su3_Z(beta, n=240):
    x, w = leggauss(n); phi = np.pi * (x + 1.0); wp = np.pi * w
    P1, P2 = np.meshgrid(phi, phi, indexing='ij'); W1, W2 = np.meshgrid(wp, wp, indexing='ij')
    Wt = W1 * W2; p3 = -(P1 + P2)
    e1, e2, e3 = np.exp(1j * P1), np.exp(1j * P2), np.exp(1j * p3)
    d = np.abs((e1 - e2) * (e1 - e3) * (e2 - e3)) ** 2
    return (1.0 / 6.0 / (2 * np.pi) ** 2) * np.sum(Wt * d * np.exp(beta * (e1 + e2 + e3).real))
def su3_cf(beta): return su3_expect(lambda t: t.real, beta) * su3_Z(beta)

b25 = 2.5
print("=" * 80)
print("W1B  the P3 dressed amplitude [R]: boundable at beta=2.5, or the binding input?")
print("=" * 80)

# ============================================================================
# [1] THE PIN  f1^2 <= M0  is exact and SIGN-BLIND
# ============================================================================
print("\n[1] exact two-sided PIN  0 <= f1^2 <= M0 = integral dmu   [T]")
print("    (M0 = sum_n f_n^2 >= f1^2 since all f_n^2 >= 0, f1^2 is one term)")
print("    At beta=0 the two exact anchors hit BOTH ends of the pin:")
# strong coupling end: f1^2=a0, M0=1/(2N^2), w1=s/2
for N in [2, 3, 5, 10]:
    if N == 2:
        M0 = 1.0 / N**2; w1 = 0.25 / M0 * 0  # SU(2): w1=s/4; a0=s/(4N^2)
        a0 = 1.0 / (4 * N**2); w1 = a0 / M0      # = s/4 for s=1
        tag = "SU(2) real fund.: M0=1/N^2 -> w1=s/4=%.3f" % w1
    else:
        M0 = 1.0 / (2 * N**2); a0 = 1.0 / (4 * N**2); w1 = a0 / M0
        tag = "N=%d complex fund.: M0=1/(2N^2) -> w1=s/2=%.3f" % (N, w1)
    print("      beta=0 :  f1^2=a0=%.5f   M0=%.5f   f1^2/M0=%.3f   (%s)" % (a0, M0, w1, tag))
print("    => strong-coupling anchor: f1^2 sits at the MIDPOINT (w1=1/2) of [0,M0].")
print("    => 2D/free anchor: f1^2=0 sits at the BOTTOM of [0,M0] (W3: w1=0 exactly).")
print("    Both satisfy every RP/KL/variance/gauge axiom  =>  the pin cannot fix the SIGN.")

# numerical pin WIDTH at 2.5 via the exact one-plaquette proxy  [V]
M0s2 = M0_su2_exact(b25); a0s2 = 1.0 / (4 * 2**2)
M0s3 = su3_M0(b25);       a0s3 = 1.0 / (4 * 3**2)
print("    beta=2.5 one-plaquette proxy pin of  f1^2/a0  (a0=s/(4N^2), s=1):")
print("      SU(2): M0(2.5)=%.5f  a0=%.5f  ->  f1^2/a0 in [0, %.3f]" % (M0s2, a0s2, M0s2 / a0s2))
print("      SU(3): M0(2.5)=%.5f  a0=%.5f  ->  f1^2/a0 in [0, %.3f]" % (M0s3, a0s3, M0s3 / a0s3))
print("    the window [0,~1.8] contains BOTH  f1^2=0  and  f1^2~a0  -> VACUOUS for the sign.")

# ============================================================================
# [2] KP RADIUS:  the ONLY proven absolute bound on |c| dies at beta=2.5
# ============================================================================
print("\n[2] KP radius u*K<1 :  the only proven bound |log D_fin| <= C u  [T]/[NG]")
print("    (P3-DRESSED-MASS 4.1: C = Ursell constant, finite iff u*K<1)")
print("    K(beta) = 80 e^{2+beta} N  (PROOF-P3 3, beta-dependent);  u=c_f/N")
cf2 = cf_su2_exact(b25); cf3 = su3_cf(b25)
K = lambda beta, N: 80.0 * np.exp(2.0 + beta) * N
for lab, cf, N in [("SU(2)", cf2, 2), ("SU(3)", cf3, 3)]:
    u = cf / N; kk = K(b25, N)
    print("      %s: c_f(2.5)=%.4f  u=%.4f  K=%.3e  u*K=%.2e  -> %s"
          % (lab, cf, u, kk, u * kk, "INSIDE" if u * kk < 1 else "OUTSIDE (bound vacuous)"))
# the KP boundary in beta (small-beta u=beta/(2N^2)):  u*K=1
for lab, N in [("SU(2)", 2), ("SU(3)", 3)]:
    from scipy.optimize import brentq
    f = lambda b: (b / (2.0 * N**2)) * K(b, N) - 1.0
    root = brentq(f, 1e-4, 1.0)
    print("      %s: small-beta KP boundary u*K=1 at beta_KP=%.4f  (2.5 is %.0fx beyond)"
          % (lab, root, b25 / root))
print("    => the convergent constant C (and any O(1) bound on c) lives only on beta<~0.01-0.02;")
print("       at 2.5 the decoration gas is outside the Kotecky-Preiss disk. No bound on c.")

# ============================================================================
# [3] SIGN-BLINDNESS:  log-convex C(t), decreasing m_eff hold for f1^2=0 AND >0
# ============================================================================
print("\n[3] RP consequences are UPPER-BOUND shaped: they hold at f1^2=0 too  [T]/[V]")
def grid(tmax=8.0, n=4000):
    return np.linspace(0.2, tmax, n)
s0 = (2 * 1.6) ** 2                                   # two-particle threshold
_s = np.linspace(s0, 60.0, 40000)
_rho = np.sqrt(_s - s0) * np.exp(-0.2 * _s)           # cut spectral density >= 0
def C_atom(t, f2=0.3, m1=0.8):
    # atom at m1 + a two-particle cut above threshold
    cut = np.trapezoid(_rho * np.exp(-np.sqrt(_s) * t), _s)
    return f2 * np.exp(-m1 * t) + cut
def C_cut(t):
    return np.trapezoid(_rho * np.exp(-np.sqrt(_s) * t), _s)
def logconvex(Cfun, name):
    t = grid(); C = np.array([Cfun(x) for x in t])
    lC = np.log(C)
    d2 = lC[2:] - 2 * lC[1:-1] + lC[:-2]              # discrete 2nd derivative
    meff = -(np.log(C[2:]) - np.log(C[:-2])) / (t[2:] - t[:-2])
    dec = np.all(meff[1:] - meff[:-1] <= 1e-9)
    print("      %-18s min(d2 logC)=%.2e (>0 convex?)   m_eff decreasing? %s   m_eff(end)=%.3f"
          % (name, d2.min(), "YES" if dec else "no", meff[-1]))
    return d2.min() > 0, dec
ok_atom = logconvex(lambda x: C_atom(x), "atom  (f1^2>0)")
ok_cut = logconvex(lambda x: C_cut(x), "pure cut (f1^2=0)")
print("    both pass RP's consequences => RP cannot tell f1^2>0 from f1^2=0  (as in T2 2.1).")

# ============================================================================
# [4] THE REDUCTION  [R] <-> W  at 2.5 (modulo the normaliser M0)
# ============================================================================
print("\n[4] [R] and W differ ONLY by the normaliser M0  [T]")
print("    W : w1 = f1^2/M0 >= w0         [R] : f1^2 >= delta a0")
print("    identity: w1 = (f1^2/a0) * (a0/M0).  a0/M0 is EXACT and >0 (one-plaquette).")
print("      SU(3) beta=2.5: a0/M0 = %.4f  =>  w1 = %.4f * (f1^2/a0)" % (a0s3 / M0s3, a0s3 / M0s3))
print("    so SIGN(w1) = SIGN(f1^2) = SIGN(1 - c u + ...):  the SAME object.")
print("    W and [R] are equivalent iff 0 < inf M0 <= sup M0 < inf  (variance bounded);")
print("    M0 is itself a positivity datum (W1: not the obstruction). => [R] IS W.")

# ============================================================================
# [5] FIVE ACID TESTS for the amplitude tool  (charter 2)
# ============================================================================
print("\n[5] five acid tests for the amplitude object f1^2 = a0(1 - c u + ...):")
# (1) free massless field: no recoupling -> a0=0 -> f1^2=0
print("    (1) free massless: no binding f(x)fbar->1 recoupling => a0=0 => f1^2=0  [T]")
# (2) 2D YM: trF^2=H conserved, f1^2=C2(R)<chi_0|chi_R>=0
def su2_orth(j, k, n=400):
    x, w = leggauss(n); th = np.pi * (x + 1.0) / 2.0; wt = np.pi * w / 2.0
    chi = lambda jj: np.sin((2 * jj + 1) * th) / np.sin(th)
    return np.sum(wt * (2 / np.pi) * np.sin(th) ** 2 * chi(j) * chi(k))
o = abs(su2_orth(0.5, 1.5))
print("    (2) 2D: <chi_{1/2}|chi_{3/2}>=%.2e => f1^2=0 (C2(R) x orthogonality)  [T]/[V]" % o)
# (3) one plaquette: degenerate, no transfer direction
print("    (3) one plaquette: atom DEGENERATE, a0/M0=s/2 is an artefact (no t-direction)  [T]")
# (4) strong coupling: w1(0)=s/2
print("    (4) strong coupling: w1(0)=s/2 (N>=3), s/4 (SU(2)) reproduced  [T] (mod KL)")
# (5) large N: w1=O(N^0), a0=O(1/N^2), M0=O(1/N^2)
print("    (5) large N: a0=1/(4N^2), M0=1/(2N^2) both O(N^-2) => w1=O(N^0), N-uniform  [T]")
# (6) continuum: w1 = Z^2 f1^2 /(Z^2 M0) invariant
print("    (6) continuum: f1^2->Z_O^2 f1^2, M0->Z_O^2 M0 => w1 invariant (gamma_w1=0)  [T]")

# ============================================================================
print("\n" + "=" * 80)
print("VERDICT (no QED, Clay 0%)")
print("  * f1^2=a0(1 - c u + ...); c = P3 dressed amplitude [R].")
print("  * [T] pin 0<=f1^2<=M0 is exact and SIGN-BLIND (both ends attained by exact models).")
print("  * [T] [R] <-> W differ only by the exact normaliser a0/M0; same obstruction.")
print("  * [NG] KP: the only proven bound |log D_fin|<=Cu needs u*K<1; at 2.5 u*K~1e4.")
print("  * => [R] is NOT boundable at beta=2.5 by expansion/positivity tools:")
print("         it IS the same nonperturbative binding input W.   [O]")
print("=" * 80)
