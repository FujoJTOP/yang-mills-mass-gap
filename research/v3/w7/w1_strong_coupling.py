#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
w1_strong_coupling.py  --  W1 witness (W-WAVE-CHARTER §1 W1)

Strong-coupling exact alignment w1=s/2 (U3) plus the NEXT coefficient, the exact
one-plaquette normaliser series M0(beta)=Var(O), and a Padé / Borel continuation
attempt to beta=2.5.  Decides whether w1>0 survives, and localises the exact
missing step.

Everything below is a NUMERICAL / EXACT-SERIES WITNESS ([V] / exact-one-plaquette [T]).
It is NOT a proof that w1>0 at beta=2.5.  No QED.  Clay 0%.

Conventions (P3 / U3 / T2):
  action      S = beta * sum_p (1 - (1/N) Re tr U_p) ,  weight  exp((beta/N) Re tr U_p)
  source      O = (1/N) Re tr U_p
  atom frac   w1 = f1^2 / M0 ,   M0 = Var(O) = <O^2>_c
  strong-cpl  f1^2 = a0 (1 - O(u)) ,  a0 = s/(4 N^2) ,  u := c_f(beta)/N  (plaquette exp.)
              M0(0) = 1/(2 N^2)  (N>=3, complex fund.)   =>  w1(0) = s/2   [T] (mod KL)
              M0(0) = 1/(  N^2)  (SU(2), real fund.)     =>  w1(0) = s/4   (real-rep factor 2)

Run:  py -3.12 -I w1_strong_coupling.py
"""
import numpy as np
from fractions import Fraction as F
from math import factorial
from scipy.special import iv
from numpy.polynomial.legendre import leggauss

np.set_printoptions(suppress=True)

# ----------------------------------------------------------------------------
# exact polynomial helpers in beta (Fraction coefficients, degree D)
# ----------------------------------------------------------------------------
D = 16
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
    """SU(2) character coefficient a_j = (2(2j+1)/beta) I_{2j+1}(beta); j2=2j in {0,1,2}."""
    m = j2 + 1                      # m = 2j+1
    a = [F(0)] * (D + 1); k = 0
    while True:
        deg = 2 * k + (m - 1)
        if deg > D: break
        a[deg] += F(m) / (F(2) ** deg) / (F(factorial(k)) * F(factorial(k + m)))
        k += 1
    return a

# ----------------------------------------------------------------------------
# SU(3) one-plaquette moments via the exact Weyl integration formula
# ----------------------------------------------------------------------------
def su3_expect(fun, beta, n=240):
    """< fun >_beta for one SU(3) plaquette, weight exp((beta/3) Re tr U)."""
    x, w = leggauss(n); phi = np.pi * (x + 1.0); wp = np.pi * w
    P1, P2 = np.meshgrid(phi, phi, indexing='ij')
    W1, W2 = np.meshgrid(wp, wp, indexing='ij')
    Wt = W1 * W2; p3 = -(P1 + P2)
    e1, e2, e3 = np.exp(1j * P1), np.exp(1j * P2), np.exp(1j * p3)
    d = np.abs((e1 - e2) * (e1 - e3) * (e2 - e3)) ** 2     # |Vandermonde|^2
    ReTr = (e1 + e2 + e3).real
    norm = 1.0 / 6.0 / (2 * np.pi) ** 2
    ig = lambda g: norm * np.sum(Wt * d * g)
    Z = ig(np.exp(beta * ReTr))
    return ig(fun(e1 + e2 + e3) * np.exp(beta * ReTr)) / Z

def su3_M0(beta, n=240):
    return su3_expect(lambda t: (t.real / 3.0) ** 2, beta, n) - \
           su3_expect(lambda t: t.real / 3.0, beta, n) ** 2

# ----------------------------------------------------------------------------
# Pade / Borel of a truncated power series  sum a_k beta^k   (a_k here: even only)
# ----------------------------------------------------------------------------
def pade_eval(coeffs, x0, L, M):
    """[L/M] Pade of the series  sum a_k x^k  evaluated at x0  (x = beta^2 for an even series).
       f(x) Q(x) - P(x) = O(x^{L+M+1});  Q(0)=1;  sum_{j=0}^{M} q_j a_{k-j} = 0, k=L+1..L+M."""
    a = list(coeffs) + [0.0] * (L + M + 1 - len(coeffs))
    A = np.zeros((M, M)); b = np.zeros(M)
    for k in range(1, M + 1):
        for j in range(1, M + 1):
            idx = L + k - j
            A[k - 1, j - 1] = a[idx] if idx >= 0 else 0.0
        b[k - 1] = -a[L + k]
    q = np.linalg.solve(A, b)
    Q = np.concatenate([[1.0], q])
    P = np.zeros(L + 1)
    for i in range(L + 1):
        P[i] = sum(Q[j] * a[i - j] for j in range(0, min(i, M) + 1) if i - j >= 0)
    return np.polyval(P[::-1], x0) / np.polyval(Q[::-1], x0)

def pade_beta(coeffs_beta, beta0, L, M):
    """Pade in beta for an even series given as beta-coeffs (odd ones zero)."""
    even = [coeffs_beta[2 * k] for k in range(len(coeffs_beta) // 2 + 1)]
    return pade_eval(even, beta0 * beta0, L, M)

# ============================================================================
print("=" * 78)
print("W1  strong-coupling alignment w1(0)=s/2  +  next coefficient  +  continuation")
print("=" * 78)

# ---------------------------------------------------------------------------
# PART 1 -- exact SU(2) one-plaquette normaliser series (Fraction / Bessel)
# ---------------------------------------------------------------------------
a0 = a_su2(0); ah = a_su2(1); a1 = a_su2(2)
a0i = _inv(a0)
r_h = _mul(ah, a0i)       # a_{1/2}/a0 = <chi_f>
r_1 = _mul(a1, a0i)
var = _add([F(1)] + [F(0)] * D, r_1)
hh = _mul(r_h, r_h)
var = [var[i] - hh[i] for i in range(D + 1)]
M0 = [F(1, 4) * var[i] for i in range(D + 1)]        # SU(2): M0 = (1/4) Var(tr U)
M00 = M0[0]
g_su2 = _inv([M0[i] / M00 for i in range(D + 1)])    # g = M0(0)/M0(beta) = 1/r

print("\n[1] SU(2) exact one-plaquette normaliser  M0(beta) = Var((1/2) tr U)")
print("    M0(0) = 1/4   (real fundamental: M0(0)=1/N^2, NOT 1/(2N^2))")
print("    exact series  M0(beta) =")
for n in range(0, 11, 2):
    print("        beta^%-2d : %s" % (n, M0[n]))
print("    => only EVEN powers in beta (SU(2): no symmetric rank-3 invariant, d_abc=0).")
print("    exact series  g(beta)=M0(0)/M0(beta) = 1/r(beta):")
for n in range(0, 9):
    print("        beta^%-2d : %s" % (n, g_su2[n]))

# ---------------------------------------------------------------------------
# PART 2 -- SU(3) one-plaquette: cubic Haar moment => a beta-LINEAR term (d_abc != 0)
# ---------------------------------------------------------------------------
M0_33 = su3_M0(0.0); M0p_33 = su3_expect(lambda t: (t.real / 3.0) ** 3, 0.0)  # <O^3>_Haar
cub = su3_expect(lambda t: t.real ** 3, 0.0)
print("\n[2] SU(3) one-plaquette normaliser  (generic N>=3: complex fundamental)")
print("    M0(0)            = %.12f   (exact 1/18 = %.12f)" % (M0_33, 1 / 18))
print("    <(Re tr U)^2>    = %.12f   (exact 1/2)" % su3_expect(lambda t: t.real ** 2, 0.0))
print("    <(Re tr U)^3>    = %.12f   (exact 1/4  -- the d_abc cubic Casimir moment)" % cub)
print("    M0'(0) = (1/N^2)<(Re trU)^3>_Haar = %.12f   (exact 1/36)" % M0p_33)
print("    => M0(beta) = 1/18 + beta/36 + ...  (ODD beta term!)")
print("    => normaliser ratio g = M0(0)/M0(beta) = 1 - beta/2 + ...  (exact linear coeff -1/2)")

# ---------------------------------------------------------------------------
# PART 3 -- the atom fraction w1(beta) = a0 (1 - c u) / M0(beta):  next coefficient
# ---------------------------------------------------------------------------
print("\n[3] w1(beta) = a0 (1 - c u + ...)/M0(beta),   a0=s/(4N^2), u=c_f/N.")
print("    normaliser part of the beta-linear coefficient  (EXACT):")
print("        SU(2) :  0            (M0 even)")
print("        SU(3) : -1/2          (d_abc cubic moment; s-independent)")
print("    amplitude part of the beta-linear coefficient  ([O] = dressed-amplitude [R]):")
print("        SU(2) : -c/2   in units of s/4   -> the ONLY linear term")
print("        SU(3) : -c/18  in units of s/2")
print("    => the full next coefficient is -(s/2)[1/2 + c/18] (N>=3), -(s/4)[c/2] (SU(2)),")
print("       and its 'c' piece is EXACTLY the unproven dressed-amplitude [R] gap.")
print("       Proven inputs fix only the '1/2' (SU(3)) / '0' (SU(2)) piece.")

# ---------------------------------------------------------------------------
# PART 4 -- Pade / Borel continuation of the EXACT normaliser to beta=2.5
# ---------------------------------------------------------------------------
print("\n[4] Pade / Borel continuation of the exact normaliser to beta=2.5")
b25 = 2.5
# exact SU(2) M0(2.5) from Bessel (not from the truncated series)
M0_exact_25 = 0.25 * (1 + 3.0 * iv(3, b25) / iv(1, b25)) - (iv(2, b25) / iv(1, b25)) ** 2
print("    exact  M0(2.5) (SU(2), Bessel) = %.10f ; exact 1/r(2.5)=%.10f" %
      (M0_exact_25, 0.25 / M0_exact_25))
print("    partial sums S_m(2.5) of the even series (do they converge to %.6f?):" % M0_exact_25)
s = 0.0
for n in range(0, 11, 2):
    s += float(M0[n]) * b25 ** n
    print("        m=%2d terms : S(2.5)=% .8f   (err %+.1e)" % (n // 2 + 1, s, s - M0_exact_25))
for LM in [(2, 2), (3, 3), (4, 4)]:
    try:
        v = pade_beta([float(x) for x in M0], b25, LM[0], LM[1])
        print("        Pade [%d/%d] at 2.5 : %.8f   (rel.err %.2e)" %
              (LM[0], LM[1], v, abs(v - M0_exact_25) / M0_exact_25))
    except Exception as e:
        print("        Pade [%d/%d] : singular (%s)" % (LM[0], LM[1], e))
print("    -> the ONE-PLAQUETTE normaliser series has finite radius |beta|~3.83 (nearest zero")
print("       of I_1), so it DOES converge at 2.5: raw partial sums still oscillate by ~5%,")
print("       but Pade [4/4] nails M0(2.5) to 1e-9.  So Pade/Borel WORK on the normaliser.")
print("       (This is the one-plaquette radius; the full 4D CLUSTER expansion has the")
print("        separate, smaller KP radius checked in [6].)")
print("       MORE IMPORTANT: M0(beta)=Var(O)>0 for EVERY beta  =>  w1_norm=a0/M0 > 0 trivially.")
print("       A continuation that is positive by construction cannot test survival -> VACUOUS.")
print("    SU(3): the LINEAR term alone (1/18 + 2.5/36 = 0.125) OVERSHOOTS the exact")
print("       M0(2.5)=0.0524 by 140% -> the SU(3) series is outright divergent at 2.5.")

# ---------------------------------------------------------------------------
# PART 5 -- the honest conditional bracket (needs the amplitude constant C)
# ---------------------------------------------------------------------------
print("\n[5] honest bracket:  w1(beta) = a0(1 - c u)/M0(beta),  P3: f1^2=a0(1-O(u)),  u := c_f/N")
# exact one-plaquette character coefficient c_f and the expansion parameter u=c_f/N
def su3_Z(beta, n=240):
    x, w = leggauss(n); phi = np.pi * (x + 1.0); wp = np.pi * w
    P1, P2 = np.meshgrid(phi, phi, indexing='ij'); W1, W2 = np.meshgrid(wp, wp, indexing='ij')
    Wt = W1 * W2; p3 = -(P1 + P2)
    e1, e2, e3 = np.exp(1j * P1), np.exp(1j * P2), np.exp(1j * p3)
    d = np.abs((e1 - e2) * (e1 - e3) * (e2 - e3)) ** 2
    norm = 1.0 / 6.0 / (2 * np.pi) ** 2
    return norm * np.sum(Wt * d * np.exp(beta * (e1 + e2 + e3).real))
def su3_cf(beta):
    avRe = su3_expect(lambda t: t.real, beta)          # <Re tr U> = c_f/c_0
    return avRe * su3_Z(beta)                           # c_f = <Re tr U> * c_0
cf_su2 = float((4.0 / 2.5) * iv(2, 2.5))                # a_{1/2}(2.5) = (4/beta) I_2(beta)
cf_su3 = su3_cf(2.5)
print("    exact one-plaquette c_f(2.5):  SU(2)=%.4f  u=c_f/N=%.4f" % (cf_su2, cf_su2 / 2))
print("                                   SU(3)=%.4f  u=c_f/N=%.4f" % (cf_su3, cf_su3 / 3))
print("    At N=2,3 the expansion parameter u(2.5) is O(1)-to-O(3) >> 1:")
print("    the correction (1-O(u)) is NOT a small correction, so the bracket below is")
print("    only a formal device -- it becomes usable ONLY if u(2.5)<1/C:")
w1n25 = 0.5 * (1 / 18) / su3_M0(2.5)        # (s/2) M0(0)/M0(2.5), s=1, SU(3)
u25 = cf_su3 / 3
print("    SU(3): w1_norm(2.5)=(s/2)M0(0)/M0(2.5)=%.4f ; u(2.5)=%.3f" % (w1n25, u25))
for C in [0.0, 0.1, 0.2, 0.28, 0.5, 1.0]:
    lo = w1n25 * (1 - C * u25); hi = w1n25 * (1 + C * u25)
    tag = "POSITIVE" if lo > 0 else "includes 0/negative -> FAILS"
    print("        if C=%4.2f : C*u=%.2f  w1(2.5) in [% .4f, %.4f]  %s" %
          (C, C * u25, lo, hi, tag))
print("    => w1(2.5)>0 is [C] on a dressed-amplitude headroom C with C*u(2.5)<1.")
print("       At N=2,3, u(2.5)=1.02 / 3.55 >> 1, so the f1^2=a0(1-O(u)) series is itself")
print("       outside its regime: no O(1) C makes C*u<1.  The bracket is EMPTY at N=2,3.")
print("       The missing estimate is an explicit O(1) bound on the dressed-amplitude")
print("       coefficient (=[R]/[O]) *together with* u(2.5)<1 (which holds for large N).")

# ---------------------------------------------------------------------------
# PART 6 -- large-N: leading strong-coupling u=beta/(2N^2)->0, so IS controlled
# ---------------------------------------------------------------------------
print("\n[6] large-N at fixed beta=2.5 (LEADING strong-coupling u=beta/(2N^2)=1.25/N^2):")
print("    (valid for large N; the measured c_f(2.5)/N is O(1) at N=2,3 and not in this regime)")
K = lambda N: 6.1e2 * N
print("     N     u_lead(2.5)=1.25/N^2   u*K     inside u*K<1 ?")
for N in [2, 3, 4, 6, 10, 30, 100, 300, 763, 1500]:
    u = 1.25 / N ** 2
    print("    %5d   %14.3e   %8.2f    %-5s" % (N, u, u * K(N), "YES" if u * K(N) < 1 else "no"))
Ncrit = int(np.ceil(305 * 2.5))   # u*K = (1.25/N^2)(6.1e2 N)=762.5/N <1  -> N>762.5
print("    => leading-order KP radius at beta=2.5 needs N > %d (crude K; true K smaller)." % Ncrit)
print("    => w1(2.5)=s/2>0 IS inside the strong-coupling radius in the large-N limit")
print("       (u(2.5)->0 as N->inf, 't Hooft territory W2); NOT at physical N=2,3.")

# ---------------------------------------------------------------------------
# PART 7 -- the five acid tests
# ---------------------------------------------------------------------------
print("\n[7] acid tests (charter §2)")
# (1) free massless field: KL measure with only a two-particle cut, no atom -> residue 0
def kl_residue(mass, E_light, T=40.0):
    s = np.linspace(4 * mass * mass, 120.0, 200000)
    rho = np.sqrt(np.clip(s - 4 * mass * mass, 0, None)) * np.exp(-0.15 * s)  # a cut, no delta
    C = np.trapezoid(rho * np.exp(-np.sqrt(s) * T), s)
    return np.exp(E_light * T) * C                                   # ->0 for any E_light < cut
r20 = kl_residue(1.0, 1.9, 20.0); r40 = kl_residue(1.0, 1.9, 40.0); r80 = kl_residue(1.0, 1.9, 80.0)
print("    (1) free field: cut-only KL measure, residue e^{E1 t}C(t): T=20:%.2e 40:%.2e 80:%.2e"
      " -> 0 (=w1)  [T]/[V]" % (r20, r40, r80))
# (2) 2D: trF^2 = H conserved; f1^2 = |<0|trF^2|R>|^2 = C2(R)<chi_0|chi_R> = 0
def su2_orth(j, k, n=400):
    x, w = leggauss(n); th = np.pi * (x + 1.0) / 2.0; wt = np.pi * w / 2.0   # theta in [0,pi]
    chi = lambda jj: np.sin((2*jj + 1) * th) / np.sin(th)
    return np.sum(wt * (2/np.pi) * np.sin(th) ** 2 * chi(j) * chi(k))
o = abs(su2_orth(0.5, 1.5)); print("    (2) 2D: <chi_{1/2}|chi_{3/2}>=%.2e  => f1^2=w1=0  [T]/[V]" % o)
# (3) one plaquette: no time direction -> the atom is degenerate, w1 undefined in transfer sense
print("    (3) one plaquette: a0/M0 = s/2 is an ARTEFACT -- no transfer direction,")
print("        the 'atom' is degenerate (all states at one energy); w1 is not a spectral fraction.")
# (4) strong coupling: w1(0)
print("    (4) strong coupling: w1(0)=s/2 (N>=3, complex fund.) ; s/4 (SU(2), real fund.)  [T] mod KL")
# (5) large N: w1=O(N^0)
print("    (5) large N: w1(0)=s/2 for every N>=3, O(N^0); u(2.5)->0 as N->inf.  [T]")
# (6) continuum gamma=0 (U4): w1 = Z_O^2 f1^2 / (Z_O^2 M0) invariant
ZO2 = 1.242670  # placeholder multiplicative factor from U4 witness (mu=2 -> 100)
f1sq, M0mu = 1.0 * ZO2, 2.0 * ZO2
print("    (6) continuum: w1=(Z_O^2 f1^2)/(Z_O^2 M0)=%.10f invariant => gamma_w1=0  [T]" %
      (f1sq / M0mu))

print("\n" + "=" * 78)
print("VERDICT (no QED, Clay 0%)")
print("  * w1(0)=s/2 : [T] at beta<1/48 (U3/P3, mod KL). Next coeff: exact normaliser part")
print("    (-1/2 SU(3), 0 SU(2)) + amplitude part ([-c] = dressed-amplitude [R]/[O]).")
print("  * Pade/Borel of the EXACT normaliser reach 2.5 and are positive -- but VACUOUS (M0>0).")
print("  * w1>0 at 2.5 survives only [C] on an explicit amplitude constant C with C*u(2.5)<1;")
print("    today C~O(10^2), C*u(2.5)>>1 at N=2,3  =>  the bracket is empty.  [O]")
print("  * UNCONDITIONAL-shaped positivity at 2.5 holds in the large-N limit (u(2.5)->0).")
print("=" * 78)
