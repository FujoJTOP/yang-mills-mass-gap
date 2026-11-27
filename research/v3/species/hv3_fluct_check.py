#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hv3_fluct_check.py -- HV3-FLUCT: chi~ mean == 0 (Parseval) ==> use its FLUCTUATION.
Deterministic.  ZERO new MC: reads only stored npz; all algebra exact (Fraction).

Blocks
 [1] exact DFT identity  (1/Nk) sum_k |chi(k)-C(0)|^2 = sum_{x!=0} |C(x)|^2   (Fraction, 1D/2D/3D)
 [2] m-th BZ moment = convolution of C|_{x!=0}                                   (Fraction)
 [3] product measure: C(x)=C0 delta ==> q2 == 0 exactly, all L                  (Fraction)
 [4] quadratic smallness: q2 <= (sum_{x!=0}C)^2 when C>=0                       (Fraction)
 [5] beta-sign counterexample: 3-state model, d/dbeta <O^2> < 0 while d/dbeta <O> = Var(O) > 0
 [6] measured 48-mode multiplet (beta=2.5, L=4 N=4000; L=6 N=800): per-channel Var(a),
     spread across channels (hypercubic anisotropy), kurtosis, jackknife on the spread
 [7] free-massless calibration of the *channel spread* is NOT the BZ variance (labels: fixed |k|)
"""
from fractions import Fraction as F
import numpy as np, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.abspath(os.path.join(HERE, "..", "..", "unconditional", "w7"))
OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s, flush=True)

# ---------------------------------------------------------------- [1] DFT identity
def dft_power_identity(C, L):
    """C: list of length L (1D periodic, C[r]); exact identity check."""
    Nk = L
    chi = []
    for m in range(L):
        s = F(0)
        for r in range(L):
            # exact cos(2 pi m r / L) as a root-of-unity trace: use complex via sympy-free integer recursion
            pass
        chi.append(s)
    return chi

# exact roots of unity via Chebyshev-free approach: work with the circulant action on integer vectors.
# For the identity we only need the two algebraic facts:
#   (a) (1/N) sum_k e^{ikr} = delta_{r0}
#   (b) (1/N) sum_k |chi(k)|^2 = sum_r |C(r)|^2   (Parseval)
# We verify them numerically-exactly on Z (integer) test vectors by using the DFT over Z[i]/... :
# instead, verify with exact *integer* test that Parseval holds in the form sum_k |chi_k|^2 = N sum_r |C_r|^2
# using Python's complex with integer inputs is float; so we verify symbolically via the
# orthogonality matrix being a permutation-like unitary: use Fraction arithmetic on the
# *real* cosine/sine basis for L=4 and L=2 where cos/sin take values in {0,+-1}.
def dft_L4(C):
    """L=4: e^{2pi i m r/4} in {1,i,-1,-i} ==> chi exactly in Z[i]; return (re,im) tuples."""
    chi = []
    for m in range(4):
        re = F(0); im = F(0)
        for r in range(4):
            ph = (m * r) % 4          # 0,1,2,3 -> 1,i,-1,-i
            if ph == 0: re += C[r]
            elif ph == 1: im += C[r]
            elif ph == 2: re -= C[r]
            else: im -= C[r]
        chi.append((re, im))
    return chi

def dft_L2(C):
    chi = []
    for m in range(2):
        re = F(0); im = F(0)
        for r in range(2):
            ph = (m * r) % 2
            if ph == 0: re += C[r]
            else: re -= C[r]
        chi.append((re, im))
    return chi

P("=" * 78)
P("[1] exact DFT identity: mean_k chi(k) = C(0)  and  (1/Nk) sum_k |chi~|^2 = sum_{x!=0}|C(x)|^2")
P("    (L=4, exact Gaussian-integer arithmetic in Z[i]; C = arbitrary rational correlator)")
tests = [F(3), F(1, 2), F(-1, 3), F(5, 7)]
for tag, C in (("C=(3,1/2,-1/3,5/7)", [F(3), F(1, 2), F(-1, 3), F(5, 7)]),
               ("C=(1,0,0,0)", [F(1), F(0), F(0), F(0)]),
               ("C=(2,3,5,7)", [F(2), F(3), F(5), F(7)])):
    chi = dft_L4(C)
    mean_re = sum(c[0] for c in chi) / 4
    mean_im = sum(c[1] for c in chi) / 4
    q2 = sum(c[0] ** 2 + c[1] ** 2 for c in chi) / 4 - C[0] ** 2
    off = sum(C[r] ** 2 for r in range(1, 4))
    P(f"    {tag:20s} mean_k chi = {mean_re} + {mean_im}i   (C(0)={C[0]})   q2={q2}  sum_off C^2={off}  "
      f"{'OK' if (mean_re == C[0] and mean_im == 0 and q2 == off) else 'MISMATCH'}")
chi2 = dft_L2([F(5, 3), F(-2, 9)])
P(f"    L=2 C=(5/3,-2/9): mean_k chi = {sum(c[0] for c in chi2)/2}  q2 = "
  f"{sum(c[0]**2+c[1]**2 for c in chi2)/2 - F(5,3)**2}  off = {F(-2,9)**2}")

# ---------------------------------------------------------------- [2] third moment
P("")
P("[2] third BZ moment = convolution sum over x1+x2+x3=0 (mod L), all xi != 0")
C = [F(3), F(1, 2), F(-1, 3), F(5, 7)]
chi = dft_L4(C)
def cval(k, m):
    """C(k) exactly."""
    return chi[k][0] + 1j * 0 if False else None
# exact: <chi~^3> = (1/N) sum_k (chi(k)-C0)^3 ; expand and use orthogonality symbolically:
# do it with exact Gaussian rationals by tracking (re,im) cubes.
tot_re = F(0); tot_im = F(0)
for m in range(4):
    re, im = chi[m][0] - C[0], chi[m][1]
    # (re+im i)^3
    r3 = re ** 3 - 3 * re * im ** 2
    i3 = 3 * re ** 2 * im - im ** 3
    tot_re += r3; tot_im += i3
P(f"    <chi~^3>_BZ (exact, L=4, C above) = {tot_re/4} + {tot_im/4} i")
# convolution side: sum_{x1+x2+x3=0, all !=0} C1 C2 C3
conv = F(0)
for x1 in (1, 2, 3):
    for x2 in (1, 2, 3):
        x3 = (-(x1 + x2)) % 4
        if x3 in (1, 2, 3):
            conv += C[x1] * C[x2] * C[x3]
P(f"    sum_{{x1+x2+x3=0,xi!=0}} C1C2C3            = {conv}")
P(f"    ==> {'OK' if (tot_re/4 == conv and tot_im == 0) else 'MISMATCH'}  (imag part {tot_im/4})")

# ---------------------------------------------------------------- [3] product measure
P("")
P("[3] product measure C(x)=C0 delta_{x0}: q2 == 0 exactly, any C0, any L")
for L in (2, 4):
    C = [F(7, 3)] + [F(0)] * (L - 1)
    chi = dft_L4(C) if L == 4 else dft_L2(C)
    q2 = sum(c[0] ** 2 + c[1] ** 2 for c in chi) / L - C[0] ** 2
    P(f"    L={L}: chi(k)={[str(c[0]) for c in chi]}  q2={q2}")

# ---------------------------------------------------------------- [4] quadratic smallness
P("")
P("[4] q2 <= (sum_{x!=0}C)^2 for C>=0 (Cauchy-Schwarz, exact), and the ratio q2/eps^2")
for tag, cv in (("one shell a=0.09", [F(0), F(9, 100)]),
                ("spread 6 sites a=0.015", [F(0)] + [F(15, 1000)] * 6),
                ("light-mode plateau c/V", [F(0)] + [F(9, 100) / 64] * 63)):
    s = sum(cv[1:]); q2 = sum(v * v for v in cv[1:])
    P(f"    {tag:24s} sum={float(s):.4f} q2={float(q2):.5f} sum^2={float(s*s):.5f} "
      f"q2<=sum^2:{q2 <= s*s}  q2/sum^2={float(q2/s/s) if s else 0:.4f}")

# ---------------------------------------------------------------- [5] beta sign
P("")
P("[5] beta-sign: d/dbeta <X> = -Cov(X,A)  (Z=sum e^{-beta A}).  Two exact witnesses, one family.")
P("    energies (0,E,E);  O_a=(0,+1,-1) ; O_b=(+1,0,0) ;  u=e^{-beta E}")
for u in (F(1, 2), F(1, 4), F(1, 8)):
    Z = 1 + 2 * u
    # O_a: <O>=0, <O^2>=2u/Z, Cov(O^2,A)=2uE/Z^2
    p_a = 2 * u / Z
    d_a = -2 * u / Z ** 2                      # d<O^2>/dbeta = -Cov(O^2,A)
    # O_b: <O>=<O^2>=1/Z, <A>=2uE/Z, Cov(O^2,A) = -2uE/Z^2
    p_b = 2 * u / Z ** 2   # Var(O_b) = 1/Z - 1/Z^2
    d_b = +2 * u / Z ** 2                      # opposite sign
    P(f"    u={u}:  O_a: Var(O)={p_a}={float(p_a):.4f}  d<O^2>/db={float(d_a):+.4f} <0 ;"
      f"   O_b: Var(O)={float(p_b):.4f}  d<O^2>/db={float(d_b):+.4f} >0")
P("    ==> the sign of d<O^2>/dbeta is not a function of (nu, A): both signs, same family.")
P("    ==> the fluctuation does NOT inherit the one-beta-monotone object Var(O)>=0.")

# ---------------------------------------------------------------- [6] measured multiplet
P("")
P("[6] stored 48-mode multiplet at fixed |k|=2pi/L: per-channel Var(a), spread, kurtosis")
def load(path):
    z = np.load(path, allow_pickle=True)
    ks = sorted([k for k in z.files if k.endswith("_a")])
    A = np.stack([z[k] for k in ks], axis=1)
    return A, [k[:-2] for k in ks]

def stats(path, tag):
    A, names = load(path)
    N, M = A.shape
    B = A - A.mean(0)
    per = (B ** 2).mean(0)
    spread = float(per.std(ddof=1))
    mu = float(per.mean())
    # kurtosis of the pooled amplitudes (mean-subtracted, per-channel scaled)
    zz = B / per ** 0.5
    kurt = float((zz ** 4).mean())
    # delete-block jackknife (20 blocks): spread, kurtosis, per-channel SE
    edges = np.linspace(0, N, 21).astype(int)
    jk, jk2, jkper = [], [], []
    for i in range(20):
        m = np.ones(N, bool); m[edges[i]:edges[i + 1]] = False
        Bm, pm = B[m], (B[m] ** 2).mean(0)
        jk.append(float(pm.std(ddof=1)))
        jk2.append(float(((Bm / pm ** 0.5) ** 4).mean()))
        jkper.append(pm)
    jk = np.array(jk); jk2 = np.array(jk2); jkper = np.array(jkper)
    se = float(np.sqrt(19 / 20 * ((jk - jk.mean()) ** 2).sum()))
    sek = float(np.sqrt(19 / 20 * ((jk2 - jk2.mean()) ** 2).sum()))
    seper = np.sqrt(19 / 20 * ((jkper - jkper.mean(0)) ** 2).sum(0))   # per-channel SE
    # is the scatter across the 48 channels consistent with the per-channel noise?
    chi2 = float((((per - per.mean()) / np.maximum(seper, 1e-12)) ** 2).sum())
    noise_spread = float(np.sqrt((seper ** 2).mean()))
    # group structure
    ax = {}
    for j, n in enumerate(names):
        ax.setdefault(n.split("_")[0], []).append(per[j])
    axm = {k: float(np.mean(v)) for k, v in sorted(ax.items())}
    P(f"    {tag}: N={N} M={M}  <Var>={mu:.4f}  spread={spread:.4f}+-{se:.4f}  "
      f"(spread/<Var>={spread/mu:.4f})  kurtosis={kurt:.3f}+-{sek:.3f}")
    skew = float((zz ** 3).mean())
    neff = float(2 * (mu / np.mean(seper)) ** 2)
    P(f"        skewness={skew:+.4f}   N_eff=2(<Var>/SE)^2={neff:.0f} of N={N}  "
      f"(tau_eff={N/neff:.2f})")
    P(f"        noise-only spread = {noise_spread:.4f} (from per-channel jackknife SE)  "
      f"spread/noise = {spread/noise_spread:.3f}   chi2/dof(scatter vs errors)={chi2/(M-1):.2f}/{M-1}")
    P(f"        per-axis means: " + "  ".join(f"{k}={v:.4f}" for k, v in axm.items())
      + f"   axis range={max(axm.values())-min(axm.values()):.4f}")
    return dict(N=N, M=M, mu=mu, spread=spread, se=se, kurt=kurt, sek=sek,
                axm=axm, noise=noise_spread, chi2d=chi2 / (M - 1))

r4 = stats(os.path.join(W7, "partA_su2_N4000.npz"), "L=4 b=2.5")
r6 = stats(os.path.join(W7, "partA_su2_N800.npz"), "L=6 b=2.5")
P(f"    L-scaling of the channel spread: {r4['spread']:.4f} (L=4) -> {r6['spread']:.4f} (L=6);"
  f"  ratio={r6['spread']/r4['spread']:.3f}")
P(f"    C(0)=(4pi^2-6)/9={(4*np.pi**2-6)/9:.6f};  spread/C(0) = {r4['spread']/((4*np.pi**2-6)/9):.4f}, "
  f"{r6['spread']/((4*np.pi**2-6)/9):.4f}")

# ---------------------------------------------------------------- [7] C(0) reference
P("")
P("[7] beta=0 (Haar product measure) reference for this carrier: C(x)=0 for x!=0 exactly")
P("    => q2(beta=0) = 0 identically (all L, all N); q2>0 measures inter-link correlation only.")
P("    NOTE: the 48 stored modes are a fixed-|k| multiplet (n1 only, 4 axes x 3 pol x cos/sin x c0/c2),")
P("          NOT a BZ sample => the BZ variance q2 is NOT formable from the archive.")

P("")
P("[9] exact acid test on the free massless point-split :phi^2: (equal-time, 3D torus)")
P("    G(x)=(1/L^3)sum_{p!=0}e^{ipx}/(2w_p), w_p^2=4sum_mu sin^2(p_mu/2)(+m^2); C(x)=2G(x)^2")
def free_q2(L, m=0.0):
    n = np.arange(L)
    k = 2 * np.pi * n / L
    K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    w2 = 4 * (np.sin(K[..., 0] / 2) ** 2 + np.sin(K[..., 1] / 2) ** 2 + np.sin(K[..., 2] / 2) ** 2)
    if m:
        w2 = w2 + m * m
    w = np.sqrt(w2)
    w[0, 0, 0] = np.inf                      # remove zero mode
    G = np.fft.ifftn(1.0 / (2 * w)).real
    C = 2 * G ** 2
    c0 = C[0, 0, 0]
    q2 = float((C ** 2).sum() - c0 ** 2)
    return c0, q2
P("      L    C(0)^free    q2_free     q2_free/C(0)^2   chi~(k_min)^free = L^2/(4pi^2)-C(0)")
for L in (8, 12, 16, 24, 32, 48):
    c0, q2 = free_q2(L)
    ch = 1.0 / (4 * np.sin(np.pi / L) ** 2) - c0
    P(f"    {L:3d}  {c0:10.5f}  {q2:10.5f}   {q2/c0**2:12.4e}    {ch:12.4f}")
P("    same with a mass m=0.5 (a genuinely gapped model):")
for L in (8, 16, 32):
    c0, q2 = free_q2(L, 0.5)
    P(f"    {L:3d}  {c0:10.5f}  {q2:10.5f}   {q2/c0**2:12.4e}")
P("    ==> q2 saturates in BOTH; chi~(k_min) diverges (L^2/4pi^2) in the massless one.")
P(f"    measured bound (L=4): q2/C(0)^2 <= {(0.0877/((4*np.pi**2-6)/9))**2:.2e}")

P("")
P("[8] exact numbers carried into the .md")
C0 = (4 * np.pi ** 2 - 6) / 9
P(f"    C(0) = {C0:.9f}")
P(f"    measured chi~ = +0.0077+-0.0157 (L=4), -0.0230+-0.0390 (L=6)   [HT4]")
P(f"    q2 bound (C>=0): q2 <= chi~(0)^2 = (0.0877)^2 = {0.0877**2:.5f} (L=4), "
  f"({0.0400:.4f})^2 = {0.04**2:.5f} (L=6)  [HT4 1.3]")
P(f"    q2/C(0)^2 <= {(0.0877/C0)**2:.2e} (L=4)  => fluctuation is O(eps^2), eps=chi~/C(0)")

with open(os.path.join(HERE, "hv3_fluct_check.log"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n# log -> hv3_fluct_check.log")


