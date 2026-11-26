"""
h08_z1.py — H08 / zero location z1 of D_O(z) = det(1 - z T_O).

Checks the five things that decide the route.  numpy only.

  B1  z1 = 1/limsup_k |C(k)|^{1/k}   (exact identity, needs k -> infinity)
  B2  Gauss/Pade [M-1/M] : max node x <= lambda_1, monotone up  =>  m1_det(M) >= m1_det
      and the f1^2 = 0 case: the lightest state is INVISIBLE => m1_det > m1 strictly
  B3  conditioning: moments are an exponentially ill-conditioned encoding of the endpoint
  B4  2D YM: C(k) == 0 identically (O = c*H conserved) => D_O == 1, no zero, yet m>0
  B5  single plaquette: T = mult by exp(-beta S_W) => spectrum is a CONTINUUM [e^-2b,1];
      ||(1-zT)^-1|| = 1/(|z| dist(1/z,[e^-2b,1]))  (no isolated zero);
      Tr T^k = Z(k beta) exactly.
  B6  the OTHER side: Z itself gives LOWER bounds m_hat(k) = -(1/ak) log(Z-1) < m1,
      converging exponentially in k -- but only at relative precision e^{-ak m1}.
"""
import numpy as np
from math import lgamma


# ---------------------------------------------------------------- moments
def Ck(lam, f2, K):
    """C(k) = sum_i f2_i lam_i^k , k = 0..K, computed with a max-shift (no overflow)."""
    lam = np.asarray(lam, float)
    f2 = np.asarray(f2, float)
    out = np.empty(K + 1)
    for k in range(K + 1):
        e = k * np.log(lam)
        top = e.max()
        out[k] = top + np.log(np.sum(f2 * np.exp(e - top)))
    return out            # returns log|C(k)|


def Cvals(lam, f2, K):
    return np.exp(Ck(lam, f2, K))


# ---------------------------------------------------------------- Pade [M-1/M]
def pade_den(m, M):
    """Q(z) = 1 + q1 z + ... + qM z^M from the [M-1/M] linear system
       sum_{j=0}^M q_j m[M-1+i-j] = 0 , i = 1..M (m[k<0] = 0)."""
    def mget(k):
        return m[k] if 0 <= k < len(m) else 0.0

    A = np.zeros((M, M))
    b = np.zeros(M)
    for i in range(1, M + 1):
        b[i - 1] = -mget(M - 1 + i)
        for j in range(1, M + 1):
            A[i - 1, j - 1] = mget(M - 1 + i - j)
    q = np.linalg.solve(A, b)
    return np.concatenate(([1.0], q))


def gauss_nodes(m, M, keep_in=(0.0, 1.0)):
    """Gauss nodes of dnu:  poles of [M-1/M] are 1/x_i.  Keep real nodes in `keep_in`."""
    q = pade_den(m, M)[::-1]
    p = np.roots(q)
    p = p[np.abs(p.imag) < 1e-8].real
    p = p[p > 1e-12]
    x = 1.0 / p
    lo, hi = keep_in
    return np.sort(x[(x > lo) & (x <= hi)])


# =================================================================== B1
print("=" * 78)
print("B1  z1 = 1/limsup_k |C(k)|^{1/k}        spectral measure with 4 atoms")
lam = np.array([0.95, 0.80, 0.65, 0.40])
f2 = np.array([1.00, 0.30, 0.20, 0.10])
z1 = 1.0 / lam.max()
print(f"    lambda_1 = {lam.max():.10f}  =>  z1 = 1/lambda_1 = {z1:.10f}")
K = 40000
lg = Ck(lam, f2, K)
kk = np.arange(1, K + 1)
r_raw = np.exp(lg[1:] / kk)                 # limsup candidate
r_norm = np.exp((lg[1:] - lg[0]) / kk)      # normalised
print(f"    {'k':>7} {'|C(k)|^(1/k)':>18} {'err':>11} {'(C(k)/C(0))^(1/k)':>19} {'err':>11}")
for k in (5, 10, 50, 100, 1000, 10000, 40000):
    print(f"    {k:7d} {r_raw[k-1]:18.12f} {abs(r_raw[k-1]-lam.max()):11.2e}"
          f" {r_norm[k-1]:19.12f} {abs(r_norm[k-1]-lam.max()):11.2e}")
delta = 1.0 - lam[1] / lam[0]
print(f"    isolation in the transfer variable: 1-lambda_2/lambda_1 = {delta:.4f}"
      f"  -> resolves once k >~ 2/delta = {2/delta:.1f}")
print("    -> limsup is exact but it is 'mass = exponential decay rate'. FINITE k: nothing.")

# =================================================================== B2
print("=" * 78)
print("B2  Gauss/Pade [M-1/M]:  x_max <= lambda_1, monotone up  =>  m1_det(M) is an UPPER bound")
a = 1.0
print(f"    {'M':>3} {'x_max':>18} {'m1_det(M)':>16} {'m1 - m1_det':>13}")
for M in range(1, 8):
    mm = Cvals(lam, f2, 2 * M - 1)
    x = gauss_nodes(mm, M)
    xmax = x.max()
    print(f"    {M:3d} {xmax:18.12f} {-np.log(xmax)/a:16.12f} "
          f"{(-np.log(lam.max())/a) + np.log(xmax)/a:13.3e}")
print(f"    exact  lambda_1 = {lam.max():.12f}   m1 = {-np.log(lam.max())/a:.12f}")
print("    -> EXACT at M = rank(nu) = 4 (Prony/Pade).  2M moments fix z1 IFF rank T_O <= M.")

print("    --- the f1^2 = 0 trap: a state at lambda=0.95 exists but has ZERO overlap ---")
lam2 = np.array([0.80, 0.65, 0.40])
f22 = np.array([0.30, 0.20, 0.10])
m1_true = -np.log(0.95) / a
mm = Cvals(lam2, f22, 19)
xmax = gauss_nodes(mm, 10).max()
print(f"    true      m1     = {m1_true:.10f}   (state at lambda=0.95 exists in spec H)")
print(f"    m1_det (M=10)    = {-np.log(xmax)/a:.10f}   x_max = {xmax:.10f}")
print(f"    m1_det - m1      = {-np.log(xmax)/a - m1_true:.10f}  > 0  : the zero JUMPS OVER")
print("    -> inf{E_n : f_n^2>0} >= inf{E_n} = m1.  Direction error; more moments cannot fix it.")

# =================================================================== B3
print("=" * 78)
print("B3  conditioning:  relative noise 1e-10 on the moments -> endpoint error")
print("    (dense nu: 10 atoms, so every M below is a genuine under-determination)")
lamD = np.array([0.990, 0.980, 0.965, 0.940, 0.900, 0.850, 0.780, 0.680, 0.550, 0.400])
f2D = np.array([1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10])
print(f"    {'M':>3} {'x_max(exact)':>16} {'m1_det(M)':>13} {'m1-m1_det':>11}"
      f" {'median err@1e-10':>17} {'max err':>11}")
rng = np.random.default_rng(0)
for M in (3, 4, 5, 6, 8, 10):
    mm = Cvals(lamD, f2D, 2 * M - 1)
    x0 = gauss_nodes(mm, M).max()
    errs = []
    for _ in range(300):
        pert = mm * (1 + 1e-10 * rng.standard_normal(len(mm)))
        x = gauss_nodes(pert, M)
        if len(x):
            errs.append(abs(x.max() - x0))
    errs = np.array(errs) if errs else np.array([np.nan])
    print(f"    {M:3d} {x0:16.12f} {-np.log(x0)/a:13.10f} {(-np.log(lamD.max())/a)+np.log(x0)/a:11.3e}"
          f" {np.median(errs):17.3e} {errs.max():11.3e}")
print("    -> exact moments at 1e-10 relative noise give a biased/bounded endpoint, and the")
print("       bias shrinks only like the DELTA of the spectrum, not like the noise.")

# =================================================================== B4
print("=" * 78)
print("B4  2D YM:  Z(t) = sum_j (2j+1)^2 exp(-t j(j+1)) ,  O = c*H (conserved)")


def Z2d(t, J=500):
    j = np.arange(J + 1)
    return float(np.sum((2 * j + 1.0) ** 2 * np.exp(-t * j * (j + 1.0))))


print(f"    {'t':>5} {'Z(t)':>14} {'Z(2t)':>14} {'Z(3t)':>14} {'m=3t/(8A),A=1':>14} {'C(k)':>6}")
for t in (0.3, 0.7, 1.5):
    print(f"    {t:5.2f} {Z2d(t):14.8e} {Z2d(2*t):14.8e} {Z2d(3*t):14.8e}"
          f" {3*t/8:14.6f} {'0':>6}")
print("    -> Tr T^k = Z(kt) exact and nonzero, but C(k) = <0|O T^k O|0> = 0 for ALL k")
print("       (O = c*H, H|0> = 0)  =>  T_O = 0, D_O(z) == 1 :  NO ZERO AT ALL,  m = 3g^2/8 > 0.")

# =================================================================== B5
print("=" * 78)
print("B5  single plaquette (d=0): T = multiplication by f(U) = exp(-beta S_W)")
beta = 2.5
Nth = 2000001
th = (np.arange(Nth) + 0.5) * np.pi / Nth
th_edges = np.linspace(0, np.pi, Nth + 1)
measure = (2.0 / np.pi) * np.sin(0.5 * (th_edges[1:] + th_edges[:-1])) ** 2 * np.diff(th_edges)
fU = np.exp(-beta * (1.0 - np.cos(th)))          # spectrum of T  = [e^{-2b}, 1]
lo, hi = np.exp(-2 * beta), 1.0

print(f"    {'k':>3} {'Tr T^k (quad)':>20} {'Z(k*beta) exact':>20} {'rel.err':>10}")
for k in (1, 2, 3, 5, 10, 20):
    trn = float(np.sum(measure * fU ** k))
    s = 0.0
    for m in range(300):
        term = np.exp((2 * m + 1) * np.log(k * beta / 2) - lgamma(m + 1) - lgamma(m + 2))
        s += term
        if term < 1e-20 * abs(s):
            break
    zexact = np.exp(-k * beta) * 2 * s / (k * beta)
    print(f"    {k:3d} {trn:20.10e} {zexact:20.10e} {abs(trn-zexact)/zexact:10.2e}")
print("    (Z(b) = e^{-b} 2 I1(b)/b is the exact SU(2) single-plaquette partition function)")
print(f"    spectrum of T = [{lo:.6f}, {hi:.6f}]  (absolutely continuous)")
print("    --- resolvent:  ||(1-zT)^-1|| = sup_U 1/|1-z f(U)|  =  1/(|z| dist(1/z,[e^-2b,1]))")
print(f"    {'z':>12} {'numerical sup':>16} {'1/(|z|dist)':>16}")
for z in (0.5, 0.99, 1.0 + 1e-6, 2.0, 50.0, hi * 1.0003, hi * 1.5, 400.0):
    num = float(np.max(1.0 / np.abs(1.0 - z * fU)))
    u = 1.0 / z
    d = 0.0 if lo <= u <= hi else min(abs(u - lo), abs(u - hi))
    pred = np.inf if d == 0 else 1.0 / (abs(z) * d)
    print(f"    {z:12.6g} {num:16.8e} {pred:16.8e}")
print(f"    => singular set = [1, e^(2beta)] = [1, {np.exp(2*beta):.4f}]: a CONTINUUM, not a point.")
print("       z=1 is only the lower EDGE; T is not trace class => no Fredholm determinant.")

# =================================================================== B6
print("=" * 78)
print("B6  the OTHER side:  Z gives LOWER bounds on m1, at exponential precision cost")
a = 1.0
E = np.array([0.0, 1.00, 1.60])       # vacuum + lightest + next  (exact 3-level model)
w = np.array([1.0, 1.0, 2.0])         # degeneracies (vacuum simple)
m1 = E[1]
print(f"    exact 3-level: E={E.tolist()}, w={w.tolist()}, a=1  =>  m1 = {m1}")
print(f"    {'k':>3} {'Z-1':>16} {'m_hat(k)':>14} {'m1 - m_hat':>12} {'rel.err needed':>15}")
for k in (1, 2, 3, 5, 8, 12, 20, 30):
    Zm1 = float(np.sum(w[1:] * np.exp(-k * a * E[1:])))
    mh = -np.log(Zm1) / (a * k)
    print(f"    {k:3d} {Zm1:16.10e} {mh:14.10f} {m1-mh:12.3e} {np.exp(-k*a*m1):15.3e}")
print("    -> m_hat(k) <= m1 for EVERY k (lower bound!), ->m1 exponentially in k.")
print("       Cost: Z-1 must be known to RELATIVE precision e^{-ak m1} (column 5).")
print("    =>  D_O zero  : cheap-ish, UPPER bound on m1.      Z : LOWER bound, exponential cost.")
print("        The mass gap needs a LOWER bound.  The zero is on the wrong side.")
print("=" * 78)
print("done.")
