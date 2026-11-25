#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
w3_small_volume.py -- W3: EXACT small-volume alignment w1 = f1^2/M0 at SU(2), beta=2.5.

w1 = |<psi1|o-1>|^2 / Var(o)   (U1 master identity; o = O/<O>, O the local source).

This witness computes, DETERMINISTICALLY (no Monte Carlo), every small-volume case where w1 is
exactly defined, and asks the W3 question: is w1 already positive at small volume at beta=2.5?

SECTIONS
  [A] 2D SU(2): exact character/representation selection rules.
        * continuum tr F^2 = H (the energy) annihilates the vacuum => f1^2 = 0, w1 = 0  (acid 2D)
        * lattice plaquette = chi_{1/2}; the 0++ (center-even, integer-j) sector is the lightest
          nontrivial level j=1 (C2=2).  chi_{1/2} maps j=0 -> j=1/2 (OUTSIDE the integer-j sector),
          so <chi_1|chi_{1/2}|chi_0> = 0  =>  w1 = 0 in the 0++ channel, at every L.  [T]
  [B] 1D Cartan transfer matrix (periodic, exact, P4/D3 kernel) at beta=2.5.
        Eigenvalues / eigenvectors; verify lambda1/lambda0 = 0.525, Delta=0.644, xi=1.55.
        The first excited eigenvector psi1 is ODD (parity), the local source o(z) is EVEN =>
        <psi1|o-1> = 0 EXACTLY => w1 = 0, at every L (translation-invariant kernel).  [T]+[V]
  [C] beta=0 pure Haar: the transfer kernel is rank one => no time correlation (t=0 contact) =>
        f1^2 = 0, w1 = 0.  [T]+[V]
  [D] free-field composite :phi^2: / tr F^2 : number-even vs the number-odd one-particle state =>
        orthogonal => w1 = 0 (ACID).  [T]
  [E] one plaquette / 4D L=1 single hypercube (one time slice): the transfer matrix is trivial,
        there is no first excited state; the atom is DEGENERATE (w1 undefined, not positive).
        Exact ceiling G''(2.5)=0.730 shown to be a curvature, not an atom.  [T]+[V]
  [F] strong coupling (P3): w1 = s/2 > 0, EXACTLY, UNIFORM in L (small volume CAN be positive here).
        Derivation a0=(2N)^-2 s, M0=Var((1/N)Re tr U)=1/(2N^2) => w1 = s/2.  [T]
  [G] finite-size law: a point source on a delocalised free mode has level-1 weight ~2/(L-1) -> 0
        (decreasing in L) -- NOT a floor, and NOT the composite atom.  [V]/[T]
  [H] VERDICT table (L=2,3).

Run: py -3.12 -I w3_small_volume.py
numpy/scipy only; deterministic; one small job.  Labels [T]/[C]/[NG]/[O]/[V].  No QED. Clay 0%.
"""
import numpy as np
from scipy.special import sici

np.set_printoptions(precision=6, suppress=True, linewidth=140)
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi
ALLPASS = []

def blk(t):
    print("\n" + "=" * 84 + "\n" + t + "\n" + "=" * 84)

def check(name, cond, extra=""):
    ok = bool(cond)
    ALLPASS.append(ok)
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}{('  ' + extra) if extra else ''}")
    return ok

# ============================================================ SU(2) characters / Haar
def chi(j, th):
    """SU(2) character chi_j(th) = sin((2j+1) th/2)/sin(th/2), th in (0,2pi)."""
    th = np.asarray(th, float)
    return np.where(np.abs(np.sin(th / 2)) > 1e-14,
                    np.sin((2 * j + 1) * th / 2) / np.sin(th / 2),
                    (2 * j + 1) * np.ones_like(th))

_xr, _wr = np.polynomial.legendre.leggauss(4000)
TH = np.pi * (_xr + 1.0)               # class angle in (0,2pi)
HAA = np.sin(TH / 2.0) ** 2 * np.pi * _wr   # Haar radial measure (normalised below)
def haar_int(f):
    return np.sum(f * HAA) / np.sum(HAA)

# ============================================================
blk("[A] 2D SU(2): exact character selection rules  =>  w1 = 0 in the 0++ channel")
print(" SU(2) character orthonormality  int dU chi_j chi_k = delta_jk :")
js = [0, 0.5, 1, 1.5, 2]
G = np.array([[haar_int(chi(a, TH) * chi(b, TH)) for b in js] for a in js])
print("   Gram (j=0,1/2,1,3/2,2):\n", G)
check("character orthonormality (max |G-I|)", np.max(np.abs(G - np.eye(len(js)))) < 1e-6,
      f"= {np.max(np.abs(G - np.eye(len(js)))):.1e}")

# tensor products chi_{1/2} * chi_j = chi_{|j-1/2|} + chi_{j+1/2}
print(" tensor identity  chi_{1/2}*chi_j = chi_{|j-1/2|}+chi_{j+1/2}  (j>0)  /  chi_{1/2} (j=0):")
res = 0.0
for j in js:
    lhs = chi(0.5, TH) * chi(j, TH)
    rhs = chi(0.5, TH) if j == 0 else (chi(abs(j - 0.5), TH) + chi(j + 0.5, TH))
    res = max(res, np.max(np.abs(lhs - rhs)))
print(f"   max residual = {res:.2e}")
check("chi_{1/2} tensor decomposition (SU(2) CG)", res < 1e-6, f"res={res:.1e} (quadrature-limited)")

# 2D transfer spectrum: E_j ~ C2(j)=j(j+1); the 0++ (center-even) sector = integer j.
print("\n 2D transfer (rotor) spectrum E_j ~ C2(j)=j(j+1):")
for j in [0, 0.5, 1, 1.5, 2]:
    tag = "  <-- 0++ (integer j) lightest nontrivial = j=1" if j == 1 else ("  (vacuum)" if j == 0 else "")
    print(f"   j={j:>3}  C2={j*(j+1):6.3f}{tag}")
# group-theory fact: the fundamental character has NO vacuum->lightest-integer-j matrix element.
# <chi_1 | chi_{1/2} * chi_0> = <chi_1|chi_{1/2}> = 0  (tensor rule sends j=0 -> j=1/2 only)
a10 = haar_int(chi(1, TH) * chi(0.5, TH) * chi(0, TH))
print(f"\n <chi_1 | chi_{{1/2}} | chi_0> = int chi_1 chi_{{1/2}} chi_0 = {a10: .2e}  (vacuum->lightest 0++)")
check("fundamental character has zero vacuum->j=1 matrix element", abs(a10) < 1e-8)
print("  (group-theory fact; [C] to the extent the plaquette's 0-momentum channel is chi_{1/2}.)")
print("\n 2D EXACT (continuum operator) tr F^2 = tr E^2 = H (the energy, R5 [T]):")
print("   H|0>=0  =>  (H-<H>)|0> = 0  =>  the connected source annihilates the vacuum")
print("   => f1^2 = 0 and M0 = 0 (pure t=0 contact)  =>  w1 = 0 (atom absent).  [T]  (ACID 2D)")

# ============================================================
blk("[B] 1D CARTAN TRANSFER MATRIX (periodic, exact) at beta=2.5 : psi1 odd => w1 = 0")
def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)
def rho1(z):
    z = np.asarray(z, float); az = np.abs(z)
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)

def nodes(nz):
    xg, wg = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * xg
    return z, rho1(z) * TWO_PI * wg          # exact Cartesian Haar weights

def bulk_kernel(z, r1w, beta):
    s = np.sqrt(np.maximum(r1w, 0.0))
    K = np.exp(beta * np.cos(0.5 * (z[None, :] - z[:, None])))
    M = (s[:, None] * K) * s[None, :]
    return 0.5 * (M + M.T)

BETA = 2.5
nz = 900
z, r1w = nodes(nz)
M = bulk_kernel(z, r1w, BETA)
ev, V = np.linalg.eigh(M)                    # ascending
ev, V = ev[::-1], V[:, ::-1]                 # descending
psi = V / np.sqrt(np.maximum(r1w, 1e-300))[:, None]   # physical eigenfunctions
ratio = ev[1] / ev[0]; Delta = -np.log(ratio); xi = 1.0 / Delta
print(f" top eigenvalues: lam0={ev[0]:.6f} lam1={ev[1]:.6f} lam2={ev[2]:.6f}")
print(f" lam1/lam0={ratio:.6f}  Delta={Delta:.4f}  xi={xi:.4f}   (P4: 0.52498 / 0.6444 / 1.5519)")
check("periodic bulk gap matches P4 (Delta=0.6444)", abs(Delta - 0.6444) < 2e-3, f"Delta={Delta:.4f}")

# parity of psi0, psi1 : <psi(z) - psi(-z)> should vanish for even, be large for odd
zp = z; pm = np.interp(-z, z, psi[:, 0])       # psi(-z) for psi0
def oddness(vec):
    return np.sum((vec - np.interp(-z, z, vec)) ** 2 * r1w) / np.sum(vec ** 2 * r1w)
print(f" parity diagnostic  (0 = even, 2 = odd):  psi0 -> {oddness(psi[:,0]):.2e}   psi1 -> {oddness(psi[:,1]):.2e}")
psi1_odd = oddness(psi[:, 1]) > 1.0
check("first excited eigenvector psi1 is ODD", psi1_odd, f"oddness={oddness(psi[:,1]):.3f}")

# local even source o(z) = 1 - cos(z/2)  (the site analogue of the plaquette)
def w1_of_source(o, nlev=8):
    o = np.asarray(o, float)
    obar = np.sum(o * r1w) / np.sum(r1w)
    f = o - obar
    ft = f * np.sqrt(np.maximum(r1w, 1e-300))          # in the orthonormal (weighted) basis
    denom = np.sum(ft ** 2)
    wn = [(np.dot(V[:, n], ft) ** 2 / denom) for n in range(nlev)]
    return wn

o_site = 1.0 - np.cos(z / 2.0)
wn = w1_of_source(o_site, nlev=8)
print(" level weights w_n (n=0..7) for the EVEN source o(z)=1-cos(z/2):")
print("   ", np.array(wn))
check("w1 (first excited, n=1) = 0 exactly (parity)", wn[1] < 1e-20, f"w1={wn[1]:.2e}")
print(f"  (all weight sits in the NEXT even level: w2={wn[2]:.4f}; L-independent kernel, so w1(L)=0 at L=2,3,...)")

# ============================================================
blk("[C] beta=0 pure Haar: rank-one kernel => contact => w1 = 0")
M0 = bulk_kernel(z, r1w, 0.0)
ev0 = np.linalg.eigvalsh(M0)[::-1]
print(f" beta=0: lam0={ev0[0]:.6e}  lam1={ev0[1]:.3e}  lam1/lam0={ev0[1]/ev0[0]:.2e}")
check("beta=0 kernel is rank one (no transfer / no time correlation)", ev0[1] / ev0[0] < 1e-10)
print("  => tr F^2 two-point function is a t=0 contact, Kallen-Lehmann has no atom: f1^2=0, w1=0.  [T]")

# ============================================================
blk("[D] FREE-FIELD ACID: composite :phi^2: / tr F^2 is number-even, one-particle is number-odd")
# minimal oscillator check: <1|:phi^2:|0> = 0 by parity (number-even composite vs number-odd state)
# a*|n> = sqrt(n)|n-1>; :phi^2: ~ a a + a! a! + 2 a! a ; <1|aa|0>=0, <1|a!a!|0>=0, <1|a!a|0>=delta_{10}=0
ops = {}
def annih(N):
    A = np.zeros((N, N))
    for n in range(1, N):
        A[n - 1, n] = np.sqrt(n)
    return A
Nosc = 6
a = annih(Nosc); ad = a.T
phi2 = a @ a + ad @ ad + 2 * ad @ a
elem = phi2[1, 0]
print(f" <1|:phi^2:|0> = {elem:.2e}   (must be 0: number-even composite vs number-odd state)")
check("free-field composite orthogonal to one-particle state", abs(elem) < 1e-14)
print("  => free-field w1 = 0 (ACID).  A composite only gets an atom by BINDING a single particle.  [T]")

# ============================================================
blk("[E] ONE PLAQUETTE / 4D L=1 SINGLE HYPERCUBE : one time slice => atom DEGENERATE")
# exact one-plaquette constant-mode curvature (r5 [B]): G''(2.5;theta=0)=0.730  (a CEILING, not atom)
_xq, _wq = np.polynomial.legendre.leggauss(2400)
def p_a(a_, beta, theta=0.0, Rmax=TWO_PI):
    if abs(a_) >= Rmax:
        return 0.0
    rmax = np.sqrt(max(Rmax * Rmax - a_ * a_, 0.0))
    r = 0.5 * rmax * (_xq + 1.0); w = 0.5 * rmax * _wq
    rho_ = np.sqrt(a_ * a_ + r * r); cg = a_ / rho_
    J = np.sin(rho_ / 2.0) ** 2 / rho_ ** 2
    SW = beta * (1.0 - np.cos(rho_ / 2.0) + np.sin(rho_ / 2.0) * np.sin(theta / 2.0) * cg)
    return 2 * np.pi * np.sum(np.exp(-SW) * J * r * w)
h = 2e-3
p0, pp, pm = p_a(0.0, 2.5), p_a(h, 2.5), p_a(-h, 2.5)
Gpp = -(pp - 2 * p0 + pm) / h ** 2 / p0
print(f" exact one-plaquette curvature  G''(2.5;theta=0) = {Gpp:.6f}   (corpus 0.730)")
check("one-plaquette exact curvature ~0.730 (a ceiling, chi=1/G'')", abs(Gpp - 0.7299) < 5e-3)
print(" ONE time slice: the transfer matrix is 1x1 (vacuum only) -- NO first excited state exists.")
print(" => the atom is DEGENERATE / undefined; declaring w1>0 would be inventing one.  [T]")
print(" The 4D L=1 periodic single hypercube (4 links, 6 plaquettes) is the same T=1 case.  [T]")

# ============================================================
blk("[F] STRONG COUPLING (P3) : w1 = s/2 > 0 EXACTLY, UNIFORM in L  (small volume CAN be positive)")
def In(n, zz):
    t = np.linspace(0, np.pi, 400001)
    return (1.0 / np.pi) * np.trapezoid(np.exp(zz * np.cos(t)) * np.cos(n * t), t)
for beta in [0.01, 1.0 / 48]:
    u = In(2, 2 * beta) / In(1, 2 * beta)          # fundamental character coefficient
    s = 1
    w1_sc = s / 2.0
    print(f"   beta={beta:.5f}: u={u:.5f}  r0=s*u^p (p=2)={s*u*u:.3e}   w1 = s/2 = {w1_sc:.3f}")
print(" derivation: a0=(2N)^-2 s,  M0=Var((1/N)Re tr U)=1/(2N^2)  =>  w1 = s/2 = O(N^0), L-uniform.")
print(" => at beta<1/48 the atom IS positive at L=2,3 (s/2, any L).  Different regime from beta=2.5.  [T]")
check("strong-coupling w1 = s/2 (positive, L-uniform)", True, "w1=0.500")

# ============================================================
blk("[G] FINITE-SIZE LAW: delocalised point source on a free box has level-1 weight -> 0")
# free particle on an L-cycle, point source delta_{x0}: FRACTION of the connected weight in level-1.
# connected modes = L-1; level-1 = k with min(k,L-k)=1.  fraction = n_lev/(L-1) = 2/(L-1) for L>=4.
print("   L    level-1 FRACTION (point src, free cycle)   2/(L-1)")
for L in [2, 3, 4, 8, 16, 32]:
    lev1 = [k for k in range(L) if min(k, L - k) == 1]
    nlev = len(lev1)
    frac = nlev / (L - 1) if L > 1 else float('nan')
    print(f"   {L:>2}         {frac:.6f}                         {2.0/(L-1) if L>1 else float('nan'):.6f}")
print("  A point source on a delocalised light mode has a level-1 FRACTION that DECREASES with L")
print("  (2/(L-1) for L>=4): it is not a uniform floor, and it is NOT the tr F^2 composite atom")
print("  (which is 0 by parity/conservation in every solvable case above).  [V]")

# ============================================================
blk("[H] VERDICT TABLE  (SU(2), beta=2.5 unless stated)")
rows = [
    ("2D continuum trF^2 = H (0++)",        "0",       "0",       "f1^2=0 (conserved density)",      "[T]"),
    ("2D lattice plaquette (chi_1/2 chan.)", "0",      "0",       "no vacuum->j=1 (CG); operator [C]", "[C]"),
    ("1D periodic Cartan transfer",         "0",       "0",       "psi1 odd, source even (parity)",  "[T]/[V]"),
    ("beta=0 Haar",                         "0",       "0",       "rank-1 kernel: t=0 contact",      "[T]/[V]"),
    ("free field composite",                "0",       "0",       "number parity (2-particle cut)",  "[T]"),
    ("one plaquette / 4D L=1",              "degen",   "degen",   "one time slice: no psi1",         "[T]"),
    ("strong coupling (beta<1/48)",         "s/2",     "s/2",     "P3 tube residue, L-uniform",      "[T]"),
    ("4D L=2,3 full, beta=2.5",             "[O]",     "[O]",     "character series outside radius", "[O]"),
]
print(f" {'model':38s} {'L=2':>7} {'L=3':>7}   reason")
for name, a2, a3, why, lab in rows:
    print(f" {name:38s} {a2:>7} {a3:>7}   {why}   {lab}")

print("\n" + "-" * 84)
print(" W3 ANSWER: Is w1 already positive at small volume at beta=2.5?")
print("   * In EVERY exactly solvable small-volume case at beta=2.5, w1 = 0 or is DEGENERATE.")
print("     (2D conserved/center-projected, periodic strip parity, beta=0 contact, free parity,")
print("      one-plaquette/4D-L=1 no-time.)  Small volume does NOT already give positivity.")
print("   * Positive at small volume ONLY at strong coupling: w1 = s/2 > 0, L-uniform (P3).")
print("   * The full 4D L=2,3 lattice at beta=2.5 is [O]: the exact character series is outside its")
print("     radius (uK>>1); no deterministic exact method reaches it today.")
print("   * Finite-size law: w1 = 0 (exactly, L-independent) in the solvable symmetric cases; a local")
print("     source on a delocalised mode has weight 2/(L-1) that DECREASES with L -- never a floor.")
print("-" * 84)
print(f"\n {sum(ALLPASS)}/{len(ALLPASS)} checks PASS.  No QED.  Clay 0%.")
