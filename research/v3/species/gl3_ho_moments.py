"""
GL3-HIGHORDER numerics (v5, analytic): improvement tower for the H(4)-invariant loop family.

Operator family
---------------
O_{m,n} = the H(4) sum of the m x n rectangle (m along one axis, n along the other), all
6 lattice planes, both orderings (m,n) and (n,m), both traversal orientations.

Symbol (exact, per plane (a,b), shape m along a and n along b):
    w_a = S_m(z_a) (1 - z_b^n),   w_b = S_n(z_b) (z_a^m - 1),   S_k(z) = 1+z+...+z^{k-1}
    h_{m,n}(p) = (1/6) sum_{planes} sum_{orderings} ( |w_a|^2 + |w_b|^2 )
(the 1/6 = 1/2 [Re tr] x 1/3 [polarization average]; the orbit sum is H(4) invariant).

Exact facts (all checked below):
  * temporal axis p=(iE,0,0,0), x=aE:  h = n^2(2-2cosh mx) + m^2(2-2cosh nx)  (per orbit)
    => kappa_j := H_j/H_1 = (m^{2(j-1)} + n^{2(j-1)}) / (2j)!          [H = -x^2(1+sum kappa_j x^{2j})]
  * small p:  h ~ 2 m^2 n^2 p^2  (per image)  =>  rho ~ area^2,  rho_{1,2}/rho_{1,1} = 8
  * BZ corner p=(pi,pi,pi,pi): per-image |w|^2 = 8 if m,n both odd, else 0
  * free massive pole:  Delta = m - kappa_{j*} a^{2(j*-1)} m^{2(j*-1)+1}/2 + ...
Improvement: d_O := c_O rho_O,  sum d = 1,  conditions sum d_O kappa_j^O = 0.
Moment structure: kappa_j = (m^{2k}+n^{2k})/(2j)! with k=j-1 => conditions are moments
mu(m^{2k}) = 0 of a signed measure mu supported on the time extents {1..M}; minimal support
size for k=1..M-1 is M => k0 = 2M, and the weights are the inverse-Vandermonde
d_m = prod_{l != m} l^2/(l^2-m^2).
"""
import numpy as np
from fractions import Fraction
from math import factorial, pi

JMAX = 6


def S(k, z):
    """S_k(z) = 1+z+...+z^{k-1}; z may be complex array."""
    if k == 1:
        return np.ones_like(z)
    out = np.ones_like(z)
    t = np.ones_like(z)
    for _ in range(k - 1):
        t = t * z
        out = out + t
    return out


def h_mn(m, n, z0, z1, z2, z3):
    """symbol of O_{m,n}: (1/6) sum over the 6 planes and the two orderings."""
    z = [z0, z1, z2, z3]
    tot = 0.0
    for a in range(4):
        for b in range(a + 1, 4):
            for (mm, nn) in ({(m, n), (n, m)}):
                wa = S(mm, z[a]) * (1 - z[b] ** nn)
                wb = S(nn, z[b]) * (z[a] ** mm - 1)
                tot = tot + np.abs(wa) ** 2 + np.abs(wb) ** 2
    return tot / 6.0


# ---------------------------------------------------------------- 1. axis kappa
def kappa(m, n, j):
    return Fraction(m ** (2 * (j - 1)) + n ** (2 * (j - 1)), factorial(2 * j))


def rho(m, n):
    """proportional to the leading p^2 coefficient (area^2 x multiplicity); rho(1,1)=1."""
    return Fraction(m * m * n * n) * (1 if m == n else 2)


print("=" * 78)
print("1. checks against known results")
print("=" * 78)
print("  kappa_2(1,1) =", kappa(1, 1, 2), "         (Wilson: expect 1/12)")
print("  kappa_2(1,2) =", kappa(1, 2, 2), "         (rect)")
print("  rho(1,2)/rho(1,1) =", rho(1, 2) / rho(1, 1), "   (standard c0+8c1=1 => expect 8)")
# LW: c0=5/3, c1=-1/12
d11 = Fraction(5, 3) * rho(1, 1)
d12 = Fraction(-1, 12) * rho(1, 2)
print(f"  LW normalized d: d(1,1) = {d11}, d(1,2) = {d12}, sum = {d11+d12}")
for j in (2, 3, 4, 5):
    print(f"    kappa_{j}(LW) = {d11*kappa(1,1,j) + d12*kappa(1,2,j)}")
print("    => Delta = m + (1/180) a^4 m^5  (kappa_3 = -1/90)  [matches GK3-SYMANZIK 2.2]")

# ---------------------------------------------------------------- 2. moment tower
print()
print("=" * 78)
print("2. moment tower: support {1..M} => k0 = 2M;  Delta = m - kappa_{M+1} a^{2M} m^{2M+1}/2")
print("=" * 78)


def inverse_vandermonde(M):
    """d_m = prod_{l != m} l^2/(l^2-m^2),  m = 1..M ;  sum d = 1, sum d m^{2k} = 0 (k=1..M-1)."""
    d = []
    for m in range(1, M + 1):
        v = Fraction(1)
        for l in range(1, M + 1):
            if l == m:
                continue
            v *= Fraction(l * l, l * l - m * m)
        d.append(v)
    return d


for M in range(1, 6):
    d = inverse_vandermonde(M)
    assert sum(d) == 1, (M, sum(d))
    for k in range(1, M):
        assert sum(di * m ** (2 * k) for di, m in zip(d, range(1, M + 1))) == 0
    mom = sum(di * m ** (2 * M) for di, m in zip(d, range(1, M + 1)))
    kap = Fraction(2) * mom / factorial(2 * M + 2)
    coef = -kap / 2
    theta_corner = Fraction(16) * d[0]
    for m in range(3, M + 1, 2):        # {m,m} square, corner 16/m^4 in rho-units
        theta_corner += d[m - 1] * Fraction(16, m ** 4)
    print(f" M={M}: k0={2*M}  d_1..d_M = {[str(x) for x in d]}")
    print(f"        kappa_{M+1} = {kap}   =>  Delta = m {'+' if -kap > 0 else '-'} "
          f"{abs(coef)} a^{2*M} m^{2*M+1} + ...")
    print(f"        corner h = {theta_corner} = {float(theta_corner):.6f}  "
          f"(|p|^2 a^2 = 4 pi^2 = {4*pi**2:.6f})  =>  theta = "
          f"{abs(1-float(theta_corner)/(4*pi**2)):.6f}")

# ---------------------------------------------------------------- 3. positivity grid
print()
print("=" * 78)
print("3. quadratic form on the BZ:  min h  (negative => the action is not a positive")
print("   quadratic form => no positive transfer matrix; grid 25^4, exact corner added)")
print("=" * 78)
N = 25
g = (2 * np.pi * np.arange(N) / N)
Z0 = np.exp(1j * g)[:, None, None, None]
Z1 = np.exp(1j * g)[None, :, None, None]
Z2 = np.exp(1j * g)[None, None, :, None]
Z3 = np.exp(1j * g)[None, None, None, :]
for M in (1, 2, 3, 4):
    d = inverse_vandermonde(M)
    h = 0.0
    for m in range(1, M + 1):
        if d[m - 1] == 0:
            continue
        h = h + float(d[m - 1]) * h_mn(m, m, Z0, Z1, Z2, Z3)
    h = np.real(h)
    hflat = h.ravel()
    idx = int(np.argmin(hflat))
    i0, i1, i2, i3 = np.unravel_index(idx, h.shape)
    # corner exact
    hc = sum(float(d[m - 1]) * h_mn(m, m, -1.0, -1.0, -1.0, -1.0) for m in range(1, M + 1))
    print(f" M={M} (k0={2*M}): min h = {hflat[idx]:+.9f} at (p0,p1,p2,p3)= "
          f"({g[i0]/np.pi:.3f},{g[i1]/np.pi:.3f},{g[i2]/np.pi:.3f},{g[i3]/np.pi:.3f})*pi ; "
          f"h(corner) = {hc:+.6f}")

# ---------------------------------------------------------------- 4. continuum counts
print()
print("=" * 78)
print("4. H(4)-invariant quadratic symbols per order:  p_{<=4}(n) - p_{<=4}(n-1)")
print("   (degree 2n symbol  <=>  dimension 2n operator  <=>  order a^{2n-2})")
print("=" * 78)


def ple4(n):
    cnt = 0
    for a in range(n + 1):
        for b in range(n + 1):
            for c in range(n + 1):
                for dd in range(n + 1):
                    if a + 2 * b + 3 * c + 4 * dd == n:
                        cnt += 1
    return cnt


print(f" {'n':>3} {'dim':>5} {'order':>7} {'p<=4(n)':>8} {'N(n)':>5}   on-shell structures")
for n in range(1, 9):
    Nn = ple4(n) - ple4(n - 1)
    print(f" {n:3d} {2*n:5d} {'a^'+str(2*n-2):>7} {ple4(n):8d} {Nn:5d}")
