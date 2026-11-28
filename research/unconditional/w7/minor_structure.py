import numpy as np
from scipy.special import iv
import itertools

# Attack (c): compute the 2x2 minors of M_f exactly for the actual c_d
#   (M_f)_{j'j} = sum_{d=|j-j'|}^{j+j'} c_d,   c_d = e^{-b}(I_{2d}(b) - I_{2d+2}(b)) > 0
# The minor for rows j1<j2 and columns j1'<j2':
#   D = M[j1,j1'] M[j2,j2'] - M[j1,j2'] M[j2,j1']
# Try to exhibit D as a POSITIVE combination of c_t (which would prove TN2).

def cvals(b, K):
    return np.array([np.exp(-b) * (iv(2 * d, b) - iv(2 * d + 2, b)) for d in range(K + 1)])

def M_of(c, jp, j):
    """(M_f)_{j'j} with j',j in half-integers (given as doubled indices)."""
    a = int(round(2 * jp)); b = int(round(2 * j))
    ja, jb = a / 2.0, b / 2.0
    lo, hi = abs(ja - jb), ja + jb
    n = int(round(hi - lo))
    idx = [int(round(lo)) + t for t in range(n + 1)]
    return sum(c[k] for k in idx if 0 <= k < len(c)), idx

def minor(c, j1, j2, k1, k2):
    A, ia = M_of(c, j1, k1); B, ib = M_of(c, j2, k2)
    C, ic = M_of(c, j1, k2); D, id_ = M_of(c, j2, k1)
    return A * B - C * D

b = 1.5
K = 30
c = cvals(b, K)

print("=== 2x2 minors of M_f for the actual coefficients (b=1.5), j in {0,1/2,1,3/2,2} ===")
js = [0.0, 0.5, 1.0, 1.5, 2.0]
rows = []
for i1, j1 in enumerate(js):
    for j2 in js[i1 + 1:]:
        for k1i, k1 in enumerate(js):
            for k2 in js[k1i + 1:]:
                v = minor(c, j1, j2, k1, k2)
                rows.append((j1, j2, k1, k2, v))
neg = [r for r in rows if r[4] < -1e-14]
print("  total minors = %d ; negative = %d" % (len(rows), len(neg)))
print("  smallest minors:")
for r in sorted(rows, key=lambda t: t[4])[:6]:
    print("    j1=%.1f j2=%.1f k1=%.1f k2=%.1f : D=%+.6e" % r)
print("  largest minors:")
for r in sorted(rows, key=lambda t: -t[4])[:4]:
    print("    j1=%.1f j2=%.1f k1=%.1f k2=%.1f : D=%+.6e" % r)

print()
print("=== structural probe: is D a positive combination of c_t? ===")
print("  For (j1,j2,k1,k2) = (0, 1/2, 0, 1/2) etc, list the index sets entering each entry:")
for (j1, j2, k1, k2) in [(0.0, 0.5, 0.0, 0.5), (0.0, 1.0, 0.5, 1.0), (0.0, 0.5, 0.5, 1.0),
                         (0.5, 1.0, 0.0, 0.5), (0.0, 1.0, 0.0, 1.0)]:
    A, ia = M_of(c, j1, k1); B, ib = M_of(c, j2, k2)
    C, ic = M_of(c, j1, k2); D, id_ = M_of(c, j2, k1)
    print("  (%3.1f,%3.1f ; %3.1f,%3.1f):  ac=%s  bd=%s  ad=%s  bc=%s  minor=%+.4e"
          % (j1, j2, k1, k2, ia, ib, ic, id_, A * B - C * D))
