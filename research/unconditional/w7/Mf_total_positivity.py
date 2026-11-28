import numpy as np
from scipy.special import iv

# In the character basis {chi_j}, j = 0, 1/2, 1, ...
#   (M_f)_{j' j} = int chi_j' f chi_j dU = sum_{j''=|j-j'|}^{j+j'} c_{j''}
#   with c_{j''} = e^{-b}(I_{2j''}(b) - I_{2j''+2}(b)) > 0.
# Test whether M_f is totally nonnegative (all contiguous 2x2 minors >= 0),
# which would give variation-diminishing / PF structure.

def ccoef(j2, b):
    j = j2 / 2.0
    return np.exp(-b) * (iv(2 * j, b) - iv(2 * j + 2, b))

def Mf(b, J):
    # indices j' = a/2, j = c/2, a,c in 0..J ; j'' = d/2 in max(|a-c|, ...) up to a+c but with |a-c|<=d<=a+c and parity
    M = np.zeros((J + 1, J + 1))
    for a in range(J + 1):
        for c in range(J + 1):
            s = 0.0
            for d in range(abs(a - c), a + c + 1):
                if (d - abs(a - c)) % 2 == 0:
                    s += ccoef(d, b)
            M[a, c] = s
    return M

print("=== is M_f entrywise positive and totally nonnegative (2x2 minors)? ===")
for b in (0.5, 1.5, 3.0):
    J = 8
    M = Mf(b, J)
    neg = (M < 0).sum()
    worst = 0.0
    for i in range(J):
        for j in range(i + 1, J + 1):
            for k in range(J):
                for l in range(k + 1, J + 1):
                    dt = M[i, k] * M[j, l] - M[i, l] * M[j, k]
                    worst = min(worst, dt)
    print("  b=%.1f : entries<0 : %d ; min contiguous 2x2 minor = %+.3e ; max entry = %.6f"
          % (b, neg, worst, M.max()))
    print("        M row0:", np.array2string(M[0], precision=6))
    print("        M row1:", np.array2string(M[1], precision=6))

print()
print("=== spectrum of T = M_f K M_f in the character basis; check r(T) <= I2/I1 ===")
def tvec(beta, J):
    return np.array([2 * iv(2 * (k / 2.0) + 1, beta) / beta for k in range(J + 1)])

bad = 0; tot = 0
for beta in (0.5, 1.0, 2.5, 4.0):
    for b in (0.5, 1.5, 3.0):
        J = 12
        M = Mf(b, J)
        K = np.diag(tvec(beta, J))
        T = M @ K @ M
        w = np.sort(np.linalg.eigvalsh(T))[::-1]
        rT = w[1] / w[0]
        rK = float(iv(2, beta) / iv(1, beta))
        tot += 1
        ok = rT <= rK + 1e-12
        if not ok: bad += 1
        print("  beta=%.1f b=%.1f : r(T)=%.8f  I2/I1=%.8f  ratio=%.6f  %s"
              % (beta, b, rT, rK, rT / rK, "OK" if ok else "FAIL"))
print("  failures = %d / %d" % (bad, tot))
