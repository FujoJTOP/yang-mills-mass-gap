import numpy as np
from scipy.special import iv

beta = 2.5
J2 = 40
t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])

def cvec(lam, n):
    return np.array([iv(k, lam) - iv(k + 2, lam) for k in range(n)])

def scvec(lam, n):
    cc = cvec(lam, n + 2)
    out = np.zeros(n)
    for k in range(n):
        if k == 0:
            out[k] = cc[1]
        else:
            out[k] = cc[k - 1] + cc[k] + cc[k + 1]
    return out

def Mmat(lam, n):
    return np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(n)] for a in range(n)])

print("=== verify explicit Bessel sums for the 2x2 block ===")
for lam in (0.5, 1.5, 3.0):
    n = J2
    cc = cvec(lam, n); sc = scvec(lam, n); tt = t[:n]
    T00 = float(np.sum(cc ** 2 * tt))
    T01 = float(np.sum(cc * sc * tt))
    T11 = float(np.sum(sc ** 2 * tt))
    Mm = Mmat(lam, n)
    T = Mm @ np.diag(tt) @ Mm
    print("  lam=%.1f" % lam)
    print("    T00: sum=%20.10f  matrix=%20.10f  diff=%.2e" % (T00, T[0, 0], abs(T00 - T[0, 0])))
    print("    T01: sum=%20.10f  matrix=%20.10f  diff=%.2e" % (T01, T[0, 1], abs(T01 - T[0, 1])))
    print("    T11: sum=%20.10f  matrix=%20.10f  diff=%.2e" % (T11, T[1, 1], abs(T11 - T[1, 1])))

print()
print("=== asymptotic peak of the profile e^{2 lam cos} sin(th) ===")
th = np.linspace(1e-5, np.pi - 1e-5, 800001)
for lam in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
    prof = np.exp(2 * lam * np.cos(th)) * np.sin(th)
    i = int(np.argmax(prof))
    print("  lam=%5.1f : argmax=%.6f   1/(2lam)=%.6f   ratio=%.4f" % (lam, th[i], 1 / (2 * lam), th[i] * 2 * lam))

print()
print("=== rank-1 quality of M_{w^2}K_op vs lambda ===")
N = 800
thg = (np.arange(N) + 0.5) * np.pi / N
TH, PH = np.meshgrid(thg, thg, indexing='ij')
Kg = (np.exp(beta * np.cos(TH - PH)) - np.exp(beta * np.cos(TH + PH))) / (2 * beta)
dth = np.pi / N
for lam in (0.5, 1.5, 3.0, 6.0):
    W = np.diag(np.exp(2 * lam * np.cos(thg)))
    A = W @ Kg * ((2 / np.pi) * dth)
    s = np.linalg.svd(A, compute_uv=False)
    s2 = s ** 2
    print("  lam=%.1f : s1^2/sum=%.8f   s2/s1=%.6f" % (lam, s2[0] / s2.sum(), s[1] / s[0]))
