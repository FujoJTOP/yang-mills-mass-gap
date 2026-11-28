import numpy as np
from scipy.special import iv

# (1) verify the explicit Bessel sums for the 2x2 block
#     c_k = I_k(lam) - I_{k+2}(lam);  sc_k = sum_{t=|k-1|}^{k+1} c_t   (with sc_0 = c_1)
#     T00 = sum c_k^2 t_k ;  T01 = sum c_k sc_k t_k ;  T11 = sum sc_k^2 t_k ;  t_k = 2 I_{k+1}(beta)/beta
# (2) asymptotic profile: M_{w^2} K_op is approximately rank-1 with profile w^2(theta) sin(theta),
#     peaked at theta* satisfying -2 lam sin(th) + cot(th) = 0, i.e. th* ~ 1/(2 lam).

beta = 2.5
J2 = 40
t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
c = np.array([iv(k, 1.5) - iv(k + 2, 1.5) for k in range(J2 + 1)])

def sc(k, lam, J2=J2):
    cc = np.array([iv(j, lam) - iv(j + 2, lam) for j in range(J2 + 1)])
    if k == 0:
        return cc[1]
    return cc[max(k - 1, 0)] + cc[k] + (cc[k + 1] if k + 1 < len(cc) else 0.0)

for lam in (0.5, 1.5, 3.0):
    cc = np.array([iv(k, lam) - iv(k + 2, lam) for k in range(J2 + 1)])
    scv = np.array([sc(k, lam, ) for k in range(J2)])
    T00 = np.sum(cc[:J2] ** 2 * t[:J2])
    T01 = np.sum(cc[:J2] * scv * t[:J2])
    T11 = np.sum(scv ** 2 * t[:J2])
    # build the actual matrix
    def Mf(a, b, lam):
        return iv(abs(a - b), lam) - iv(a + b + 2, lam)
    Mm = np.array([[Mf(a, b, lam) for b in range(J2 + 1)] for a in range(J2 + 1)])
    T = Mm @ np.diag(t) @ Mm
    print("lam=%.1f : T00 sum=%.10f vs matrix=%.10f" % (lam, T00, T[0, 0]))
    print("           T01 sum=%.10f vs matrix=%.10f" % (lam, T01, T[0, 1]))
    print("           T11 sum=%.10f vs matrix=%.10f" % (lam, T11, T[1, 1]))

print()
print("=== asymptotic peak of the profile w^2 sin(theta) = e^{2lam cos} sin(th) ===")
th = np.linspace(1e-4, np.pi - 1e-4, 400001)
for lam in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
    prof = np.exp(2 * lam * np.cos(th)) * np.sin(th)
    i = np.argmax(prof)
    pred = 1.0 / (2 * lam)
    print("  lam=%5.1f : argmax th = %.6f   1/(2lam) = %.6f   ratio = %.4f"
          % (lam, th[i], pred, th[i] / pred))

print()
print("=== rank-1 quality of M_{w^2}K_op as a function of lambda ===")
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
    print("  lam=%.1f : s1^2/sum = %.8f  ; second/first = %.6f" % (lam, s2[0] / s2.sum(), s[1] / s[0]))

