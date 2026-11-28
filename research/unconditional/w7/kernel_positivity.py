import numpy as np
from scipy.special import iv

# (1) Is the transfer-matrix kernel  K(th,ph) = sum_k t_k u_k(th) u_k(ph)  POSITIVE?
#     (t_k = 2 I_{2k+1}(beta)/beta,  u_k = sin((2k+1)*))
# (2) Is the chessboard sign of M_f^{-1} robust in lambda and J2?

beta = 2.5
KMAX = 40
t = np.array([2 * iv(2 * k + 1, beta) / beta for k in range(KMAX + 1)])

def Kern(th, ph, t):
    s = 0.0
    for k in range(len(t)):
        s = s + t[k] * np.sin((2 * k + 1) * th) * np.sin((2 * k + 1) * ph)
    return s

th = np.linspace(0.001, np.pi - 0.001, 201)
TH, PH = np.meshgrid(th, th, indexing='ij')
KK = Kern(TH, PH, t)
print("=== K(th,ph) = sum_k t_k u_k(th) u_k(ph) at beta=%.1f, KMAX=%d ===" % (beta, KMAX))
print("  min = %+.6e  max = %+.6e  negative points = %d / %d"
      % (KK.min(), KK.max(), (KK < -1e-12).sum(), KK.size))

print()
print("  compare with the WKB/saddle estimate e^{beta*cos(th-ph)} etc:")
for (a, b) in [(0.5, 0.5), (1.0, 2.0), (0.5, 2.8), (2.5, 0.3)]:
    print("   K(%.1f,%.1f) = %+.6e ;  e^{b cos(a-b)} = %+.6e ; e^{b cos(a+b)} = %+.6e"
          % (a, b, Kern(a, b, t), np.exp(beta * np.cos(a - b)), np.exp(beta * np.cos(a + b))))

print()
print("=== (2) chessboard sign of M_f^{-1}: sweep lambda and J2 ===")
def Mf(J2, lam):
    return np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                     for a in range(J2 + 1)])
for lam in (0.2, 0.5, 1.0, 2.0, 4.0, 8.0):
    for J2 in (8, 14, 20):
        M = Mf(J2, lam)
        Mi = np.linalg.inv(M)
        D = np.diag([(-1.0) ** a for a in range(J2 + 1)])
        Mt = D @ Mi @ D
        print("  lam=%4.1f J2=%2d : min(D M^-1 D) = %+.3e   max = %+.3e   neg=%d"
              % (lam, J2, Mt.min(), Mt.max(), (Mt < -1e-14).sum()))
