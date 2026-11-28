import numpy as np
from scipy.special import iv

# CORRECT indexing: j = k/2, k=0,1,2,...
#   chi_j(th) = sin((2j+1)th)/sin(th) = sin((k+1)th)/sin(th)
#   u_k(th)   = sin((k+1)th)          <-- orthonormal in L^2([0,pi], (2/pi) dth)
#   t_k       = 2 I_{k+1}(beta)/beta
# Claim (exact closed form):
#   K(th,ph) = sum_{k>=0} t_k u_k(th) u_k(ph)
#            = (1/beta) [ e^{beta cos(th-ph)} - e^{beta cos(th+ph)} ]  >= 0 on [0,pi]^2

beta = 2.5

def K_series(th, ph, KMAX=60):
    s = 0.0
    for k in range(KMAX + 1):
        s += 2 * iv(k + 1, beta) / beta * np.sin((k + 1) * th) * np.sin((k + 1) * ph)
    return s

def K_closed(th, ph):
    return (np.exp(beta * np.cos(th - ph)) - np.exp(beta * np.cos(th + ph))) / beta

print("=== verify the closed form ===")
for (a, b) in [(0.3, 0.3), (0.5, 1.2), (0.5, 2.8), (1.0, 2.0), (0.1, 3.0)]:
    print("  K(%.1f,%.1f): series=%.10f  closed=%.10f  diff=%.2e"
          % (a, b, K_series(a, b), K_closed(a, b), abs(K_series(a, b) - K_closed(a, b))))

print()
print("=== positivity of the closed form on a grid ===")
th = np.linspace(0, np.pi, 401)
TH, PH = np.meshgrid(th, th, indexing='ij')
KK = K_closed(TH, PH)
print("  min = %+.6e  max = %+.6e  negative points = %d / %d"
      % (KK.min(), KK.max(), (KK < 0).sum(), KK.size))
print("  (cos(th-ph) >= cos(th+ph) on [0,pi]^2  =>  K >= 0 by monotonicity of exp)")

print()
print("=== spectrum of T = M_w K_op M_w in the orthonormal sine basis {u_k} ===")
def T_op(J, lam, N=2000):
    # grid: theta_i, orthonormal measure (2/pi) dth
    th_i = (np.arange(N) + 0.5) * np.pi / N
    dth = np.pi / N
    w = np.exp(lam * np.cos(th_i))
    U = np.array([np.sin((k + 1) * th_i) for k in range(J + 1)])       # (J+1,N)
    # kernel acts: (K_op f)(th) = (2/pi) int K(th,ph) f(ph) dph
    # so T[j,j'] = <u_j, M_w K_op M_w u_j'>_nu
    #           = (2/pi)^2 int int u_j(th) w(th) K(th,ph) w(ph) u_j'(ph) dth dph
    THi, PHi = np.meshgrid(th_i, th_i, indexing='ij')
    KKg = K_closed(THi, PHi)
    wu = w[None, :] * U                                                # (J+1,N)
    Tm = (wu @ KKg @ wu.T) * (2 / np.pi) ** 2 * dth ** 2
    return Tm

for lam in (0.0, 0.5, 1.5, 3.0):
    J = 10
    Tm = T_op(J, lam)
    w = np.sort(np.linalg.eigvalsh((Tm + Tm.T) / 2))[::-1]
    rT = w[1] / w[0]
    rK = float(iv(2, beta) / iv(1, beta))
    print("  lam=%.1f : top3 = %s ; r(T)=%.10f  I2/I1=%.10f  ratio=%.6f  %s"
          % (lam, np.array2string(w[:3], precision=8), rT, rK, rT / rK,
             "OK" if rT <= rK + 1e-10 else "FAIL"))
