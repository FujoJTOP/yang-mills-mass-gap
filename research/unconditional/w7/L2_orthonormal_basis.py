import numpy as np
from scipy.special import iv

# KEY STRUCTURAL CLAIM (round 19):
#   With the measure  dnu = (2/pi) dtheta  on [0,pi],  the functions u_k = sin((2k+1)*)
#   are ORTHONORMAL:  (2/pi) int_0^pi u_j u_k dtheta = delta_{jk}.
#   Hence in L^2([0,pi], dnu):
#       K_op = sum_k t_k u_k (x) u_k        (diagonal, eigenvalues t_k)
#       M_w  = multiplication by w = exp(lam cos .)
#       T[j,j'] = <u_j, M_w K_op M_w u_{j'}>_nu
#   => T is EXACTLY the matrix of  M_w K_op M_w  in the orthonormal basis {u_k}.
#   => spec(T) = spec(M_w K_op M_w) = spec(K_op M_{w^2})  (similarity by M_w).
# Verify numerically on a fine theta grid with a truncation, and check r(T) <= I2/I1.

beta = 2.5

def build_grid(N):
    th = (np.arange(N) + 0.5) * np.pi / N
    dth = np.pi / N
    return th, dth

def u(k, th):
    return np.sin((2 * k + 1) * th)

def T_matrix(J, lam, beta, N=1200, KMAX=40):
    th, dth = build_grid(N)
    w = np.exp(lam * np.cos(th))
    U = np.array([u(k, th) for k in range(J + 1)])            # (J+1, N)
    t = np.array([2 * iv(2 * k + 1, beta) / beta for k in range(KMAX + 1)])
    # K_op action on grid: Kf(th) = sum_k t_k u_k(th) * <u_k, f>_nu ;  <f,g>_nu = (2/pi) int f g dth
    Kmat = (U.T * t[:J + 1]) @ U * (2 / np.pi) * dth         # (N,N) kernel matrix for truncation j<=J
    # T[j,j'] = (2/pi)^2 int int w w u_j K u_j'  = (2/pi)^2 * (wU) Kmat (wU)^T * dth^2
    wU = w[None, :] * U
    Tm = (wU @ Kmat @ wU.T) * (2 / np.pi) ** 2 * dth ** 2
    return Tm, np.diag(Tm), t

print("=== compare spec(T) from the matrix construction with spec(M_w K_op M_w) on the grid ===")
for lam in (0.0, 0.5, 1.5):
    for J in (8, 12):
        Tm, _, t = T_matrix(J, lam, beta)
        w = np.linalg.eigvalsh((Tm + Tm.T) / 2)[::-1]
        rT = w[1] / w[0]
        rK = float(iv(2, beta) / iv(1, beta))
        print("  lam=%.1f J=%2d : top eigenvalues = %s"
              % (lam, J, np.array2string(w[:3], precision=6)))
        print("               r(T)=%.8f   I2/I1=%.8f   ratio=%.6f  %s"
              % (rT, rK, rT / rK, "OK" if rT <= rK + 1e-12 else "FAIL"))
    print()

print("=== sanity: at lam=0 T should equal K exactly, i.e. eigenvalues = t_k ===")
Tm, _, t = T_matrix(8, 0.0, beta)
w = np.linalg.eigvalsh((Tm + Tm.T) / 2)[::-1]
print("  spec(T) at lam=0 :", np.array2string(w[:5], precision=8))
print("  t_k              :", np.array2string(t[:5], precision=8))
