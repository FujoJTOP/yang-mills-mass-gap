import numpy as np
from scipy.special import iv

# CORRECT reduction (round 21):
#   A := M_w K_op^{1/2}  on L^2([0,pi], (2/pi)dth) ;  T = A A^T
#   t_1(T) = ||A||_op^2
#   For ANY unit g:  t_2(T) <= ||(I - |g><g|) T (I - |g><g|)||_op = ||(I-|g><g|) A||_op^2
#   Choose g = u_0 = sin(theta)  (the top eigenfunction of K_op):
#        t_2(T) <= ||P_perp A||_op^2 ,   P_perp = I - |u_0><u_0|
#   => (R) follows if   ||P_perp A||^2 / ||A||^2  <=  I_2/I_1
#   At lam=0: A = K^{1/2}, LHS ratio = t_1/t_0 = I_2/I_1  (EQUALITY).
# Test this.

beta = 2.5
N = 1200
th = (np.arange(N) + 0.5) * np.pi / N
dth = np.pi / N
mu = (2 / np.pi) * dth

def K_closed(a, b, beta):
    return (np.exp(beta * np.cos(a - b)) - np.exp(beta * np.cos(a + b))) / (2 * beta)

TH, PH = np.meshgrid(th, th, indexing='ij')
Kg = K_closed(TH, PH, beta)
Kmat = Kg * mu                                   # discretized operator
Kmat = (Kmat + Kmat.T) / 2
ev, U = np.linalg.eigh(Kmat)
ev = np.maximum(ev, 0)
Khalf = U @ np.diag(np.sqrt(ev)) @ U.T           # K_op^{1/2}

u0 = np.sin(th)
u0 = u0 / np.sqrt((2 / np.pi) * dth * np.sum(u0 ** 2))

def ratios(beta, lam):
    w = np.exp(lam * np.cos(th))
    A = np.diag(w) @ Khalf
    # operator norms via the Gram matrix  A A^T = T ;  use eigenvalues of A^T A (same spectrum)
    G = A.T @ A * dth * (2 / np.pi)              # this is A^T A with the dnu measure
    evA2 = np.sort(np.linalg.eigvalsh((G + G.T) / 2))[::-1]
    P = np.eye(N) - np.outer(u0, u0) * ((2 / np.pi) * dth)
    Ap = P @ A
    Gp = Ap.T @ Ap * dth * (2 / np.pi)
    evP = np.sort(np.linalg.eigvalsh((Gp + Gp.T) / 2))[::-1]
    return evA2[0], evP[0], evA2[1] if len(evA2) > 1 else np.nan

print("=== ||P_perp A||^2 / ||A||^2  vs  I2/I1  (corrected reduction) ===")
rK = float(iv(2, beta) / iv(1, beta))
print("  I2/I1 = %.10f" % rK)
for lam in (0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0):
    t1A, tP, t2 = ratios(beta, lam)
    ratio = tP / t1A
    print("  lam=%.2f : ||A||^2=%.8f  ||P_perp A||^2=%.8f  ratio=%.10f  t_2/t_1=%.8f  %s"
          % (lam, t1A, tP, ratio, t2 / t1A, "OK" if ratio <= rK + 1e-9 else "FAIL"))

print()
print("=== sweep beta ===")
for beta in (1.0, 2.0, 4.0):
    Kgb = K_closed(TH, PH, beta)
    Kmb = (Kgb * mu); Kmb = (Kmb + Kmb.T) / 2
    evb, Ub = np.linalg.eigh(Kmb); evb = np.maximum(evb, 0)
    Khb = Ub @ np.diag(np.sqrt(evb)) @ Ub.T
    rKb = float(iv(2, beta) / iv(1, beta))
    globals()['Khalf'] = Khb
    print("  beta=%.1f : I2/I1 = %.10f" % (beta, rKb))
    for lam in (0.25, 1.0, 2.5):
        t1A, tP, t2 = ratios(beta, lam)
        print("     lam=%.2f : ratio=%.10f  %s" % (lam, tP / t1A, "OK" if tP / t1A <= rKb + 1e-9 else "FAIL"))
