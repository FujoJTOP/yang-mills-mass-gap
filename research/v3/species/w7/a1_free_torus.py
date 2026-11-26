"""A1-FREE-TORUS witness: 4D free massless field on L^4 torus, transfer-matrix lowest
nonzero eigenvalue lambda_1(a,L), its scaling, and the comparison with gamma_BD.

Convention (declared, aligned with R1-BRIDGE Sec.1.1(1b)):
    omega(k) = (2/a) sin(k a/2)      k = 2 pi /(L a)        (lattice dispersion, one direction)
    m_1      = omega(k_min)                                  (transfer-matrix gap / one-particle energy)
    lambda_1 = 2 m_1 = (4/a) sin(pi/L)                       (physical normalisation: kinetic 1/2)
    r_1      = e^{-a m_1}  = lowest nonzero eigenvalue ratio of T = e^{-aH}
    =>  lambda_1 = -(1/a) ln(r_1^2) = -(2/a) ln r_1          (NOT e^{-2 a omega})
    lambda_conf = lowest nonzero eigenvalue of -Delta_{L^4}   (config-space Poincare constant)
    lambda_slice = Poincare constant of the equal-time slice measure
Exact 4D isotropic-lattice identities checked here:
    lambda_conf = 4 sin^2(pi/L)  (a=1);  C_slice(k)^-1 -> sqrt(z(z+4)), z = 4 sin^2(pi/L)
    => lambda_slice/lambda_1 = sqrt(z+4)/2   (R1-BRIDGE Sec.0(c) closed form)
Deterministic, single process, no Monte Carlo.  py -3.12 -I a1_free_torus.py
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

N_GROUP = 3          # SU(N) label in Ric = (N/4) g  -> N/4 = 0.75 (repo normalisation g = -2tr)
G_BD = N_GROUP / 4.0  # gamma_BD at Q(inf) = 0  (Gaussian end), lattice units


def omega(a, L):
    return 2.0 * np.sin(np.pi / L) / a


def lam1(a, L):
    return 2.0 * omega(a, L)          # = (4/a) sin(pi/L)


def zfac(L):
    return 4.0 * np.sin(np.pi / L) ** 2   # spatial k-hat^2 of the lightest mode, a=1


# ---------------------------------------------------------------- [1] L scan, a=1
print("[1] lambda_1(a=1,L) = 2 m_1 = 4 sin(pi/L);  rate: L*lambda_1 -> 4 pi")
for L in [4, 8, 16, 32, 64, 128, 1024, 16384]:
    x = lam1(1.0, L)
    print(f"  L={L:6d}  lam1={x:.12f}   L*lam1={L*x:.10f}   lam1/(4pi/L)={x/(4*np.pi/L):.10f}")
Ls = np.array([4., 8., 16., 32., 64., 128., 256., 512.])
slope = np.polyfit(np.log(Ls), np.log(lam1(1.0, Ls)), 1)[0]
print(f"  fitted  d ln lam1 / d ln L = {slope:+.6f}      (analytic: -1)")

# ---------------------------------------------------------------- [2] a -> 0 at fixed L a = 8
print("\n[2] a->0 with L a = 8:  lam1 = (4/a) sin(a pi/8) -> pi/2 ; O(a^2) coeff = -pi^3/768")
for a in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125]:
    L = 8.0 / a
    x = lam1(a, L)
    print(f"  a={a:<9g} L={L:9.1f}  lam1={x:.12f}   (lam1-pi/2)/a^2 = {(x-np.pi/2)/a**2:+.9f}"
          f"   [-pi^3/768 = {-np.pi**3/768:+.9f}]")
print(f"  d ln lam1 / d ln a  at fixed L = -1 exactly (lam1 = (4/a) sin(pi/L))")

# ---------------------------------------------------------------- [3] exact 4D lattice Laplacian
def lap4(L):
    """-Delta on L^4, unit spacing, periodic. Dim L^4."""
    I1 = sparse.identity(L, format="csr")
    T = sparse.csr_matrix((np.ones(L), (np.arange(L), (np.arange(L) + 1) % L)), shape=(L, L))
    T = T + sparse.csr_matrix((np.ones(L), (np.arange(L), (np.arange(L) - 1) % L)), shape=(L, L))
    s = 2.0 * I1 - T
    K = sparse.csr_matrix((L ** 4, L ** 4))
    for mu in range(4):
        M = sparse.identity(1, format="csr")
        for nu in range(4):
            M = sparse.kron(M, s if nu == mu else I1, format="csr")
        K = K + M
    return K


print("\n[3] exact -Delta on L^4 (sparse, no MC): smallest nonzero eigenvalue vs 4 sin^2(pi/L)")
for L in [4, 6, 8]:
    K = lap4(L)
    ev = np.sort(eigsh(K, k=4, which="SM", return_eigenvectors=False))
    print(f"  L={L}:  ev[0]={ev[0]:+.3e}  ev_min_nonzero={ev[1]:.12f}   4 sin^2(pi/L)={zfac(L):.12f}"
          f"   [lambda_conf]     2 m_1 = {lam1(1.0, L):.12f}")

# ---------------------------------------------------------------- [4] slice measure, exact mode sum
print("\n[4] equal-time slice measure (4D isotropic lattice, unit spacing, L_t = L):")
print("    C(k) = (1/L) sum_n 1/(z + 2 - 2 cos(2 pi n / L));  lambda_slice = 1/C ;  C^-1 -> sqrt(z(z+4))")
for L in [4, 8, 16, 32, 64]:
    z = zfac(L)
    n = np.arange(L)
    C = np.mean(1.0 / (z + 2.0 - 2.0 * np.cos(2.0 * np.pi * n / L)))
    ls = 1.0 / C
    zz = np.sqrt(z * (z + 4.0))
    print(f"  L={L:3d}  lambda_slice={ls:.10f}   sqrt(z(z+4))={zz:.10f}   ratio={ls/zz:.6f}"
          f"   lambda_slice/lambda_1={ls/lam1(1.0, L):.6f}   sqrt(z+4)/2={np.sqrt(z+4)/2:.6f}")

# ---------------------------------------------------------------- [5] gamma_BD vs lambda_1
print("\n[5] gamma_BD = (N/4) e^{-2Q} = 0.75 (Q=0, L- and a-independent) vs lambda_1 = (4/a) sin(pi/L)")
print("    L      a      lam1        lam_conf    gamma_BD    gamma/lam1   lam1/gamma")
rows = [(4, 1.0), (8, 1.0), (16, 1.0), (17, 1.0), (32, 1.0), (64, 1.0), (10 ** 4, 1.0),
        (8, 1.0), (16, 0.5), (32, 0.25), (64, 0.125), (10 ** 6, 8e-6)]
for L, a in rows:
    x = lam1(a, L)
    conf = zfac(L) / a ** 2
    print(f"  {L:<8.0f} {a:<7g} {x:<11.6f} {conf:<11.6f} {G_BD:<11.6f} {G_BD/x:<12.6f} {x/G_BD:.6f}")

Lstar = np.pi / np.arcsin(G_BD / 4.0)      # solve (4/a) sin(pi/L) = N/4 at a=1
print(f"\n  crossing lambda_1 = gamma_BD at a=1:  L* = pi/arcsin(N/16) = {Lstar:.6f}")
for a in [1.0, 0.5, 0.25, 0.1, 0.01]:
    print(f"    a={a:<7g} L*(a) = pi/arcsin(N a/16) = {np.pi/np.arcsin(G_BD*a/4.0):.6f}"
          f"   (physical size a L* = {a*np.pi/np.arcsin(G_BD*a/4.0):.6f})")

print("\n[6] exponent pairs  (d ln . / d ln a ,  d ln . / d ln L)")
print("    lambda_1 = 2 m_1            (-1, -1)")
print("    lambda_conf                 (-2, -2)")
print("    gamma_BD  (lattice units)   ( 0,  0);   in A-units it carries a^{-2} (R1-BRIDGE Sec.2.3)")
print("\n[done]")
