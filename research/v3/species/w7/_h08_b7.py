import numpy as np
lam=np.array([0.95,0.80,0.65,0.40]); f2=np.array([1.0,0.30,0.20,0.10])
C=lambda k: float(np.sum(f2*lam**k))
def nodes(N):
    S=np.array([[C(i+j) for j in range(N)] for i in range(N)])
    A=np.array([[C(i+j+1) for j in range(N)] for i in range(N)])
    ev=np.sort(np.linalg.eigvals(np.linalg.solve(S,A)).real)
    return ev
print("  eigenvalues of S^{-1}A  (should be the Gauss nodes of dnu)")
for N in (2,3,4,6,10):
    print(f"   N={N:2d}  nodes={np.round(nodes(N),10)}")
print("  exact support:", np.sort(lam))
print()
print("  det(S_N - z A_N)/det(S_N) at rank N=4 :  should equal exact D_O")
S=np.array([[C(i+j) for j in range(4)] for i in range(4)])
A=np.array([[C(i+j+1) for j in range(4)] for i in range(4)])
for z in (0.3,0.8,1.0,1.03,1.5,2.5):
    print(f"   z={z:6.2f}  ratio={np.linalg.det(S-z*A)/np.linalg.det(S):20.12e}   exact={np.prod(1-z*lam):20.12e}")
