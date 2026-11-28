import numpy as np
from scipy.linalg import expm

# (a) step 1: Savvidy 点单 plaquette Hessian 的【横截面谱】
#   S_P = beta (1 - (1/2) Re tr( U1 U2 U3^-1 U4^-1 ))
#   Savvidy: q_P = 0, |a_P| = 1  =>  U_P = i sigma_3  =>  取 U1 = i sigma_3, U2=U3=U4 = I
#   Hessian: 12x12 in eps_{i,a} (i=1..4 links, a=1..3 Lie algebra)
#   gauge directions: U_i -> g U_i g^-1,  g = exp(i phi.sigma/2), phi in R^3  (3 params)
#   project H onto the orthogonal complement of the gauge span (9-dim), report the spectrum.

s1 = np.array([[0,1],[1,0]], dtype=complex)
s2 = np.array([[0,-1j],[1j,0]], dtype=complex)
s3 = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
sig = [s1, s2, s3]

def su2(v):
    """exp(i v.sigma/2) with v in R^3."""
    return expm(1j * (v[0]*s1 + v[1]*s2 + v[2]*s3) / 2.0)

U0 = [1j*s3, I2.copy(), I2.copy(), I2.copy()]

def Sp(eps, beta):
    Us = []
    for i in range(4):
        Us.append(su2(eps[3*i:3*i+3]) @ U0[i])
    P = Us[0] @ Us[1] @ np.linalg.inv(Us[2]) @ np.linalg.inv(Us[3])
    return beta * (1.0 - 0.5*np.real(np.trace(P)))

beta = 1.0
h = 1e-5
H = np.zeros((12,12))
e0 = np.zeros(12)
for a in range(12):
    for b in range(12):
        e = np.zeros(12)
        ea = e.copy(); ea[a] = h
        eb = e.copy(); eb[b] = h
        eab = e.copy(); eab[a] = h; eab[b] = h
        ema = e.copy(); ema[a] = -h
        emb = e.copy(); emb[b] = -h
        f = lambda v: Sp(v, beta)
        H[a,b] = (f(eab) - f(ea) - f(eb) + f(e0)
                  - f(ema+eb) + f(ema) + f(eb) - f(e0)) / (2*h*h) if False else \
                 (f(eab) - f(ea) - f(eb) + f(e0)) / (h*h)
H = (H + H.T) / 2

print("=== unprojected 12x12 Hessian at Savvidy (beta=1) ===")
ev = np.sort(np.linalg.eigvalsh(H))
print("  eigenvalues:", np.array2string(ev, precision=6, max_line_width=200))

# gauge tangent vectors
G = np.zeros((12,3))
hg = 1e-6
for k in range(3):
    phi = np.zeros(3); phi[k] = hg
    g = su2(phi)
    gi = np.linalg.inv(g)
    for i in range(4):
        dU = g @ U0[i] @ gi - U0[i]
        # convert dU = i (v.sigma/2) U_i  =>  v
        M = dU @ np.linalg.inv(U0[i])
        v = np.array([-1j*np.real(np.trace(M @ sig[a])) for a in range(3)])
        G[3*i:3*i+3, k] = v / hg

print()
print("=== gauge tangent vectors (3 columns), norms ===")
print("  ", np.array2string(np.linalg.norm(G, axis=0), precision=8))

# orthonormalize gauge span and project
Q, _ = np.linalg.qr(G)
P = np.eye(12) - Q @ Q.T
Hp = P @ H @ P
evp = np.sort(np.linalg.eigvalsh(Hp))
print()
print("=== projected Hessian (gauge removed): 12 eigenvalues, 3 should be ~0 ===")
print("  ", np.array2string(evp, precision=8, max_line_width=220))
w, V = np.linalg.eigh(Hp)
# restrict to the 9-dim transverse space: drop the 3 gauge zero directions
idx = np.argsort(np.abs(evp))[:3]
mask = np.ones(12, bool); mask[idx] = False
print()
print("=== transverse (9-dim) spectrum ===")
print("  ", np.array2string(evp[mask], precision=8, max_line_width=220))
print("  min = %.8f   (need >= 0 for transverse Hessian control)" % evp[mask].min())
