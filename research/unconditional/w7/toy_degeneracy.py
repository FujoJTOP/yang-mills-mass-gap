import numpy as np

def build(N):
    return np.array([[1 if (k >> i) & 1 else -1 for i in range(N)] for k in range(2 ** N)], float)

def spectrum(N, KI, beta, k=6):
    S = build(N); n = S.shape[0]
    K = np.exp(KI * (S @ S.T))
    bond = np.array([np.sum([S[j, i] * S[j, (i + 1) % N] for i in range(N)]) for j in range(n)])
    f = np.exp(-(beta / 2.0) * (N - bond))
    T = np.diag(f) @ K @ np.diag(f)
    wK, VK = np.linalg.eigh(K); iK = np.argsort(-wK); wK = wK[iK]; VK = VK[:, iK]
    wT, VT = np.linalg.eigh(T); iT = np.argsort(-wT); wT = wT[iT]; VT = VT[:, iT]
    return wK[:k], wT[:k], VT, f, S

print("=== degeneracy check of the toy (N=5) ===")
for KI in (0.3, 1.0, 2.5):
    wK, wT, VT, f, S = spectrum(5, KI, 2.5, k=6)
    print("KI=%.2f  beta=2.5" % KI)
    print("   top K eigenvalues:", np.array2string(wK, precision=8))
    print("   top T eigenvalues:", np.array2string(wT, precision=8))
    print("   t1-t2 = %.3e   t2-t3 = %.3e" % (wT[0] - wT[1], wT[1] - wT[2]))
    # symmetry: is the ground state symmetric under sigma -> -sigma ?
    # index map: flipping all spins = index n-1-k
    n = S.shape[0]
    perm = np.array([n - 1 - k for k in range(n)])
    for j in range(3):
        v = VT[:, j]
        par = np.sum(v * v[perm])
        print("   phi_%d symmetry weight sum v*v[perm] = %+0.6f" % (j + 1, par))
    print()

print("=== also: is V trivial (f == 1) when beta=0 ? check r(T)=r(K) then ===")
wK, wT, VT, f, S = spectrum(5, 1.0, 0.0, k=3)
print("  beta=0 : r(K)=%.8f  r(T)=%.8f" % (wK[1] / wK[0], wT[1] / wT[0]))
