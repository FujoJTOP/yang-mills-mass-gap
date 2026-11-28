import numpy as np
from scipy.special import iv

# Faithful 2-link SU(2) test.
#   Represent U in SU(2) by its quaternion q in S^3 (unit 4-vector).
#   Then  (1/2) tr(U V^{-1}) = <q_U, q_V>_{R^4}.
#   "temporal" kernel  K = K1 (x) K2 ,  K(u,v) = exp(beta * <q_u,q_v>)      (constant ground vector)
#   "spatial"  weight  f(u1,u2) = exp(-beta' * (1 - <q1,q2>))  <= 1         (couples the two links)
#   T = M_f K M_f  on L^2(SU(2)^2)
# Test claim (c) and r(T) <= I2/I1.

rng = np.random.default_rng(11)

def sobol_s3(n):
    # deterministic low-discrepancy points on S^3 via normalising a scrambled Halton-like set
    # (simple: use a fixed random seed for reproducibility; adequate for a spectral check)
    g = rng.normal(size=(n, 4))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    return g

def build(n):
    Q = sobol_s3(n)
    G = Q @ Q.T                      # <q_i, q_j>
    return Q, G

def run(n, beta, betap):
    Q, G = build(n)
    K1 = np.exp(beta * G)            # single-link kernel
    K = np.kron(K1, K1)              # tensor product over the two links
    # f(u1,u2) = exp(-betap*(1 - <q_u1, q_u2>))
    fpair = np.exp(-betap * (1.0 - G))      # [i,j] for link1=i, link2=j
    fvec = fpair.reshape(-1)                # matches kron ordering (i outer, j inner)
    assert fvec.max() <= 1 + 1e-12
    Mf = np.diag(fvec)
    T = Mf @ K @ Mf
    w, V = np.linalg.eigh(T)
    idx = np.argsort(-w); w = w[idx]; V = V[:, idx]
    p1, p2 = V[:, 0], V[:, 1]
    n1 = (fvec * p1) @ (fvec * p1)
    n2 = (fvec * p2) @ (fvec * p2)
    return n1, n2, w[0], w[1]

print("=== 2-link SU(2) (quaternion grid, exact <q1,q2> structure) ===")
print("%4s %6s %7s %14s %14s %10s %10s %10s %6s"
      % ("n", "beta", "beta'", "||f p1||^2", "||f p2||^2", "n2/n1", "I2/I1", "r(T)", "(c)"))
bad = 0; tot = 0
for n in (12, 18, 24):
    for beta in (0.5, 1.0, 2.5):
        for betap in (0.5, 1.5, 3.0):
            n1, n2, t1, t2 = run(n, beta, betap)
            rK = float(iv(2, beta) / iv(1, beta))
            rT = t2 / t1
            ratio = n2 / n1
            ok = ratio <= 1 + 1e-9
            tot += 1
            if not ok: bad += 1
            flag = "OK" if ok else "FAIL"
            print("%4d %6.1f %7.1f %14.8f %14.8f %10.6f %10.6f %10.6f %6s"
                  % (n, beta, betap, n1, n2, ratio, rK, rT, flag))
    print()
print("  (c)-failures = %d / %d" % (bad, tot))
