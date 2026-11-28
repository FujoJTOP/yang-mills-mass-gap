import numpy as np
from scipy.special import iv

# Consistent comparison: on the SAME discretisation, compute r(K) and r(T).
# The valid question is  r(T) <= r(K)  (does the spatial weight shrink the ratio?),
# NOT r(T) <= I2/I1 (which is only meaningful for the exact SU(2) single-link kernel).

rng = np.random.default_rng(11)

def run(n, beta, betap, seed=11):
    r = np.random.default_rng(seed)
    Q = r.normal(size=(n, 4)); Q /= np.linalg.norm(Q, axis=1, keepdims=True)
    G = Q @ Q.T
    K1 = np.exp(beta * G)
    K = np.kron(K1, K1)
    fpair = np.exp(-betap * (1.0 - G))
    assert fpair.max() <= 1 + 1e-12
    fvec = fpair.reshape(-1)
    Mf = np.diag(fvec)
    T = Mf @ K @ Mf
    wK = np.sort(np.linalg.eigvalsh(K))[::-1]
    wT = np.sort(np.linalg.eigvalsh(T))[::-1]
    p1 = None
    wTv, VT = np.linalg.eigh(T); iT = np.argsort(-wTv); VT = VT[:, iT]; wTv = wTv[iT]
    n1 = (fvec * VT[:, 0]) @ (fvec * VT[:, 0])
    n2 = (fvec * VT[:, 1]) @ (fvec * VT[:, 1])
    return wK[1] / wK[0], wTv[1] / wTv[0], n2 / n1, wTv[0], wTv[1]

print("=== 2-link: compare r(T) against r(K) on the SAME discretisation ===")
print("%4s %6s %7s %6s %10s %10s %10s %10s %8s %6s"
      % ("n", "beta", "beta'", "seed", "r(K)", "r(T)", "rT/rK", "n2/n1", "t1-t2", "rT<=rK"))
bad = 0; tot = 0
for n in (14, 20):
    for beta in (0.5, 1.0, 2.5):
        for betap in (0.5, 1.5, 3.0):
            for seed in (11, 23):
                rK, rT, ratio, t1, t2 = run(n, beta, betap, seed)
                tot += 1
                ok = rT <= rK + 1e-9
                if not ok: bad += 1
                print("%4d %6.1f %7.1f %6d %10.6f %10.6f %10.6f %10.6f %8.2e %6s"
                      % (n, beta, betap, seed, rK, rT, rT / rK, ratio, t1 - t2, "OK" if ok else "FAIL"))
    print()
print("  r(T)<=r(K) failures = %d / %d" % (bad, tot))
