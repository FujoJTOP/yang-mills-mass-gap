import numpy as np

# In the faithful structural toy (round 8), map the actual spectral ratio t1/t2 of T
# against the bare ratio t1(K)/t2(K). Goal: see whether V INCREASES or DECREASES the gap,
# and whether any uniform bound r < 1 survives as KI -> infinity.

def build(N):
    return np.array([[1 if (k >> i) & 1 else -1 for i in range(N)] for k in range(2 ** N)], float)

def analyse(N, KI, beta):
    S = build(N); n = S.shape[0]
    ov = S @ S.T
    K = np.exp(KI * ov)
    bond = np.array([np.sum([S[k, i] * S[k, (i + 1) % N] for i in range(N)]) for k in range(n)])
    f = np.exp(-(beta / 2.0) * (N - bond))
    Mf = np.diag(f)
    T = Mf @ K @ Mf
    wK = np.sort(np.linalg.eigvalsh(K))[::-1]
    wT = np.sort(np.linalg.eigvalsh(T))[::-1]
    rK = wK[1] / wK[0]
    rT = wT[1] / wT[0]
    return rK, rT

print("=== toy: does the V-weight increase or decrease the spectral ratio r = t2/t1 ? ===")
print("%4s %8s %6s %12s %12s %10s" % ("N", "KI", "beta", "r(K)", "r(T)", "rT/rK"))
for N in (4, 5):
    for KI in (0.3, 1.0, 2.5, 5.0):
        for beta in (0.5, 1.0, 2.5):
            rK, rT = analyse(N, KI, beta)
            print("%4d %8.2f %6.1f %12.6f %12.6f %10.4f" % (N, KI, beta, rK, rT, rT / rK))
    print()

print("=== large-KI limit: r -> 1 ?  (N=5) ===")
for KI in (2.5, 5.0, 10.0, 20.0):
    for beta in (2.5,):
        rK, rT = analyse(5, KI, beta)
        print("  KI=%6.2f beta=%.1f : r(K)=%.6f  r(T)=%.6f" % (KI, beta, rK, rT))
