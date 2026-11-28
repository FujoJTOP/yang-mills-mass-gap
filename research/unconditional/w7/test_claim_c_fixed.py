import numpy as np

# CORRECTED toy: break the sigma -> -sigma symmetry (which caused a spurious
# cross-parity comparison in the round-8 test).
#   K(s,s') = exp(KI * sum_i s_i s'_i + h * sum_i (s_i + s'_i)/2)
#   f(s)    = exp(-(beta/2) * sum_i (1 - s_i s_{i+1}) - h * sum_i s_i)
# With h != 0 both the ground and first excited states lie in the same sector.

def build(N):
    return np.array([[1 if (k >> i) & 1 else -1 for i in range(N)] for k in range(2 ** N)], float)

def analyse(N, KI, beta, h):
    S = build(N); n = S.shape[0]
    ov = S @ S.T
    mag = S.sum(axis=1)
    K = np.exp(KI * ov + h * (mag[:, None] + mag[None, :]) / 2.0)
    bond = np.array([np.sum([S[j, i] * S[j, (i + 1) % N] for i in range(N)]) for j in range(n)])
    f = np.exp(-(beta / 2.0) * (N - bond) - h * mag)
    T = np.diag(f) @ K @ np.diag(f)
    wK = np.sort(np.linalg.eigvalsh(K))[::-1]
    wT, VT = np.linalg.eigh(T)
    iT = np.argsort(-wT); wT = wT[iT]; VT = VT[:, iT]
    p1, p2 = VT[:, 0], VT[:, 1]
    n1 = (f * p1) @ (f * p1); n2 = (f * p2) @ (f * p2)
    return n1, n2, wK[1] / wK[0], wT[1] / wT[0]

print("=== CORRECTED toy (symmetry broken by h) ===")
print("%3s %6s %6s %5s %13s %13s %9s %9s %9s %6s"
      % ("N", "KI", "beta", "h", "||f p1||^2", "||f p2||^2", "n2/n1", "r(K)", "r(T)", "(c)"))
fails = 0; tot = 0
for N in (4, 5):
    for KI in (0.3, 1.0, 2.5):
        for beta in (0.5, 1.0, 2.5):
            for h in (0.4,):
                n1, n2, rK, rT = analyse(N, KI, beta, h)
                ratio = n2 / n1
                tot += 1
                ok = ratio <= 1 + 1e-10
                if not ok: fails += 1
                print("%3d %6.2f %6.1f %5.1f %13.8f %13.8f %9.6f %9.6f %9.6f %6s"
                      % (N, KI, beta, h, n1, n2, ratio, rK, rT, "OK" if ok else "FAIL"))
print()
print("  (c)-failures = %d / %d" % (fails, tot))
print()
print("=== with symmetry, confirm the spurious degeneracy is gone ===")
for N in (5,):
    for KI in (1.0, 2.5):
        n1, n2, rK, rT = analyse(N, KI, 2.5, 0.4)
        print("  N=%d KI=%.2f : r(K)=%.6f  r(T)=%.6f  (1-r(T) = %.3e)" % (N, KI, rK, rT, 1 - rT))
