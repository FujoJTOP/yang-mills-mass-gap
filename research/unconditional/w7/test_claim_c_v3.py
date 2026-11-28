import numpy as np

# FIX: keep f <= 1 (as in the gauge setting) while breaking the Z_2 symmetry.
#   break symmetry in K only:
#       K(s,s') = exp(KI * sum_i s_i s'_i + (h/2) * sum_i (s_i + s'_i))
#   f(s) = exp(-(beta/2) * sum_i (1 - s_i s_{i+1}))   <= 1   ALWAYS
# T = M_f K M_f is then not Z_2-invariant, so phi_1 and phi_2 share a sector,
# and f <= 1 holds as in T = M_{V^{1/2}} K M_{V^{1/2}} with V = e^{-beta S_3D} <= 1.

def build(N):
    return np.array([[1 if (k >> i) & 1 else -1 for i in range(N)] for k in range(2 ** N)], float)

def analyse(N, KI, beta, h):
    S = build(N); n = S.shape[0]
    mag = S.sum(axis=1)
    K = np.exp(KI * (S @ S.T) + (h / 2.0) * (mag[:, None] + mag[None, :]))
    bond = np.array([np.sum([S[j, i] * S[j, (i + 1) % N] for i in range(N)]) for j in range(n)])
    f = np.exp(-(beta / 2.0) * (N - bond))
    assert f.max() <= 1 + 1e-12, "f must be <= 1"
    T = np.diag(f) @ K @ np.diag(f)
    wK = np.sort(np.linalg.eigvalsh(K))[::-1]
    wT, VT = np.linalg.eigh(T)
    iT = np.argsort(-wT); wT = wT[iT]; VT = VT[:, iT]
    p1, p2 = VT[:, 0], VT[:, 1]
    n1 = (f * p1) @ (f * p1); n2 = (f * p2) @ (f * p2)
    return n1, n2, wK[1] / wK[0], wT[1] / wT[0], wT[0], wT[1]

print("=== corrected test: f <= 1 enforced, symmetry broken in K only ===")
print("%3s %6s %6s %5s %12s %12s %10s %9s %9s %6s"
      % ("N", "KI", "beta", "h", "||f p1||^2", "||f p2||^2", "n2/n1", "r(K)", "r(T)", "(c)"))
fails = 0; tot = 0
for N in (4, 5):
    for KI in (0.3, 1.0, 2.5, 5.0):
        for beta in (0.5, 1.0, 2.5):
            h = 0.3
            n1, n2, rK, rT, t1, t2 = analyse(N, KI, beta, h)
            ratio = n2 / n1
            tot += 1
            ok = ratio <= 1 + 1e-10
            if not ok: fails += 1
            print("%3d %6.2f %6.1f %5.1f %12.6f %12.6f %10.6f %9.6f %9.6f %6s"
                  % (N, KI, beta, h, n1, n2, ratio, rK, rT, "OK" if ok else "FAIL"))
print()
print("  (c)-failures = %d / %d" % (fails, tot))
print()
print("=== sanity: at beta=0, f==1 so T==K and r(T)=r(K) exactly; (c) becomes equality ===")
for N in (4, 5):
    n1, n2, rK, rT, t1, t2 = analyse(N, 1.0, 0.0, 0.3)
    print("  N=%d beta=0 : n1=%.10f n2=%.10f  rK=%.8f rT=%.8f" % (N, n1, n2, rK, rT))
