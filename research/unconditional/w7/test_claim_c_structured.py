import numpy as np

# Faithful structural toy for the gauge setting:
#   degrees of freedom: sigma in {+-1}^N  (analogue of the spatial links)
#   "temporal" kernel K(sigma,sigma') = exp(KI * sum_i sigma_i sigma'_i)   (product over sites)
#   "spatial" weight  f(sigma)       = exp(-(beta/2) * sum_i (1 - sigma_i sigma_{i+1}))  (product over bonds)
#   T = M_f K M_f  ;  check  ||f phi_2||^2 <= ||f phi_1||^2
# This mirrors the gauge case: K couples "across time" at the same site, f is a
# spatial-plaquette weight on the same degrees of freedom.

def build(N):
    states = np.array([[1 if (k >> i) & 1 else -1 for i in range(N)] for k in range(2 ** N)], float)
    return states

def test(N, KI, beta):
    S = build(N)
    n = S.shape[0]
    # K
    ov = S @ S.T                       # sum_i sigma_i sigma'_i
    K = np.exp(KI * ov)
    # f
    bond = np.array([np.sum([S[k, i] * S[k, (i + 1) % N] for i in range(N)]) for k in range(n)])
    f = np.exp(-(beta / 2.0) * (N - bond))
    Mf = np.diag(f)
    T = Mf @ K @ Mf
    w, V = np.linalg.eigh(T)
    idx = np.argsort(-w); w = w[idx]; V = V[:, idx]
    p1, p2 = V[:, 0], V[:, 1]
    n1 = (f * p1) @ (f * p1)
    n2 = (f * p2) @ (f * p2)
    return n1, n2, w[0] / w[1] if w[1] > 0 else np.inf

print("=== structured toy: K = temporal coupling, f = spatial weight (same d.o.f.) ===")
print("%5s %7s %7s %14s %14s %10s %8s %6s" % ("N", "KI", "beta", "||f phi1||^2", "||f phi2||^2", "ratio", "t1/t2", "claim"))
bad = 0; tot = 0
for N in (3, 4, 5):
    for KI in (0.3, 1.0, 2.5):
        for beta in (0.5, 1.0, 2.5, 4.0):
            n1, n2, r = test(N, KI, beta)
            ratio = n2 / n1
            ok = ratio <= 1 + 1e-10
            tot += 1
            if not ok: bad += 1
            print("%5d %7.2f %7.2f %14.8f %14.8f %10.6f %8.4f %6s"
                  % (N, KI, beta, n1, n2, ratio, r, "OK" if ok else "FAIL"))
print()
print("  failures = %d / %d" % (bad, tot))
