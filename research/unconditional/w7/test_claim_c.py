import numpy as np

# Claim (c): for T = M_f K M_f with f <= 1 and K >= 0 whose ground vector is CONSTANT,
#            ||f phi_2|| <= ||f phi_1||   (phi_i orthonormal eigenvectors of T)
# Reformulation: ||f phi||^2 = <phi, f^2 phi>.
# Test on random K with constant ground state and random f in (0,1].

rng = np.random.default_rng(2026)

def random_K_with_constant_ground(n, rng):
    # build a symmetric positive matrix K with K*1 = t0*1  (constant is an eigenvector)
    # construct via conjugation: K = P diag(t) P^T with P orthogonal and P[:,0] = 1/sqrt(n)
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    if Q[0, 0] < 0:
        Q[:, 0] *= -1
    P = Q.copy()
    P[:, 0] = 1.0 / np.sqrt(n)
    Q2, _ = np.linalg.qr(P)
    P = Q2.copy()
    P[:, 0] = 1.0 / np.sqrt(n)
    t = np.sort(rng.uniform(0.2, 1.0, size=n))[::-1]
    K = P @ np.diag(t) @ P.T
    K = (K + K.T) / 2
    return K

print("=== test claim (c) on random K (constant ground state) + random f in (0,1] ===")
fails = 0
trials = 3000
worst = -1e9
for _ in range(trials):
    n = 5
    K = random_K_with_constant_ground(n, rng)
    f = rng.uniform(0.05, 1.0, size=n)
    Mf = np.diag(f)
    T = Mf @ K @ Mf
    w, V = np.linalg.eigh(T)
    # sort descending
    idx = np.argsort(-w)
    w = w[idx]; V = V[:, idx]
    p1 = V[:, 0]; p2 = V[:, 1]
    n1 = (f * p1) @ (f * p1)
    n2 = (f * p2) @ (f * p2)
    ratio = n2 / n1 if n1 > 0 else np.inf
    worst = max(worst, ratio)
    if ratio > 1 + 1e-10:
        fails += 1

print("  trials = %d ; failures (||f phi_2||^2 > ||f phi_1||^2) = %d" % (trials, fails))
print("  worst ratio = %.6f" % worst)
print()
print("  => if failures > 0, claim (c) is FALSE for generic (K,f);")
print("     it would then need the specific gauge structure (K = temporal convolution,")
print("     f = spatial-plaquette weight) to hold.")
