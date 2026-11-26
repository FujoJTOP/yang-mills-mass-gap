"""
O09b - E[Phi''] on the pure-Haar slice for N = 2,3,4 (Haar-at-centre conditioning slice).
Haar sample: U = QR of a complex Ginibre matrix (exact Haar).  Y = i log U (Hermitian),
eigenvalues mu.  Phi = sum_{i<j} g(mu_i - mu_j);  Phi'' along v = T_3 by 5-point FD on Y + eps*i*T_3.
This is the same slice H1A-EXACT.md §2.1 uses (there S is constant on it, so the measure is Haar).
"""
import numpy as np

LOG = open(r"D:\Dev\math\w7\o09_suN.log", "w", encoding="utf-8")


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def g(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-4
    safe = np.where(small, 1.0, d)
    out = -2.0 * np.log(np.sin(safe / 2) / (safe / 2))
    return np.where(small, d ** 2 / 12 + d ** 4 / 720, out)


def phi_mu(mu):
    n = mu.shape[-1]
    i, j = np.triu_indices(n, 1)
    return g(mu[..., i] - mu[..., j]).sum(axis=-1)


def haar(N, n, rng):
    A = (rng.normal(size=(n, N, N)) + 1j * rng.normal(size=(n, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(A)
    d = np.diagonal(R, axis1=1, axis2=2)
    ph = d / np.abs(d)
    Q = Q * ph[:, None, :]
    det = np.linalg.det(Q)
    Q = Q / (det ** (1.0 / N))[:, None, None]          # force det = 1 -> SU(N)
    return Q


rng = np.random.default_rng(11)
say("=" * 96)
say("[O09b] E[Phi''] over the pure-Haar slice, v = T_3, 5-point FD, numpy only")
say("  model: U ~ Haar(SU(N));  Y = i log U (Hermitian, eigenvalues mu);")
say("         Phi = sum_{i<j} g(mu_i-mu_j), g(d) = -2 log[sin(d/2)/(d/2)] = -log J  (+const)")
say("  conditioning slice: S_W constant there (H1A-EXACT.md §2.1) => measure = Haar, beta-independent")
say("-" * 96)
h = 1e-5
for N, ns in ((2, 100000), (3, 60000), (4, 30000), (5, 12000)):
    U = haar(N, ns, rng)
    # Y = i log U : log U = V diag(i theta) V^dagger with theta = angle(eig U)
    w, Vv = np.linalg.eigh(U)
    theta = np.angle(w)                              # (-pi, pi]
    Y = np.einsum("nij,nj,nkj->nik", Vv, -theta, Vv.conj())   # Y = V diag(-theta) V^dag
    T = np.zeros((N, N), complex)
    T[0, 0], T[1, 1] = -0.5, 0.5                     # Hermitian lift of v = i*T_3 in su(N)
    F = {}
    for k in (-2, -1, 0, 1, 2):
        ev = np.linalg.eigvalsh(Y + k * h * T)
        F[k] = phi_mu(ev)
    d2 = (-F[2] + 16 * F[1] - 30 * F[0] + 16 * F[-1] - F[-2]) / (12 * h ** 2)
    d1 = (F[1] - F[-1]) / (2 * h)
    m2, m1 = d2.mean(), d1.mean()
    se = d2.std(ddof=1) / np.sqrt(ns)
    say(f"  N={N}: E[Phi''] = {m2:.5f} +- {se:.5f}   (med {np.median(d2):.4f})   "
        f"floor N/12 = {N/12:.5f}   margin = {100*(m2/(N/12)-1):.1f}%")
    say(f"        E[Phi'] = {m1:+.5f}  (0 by symmetry)   Var(Phi') = {d1.var(ddof=1):.5f}   "
        f"min Phi'' = {d2.min():.4f}")
say("-" * 96)
say("  paper (ch05-theorems.tex:208-209) reports 'measured 0.34, 0.46' for SU(2), SU(3);")
say("  the SU(2) exact value here is 0.20512 (= 1/(2 Cin(2pi)) , betastar_scaling.log:5) ->")
say("  the 0.34 is NOT the single-plaquette Haar slice; its provenance is not in the repo [E].")
LOG.close()
