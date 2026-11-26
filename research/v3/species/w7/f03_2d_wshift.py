"""
f03_2d_wshift.py -- D1 W-SHIFT sanity check on 2D YM (d06 discipline).

H_munu^{ab} = -D^2 delta_munu delta^{ab} + 2 f^{abc} F_munu^c   (Savvidy form).
2D: F_01 is the only component => F_B = F_01 = B (a Lie-algebra-valued scalar).

Real-space grid (no basis truncation), SU(2), background A_1 = B x0 T3.
(A) exact 2nd variation: Q = 1/2 int A_01^2 + int F_01^a[a_0,a_1]^a
    => the operator is -D^2 + c*(adF) with c fixed by the action.
(B) spectrum: PSD? negative mode? flat mode?
(C) W-SHIFT: W^1/2 (H-E0) W^-1/2 >= Delta^2 -- similarity theorem, numbers.
(D) Delta^2 of the Hessian vs the exact 2D gap (Migdal/Casimir) => D06 test (ii).
"""
import numpy as np

LOG = []
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)

B = 1.0
X, nx = 5.0, 200
x = np.linspace(-X, X, nx)
h = x[1] - x[0]
rot = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])   # ad T3

# d/dx : symmetric, Dirichlet
D1 = np.zeros((nx, nx))
D1[range(nx - 1), range(1, nx)] = 0.5 / h
D1[range(1, nx), range(nx - 1)] = -0.5 / h


def hessian(c=1.0):
    """Q = 1/2 a^T M a ,  a = (a_mu^colour)(x) ,  layout (mu, colour, x)."""
    dim = 6 * nx
    M = np.zeros((dim, dim))
    # kinetic: A_01 = d0 a_1 - B x Rot a_0    -> M_kin = T^T T
    T = np.zeros((3 * nx, dim))
    for a in range(3):
        for al in range(3):
            if rot[a, al] == 0.0:
                continue
            for i in range(nx):
                T[a * nx + i, (0 * 3 + al) * nx + i] = -B * x[i] * rot[a, al]
        for i in range(nx):
            for j in range(nx):
                if D1[i, j]:
                    T[a * nx + i, (1 * 3 + a) * nx + j] = D1[i, j]
    M += T.T @ T
    # mass: Q_mass = c * B int (a_0^1 a_1^2 - a_0^2 a_1^1)
    for i in range(nx):
        M[(0 * 3 + 0) * nx + i, (1 * 3 + 1) * nx + i] += c * B
        M[(1 * 3 + 1) * nx + i, (0 * 3 + 0) * nx + i] += c * B
        M[(0 * 3 + 1) * nx + i, (1 * 3 + 0) * nx + i] -= c * B
        M[(1 * 3 + 0) * nx + i, (0 * 3 + 1) * nx + i] -= c * B
    return M, T


out("=" * 74)
out("(A/B) 2D Savvidy Hessian  H = -D^2 + 2 adF_B ,  F_B = F_01 = B  (a scalar)")
out("=" * 74)
out("  grid: x in [-%.0f,%.0f], nx=%d, h=%.4f, Dirichlet, B=%.1f" % (X, X, nx, h, B))
M, T = hessian(c=1.0)
ev = np.linalg.eigvalsh(M)
sc = ev[-1]
out("  coefficient c from the action = 1  (Q = 1/2 int A_01^2 + int F_01[a_0,a_1])")
out("  spectrum: min = %+.3e  next = %+.6e  max = %.1f   #neg = %d"
    % (ev[0], ev[1], sc, int(np.sum(ev < -1e-9 * sc))))
M2, _ = hessian(c=2.0)
ev2 = np.linalg.eigvalsh(M2)
out("  control (c=2, the literal '2adF' of the literature): min=%+.6f  #neg=%d"
    % (ev2[0], int(np.sum(ev2 < -1e-9 * ev2[-1]))))
out("  analytic spectrum of c=1 charged sector: 2B(n+1/2) -/+ B = {0,2B,2B,2B,4B,...}")
out("  => min = -cB (infimum, NOT attained): the 2D Hessian IS indefinite.")

# ---------------------------------------------------------------- (C)
out("")
out("=" * 74)
out("(C) W-SHIFT :  W^{1/2}(H-E0)W^{-1/2} >= Delta^2  ?")
out("=" * 74)
n = M.shape[0]
rng = np.random.default_rng(5)


def randW(n, cond):
    A = rng.normal(size=(n, n)); A = (A + A.T) / 2
    w, V = np.linalg.eigh(A)
    w = (w - w.min()) / max(w.max() - w.min(), 1e-30) * np.log(cond)
    return V @ np.diag(np.exp(w / 2)) @ V.T


for name, E0 in (("E0 = 0 (background value)", 0.0),
                 ("E0 = min spec(H) (negative direction lifted by hand)", ev[0])):
    A = M - E0 * np.eye(n)
    lam = np.linalg.eigvalsh(A)[0]
    nrm = max(abs(np.linalg.eigvalsh(A)).max(), 1e-30)
    err = 0.0
    for k in range(24):
        S = randW(n, 10 ** (k % 3))
        err = max(err, abs(np.real(np.linalg.eigvals(S @ A @ np.linalg.inv(S))).min() - lam) / nrm)
    out("  %s : min spec(H-E0) = %+.6e" % (name, lam))
    out("     24 random W (cond 1..1e2): max RELATIVE shift of min spec = %.2e" % err)
    out("     => sup_W Delta^2 = min spec(H-E0) = %+.6e  [T: conjugation is isospectral]"
        % lam)

w_, V_ = np.linalg.eigh(M)
v = V_[:, 0]
norms = np.array([np.linalg.norm(v[k * nx:(k + 1) * nx]) for k in range(6)])
labs = ["a0^1", "a0^2", "a0^3", "a1^1", "a1^2", "a1^3"]
out("  most negative eigenvector, |v| per component: " +
    "  ".join("%s=%.3f" % (labs[k], norms[k]) for k in range(6)))
out("  -> channel (a_0^1, a_1^2):  Q = 1/2 (B x u)^2 - B u^2  (<0, inf = -B)")

# ---------------------------------------------------------------- (D)
out("")
out("=" * 74)
out("(D) Hessian Delta^2 vs exact 2D gap  =>  D06 acceptance (ii)")
out("=" * 74)
out("  torus area A, flux m:  B = 2 pi m / A ;  Hessian gap Delta^2 = 2B = 4 pi m / A")
out("  2D exact gap (Migdal/Casimir): m_true^2 = g^2 C2(adj)/2 = g^2 N/2 ,  g^2=1/beta")
for A_ in (16, 64, 256, 1024, 10 ** 5, 10 ** 7):
    D2 = 2 * (2 * np.pi / A_)
    m2 = 0.4 * 2 / 2.0
    out("     A=%8d  Delta^2=%.4e   Delta^2/m_true^2=%.3e" % (A_, D2, D2 / m2))
out("  Delta^2 -> 0 (like 1/A) while m_true is A-independent => acceptance (ii) fails.")

open("research/v3/species/w7/f03_2d_wshift.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")
print("OK log written")


