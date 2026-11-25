"""mc_tensor.py -- the exact structure of the single-plaquette Hessian.

TWO DIFFERENT OBJECTS, and only one of them is the Hessian the BD criterion needs.

 (a) CHART second derivative (what w7_pl_hess.py computes):
       D2_ab = Re tr(lam_a lam_b W)          [second derivative in the exponential chart
                                              U = exp(x^a lam_a), at x = 0]
     D2 is NOT a tensor: it is not Ad-equivariant, so its eigenvalues are not conjugation
     invariants, and it does not transform by an orthogonal similarity under a change of
     orthonormal frame.

 (b) RIEMANNIAN (covariant) Hessian, the one in  Hess_G V >= lam_dot G:
       Hess_G phi(e_a,e_b) = e_a e_b phi - Gamma^c_ab(G) e_c phi
     with e_a the left-invariant field and Gamma the Koszul connection of the left-invariant
     metric G.  Hess is Ad-equivariant, so at U = 1 its spectrum depends only on the
     conjugacy class of W; and for the bi-invariant metric (G = delta) it is ISOTROPIC.

Consequence (verified below):  sup_W lam_max(Hess Re tr(UW)|_{U=1}) = 1/2  for every N,
attained iff tr W = -N (centre element).  Hence
       Hess S_W  >=  -(3 beta/N) delta     (sharp, and the tensor is a multiple of delta at
                                            the extremal configuration)
and NOT  -(6 c_N/N) beta delta  with c_N ~ 0.78 : that constant came from the chart object.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=170)
rng = np.random.default_rng(20261125)


def su_generators(N):
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = T[k, j] = 0.5; Ts.append(T)
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = -0.5j; T[k, j] = 0.5j; Ts.append(T)
    for l in range(1, N):
        T = np.zeros((N, N), dtype=complex)
        d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d
        Ts.append(T)
    return Ts


def frame(N):
    lam = [1j * T for T in su_generators(N)]
    n = len(lam)
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            C = lam[a] @ lam[b] - lam[b] @ lam[a]
            f[:, a, b] = [-2.0 * np.trace(C @ lam[c]).real for c in range(n)]
    return lam, f


def gamma_of(G, f):
    Gi = np.linalg.inv(G)
    T = np.einsum('cd,dab->abc', G, f)
    Glow = 0.5 * (np.transpose(T, (1, 2, 0)) - np.transpose(T, (1, 0, 2)) - np.transpose(T, (2, 0, 1)))
    return np.einsum('cd,dab->cab', Gi, Glow)


def haar_su(N):
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R)
    Q = Q * (d / np.abs(d)).conj()[None, :]
    return Q / np.linalg.det(Q) ** (1.0 / N)


print("=" * 100)
print("STEP A  at U=1, W = -1 (centre):  D2 vs covariant Hess, for several left-invariant G")
for N in (2, 3):
    lam, f = frame(N)
    n = len(lam)
    Wm = -np.eye(N, dtype=complex)
    D2 = np.array([[np.trace(lam[a] @ lam[b] @ Wm).real for b in range(n)] for a in range(n)])
    d = np.array([np.trace(l @ Wm).real for l in lam])
    print(f"  SU({N}):  D2 = {0.5*np.eye(n)[0]} ...  (D2 = delta/2),  d = {np.abs(d).max():.1e}")
    for nm, g in [("delta", [1.0] * n), ("stretch e1", [np.e, 1.0, 1 / np.e] + [1.0] * (n - 3))]:
        G = np.diag(np.array(g[:n]))
        Hc = D2 - np.einsum('cab,c->ab', gamma_of(G, f), d)
        print(f"     G={nm:12s}: eig(Hess Re tr(UW)) = {np.linalg.eigvalsh(Hc)}   "
              f"(isotropic: {np.ptp(np.linalg.eigvalsh(Hc)):.2e})")

print()
print("=" * 100)
print("STEP B  sup_W lam_max(Hess phi_W):  is it exactly 1/2 ?")
for N in (2, 3):
    lam, f = frame(N)
    n = len(lam)
    Gm = np.eye(n)
    best_chart, best_cov = -9, -9
    for _ in range(4000):
        W = haar_su(N)
        D2 = np.array([[np.trace(lam[a] @ lam[b] @ W).real for b in range(n)] for a in range(n)])
        d = np.array([np.trace(l @ W).real for l in lam])
        Hc = D2 - np.einsum('cab,c->ab', gamma_of(Gm, f), d)
        best_chart = max(best_chart, np.linalg.eigvalsh(D2)[-1])
        best_cov = max(best_cov, np.linalg.eigvalsh(Hc)[-1])
    print(f"  SU({N}) over 4000 Haar W:  sup lam_max(CHART D2) = {best_chart:.6f}   "
          f"sup lam_max(COVARIANT Hess) = {best_cov:.6f}   (exact value: 0.5)")

print()
print("=" * 100)
print("STEP C  is the covariant Hessian isotropic at the centre?  (SU(3), off-torus W)")
N = 3
lam, f = frame(3)
n = 8
for _ in range(4):
    W = haar_su(3)
    D2 = np.array([[np.trace(lam[a] @ lam[b] @ W).real for b in range(n)] for a in range(n)])
    d = np.array([np.trace(l @ W).real for l in lam])
    Hc = D2 - np.einsum('cab,c->ab', gamma_of(np.eye(n), f), d)
    ev = np.linalg.eigvalsh(Hc)
    print(f"  Re tr W = {np.trace(W).real:+.4f}  eig(Hess) = {np.round(ev,6)}  "
          f"spread = {np.ptp(ev):.4f}")

print()
print("=" * 100)
print("STEP D  consequence at the centre-element configuration (exact, all N)")
print("  config:  U_1(x)=(-1)^(x2+x3+x4), U_mu(x)=(-1)^x1 (mu>=2) on an even torus")
print("  => every plaque in the (1,mu) planes equals -1, so a direction-1 link sits in 6")
print("     plaques all equal to -1 and Hess S_W = -(beta/N)*6*(delta/2) = -(3 beta/N) delta.")
print("  (star) at t=0 then requires  Ric_G > (6 beta/N) delta ;  sup_G lam_min(Ric_G) = N/4")
print("  => feasible iff  6 beta/N < N/4  iff  beta < N^2/24.")
for N in (2, 3, 4, 5):
    print(f"    SU({N}): threshold beta* = N^2/24 = {N*N/24:.4f} ; deficit at beta=2.5: "
          f"{(6*2.5/N)/(N/4):.2f}x")
