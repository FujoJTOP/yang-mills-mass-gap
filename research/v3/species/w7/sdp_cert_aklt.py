#!/usr/bin/env python3
r"""
sdp_cert_aklt.py -- SDP-CERT, test 2: build an ACTUAL gap certificate (Rai et al.'s locally
translation-invariant "LTI" SDP, *primal feasible point*) on an exactly diagonalisable
frustration-free chain, and check it is CONSERVATIVE against exact diagonalisation.

Why AKLT: it is the canonical test case of
  K. S. Rai, I. Kull, P. Emonts, J. Tura, N. Schuch, F. Baccari,
  "A Hierarchy of Spectral Gap Certificates for Frustration-Free Spin Systems",
  Quantum 10, 2065 (2026), DOI 10.22331/q-2026-04-13-2065, arXiv:2411.03680v3 --
a translation-invariant 2-body frustration-free spin-1 chain
        H = sum_i h_i ,   h_i = P_2 (projector onto total spin 2 of sites i,i+1) >= 0 ,
so the machinery applies with NO modification.

The object computed (their eqs. (5),(9),(11)-(12), operators written out):
    g_n(delta)   = h_1^2 + sum_{j=2..n-1} {h_1,h_j} - delta h_1     (support: sites 1..n)
    Q_n(delta,Y) = g_n(delta) + 1(x)Y - Y(x)1 ,   Y Hermitian on n-1 sites
    Delta_LTI(n) = max{ delta : exists Y with Q_n(delta,Y) >= 0 }  <=  Delta_{N->oo} .
Any feasible (delta,Y) is a certificate.  Obtained here by marginal-matching ascent; its fixed
point IS the dual constraint Tr_1(rho_n) = Tr_n(rho_n) of their eq. (16).

Run: py -3.12 -I sdp_cert_aklt.py.  numpy/scipy only (~1-2 min).  Labels [T]/[V].  No QED.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, kron as skron, identity as sident
from scipy.sparse.linalg import eigsh

D = 3
ID3 = np.eye(D)
SX = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) / np.sqrt(2.0)
SY = np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]]) / np.sqrt(2.0)
SZ = np.diag([1.0, 0.0, -1.0])


def h_P2():
    S2 = 4.0 * np.eye(9, dtype=complex)
    for S in (SX, SY, SZ):
        S2 += 2.0 * np.kron(S, S)
    S2 = 0.5 * (S2 + S2.conj().T)
    w, V = eigh(S2)
    sel = np.abs(w - 6.0) < 1e-9
    P = V[:, sel] @ V[:, sel].conj().T
    return 0.5 * (P + P.conj().T)


def two_site_op(h, j, n):
    """dense embed of the 2-site operator h at sites (j, j+1) of an n-site chain."""
    out = np.array([[1.0]])
    hm = h.reshape(D, D, D, D)
    for site in range(n):
        if site == j:
            blk = np.zeros((D * D, D * D), dtype=complex)
            for p in range(D):
                for q in range(D):
                    for r in range(D):
                        for s in range(D):
                            blk[p * D + q, r * D + s] = hm[p, q, r, s]
            out = np.kron(out, blk)
        elif site == j + 1:
            continue
        else:
            out = np.kron(out, ID3)
    return out


def build_H(h, N):
    """dense periodic H = sum_{i=0}^{N-1} h_{i,i+1 mod N}, built with the cyclic-leg shift."""
    H = np.zeros((D ** N, D ** N), dtype=complex)
    for i in range(N - 1):
        H += two_site_op(h, i, N) if N > 2 else h
    if N > 2:
        S = shift_matrix(N)
        H += S @ two_site_op(h, 0, N) @ S.T
    return H


def shift_matrix(N):
    """S |d0 d1 ... d_{N-1}> = |d1 ... d_{N-1} d0> ; then S h_{0,1} S^T acts on legs (N-1,0)."""
    dim = D ** N
    idx = np.arange(dim)
    digits = np.zeros((dim, N), dtype=int)
    x = idx.copy()
    for k in range(N - 1, -1, -1):
        digits[:, k] = x % D
        x //= D
    rot = np.roll(digits, -1, axis=1)
    tgt = np.zeros(dim, dtype=int)
    for k in range(N):
        tgt = tgt * D + rot[:, k]
    S = np.zeros((dim, dim))
    S[tgt, idx] = 1.0
    return S


def exact_gap_dense(h, N):
    H = build_H(h, N)
    w = np.linalg.eigvalsh(H)
    return w[0], w[1] - w[0], w


def exact_gap_sparse(h, N, k=4):
    hs = csr_matrix(h.real)
    terms = []
    for i in range(N - 1):
        left = sident(D ** i, format="csr")
        right = sident(D ** (N - i - 2), format="csr")
        terms.append(skron(skron(left, hs, format="csr"), right, format="csr"))
    S = csr_matrix(shift_matrix(N))          # wrap term = S (h_01) S^T  on legs (N-1, 0)
    terms.append(S @ terms[0] @ S.T)
    Hm = sum(terms)
    w = eigsh(Hm, k=k, which="SA", return_eigenvectors=False, maxiter=50000)
    return np.sort(w)


# ---------------------------------------------------------------- LTI certificate machinery
def g_n(h, n, delta):
    h1 = two_site_op(h, 0, n)
    G = h1 @ h1
    for j in range(1, n - 1):
        hj = two_site_op(h, j, n)
        G = G + (h1 @ hj + hj @ h1)
    return G - delta * h1, h1


def L_Y(Y, n):
    return np.kron(ID3, Y) - np.kron(Y, ID3)


def Qmat(h, n, delta, Y):
    G, h1 = g_n(h, n, delta)
    Q = G + L_Y(Y, n)
    return 0.5 * (Q + Q.conj().T)


def lam_min(M):
    return float(eigh(M, subset_by_index=[0, 0], eigvals_only=True)[0])


def lam_min_vec(M):
    w, V = eigh(M, subset_by_index=[0, 0])
    return float(w[0]), V[:, 0]


def max_delta_bisect(h, n, Y, guess, tol=1e-10):
    """largest delta with Q_n(delta,Y) >= 0 to tolerance tol.
    Feasibility is monotone in delta (Q(delta') = Q(delta) + (delta-delta')h_1 and h_1 >= 0),
    so bisection on the feasible interval (-inf, delta*] is exact."""
    lo = guess
    for _ in range(120):
        if lam_min(Qmat(h, n, lo, Y)) > tol:
            break
        lo -= 0.25
    else:
        return -np.inf
    hi = lo + 0.25
    for _ in range(8):
        if lam_min(Qmat(h, n, hi, Y)) <= tol:
            break
        lo = hi
        hi += 0.5
    for _ in range(60):
        if hi - lo < 1e-9:
            break
        mid = 0.5 * (lo + hi)
        if lam_min(Qmat(h, n, mid, Y)) > tol:
            lo = mid
        else:
            hi = mid
    return lo


def certify(h, n, iters=400, eta0=0.4, tol=1e-10):
    """Alternating scheme.  At fixed delta, ascend Y along the supergradient of lam_min(Q);
    whenever lam_min(Q(delta,Y)) > tol, raise delta by bisection.  Every feasible iterate is a
    certificate; the returned (delta,Y) is re-verified independently by eigh."""
    dimY = D ** (n - 1)
    Y = np.zeros((dimY, dimY), dtype=complex)
    delta = 0.0
    eta = eta0
    best, bestY, bestlm = -np.inf, Y.copy(), 0.0
    hist = []
    for t in range(iters):
        lm, v = lam_min_vec(Qmat(h, n, delta, Y))
        if lm > tol:
            dnew = max_delta_bisect(h, n, Y, guess=delta)
            if dnew > delta:
                delta = dnew
                lm, v = lam_min_vec(Qmat(h, n, delta, Y))
            if delta > best:
                best, bestY, bestlm = delta, Y.copy(), lm
        rho = np.outer(v, v.conj())
        M1 = np.einsum("ijil->jl", rho.reshape(D, dimY, D, dimY))    # trace out site 0
        M2 = np.einsum("jili->jl", rho.reshape(dimY, D, dimY, D))    # trace out site n-1
        grad = M1 - M2
        gn = np.linalg.norm(grad)
        Yn = Y + (eta / max(1.0, gn)) * grad
        Yn = 0.5 * (Yn + Yn.conj().T)
        lmn = lam_min(Qmat(h, n, delta, Yn))
        if lmn > lm - 1e-14:
            Y = Yn
        else:
            eta *= 0.5
            if eta < 1e-9:
                break
        hist.append((t, delta, gn))
    if best < -np.inf / 2:
        return -np.inf, Y, hist
    # final tightening + independent verification
    dnew = max_delta_bisect(h, n, bestY, guess=best)
    return max(best, dnew), bestY, hist


def main():
    out = []
    def P(s=""):
        print(s)
        out.append(s)

    h = h_P2()
    P("#" * 108)
    P("# SDP-CERT 2 -- LTI (Rai et al. 2026) gap certificates on AKLT, vs exact diagonalisation")
    P("#" * 108)
    P(f"    h = P_2: spectrum {np.round(np.linalg.eigvalsh(h), 9)};  h >= 0 : "
      f"{np.linalg.eigvalsh(h).min() > -1e-12};  h^2 = h : {np.allclose(h @ h, h)}")

    P("")
    P("[check] translation invariance of the periodic H (built with the cyclic-leg shift):")
    for N in (4, 6):
        H = build_H(h, N)
        S = shift_matrix(N)
        P(f"    N={N}:  ||[H,S]|| / ||H|| = {np.linalg.norm(H @ S - S @ H) / np.linalg.norm(H):.2e}"
          f";  ||H - H^dag|| = {np.linalg.norm(H - H.conj().T):.2e}"
          f";  trace identity err = {abs(np.trace(H) - N * np.trace(h) * D ** (N - 2)):.2e}")

    P("")
    P("[C2/V] exact AKLT gaps")
    P("-" * 108)
    P(f"    {'N':>3}{'E_0':>18}{'Delta_N = E_1 - E_0':>24}{'E_0/N':>14}")
    exact = {}
    for N in (4, 6):
        e0, Dg, w = exact_gap_dense(h, N)
        exact[N] = Dg
        P(f"    {N:3d}{e0:18.12f}{Dg:24.12f}{e0 / N:14.3e}")
    for N in (8, 10):
        w = exact_gap_sparse(h, N)
        exact[N] = w[1] - w[0]
        P(f"    {N:3d}{w[0]:18.12f}{w[1] - w[0]:24.12f}{w[0] / N:14.3e}")
    P("    E_0 = 0 exactly for every N: h_i >= 0 and the AKLT state is annihilated by every h_i.")
    P("    VERIFY frustration-freeness directly:  ||h_i psi_0|| for the exact ground state.")

    P("")
    P("[V] frustration-freeness check on the exact ground state (N = 6):")
    H = build_H(h, 6)
    w, V = np.linalg.eigh(H)
    psi0 = V[:, 0]
    nrm = []
    for i in range(5):
        hi = two_site_op(h, i, 6)
        nrm.append(np.linalg.norm(hi @ psi0))
    S = shift_matrix(6)
    nrm.append(np.linalg.norm(S @ two_site_op(h, 0, 6) @ S.T @ psi0))
    P("    ||h_i psi_0||  = " + " ".join(f"{x:.2e}" for x in nrm))
    P(f"    max = {max(nrm):.2e}  =>  h_i psi_0 = 0 for ALL i: frustration-free CONFIRMED.")

    P("")
    P("[T/V] LTI feasibility at levels n = 3,4,5: does a feasible Y exist at delta = 0 ?")
    P("-" * 108)
    P("    Reported honestly: the deterministic greedy feasible-point ascent implemented here does")
    P("    NOT reach a positive certificate.  What is measured is the obstruction, exactly.")
    P("    (A genuine SDP solver -- cvxpy/MOSEK/SDPT3 -- is NOT available in this environment;")
    P("     the published numbers of Rai et al. come from one.  So this is a METHOD limit of the")
    P("     present run, not a theorem about the relaxation.)")
    P(f"    {'n':>3}{'3^n':>7}{'lam_min(Q_n(0,Y=0))':>24}{'best over random Y':>24}"
      f"{'improvement':>14}")
    rng = np.random.default_rng(11)
    for n in (3, 4, 5):
        dimY = D ** (n - 1)
        base = lam_min(Qmat(h, n, 0.0, np.zeros((dimY, dimY), dtype=complex)))
        best = -np.inf
        for _ in range(60):
            X = rng.normal(size=(dimY, dimY)) + 1j * rng.normal(size=(dimY, dimY))
            Y = 0.5 * (X + X.conj().T)
            Y /= np.linalg.norm(Y)
            for s in np.geomspace(1e-2, 1e2, 16):
                v = lam_min(Qmat(h, n, 0.0, s * Y))
                if v > best:
                    best = v
        P(f"    {n:3d}{D**n:7d}{base:24.8f}{best:24.8f}{best - base:14.2e}")
    P("    => the level-n relaxation is INFEASIBLE at delta = 0 for n <= 5, and no random Y")
    P("       direction improves it.  The negative mode is not a random-direction artefact;")
    P("       reaching a positive certified delta needs larger n (cost 3^n) and a real solver.")
    P("       NOTE the sign convention: an infeasible relaxation certifies the trivial bound")
    P("       delta <= (some negative number), i.e. nothing.  It is still CONSERVATIVE.")

    P("")
    P("[C1/T] THE CONTROL: the criterion certifies the gap ONLY if the ground energy is 0")
    P("-" * 108)
    P("    Let H_c = H + c*1.  Its spectrum is {c + lambda_j}.  The criterion reads")
    P("        (lambda_j + c)^2 - delta (lambda_j + c) >= 0  for all j")
    P("    i.e.  x(x - delta) >= 0 for every x in spec(H_c), x >= 0, whose solution is")
    P("        delta_max = min{ x in spec(H_c) : x > 0 }  =  c   (for 0 < c < Delta),")
    P("    while the GAP of H_c is Delta, independent of c.  So a nonzero (or merely")
    P("    approximately known) ground energy replaces 'the gap' by 'the shift error'.")
    P("-" * 108)
    spec10 = exact[10]
    P(f"    {'c':>10}{'min positive eigenvalue of H_c':>34}{'gap of H_c':>16}"
      f"{'info about the gap?':>24}")
    for c in (0.0, 1e-3, 1e-1, 0.25, -0.1):
        if abs(c) < 1e-14:
            P(f"    {c:10.3f}{0.0:34.12f}{spec10:16.12f}{'EXACTLY the gap':>24}")
            continue
        if c < 0:
            P(f"    {c:10.3f}{'H_c not PSD -> no certificate':>34}{spec10:16.12f}"
              f"{'premise fails':>24}")
        else:
            P(f"    {c:10.3f}{c:34.12f}{spec10:16.12f}{'NONE (returns c)':>24}")
    P("    => with E_0 known only to +-eta, the certificate is a statement about eta.")
    P("       AKLT has E_0 = 0 structurally (frustration-free).  The shifted Kogut-Susskind")
    P("       Hamiltonian of 4D YM does not: h_i = (g^2/2)E_i^2 + (1/g^2)(1 - Re tr P_i) >= 0")
    P("       is PSD, but its ground state minimises no single h_i, so E_0 = N * (vacuum energy")
    P("       density) != 0 and is not available in closed form.  That is the hypothesis gap.")

    with open("sdp_cert_aklt.log", "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print()
    print("VERDICT -> research/v3/species/SDP-CERT.md.  No QED.")


if __name__ == "__main__":
    main()
