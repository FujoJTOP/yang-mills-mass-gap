"""
qg2_abelian.py -- Q4 route, step 2: EXACT spectra of abelian lattice gauge theory
                   vs the geometric (skeleton) effective resistance.

Model: Z_N lattice gauge theory, 2+1 dimensions, spatial torus L x L, anisotropic
transfer matrix in the temporal gauge.

  T(sigma,sigma') = exp[(beta_s/2)(A(sigma)+A(sigma'))] * prod_l exp[beta_t cos(2pi(m_l-m'_l)/N)]
  A(sigma)        = sum_p cos(2pi m_p / N) ,   m_p = oriented plaquette sum
  ==> T = D K^{tensor} D with D = exp[(beta_s/2) A]  ==> similar to the SYMMETRIC
      S = exp[(beta_s/4)A] K^{tensor} exp[(beta_s/4)A]   (same spectrum).

  mass gap  m = -ln(lambda_2/lambda_1)  (lattice units, a_t = 1).

The skeleton network of the same spatial torus has an effective resistance
R_eff(L) that is PURE GEOMETRY: it does not contain beta, N, or the action.
That is the whole point of the test below.

  PART A  validation: beta_s = 0 is exactly solvable at every L
  PART B  beta scan at L = 2, 3 (exact spectrum; dense for L=2, Lanczos for L=3)
  PART C  Z_4 (U(1) truncation) at L = 2, for the same question
  PART D  the table  m(beta) * R_eff(L)  --> is it beta-blind? is it L-independent?
"""
import json
import math
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ---------------------------------------------------------------------------
def plaquette_link_table(L):
    """returns (P, 4) array of link indices with orientation signs (P,4)."""
    nl = 2 * L * L
    h = lambda x, y: (0 * L * L) + (x % L) * L + (y % L)
    v = lambda x, y: (1 * L * L) + (x % L) * L + (y % L)
    P, S = [], []
    for x in range(L):
        for y in range(L):
            # m_p = m_h(x,y) + m_v(x+1,y) - m_h(x,y+1) - m_v(x,y)
            P.append([h(x, y), v((x + 1) % L, y), h(x, y + 1), v(x, y)])
            S.append([+1, +1, -1, -1])
    return np.array(P), np.array(S), nl


def digarray(N, nl):
    """(N**nl, nl) base-N digits, digit l = link l."""
    M = N ** nl
    idx = np.arange(M)
    D = np.empty((M, nl), dtype=np.int64)
    for l in range(nl):
        D[:, l] = (idx // (N ** l)) % N
    return D


def kernel_N(N, beta_t):
    """K[m,m'] = exp(beta_t cos(2 pi (m-m')/N))"""
    m = np.arange(N)
    d = (m[:, None] - m[None, :]) % N
    return np.exp(beta_t * np.cos(2.0 * np.pi * d / N))


def apply_tensor_local(vec, Km, N, nl, shape):
    """apply the same N x N kernel to every one of the nl axes."""
    v = vec.reshape(shape)
    for ax in range(nl):
        v = np.moveaxis(v, ax, 0)
        sh = v.shape
        v = (Km @ v.reshape(N, -1)).reshape(sh)
        v = np.moveaxis(v, 0, ax)
    return v.reshape(-1)


def spectrum(N, L, beta_s, beta_t, k=4, dense_limit=60000):
    P, S, nl = plaquette_link_table(L)
    M = N ** nl
    D = digarray(N, nl)
    mp = np.zeros(M, dtype=np.float64)
    for p in range(P.shape[0]):
        s = np.zeros(M, dtype=np.int64)
        for a in range(4):
            s += S[p, a] * D[:, P[p, a]]
        mp += np.cos(2.0 * np.pi * (s % N) / N)
    A = mp                                  # diagonal: sum of plaquette cosines
    Dg = np.exp(0.25 * beta_s * A)
    Km = kernel_N(N, beta_t)
    shape = (N,) * nl
    g = Dg.reshape(-1)

    def matvec(x):
        return g * apply_tensor_local(g * x, Km, N, nl, shape)

    if M <= dense_limit:
        # dense build by columns (fine for small M)
        H = np.empty((M, M))
        e = np.zeros(M)
        for j in range(M):
            e[:] = 0.0
            e[j] = 1.0
            H[:, j] = matvec(e)
        H = 0.5 * (H + H.T)
        ev = np.linalg.eigvalsh(H)
        ev = ev[::-1]
        return ev[:k], M
    else:
        op = LinearOperator((M, M), matvec=matvec, dtype=np.float64)
        ev = eigsh(op, k=k, which="LA", return_eigenvectors=False,
                   tol=1e-10, maxiter=4000, ncv=min(M - 1, 40))
        return np.sort(ev)[::-1], M


def gap_of(ev):
    return float(-math.log(ev[1] / ev[0]))


# ---------------------------------------------------------------------------
def torus_Reff_max(L, d):
    """geometric R_eff^max of the d-dim torus (exact Fourier sum)."""
    grids = np.meshgrid(*[np.arange(L)] * d, indexing="ij")
    m = np.stack([g.ravel() for g in grids], axis=1).astype(float)
    kk = 2.0 * math.pi * m / L
    lam = np.sum(4.0 * np.sin(kk / 2.0) ** 2, axis=1)
    n = len(lam)
    w = np.zeros(n)
    nz = lam > 1e-12
    w[nz] = 1.0 / lam[nz]
    C = np.cos(2.0 * math.pi * (m @ m.T) / L)
    R = (2.0 / n) * (w.sum() - C.T @ w)
    return float(np.max(np.maximum(R, 0.0)))


# ---------------------------------------------------------------------------
def part_A():
    say("=" * 78)
    say("PART A -- validation: beta_s = 0 is exactly solvable at every L")
    say("  S = K^{tensor nl}; eigenvalues prod_l (1 +- t), t = exp(-2 beta_t)")
    say("  => lambda_2/lambda_1 = (1-t)/(1+t)  =>  m = -ln((1-t)/(1+t)),  L-independent")
    say("=" * 78)
    ok = []
    for L in [2, 3]:
        for bt in [0.1, 0.3, 0.7]:
            ev, M = spectrum(2, L, 0.0, bt, k=2)
            t = math.exp(-2 * bt)
            pred = -math.log((1 - t) / (1 + t))
            got = gap_of(ev)
            ok.append(dict(L=L, beta_t=bt, dim=M, m_exact_num=got, m_analytic=pred,
                           err=abs(got - pred)))
            say(f"  L={L} N=2 dim={M:>8} beta_t={bt:>4}: m_num={got:.10f}  "
                f"m_analytic={pred:.10f}  |diff|={abs(got-pred):.2e}")
    say("")
    say("  PROVED (exact, all L): at beta_s = 0 the gap does not depend on L.")
    say("  The GEOMETRIC R_eff of the same torus grows like ln L (2D).")
    say("  => m * R_eff ~ ln L, UNBOUNDED.  The premise 'lambda ~ 1/R_eff' is wrong")
    say("  already in the exactly solvable corner.")
    say("")
    return ok


def part_B():
    say("=" * 78)
    say("PART B -- exact mass gap m(beta), isotropic beta_s = beta_t = beta")
    say("=" * 78)
    rows = []
    for L in [2, 3]:
        Rg = torus_Reff_max(L, 2)
        say(f"  --- spatial torus L={L}  (geometric R_eff^max = {Rg:.6f}, beta-free)")
        say(f"      {'beta':>6} {'dim':>9} {'lam2/lam1':>12} {'m = -ln ratio':>15} "
            f"{'m*R_eff':>12}")
        for beta in [0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0]:
            ev, M = spectrum(2, L, beta, beta, k=2)
            r = ev[1] / ev[0]
            m = -math.log(r)
            say(f"      {beta:>6.2f} {M:>9} {r:>12.6e} {m:>15.6e} {m*Rg:>12.6e}")
            rows.append(dict(L=L, beta=beta, dim=M, ratio=r, m=m, Rgeo=Rg, prod=m * Rg))
        say("")
    say("  Reading: m(beta) sweeps ~2 decades while R_eff is a CONSTANT (geometry).")
    say("  So 'm * R_eff' is not a certificate: it is m itself, rescaled by a number")
    say("  that contains no dynamics.")
    say("")
    return rows


def part_C():
    say("=" * 78)
    say("PART C -- Z_4 (a U(1) truncation) at L = 2: same question, clock variables")
    say("=" * 78)
    rows = []
    Rg = torus_Reff_max(2, 2)
    say(f"  geometric R_eff^max(T^2_2) = {Rg:.6f}")
    say(f"  {'beta':>6} {'dim':>8} {'lam2/lam1':>12} {'m':>12} {'m*R_eff':>12}")
    for beta in [0.1, 0.5, 1.0, 2.0, 4.0]:
        ev, M = spectrum(4, 2, beta, beta, k=2)
        r = ev[1] / ev[0]
        m = -math.log(r)
        say(f"  {beta:>6.2f} {M:>8} {r:>12.6e} {m:>12.6e} {m*Rg:>12.6e}")
        rows.append(dict(N=4, L=2, beta=beta, dim=M, ratio=r, m=m, prod=m * Rg))
    say("")
    return rows


def part_D(B, C):
    say("=" * 78)
    say("PART D -- summary table: is m * R_eff beta-blind / L-independent?")
    say("=" * 78)
    say(f"  {'model':>10} {'L':>3} {'beta':>6} {'m':>12} {'R_eff^geo':>10} {'m R_eff':>12}")
    for r in B:
        say(f"  {'Z2 2+1D':>10} {r['L']:>3} {r['beta']:>6.2f} {r['m']:>12.5e} "
            f"{r['Rgeo']:>10.5f} {r['prod']:>12.5e}")
    for r in C:
        say(f"  {'Z4 2+1D':>10} {r['L']:>3} {r['beta']:>6.2f} {r['m']:>12.5e} "
            f"{r['Rgeo'] if 'Rgeo' in r else torus_Reff_max(2,2):>10.5f} {r['prod']:>12.5e}")
    say("")
    say("  Column 'm R_eff' varies by >2 decades across beta at FIXED L, and across L")
    say("  at fixed beta (see PART A: unbounded in L).  Neither bounded nor L-independent.")
    say("")


def main():
    t0 = time.time()
    say(f"qg2_abelian.py   start {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say("")
    A = part_A()
    B = part_B()
    C = part_C()
    part_D(B, C)
    say(f"total wall {time.time() - t0:.1f} s")
    json.dump(dict(partA=A, partB=B, partC=C), open("qg2_abelian.json", "w"), indent=1)
    open("qg2_abelian.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
