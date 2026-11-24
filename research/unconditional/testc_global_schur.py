#!/usr/bin/env python
"""
Test C -- FEASIBLE version of the decisive SRF-2.1 falsification probe.

Object (SRF 2.1): block-average decomposition of the quadratic (Hessian) form
    K_L = [[ A, B^T ],[ B, D ]]      V = coarse (block-constant),  eta = fiber
Schur / Feshbach complement
    S_L = A - B^T D^{-1} B ,         lambda_soft(L) = inf spec( S_L | coarse modes
                                                       orthogonal to the vacuum ).
Question: as L -> inf, does lambda_soft(L) decay like L^{-2} (no floor, the
free/kinetic behaviour -> SRF mechanism dead), or approach a positive value
(a mass floor, i.e. mass generation)?

This is NOT a true global transfer-operator Schur complement (infeasible at
L=12..96).  It is the feasible, honest version:

  free : FREE SCALAR calibration (exact, cheap).  Block average, form S_L,
         verify lambda_soft ~ L^{-2} for massless and -> m^2 for massive.
         This is the correctness gate: the machinery must reproduce the
         known free-field answer, non-vacuously.
  u1   : FREE U(1) MAXWELL calibration (3D, gauge theory): project out the
         F=0 kernel (pure gauge + torons), block-average the transverse photon,
         form S_L, verify the massless photon lambda_soft ~ L^{-2} and a
         (Proca-)massive photon -> m^2.  This validates the *gauge* machinery
         (kernel removal + transverse projection + block Schur) that SU(2) uses.
  su2  : SU(2) static-action Hessian proxy at L=4 (thermalised, beta=2.2-2.4).
         The faithful static spatial Hessian H=d^2 S, gauge-transverse H_t,
         block-average, global Schur complement, lambda_soft over gauge-
         invariant non-uniform coarse channels.  *** PROXY ONLY ***: the static
         spatial Hessian is NOT the time-transfer operator; SRF 2.1 itself
         warns static Hessian != physical time gap.  It can nevertheless tell
         whether the softest coarse Schur eigenvalue shows an L^{-2}-type
         kinetic tail (like free) or a positive floor.

Memory discipline (charter section 4): ONE job at a time; the SU(2) part is
dense but only at L=4 (H: 3072^2 = 75 MB, the rest < ~150 MB).  L=6 is NOT run
as a global Schur complement: its transverse basis alone would be ~1.5 GB and
the fiber block ~0.9 GB (multi-GB temporaries) => OOM risk; that is reported
honestly as a feasibility limit, not faked.

Usage:
  python testc_global_schur.py --mode free
  python testc_global_schur.py --mode u1
  python testc_global_schur.py --mode su2 --L 4 --betas 2.2,2.3,2.4 --nkeep 5
  python testc_global_schur.py --mode all
"""
import sys
import os
import time
import argparse
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LATTICE_DIR = os.path.join(REPO, "research", "lattice-cpb")

# ---------------------------------------------------------------- memory probe
try:
    import psutil
    _PROC = psutil.Process(os.getpid())
    def rss_mb():
        return _PROC.memory_info().rss / 1e6
    def peak_mb():
        # Windows exposes the true peak working set; fall back to sampled rss.
        mi = _PROC.memory_info()
        return (getattr(mi, "peak_wset", 0) or mi.rss) / 1e6
except Exception:                                    # pragma: no cover
    def rss_mb():
        return float("nan")
    def peak_mb():
        return float("nan")

_PEAK = [0.0]
def note_peak(tag=""):
    r = rss_mb()
    if r == r and r > _PEAK[0]:
        _PEAK[0] = r
    pk = peak_mb()
    if pk == pk and pk > _PEAK[0]:
        _PEAK[0] = pk
    return r

def report_peak(tag=""):
    note_peak()
    print(f"  [peak RSS so far: {_PEAK[0]:.0f} MB {tag}]", flush=True)


# ======================================================================
# generic Schur-complement driver on an abstract coarse/fiber split
# ======================================================================
def schur_softest(K, Cb, vacuum_cols=None, extra_mass2=0.0, label=""):
    """K : (D,D) symmetric quadratic-form matrix (Hessian-like).
    Cb : (D,nc) orthonormal coarse basis.  Fiber = orthonormal complement.
    vacuum_cols : (D,k) *full-space* vectors spanning the 'vacuum' directions to
                  be removed from the coarse space (e.g. the constant mode).
    Returns dict with spectra of A, D, S and lambda_soft (min S-eigenvalue on the
    coarse space orthogonal to the coarse image of vacuum_cols).
    """
    D_ = K.shape[0]
    nc = Cb.shape[1]
    # fiber basis = orthonormal complement of Cb
    Uf = _orth_complement(Cb)
    # blocks
    A = Cb.T @ K @ Cb
    Dm = Uf.T @ K @ Uf
    W = Cb.T @ K @ Uf
    A = 0.5 * (A + A.T)
    Dm = 0.5 * (Dm + Dm.T)
    if extra_mass2:
        A = A + extra_mass2 * np.eye(nc)
        Dm = Dm + extra_mass2 * np.eye(Dm.shape[0])
    # Schur complement
    if Dm.shape[0] == 0:
        S = A.copy()
    else:
        try:
            Sol = np.linalg.solve(Dm, W.T)          # (f, nc)
        except np.linalg.LinAlgError:
            Sol = np.linalg.pinv(Dm) @ W.T
        S = A - W @ Sol
        S = 0.5 * (S + S.T)
    evA = np.sort(np.linalg.eigvalsh(A))
    evD = np.sort(np.linalg.eigvalsh(Dm)) if Dm.shape[0] else np.array([])
    evS = np.sort(np.linalg.eigvalsh(S))
    # lambda_soft: min S-eigenvalue on coarse space orthogonal to vacuum image
    if vacuum_cols is not None and vacuum_cols.size:
        Vc = Cb.T @ vacuum_cols                       # coarse coords of vacuum
        Qv = np.linalg.qr(Vc)[0]
        rk = _rank(Vc)
        Qv = Qv[:, :rk]
        if rk > 0 and rk < nc:
            Qperp = _orth_complement(Qv)              # nc x (nc-rk)
            Sred = Qperp.T @ S @ Qperp
            Sred = 0.5 * (Sred + Sred.T)
            lam = np.sort(np.linalg.eigvalsh(Sred))
            lam_soft = float(lam[0]) if lam.size else float(evS[0])
        else:
            lam_soft = float(evS[0])
    else:
        lam_soft = float(evS[0])
    return dict(evA=evA, evD=evD, evS=evS, lam_soft=lam_soft,
                nc=nc, nf=Dm.shape[0], D=D_, S=S)


def _rank(M, tol=1e-9):
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > tol * max(1.0, s[0])))


def _orth_complement(B, tol=1e-9):
    """Orthonormal basis (D x (D-r)) of the orthogonal complement of span(B)."""
    D_ = B.shape[0]
    if B.shape[1] == 0:
        return np.eye(D_)
    from scipy.linalg import null_space
    return null_space(B.T, rcond=tol)


# ======================================================================
# FREE SCALAR calibration
# ======================================================================
def lap1d(L, m2=0.0):
    M = np.zeros((L, L))
    for x in range(L):
        M[x, x] += 2.0 + m2
        M[x, (x - 1) % L] += -1.0
        M[x, (x + 1) % L] += -1.0
    return M

def lap2d(L, m2=0.0):
    N = L * L
    M = np.zeros((N, N))
    def ix(x, y): return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = ix(x, y)
            M[i, i] += 4.0 + m2
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                M[i, ix(x + dx, y + dy)] += -1.0
    return M

def block_basis_1d(L, b):
    n = L // b
    C = np.zeros((L, n))
    for I in range(n):
        C[I*b:(I+1)*b, I] = 1.0 / np.sqrt(b)
    return C

def block_basis_2d(L, bx, by):
    nx, ny = L // bx, L // by
    N = L * L
    C = np.zeros((N, nx * ny))
    def sidx(x, y): return x * L + y
    for I in range(nx):
        for J in range(ny):
            col = I * ny + J
            for a in range(bx):
                for b in range(by):
                    C[sidx(I*bx + a, J*by + b), col] = 1.0 / np.sqrt(bx * by)
    return C

def run_free_scalar():
    print("=" * 78)
    print("FREE SCALAR calibration:  lambda_soft(L) on coarse modes orthogonal to")
    print("the uniform (vacuum) mode.  Expect massless ~ c/L^2 -> 0 ; massive -> m^2.")
    print("=" * 78)
    rows1 = []
    print("\n-- 1D periodic, block b=2 --")
    print(f"{'L':>5} {'N':>6} {'n':>5} {'f':>6} {'lam_soft':>13} "
          f"{'L^2*lam':>11} {'infM':>10}  (m=0)")
    for L in (24, 48, 96, 192, 384):
        M = lap1d(L, 0.0)
        C = block_basis_1d(L, 2)
        e0 = np.ones((L, 1)) / np.sqrt(L)
        r = schur_softest(M, C, vacuum_cols=e0)
        evM = np.linalg.eigvalsh(M)
        rows1.append((L, r['lam_soft'], evM[0]))
        print(f"{L:>5} {L:>6} {r['nc']:>5} {r['nf']:>6} {r['lam_soft']:>13.6e} "
              f"{L*L*r['lam_soft']:>11.4f} {evM[0]:>10.2e}")
    print("\n-- 1D massive m^2=0.5, block b=2 (expect lam_soft -> 0.5) --")
    m2 = 0.5
    print(f"{'L':>5} {'lam_soft':>13} {'lam-m^2':>11} {'infM':>10}")
    for L in (24, 48, 96, 192, 384):
        M = lap1d(L, m2)
        C = block_basis_1d(L, 2)
        e0 = np.ones((L, 1)) / np.sqrt(L)
        r = schur_softest(M, C, vacuum_cols=e0)
        evM = np.linalg.eigvalsh(M)
        print(f"{L:>5} {r['lam_soft']:>13.6e} {r['lam_soft']-m2:>11.3e} {evM[0]:>10.6f}")
    report_peak("after 1D")
    rows2 = []
    print("\n-- 2D periodic, block 2x2 --")
    print(f"{'L':>4} {'N':>6} {'n':>5} {'f':>6} {'lam_soft(m0)':>13} "
          f"{'L^2*lam':>10} {'lam_soft(m=.5)':>14}")
    for L in (8, 12, 16, 24):
        M0 = lap2d(L, 0.0); M1 = lap2d(L, 0.5)
        C = block_basis_2d(L, 2, 2)
        e0 = np.ones((L * L, 1)) / L
        r0 = schur_softest(M0, C, vacuum_cols=e0)
        r1 = schur_softest(M1, C, vacuum_cols=e0)
        rows2.append((L, r0['lam_soft']))
        print(f"{L:>4} {L*L:>6} {r0['nc']:>5} {r0['nf']:>6} {r0['lam_soft']:>13.6e} "
              f"{L*L*r0['lam_soft']:>10.4f} {r1['lam_soft']:>14.6f}")
    report_peak("after 2D")
    # power-law fit (massless): log lam vs log L
    if len(rows1) >= 3:
        Ls = np.array([r[0] for r in rows1], float)
        lam = np.array([max(r[1], 1e-30) for r in rows1])
        p = np.polyfit(np.log(Ls), np.log(lam), 1)[0]
        print(f"\n  1D massless: d log(lambda_soft)/d log(L) = {p:.3f}   (expect ~ -2)")
    if len(rows2) >= 3:
        Ls = np.array([r[0] for r in rows2], float)
        lam = np.array([max(r[1], 1e-30) for r in rows2])
        p = np.polyfit(np.log(Ls), np.log(lam), 1)[0]
        print(f"  2D massless: d log(lambda_soft)/d log(L) = {p:.3f}   (expect ~ -2)")
    return rows1, rows2


# ======================================================================
# FREE U(1) MAXWELL calibration (gauge theory; d=3): transverse photon
# ======================================================================
def u1_fieldstrength_op(L):
    """B: maps link field A[(mu,x)] (3*L^3) -> F_{mu nu}(x)  (3*L^3 rows)."""
    ns = L ** 3
    d = 3
    planes = [(0, 1), (0, 2), (1, 2)]
    def sidx(x, y, z): return (x * L + y) * L + z
    def aidx(mu, s): return mu * ns + s
    rows, cols, vals = [], [], []
    row = 0
    for (mu, nu) in planes:
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    s = sidx(x, y, z)
                    c = [x, y, z]
                    cp = c.copy(); cp[mu] = (cp[mu] + 1) % L
                    sp = sidx(*cp)
                    cq = c.copy(); cq[nu] = (cq[nu] + 1) % L
                    sq = sidx(*cq)
                    # F = A_nu(x+mu) - A_nu(x) - A_mu(x+nu) + A_mu(x)
                    rows += [row, row, row, row]
                    cols += [aidx(nu, sp), aidx(nu, s), aidx(mu, sq), aidx(mu, s)]
                    vals += [1.0, -1.0, -1.0, 1.0]
                    row += 1
    from scipy.sparse import csr_matrix
    B = csr_matrix((vals, (rows, cols)), shape=(row, d * ns))
    return B, d

def u1_block_basis(L, b):
    """block-constant basis on link fields (block factor b in all 3 dirs)."""
    ns = L ** 3
    Lc = L // b
    def sidx(x, y, z): return (x * L + y) * L + z
    nblk = Lc ** 3
    D = 3 * ns
    cols = []
    col_index = {}
    ncol = 0
    # each (mu, block) -> one basis vector
    C = np.zeros((D, 3 * nblk))
    k = 0
    for mu in range(3):
        for I in range(Lc):
            for J in range(Lc):
                for K in range(Lc):
                    for a in range(b):
                        for bb in range(b):
                            for cc in range(b):
                                s = sidx(I*b + a, J*b + bb, K*b + cc)
                                C[mu * ns + s, k] = 1.0 / np.sqrt(b ** 3)
                    k += 1
    return C

def run_u1():
    print("=" * 78)
    print("FREE U(1) MAXWELL calibration (d=3): remove F=0 kernel (pure gauge +")
    print("torons) -> transverse photon, block-average, Schur S_L, lambda_soft.")
    print("Expect massless ~ c/L^2 -> 0 ; Proca massive -> m^2.")
    print("=" * 78)
    print(f"{'L':>4} {'D':>6} {'ker':>5} {'trans':>6} {'ncoarse':>8} "
          f"{'lam_soft(m0)':>13} {'L^2*lam':>10} {'lam_soft(m=.5)':>14}")
    res = []
    for L in (6, 8, 10):
        B, d = u1_fieldstrength_op(L)
        ns = L ** 3
        D = d * ns
        M = (B.T @ B).toarray()
        M = 0.5 * (M + M.T)
        ev, V = np.linalg.eigh(M)
        kmax = max(1.0, abs(ev).max())
        ker = int(np.sum(np.abs(ev) < 1e-9 * kmax))
        Pk = V[:, ker:]                                   # transverse/coexact basis
        # coarse block basis, projected to transverse space
        Cb = u1_block_basis(L, 2)
        Ct = Pk.T @ Cb
        Uc, _ = np.linalg.qr(Ct)
        rk = _rank(Ct)
        Uc = Uc[:, :rk]
        # fibre = complement of Uc within transverse space Pk
        Uf = _orth_complement(Uc)                         # in R^{trans}, expressed in tx coords
        # transverse Maxwell operator in tx coords
        Mt = Pk.T @ M @ Pk
        Mt = 0.5 * (Mt + Mt.T)
        A = Uc.T @ Mt @ Uc
        Dm = Uf.T @ Mt @ Uf
        W = Uc.T @ Mt @ Uf
        A = 0.5 * (A + A.T); Dm = 0.5 * (Dm + Dm.T)
        S = A - W @ np.linalg.solve(Dm, W.T) if Dm.shape[0] else A.copy()
        S = 0.5 * (S + S.T)
        evS = np.sort(np.linalg.eigvalsh(S))
        lam0 = float(evS[0])
        # Proca massive: add m^2 to transverse operator
        m2 = 0.5
        At = A + m2 * np.eye(A.shape[0]); Dt = Dm + m2 * np.eye(Dm.shape[0])
        St = At - W @ np.linalg.solve(Dt, W.T) if Dm.shape[0] else At.copy()
        lam1 = float(np.sort(np.linalg.eigvalsh(0.5 * (St + St.T)))[0])
        res.append((L, lam0))
        print(f"{L:>4} {D:>6} {ker:>5} {Pk.shape[1]:>6} {rk:>8} "
              f"{lam0:>13.6e} {L*L*lam0:>10.4f} {lam1:>14.6f}")
        note_peak(); report_peak(f"L={L} u1")
    if len(res) >= 3:
        Ls = np.array([r[0] for r in res], float)
        lam = np.array([max(r[1], 1e-30) for r in res])
        p = np.polyfit(np.log(Ls), np.log(lam), 1)[0]
        print(f"\n  U(1) massless photon: d log(lambda_soft)/d log(L) = {p:.3f}  (expect ~ -2)")
    return res


# ======================================================================
# SU(2) static-Hessian proxy at L=4
# ======================================================================
def su2_block_basis_full(L, b=2):
    """(D, 4*(L/b)^4*3) block-constant basis in the full link-algebra space,
    fidx(mu,s,a)=(mu*L^4+s)*3+a.  Block factor b in all four directions."""
    ns = L ** 4
    D = 4 * ns * 3
    Lc = L // b
    def sidx(c): return ((c[0] * L + c[1]) * L + c[2]) * L + c[3]
    ncol = 4 * (Lc ** 4) * 3
    C = np.zeros((D, ncol))
    k = 0
    import itertools
    for mu in range(4):
        for I in itertools.product(range(Lc), repeat=4):
            for a in range(3):
                # fill the b^4 sites of this block
                for off in itertools.product(range(b), repeat=4):
                    c = [I[j] * b + off[j] for j in range(4)]
                    s = sidx(c)
                    C[(mu * ns + s) * 3 + a, k] = 1.0 / np.sqrt(b ** 4)
                k += 1
    return C

def su2_uniform_basis_full(L):
    """(D, 12) spatially-uniform full vectors for each (mu,color)."""
    ns = L ** 4
    D = 4 * ns * 3
    C = np.zeros((D, 12))
    for mu in range(4):
        for a in range(3):
            k = mu * 3 + a
            for s in range(ns):
                C[(mu * ns + s) * 3 + a, k] = 1.0 / np.sqrt(ns)
    return C

def run_su2(L, betas, nkeep, ntherm, gapevery, seed):
    sys.path.insert(0, LATTICE_DIR)
    import w4_blocked as W
    print("=" * 78)
    print(f"SU(2) STATIC-ACTION HESSIAN PROXY  L={L}  betas={betas}  "
          f"nkeep={nkeep}")
    print("*** PROXY ONLY: static spatial Hessian != time-transfer gap ***")
    print("=" * 78)
    ns = L ** 4
    Cb = su2_block_basis_full(L, 2)
    Vun = su2_uniform_basis_full(L)                     # full-space uniform modes
    allres = {}
    for beta in betas:
        rng = np.random.default_rng(seed)
        t0 = time.time()
        cfgs = W.sample2(L, beta, rng, ntherm, nkeep, gapevery)
        print(f"\n[beta={beta}] thermalised {ntherm}+{nkeep*gapevery} sweeps "
              f"in {time.time()-t0:.1f}s", flush=True)
        rows = []
        for ci, U in enumerate(cfgs):
            t1 = time.time()
            H = W.hessian2(U, beta)
            Q = W.gauge_basis2(U)
            Z = W.transverse_basis(Q)
            Ht = Z.T @ H @ Z
            Ht = 0.5 * (Ht + Ht.T)
            evHt = np.linalg.eigvalsh(Ht)
            # coarse-transverse basis
            Ct = Z.T @ Cb
            Uc, _ = np.linalg.qr(Ct)
            rk = _rank(Ct)
            Uc = Uc[:, :rk]
            Uf = _orth_complement(Uc)                    # complement in R^{dim Z}
            A = Uc.T @ Ht @ Uc
            Dm = Uf.T @ Ht @ Uf
            Wc = Uc.T @ Ht @ Uf
            A = 0.5 * (A + A.T); Dm = 0.5 * (Dm + Dm.T)
            try:
                Sol = np.linalg.solve(Dm, Wc.T)
            except np.linalg.LinAlgError:
                Sol = np.linalg.pinv(Dm) @ Wc.T
            S = A - Wc @ Sol
            S = 0.5 * (S + S.T)
            evA = np.sort(np.linalg.eigvalsh(A))
            evS = np.sort(np.linalg.eigvalsh(S))
            # vacuum = coarse image of the 12 uniform full vectors, then removed
            Vc = Uc.T @ (Z.T @ Vun)
            Qv = np.linalg.qr(Vc)[0]
            rkv = _rank(Vc)
            Qv = Qv[:, :rkv]
            if rkv > 0 and rkv < rk:
                Qp = _orth_complement(Qv)
                Sred = Qp.T @ S @ Qp
                Sred = 0.5 * (Sred + Sred.T)
                lamS_nu = float(np.sort(np.linalg.eigvalsh(Sred))[0])
            else:
                lamS_nu = float(evS[0])
            # D-fibre spectrum (indefinite? Savvidy), conditioning, inertia check
            evD = np.linalg.eigvalsh(Dm)
            nHt = int(np.sum(evHt < 0)); nD = int(np.sum(evD < 0)); nS = int(np.sum(evS < 0))
            ld_res = (np.linalg.slogdet(Ht)[1] - np.linalg.slogdet(Dm)[1]
                      - np.linalg.slogdet(S)[1])
            rows.append(dict(
                trans_negfrac=float(np.mean(evHt < 0)),
                trans_min=float(evHt[0]),
                A_min=float(evA[0]), A_low=evA[:6],
                S_min=float(evS[0]), S_low=evS[:8],
                lamS_nonuniform=lamS_nu,
                Dmin=float(evD[0]) if evD.size else np.nan,
                Dminabs=float(np.abs(evD).min()) if evD.size else np.nan,
                Dcond=float(np.abs(evD).max()/max(np.abs(evD).min(),1e-30)) if evD.size else np.nan,
                nHt=nHt, nD=nD, nS=nS, det_res=float(abs(ld_res)),
                nc=rk, ntrans=Ht.shape[0]))
            print(f"  cfg{ci}: trans[neg {rows[-1]['trans_negfrac']*100:.1f}% "
                  f"min {rows[-1]['trans_min']/beta:+.3f}b] | "
                  f"A_min/b {rows[-1]['A_min']/beta:+.3f} | "
                  f"S_min/b {rows[-1]['S_min']/beta:+.3f} | "
                  f"lamS(nonunif)/b {rows[-1]['lamS_nonuniform']/beta:+.3f} | "
                  f"Dmin/b {rows[-1]['Dmin']/beta:+.3f} | "
                  f"inertia neg Ht/Df/S {nHt}/{nD}/{nS} (=? {nHt==nD+nS}) | "
                  f"cond(Df) {rows[-1]['Dcond']:.1e} | det_res {rows[-1]['det_res']:.1e} | "
                  f"nc={rk} | ({time.time()-t1:.1f}s)", flush=True)
            note_peak()
        report_peak(f"beta={beta}")
        allres[beta] = rows
        # summary
        lamS = np.array([r['lamS_nonuniform'] for r in rows]) / beta
        Smin = np.array([r['S_min'] for r in rows]) / beta
        Amin = np.array([r['A_min'] for r in rows]) / beta
        print(f"  [beta={beta}] lam_soft(non-uniform)/beta: "
              f"median {np.median(lamS):+.3f}  mean {lamS.mean():+.3f}  "
              f"[min {lamS.min():+.3f}, max {lamS.max():+.3f}]")
        print(f"  [beta={beta}] A_min/beta mean {Amin.mean():+.3f}; "
              f"S_min/beta (incl uniform) mean {Smin.mean():+.3f}; "
              f"trans neg-frac mean {np.mean([r['trans_negfrac'] for r in rows])*100:.1f}%")
    return allres


# ======================================================================
# SU(2) coarse blocked-Hessian A -- L-scaling only (matrix-free, memory-lean).
# A = Vc^T H Vc, Vc = gauge-transverse projection of the block-constant basis.
# This is the leading (well-defined, positive) piece of S_L.  The full Schur
# complement is ill-posed already at L=4 (indefinite near-singular fibre), so
# only A's L-scaling is reported; it is the honest "kinetic tail" diagnostic.
# ======================================================================
def _su2_expm_alg(a):
    """Vectorised exp(i t^a a_a) for SU(2): a has shape (...,3). Returns (...,2,2)."""
    th = np.sqrt(np.einsum('...a,...a->...', a, a))
    half = 0.5 * th
    # exp(i th n.sigma/2) = cos(half) I + i sin(half) n.sigma
    n = np.zeros_like(a)
    nz = th > 1e-12
    n[nz] = a[nz] / th[nz][..., None]
    s = np.sin(half)[..., None, None]
    c = np.cos(half)[..., None, None]
    out = c * np.eye(2, dtype=np.complex128)
    # i n.sigma : sigma_1=[[0,1],[1,0]], sigma_2=[[0,-i],[i,0]], sigma_3=[[1,0],[0,-1]]
    i = 1j
    out = out + s * i * (n[..., 0, None, None] * np.array([[0, 1], [1, 0]]) +
                         n[..., 1, None, None] * np.array([[0, -i], [i, 0]]) +
                         n[..., 2, None, None] * np.array([[1, 0], [0, -1]]))
    return out

def su2_coarse_A_scaling(L, beta, nkeep, ntherm, gapevery, seed, eps=1e-4):
    sys.path.insert(0, LATTICE_DIR)
    import w4_blocked as W
    ns = L ** 4
    D = 4 * ns * 3
    print(f"\n--- SU(2) coarse blocked Hessian A, L={L}, beta={beta} (matrix-free) ---")
    rng = np.random.default_rng(seed)
    cfgs = W.sample2(L, beta, rng, ntherm, nkeep, gapevery)
    Cb = su2_block_basis_full(L, 2)
    res = []
    for ci, U in enumerate(cfgs):
        t1 = time.time()
        Q = W.gauge_basis2(U)                      # gauge orbit (dense) -- memory heavy
        G = Q.T @ Q
        # transverse projection of block basis: Y = Cb - Q G^{-1} Q^T Cb
        Y = Cb - Q @ np.linalg.solve(G, Q.T @ Cb)
        Vc = np.linalg.qr(Y)[0]
        rc = _rank(Y); Vc = Vc[:, :rc]
        # matrix-free A_ij = v_i . grad(U(+eps v_j)) via directional FD
        A = np.empty((rc, rc))
        base = U.copy()
        Uf = base.reshape(4, ns, 2, 2)
        Vc4 = Vc.reshape(4, ns, 3, rc)             # (4,ns,3,rc)
        for j in range(rc):
            aj = Vc4[:, :, :, j]                    # (4,ns,3)
            Ep = _su2_expm_alg(+eps * aj)           # (4,ns,2,2)
            Em = _su2_expm_alg(-eps * aj)
            Up = base.copy().reshape(4, ns, 2, 2)
            Um = base.copy().reshape(4, ns, 2, 2)
            Up = Ep @ Uf
            Um = Em @ Uf
            gp = W.grad2(Up.reshape(U.shape), beta).reshape(-1)
            gm = W.grad2(Um.reshape(U.shape), beta).reshape(-1)
            colvec = Vc.T @ ((gp - gm) / (2.0 * eps))
            A[:, j] = colvec
        A = 0.5 * (A + A.T)
        evA = np.sort(np.linalg.eigvalsh(A))
        # split off the 12 uniform (toron) transverse directions: A on the rest
        Vun = su2_uniform_basis_full(L)
        y = Vc.T @ (Vun - Q @ np.linalg.solve(G, Q.T @ Vun))
        Qv = np.linalg.qr(y)[0]; rkv = _rank(y); Qv = Qv[:, :rkv]
        if 0 < rkv < rc:
            Qp = _orth_complement(Qv)
            Ar = Qp.T @ A @ Qp; Ar = 0.5 * (Ar + Ar.T)
            evAnu = np.sort(np.linalg.eigvalsh(Ar))
            toron_ov = float(np.linalg.norm(Qv.T @ np.linalg.eigh(A)[1][:, 0]))
        else:
            evAnu = evA; toron_ov = np.nan
        res.append(evA)
        print(f"  cfg{ci}: A_min/beta {evA[0]/beta:+.4f}  "
              f"A_min(non-uniform)/beta {evAnu[0]/beta:+.4f}  "
              f"toron-overlap {toron_ov:.2f}  A_max/beta {evA[-1]/beta:+.3f} "
              f" rc={rc}  ({time.time()-t1:.1f}s)", flush=True)
        del Q, G, Y
        note_peak(); report_peak(f"L={L} coarseA")
    mins = np.array([e[0] for e in res])
    med = np.median(mins)
    print(f"  L={L} beta={beta}: A_min/beta median {med/beta:+.4f} "
          f"(mean {mins.mean()/beta:+.4f})  [coarse dim {res[0].size}]")
    return med, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='all',
                    choices=['free', 'u1', 'su2', 'su2A', 'all'])
    ap.add_argument('--Ls', type=str, default='4,6', help='L list for su2A mode')
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--betas', type=str, default='2.2,2.3,2.4')
    ap.add_argument('--nkeep', type=int, default=5)
    ap.add_argument('--ntherm', type=int, default=200)
    ap.add_argument('--gapevery', type=int, default=20)
    ap.add_argument('--seed', type=int, default=1234)
    a = ap.parse_args()
    betas = [float(x) for x in a.betas.split(',')]
    t0 = time.time()
    if a.mode in ('free', 'all'):
        run_free_scalar()
    if a.mode in ('u1', 'all'):
        run_u1()
    if a.mode in ('su2', 'all'):
        if a.L > 4:
            print(f"\n[refuse] SU(2) global Schur at L={a.L} is NOT run: the "
                  f"transverse basis at L={a.L} would be >~1.5 GB and the fiber "
                  f"block >~0.9 GB (multi-GB temporaries) -> OOM per charter 4. "
                  f"Use L=4; see TEST-C-GLOBAL-SCHUR.md feasibility section.")
        else:
            run_su2(a.L, betas, a.nkeep, a.ntherm, a.gapevery, a.seed)
    if a.mode == 'su2A':
        Ls = [int(x) for x in a.Ls.split(',')]
        for L in Ls:
            su2_coarse_A_scaling(L, betas[0], a.nkeep, a.ntherm, a.gapevery, a.seed)
    print(f"\n[done] total wall {time.time()-t0:.1f}s, peak RSS {_PEAK[0]:.0f} MB")


if __name__ == '__main__':
    main()
