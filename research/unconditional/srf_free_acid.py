"""
SRF 2.0 free-field acid test -- exact, tiny, single-process.
Free Gaussian scalar on a lattice, block-average coarse map C, orthogonal split
    R^N = C_space (block-constant, the 'coarse V' that RETAINS soft modes)
          (+) F_space (block-zero-sum, the 'fiber eta')
Hessian M = -Delta + m^2.  Block form  [[A, W],[W^T, D]].  Feshbach/Schur:
    S(z) = A - z - Sigma(z),   Sigma(z) = W (D - z)^{-1} W^T  (>=0 for z<gamma_fib).
Checks:
  (1) free massless m=0:  inf spec(coarse S) must be 0 EXACTLY; k=0 constant mode
      lives in C_space, survives (W^T const = 0, Sigma const = 0).
  (2) free massive m>0:   physical mass m^2 reproduced, no shift.
  (3) self-energy sign Sigma >= 0;  A - Sigma lowers, never gaps.
  (4) fiber gap gamma_fib > 0 already at m=0  => fiber-control is free-field-blind.
  (5) the LOCAL single-block 'SRF-CUBE' reading -> spurious gap 2-2cos(pi/b) > 0.
  (6) 1D and 2D; L-scaling of first nonzero mode -> 0 (no L-uniform free gap).
Memory: arrays are N x N with N <= a few hundred.  One process. Safe.
"""
import numpy as np
np.set_printoptions(precision=6, suppress=False, linewidth=140)

# ---------- lattice builders ----------
def lap1d(L, periodic=True):
    M = np.zeros((L, L))
    for x in range(L):
        M[x, x] += 2.0
        xm = (x - 1) % L if periodic else x - 1
        xp = (x + 1) % L if periodic else x + 1
        if 0 <= xm < L:
            M[x, xm] += -1.0
        if 0 <= xp < L:
            M[x, xp] += -1.0
    return M

def lap1d_open(b):
    """True OPEN path Laplacian on b sites (degree 1 at the ends): the honest
    'isolated free block'.  Lowest non-constant eigenvalue = 2 - 2cos(pi/b)."""
    M = -np.eye(b, k=1) - np.eye(b, k=-1)     # off-diagonals -1
    for x in range(b):
        M[x, x] = 0.0
        if x > 0:    M[x, x] += 1.0
        if x < b-1:  M[x, x] += 1.0
    return M

def lap1d_dirichlet(b):
    """DIRICHLET path Laplacian (neighbours pinned to 0): 'block with neighbours held
    at the block mean'.  Identical to the lap1d(periodic=False) above.
    Lowest eigenvalue = 2 - 2cos(pi/(b+1))."""
    return lap1d(b, periodic=False)

def lap2d(Lx, Ly, periodic=True):
    N = Lx * Ly
    def idx(x, y): return x * Ly + y
    M = np.zeros((N, N))
    for x in range(Lx):
        for y in range(Ly):
            i = idx(x, y); M[i, i] += 4.0
            for dx, dy in ((1,0),(-1,0),(0,1),(0,-1)):
                xn, yn = x+dx, y+dy
                if periodic:
                    xn %= Lx; yn %= Ly
                elif not (0 <= xn < Lx and 0 <= yn < Ly):
                    continue
                M[i, idx(xn, yn)] += -1.0
    return M

def block_map_1d(L, b):
    """C: R^L -> R^n, block average (rows sum to 1 over each block of b sites)."""
    assert L % b == 0
    n = L // b
    C = np.zeros((n, L))
    for I in range(n):
        for k in range(b):
            C[I, I*b + k] = 1.0 / b
    return C

def block_map_2d(Lx, Ly, bx, by):
    assert Lx % bx == 0 and Ly % by == 0
    nx, ny = Lx//bx, Ly//by
    N = Lx*Ly; n = nx*ny
    def sidx(x, y): return x*Ly + y
    def cidx(I, J): return I*ny + J
    C = np.zeros((n, N))
    for I in range(nx):
        for J in range(ny):
            for a in range(bx):
                for bb in range(by):
                    C[cidx(I,J), sidx(I*bx+a, J*by+bb)] = 1.0/(bx*by)
    return C

# ---------- orthogonal split + blocks ----------
def split_and_blocks(M, C):
    """Return A_c (n x n), D_f (f x f), W (n x f), Uc (N x n), Uf (N x f),
    and the constant vector e0 lifted, for the orthogonal split R^N=C^perp F."""
    N = M.shape[0]; n = C.shape[0]
    # orthonormal basis of C_space = row space of C:
    #   Uc = C^T (C C^T)^{-1/2}, columns orthonormal.
    CCt = C @ C.T
    # symmetric inverse sqrt
    w, V = np.linalg.eigh(CCt)
    CCt_isqrt = V @ np.diag(1.0/np.sqrt(w)) @ V.T
    Uc = C.T @ CCt_isqrt                      # N x n, orthonormal cols
    # fiber = null space of C (block-zero-sum). Orthonormal complement:
    #   columns of Uf span ker(C).
    # use SVD of C
    Uu, S, Vt = np.linalg.svd(C)
    f = N - n
    Uf = Vt[n:].T                             # N x f, orthonormal, C Uf = 0
    # sanity
    assert Uf.shape[1] == f
    # blocks
    A = Uc.T @ M @ Uc
    D = Uf.T @ M @ Uf
    W = Uc.T @ M @ Uf
    # symmetrize (numerical)
    A = 0.5*(A + A.T); D = 0.5*(D + D.T)
    return A, D, W, Uc, Uf

def sigma(W, D, z):
    return W @ np.linalg.solve(D - z*np.eye(D.shape[0]), W.T)

def rep(name, x): print(f"    {name:<46s} {x}")

def run_1d(L, b, m2, tag, do_det_check=True, do_spurious=True):
    print(f"\n====================  {tag}: 1D L={L}, b={b}, m^2={m2}  ====================")
    M = lap1d(L, periodic=True) + m2*np.eye(L)
    C = block_map_1d(L, b)
    A, D, W, Uc, Uf = split_and_blocks(M, C)
    N, n = L, C.shape[0]; f = N - n
    print(f"  dims: N(sites)={N}  n(coarse C)={n}  f(fiber)={f}")

    evM = np.sort(np.linalg.eigvalsh(M))
    evA = np.sort(np.linalg.eigvalsh(A))
    evD = np.sort(np.linalg.eigvalsh(D))
    rep("inf spec FULL Hessian M", evM[0])
    rep("first nonzero of M (mode>=1)", evM[1])
    rep("inf spec coarse block A", evA[0])
    rep("fiber gap gamma_fib = inf spec D", evD[0])

    # ---- constant / soft mode lives in C?  and survives Feshbach ----
    e0 = np.ones(N)/np.sqrt(N)
    inC = np.linalg.norm(Uc @ (Uc.T @ e0) - e0)
    rep("||e0 - P_C e0||  (0 => k=0 mode in V)", inC)
    # coarse coords of e0
    e0c = Uc.T @ e0
    # A e0c should be m2 * e0c (e0 is eigenvector of M, in C)
    Ae0c = A @ e0c
    rep("||A e0c - m2 e0c||", np.linalg.norm(Ae0c - m2*e0c))
    WT_e0 = Uf.T @ (M @ e0)                    # = W^T e0c  (fiber coupling of soft mode)
    rep("||W^T e0|| = ||P_F M e0||  (0 => soft decoupled)", np.linalg.norm(WT_e0))
    S0 = sigma(W, D, 0.0)
    rep("||Sigma(z) e0|| at z=0", np.linalg.norm(S0 @ e0c))

    # ---- Feshbach coarse operator spectrum, evaluated properly ----
    # The physical soft eigenvalue is the smallest z with S(z) singular on the
    # coarse space.  We confirm it = inf spec(M) and that the minimiser is e0c.
    # For the massless case the answer must be EXACTLY the full inf spec.
    # (a) static z=0 coarse operator for massless:
    if abs(m2) < 1e-14:
        scorer = A - sigma(W, D, 0.0)         # = A - Sigma(0), the z=0 coarse operator
        evSc = np.sort(np.linalg.eigvalsh(scorer))
        rep("inf spec (A - Sigma(0))  [massless]", evSc[0])
        rep("  its minimiser overlap^2 with e0c",
            (scorer @ e0c) @ e0c / (e0c@e0c) / np.linalg.norm(e0c)**2 if False else
            float(np.dot(e0c, scorer@e0c)/np.dot(e0c,e0c)))
    # (b) physical: smallest z where det S(z)=0  (scan + polish, avoid fiber band)
    #     we instead certify via Schur identity and the exact eigenvector argument.
    # Schur identity check: det(M-z) = det(D-z) det(S(z)) at several generic z.
    if do_det_check:
        gl = []
        for z in [-0.37, 0.123, 0.777]:
            lhs = np.linalg.slogdet(M - z*np.eye(N))[1]
            rhs = np.linalg.slogdet(D - z*np.eye(f))[1] + np.linalg.slogdet(A - z*np.eye(n) - sigma(W,D,z))[1]
            gl.append(abs(lhs-rhs))
        rep("max |logdet(M-z) - logdet(D-z) - logdet S(z)|", max(gl))

    # ---- self-energy sign + 'lowering' ----
    Sig0 = sigma(W, D, 0.0)
    evS = np.sort(np.linalg.eigvalsh(Sig0))
    rep("min eig Sigma(0)  (>=0 => self-energy PSD)", evS[0])
    rep("max eig Sigma(0)", evS[-1])
    rep("min eig A  vs  min eig (A - Sigma(0))",
        f"{evA[0]:.6g}  vs  {np.sort(np.linalg.eigvalsh(A - Sig0))[0]:.6g}")

    # ---- local single-block 'SRF-CUBE' reading -> spurious gap ----
    if do_spurious:
        ev_open = np.sort(np.linalg.eigvalsh(lap1d_open(b) + m2*np.eye(b)))
        ev_dir  = np.sort(np.linalg.eigvalsh(lap1d_dirichlet(b) + m2*np.eye(b)))
        # lowest mode of an isolated block WITH its own constant allowed = 0;
        # the spurious gap is the first *relative* (gradient) mode:
        rep(f"LOCAL block OPEN   first relative mode (2-2cos(pi/b)={2-2*np.cos(np.pi/b):.6f})",
            ev_open[1])
        rep(f"LOCAL block DIRICHLET lowest mode (2-2cos(pi/(b+1))={2-2*np.cos(np.pi/(b+1)):.6f})",
            ev_dir[0])

    return dict(infM=evM[0], firstNZ=evM[1], infA=evA[0], gamma=evD[0], evD=evD, evM=evM)

def run_2d(L, b, m2, tag):
    print(f"\n====================  {tag}: 2D {L}x{L}, block {b}x{b}, m^2={m2}  ====================")
    M = lap2d(L, L, periodic=True) + m2*np.eye(L*L)
    C = block_map_2d(L, L, b, b)
    A, D, W, Uc, Uf = split_and_blocks(M, C)
    N = L*L; n = C.shape[0]; f = N-n
    print(f"  dims: N={N}  n(coarse)={n}  f(fiber)={f}")
    evM = np.sort(np.linalg.eigvalsh(M))
    evD = np.sort(np.linalg.eigvalsh(D))
    rep("inf spec FULL Hessian M", evM[0])
    rep("first nonzero of M (mode>=1)", evM[1])
    rep("fiber gap gamma_fib", evD[0])
    e0 = np.ones(N)/np.sqrt(N)
    rep("||e0 - P_C e0||", np.linalg.norm(Uc@(Uc.T@e0)-e0))
    rep("||P_F M e0||", np.linalg.norm(Uf.T@(M@e0)))
    if abs(m2) < 1e-14:
        scorer = A - sigma(W, D, 0.0)
        rep("inf spec (A - Sigma(0))", np.sort(np.linalg.eigvalsh(scorer))[0])
    return dict(infM=evM[0], firstNZ=evM[1], gamma=evD[0])

if __name__ == "__main__":
    print("############ SRF 2.0 FREE-FIELD ACID TEST -- EXACT ############")
    # ---- (1) FREE MASSLESS: the decisive case ----
    r_massless = run_1d(L=12, b=3, m2=0.0, tag="FREE MASSLESS")
    # ---- (3) FREE MASSIVE ----
    m2 = 0.5
    r_massive = run_1d(L=12, b=3, m2=m2, tag="FREE MASSIVE")
    # extra block sizes for gamma_fib dependence
    print("\n--- fiber gap vs block size (massless, L=12) ---")
    for b in (2,3,4):
        M = lap1d(12, True); C = block_map_1d(12, b)
        A,D,W,_,_ = split_and_blocks(M, C)
        print(f"    b={b}: gamma_fib={np.sort(np.linalg.eigvalsh(D))[0]:.6f}   "
              f"sharp Poincare const 1/gamma={1.0/np.sort(np.linalg.eigvalsh(D))[0]:.6f}   "
              f"single-block spurious 2-2cos(pi/b)={2-2*np.cos(np.pi/b):.6f}")
    # ---- (6) L-scaling: no L-uniform free gap ----
    print("\n--- first nonzero mode of free massless M: L-scaling (1D periodic) ---")
    for L in (12,24,48,96):
        ev = np.sort(np.linalg.eigvalsh(lap1d(L, True)))
        print(f"    L={L:3d}: first nonzero={ev[1]:.8f}   4 sin^2(pi/L)={4*np.sin(np.pi/L)**2:.8f}   "
              f"4pi^2/L^2={4*np.pi**2/L**2:.8f}")
    print("\n--- first nonzero mode of free massless M: L-scaling (2D LxL) ---")
    for L in (4,6,8):
        ev = np.sort(np.linalg.eigvalsh(lap2d(L,L,True)))
        print(f"    L={L:2d}: first nonzero={ev[1]:.8f}   4 sin^2(pi/L)={4*np.sin(np.pi/L)**2:.8f}")
    # ---- (2D) massless decisive + massive ----
    r2d0 = run_2d(4, 2, 0.0, "2D FREE MASSLESS")
    r2d1 = run_2d(4, 2, 0.5, "2D FREE MASSIVE")

    # ---- summary ----
    print("\n############ SUMMARY (exact inf-spec) ############")
    print(f"  FREE MASSLESS : inf spec FULL M = {r_massless['infM']:.3e} (truth 0); "
          f"gamma_fib = {r_massless['gamma']:.6f} (>0, UV); first nonzero M = {r_massless['firstNZ']:.6f} -> 0 as L->inf")
    print(f"  FREE MASSIVE  : inf spec FULL M = {r_massive['infM']:.6f} (truth m^2={m2}); "
          f"gamma_fib = {r_massive['gamma']:.6f}")
    print(f"  2D MASSLESS   : inf spec FULL M = {r2d0['infM']:.3e}; gamma_fib={r2d0['gamma']:.6f}")
    print(f"  2D MASSIVE    : inf spec FULL M = {r2d1['infM']:.6f}; gamma_fib={r2d1['gamma']:.6f}")
