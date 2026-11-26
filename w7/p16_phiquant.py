"""
P16 - quantify  Phi'' >= N/12 + c a0^2.   py -3.12, numpy only, deterministic.
Objects verbatim from w7/o09_phicurv.py (= O09-PHICURV.md S1) and O06-CURV.md S4.4:
    Phi = sum_{i<j} g(lam_i-lam_j),   g(d) = -2 log( sin(d/2)/(d/2) )
    g'(d)  = 2/d - cot(d/2)  = d/6 + d^3/360 + ...            >= 0      on (0,2pi)
    g''(d) = csc^2(d/2)/2 - 2/d^2 = 1/6 + d^2/120 + ...       >= 1/6    on (0,2pi)
    flat chart |lam_i-lam_j| < 2pi,  X = diag(lam) WLOG (Phi is Ad-invariant)
    Phi''(X;v) = sum_{i<j}[ g''(d_ij)(W_ii-W_jj)^2 + g'(d_ij)(R_ii-R_jj) ],
        W = U^dag v U,  R_ii = delta^2 lam_i = 2 sum_{k!=i} |W_ik|^2/(lam_i-lam_k)
    normalization <v,v> := 2 tr v^2,  v=T_3=diag(1/2,-1/2,..) -> Phi''(0)=N/12
    phi := g - d^2/12  (>=0, even, phi''=g''-1/6>=0) ;  h := Phi - (N/12) tr X^2 = sum_{i<j} phi(lam_i-lam_j)
    CLAIM (this file): Phi'' >= (N/12)<v,v>  <=>  h convex  <=> Hess Phi >= (N/6) tr(.)^2
"""
import numpy as np

TWO_PI = 2.0 * np.pi
LOG = open(r"D:\Dev\math\w7\p16_phiquant.log", "w", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


np.set_printoptions(precision=8, suppress=True)


# ---------------- g, g', g'' (series-safe |d|<1e-2) ----------------
def _prep(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-2
    return d, small, np.where(small, 1.0, d)


def g(d):
    d, small, safe = _prep(d)
    with np.errstate(all="ignore"):
        out = -2.0 * np.log(np.sin(safe / 2) / (safe / 2))
    return np.where(small, d ** 2 / 12 + d ** 4 / 720 + d ** 6 / 45360 + d ** 8 / 2419200, out)


def gp(d):
    d, small, safe = _prep(d)
    with np.errstate(all="ignore"):
        out = 2.0 / safe - 1.0 / np.tan(safe / 2)
    return np.where(small, d / 6 + d ** 3 / 360 + d ** 5 / 15120, out)


def gpp(d):
    d, small, safe = _prep(d)
    with np.errstate(all="ignore"):
        out = 0.5 / np.sin(safe / 2) ** 2 - 2.0 / safe ** 2
    return np.where(small, 1.0 / 6 + d ** 2 / 120 + d ** 4 / 3024 + d ** 6 / 86400, out)


def phi(d):
    return g(d) - d ** 2 / 12.0


# ---------------- exact Phi'' for X = diag(lam), direction W (same basis) ----------------
def phi2(lam, W):
    n = lam.size
    M = lam[:, None] - lam[None, :]
    off = ~np.eye(n, dtype=bool)
    Mo = np.where(off, M, 1.0)
    R = 2.0 * np.sum((np.abs(W) ** 2) * (off / Mo), axis=1)          # R_ii
    w = np.where(off, gp(M), 0.0).sum(axis=1)                         # w_i = sum_j g'(lam_i-lam_j)
    iu = np.triu_indices(n, 1)
    dg = np.real(np.diag(W))
    t1 = np.sum(gpp(M[iu]) * (dg[iu[0]] - dg[iu[1]]) ** 2)
    return float(t1 + np.sum(w * R))


def su_gens(N):
    G = []
    for i in range(N):
        for j in range(i + 1, N):
            A = np.zeros((N, N), complex); A[i, j] = 1; A[j, i] = 1; G.append(A / 2)
            B = np.zeros((N, N), complex); B[i, j] = -1j; B[j, i] = 1j; G.append(B / 2)
    for k in range(1, N):
        D = np.zeros((N, N), complex)
        for i in range(k):
            D[i, i] = 1.0
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        G.append(D)
    return G


def hess_min(lam, Gb):
    """min/max eigenvalue of Hess Phi at X=diag(lam), orthonormal basis Gb (<A,B>=2tr AB)."""
    m = len(Gb)
    Q = np.array([phi2(lam, B) for B in Gb])
    H = np.diag(Q).astype(float)
    for a in range(m):
        for b in range(a + 1, m):
            v = 0.5 * (phi2(lam, Gb[a] + Gb[b]) - Q[a] - Q[b])
            H[a, b] = H[b, a] = v
    ev = np.linalg.eigvalsh(H)
    return ev[0], ev[-1]


def fd_phi2(lam, W, h=1e-4):
    def F(t):
        return float(g((lam[:, None] - lam[None, :] + t * (np.diag(W)[:, None] - np.diag(W)[None, :]))[np.triu_indices(lam.size, 1)]).sum())
    # 2nd-order FD along the full matrix direction is only valid for diagonal W; use exact for the rest
    return (-F(2 * h) + 16 * F(h) - 30 * F(0) + 16 * F(-h) - F(-2 * h)) / (12 * h * h)


say("=" * 104)
say("[1] derivative sanity:  g'' >= 1/6,  g'/d >= 1/6,  phi = g - d^2/12 >= 0,  phi'' = g''-1/6 >= 0")
dd = np.linspace(1e-3, TWO_PI - 1e-3, 300001)
g2, g1, ph, ph2 = gpp(dd), gp(dd), phi(dd), gpp(dd) - 1.0 / 6
say(f"  min g''            = {g2.min():.12f} at d={dd[g2.argmin()]:.6f}   (1/6 = {1/6:.12f})")
say(f"  min g'/d           = {(g1/dd).min():.12f} at d={dd[(g1/dd).argmin()]:.6f}")
say(f"  min phi = g-d^2/12 = {ph.min():.3e}   min phi'' = g''-1/6 = {ph2.min():.3e}  (>=0)")
say(f"  monotone steps of g'' : up={int((np.diff(g2)>0).sum())}/{dd.size-1}  (increasing on (0,2pi))")
say(f"  phi ~ d^4/720 + ...  : 720*phi(d)/d^4 at d=1e-2 = {720*phi(1e-2)/1e-8:.9f}, d=1e-1 = {720*phi(1e-1)/1e-4:.9f}")

say("=" * 104)
say("[2] exact Phi'' formula vs 5-pt FD (diagonal directions), N=2,3,4, random X")
rng = np.random.default_rng(7)
for N in (2, 3, 4):
    err = []
    for _ in range(5):
        R = rng.uniform(0.2, 4.0)
        lam = rng.uniform(-1, 1, N); lam = (lam - lam.mean()) / max(abs(lam.max() - lam.min()), 1e-9) * R
        Wd = np.diag(rng.normal(size=N)); Wd = Wd - np.trace(Wd) / N * np.eye(N)
        v = phi2(lam, Wd.astype(complex))
        f = fd_phi2(lam, Wd)
        err.append(abs(v - f) / max(abs(v), 1e-12))
    say(f"  N={N}: max rel err exact vs FD = {max(err):.2e}")

say("=" * 104)
say("[3] DECISIVE: min over ALL directions W of Phi''/ <v,v>/(N/12)  (X = diag(lam), Ad-invariance -> WLOG)")
say("    spectra battery: equally spaced (R = total spread), two-cluster (near-degenerate), random, tiny")
rows = []
for N in (2, 3, 4, 5, 6):
    Gb = su_gens(N)
    spec = []
    for R in (1e-6, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 3.5, 5.0, 6.0):
        lam = np.arange(N)[::-1].astype(float); lam = (lam - lam.mean()) / max(abs(lam - lam.mean()).max(), 1e-12) * (R / 2)
        spec.append((f"equi R={R:g}", lam))
    for R in (0.5, 2.0, 5.5):
        lam = np.zeros(N)
        lam[: N // 2] = R / 2; lam[N // 2:] = -R / 2
        spec.append((f"2cluster R={R:g}", lam))
    for _ in range(25):
        R = rng.uniform(0.05, 6.1)
        lam = rng.uniform(-1, 1, N)
        lam = (lam - lam.mean()) / max(abs(lam.max() - lam.min()), 1e-12) * R
        spec.append((f"rand R={R:.2f}", lam))
    best = (1e9, None, None)
    for tag, lam in spec:
        lam = lam + 1e-9 * (np.arange(N) - 0.5 * (N - 1))    # break exact degeneracies (limit is continuous)
        lmin, lmax = hess_min(lam, Gb)
        if lmin < best[0]:
            best = (lmin, tag, lam.copy())
    ratio = best[0] / (N / 12.0)
    rows.append((N, ratio, best))
    say(f"  N={N}: inf over battery  min eig(Hess Phi) = {best[0]:.10f}   N/12 = {N/12:.10f}"
        f"   ratio = {ratio:.6f}   {'OK >=1' if ratio >= 1 - 1e-9 else 'BELOW FLOOR'}   at {best[1]}")
    say(f"        argmin spectrum = {np.array2string(np.sort(best[2])[::-1], precision=6)}")

say("=" * 104)
say("[4] delta^2 lam : exact formula, UNBOUNDED below, but the g'*delta^2 lam ledger is finite")
for D in (1.0, 1e-1, 1e-2, 1e-3, 1e-4):
    lam = np.array([D, 0.0, -D])
    Gb = su_gens(3)
    r = []
    for B in Gb:
        W = B
        M = lam[:, None] - lam[None, :]
        off = ~np.eye(3, dtype=bool)
        Mo = np.where(off, M, 1.0)
        R = 2.0 * np.sum((np.abs(W) ** 2) * (off / Mo), axis=1)
        nrm = 2 * np.real(np.trace(W @ W))
        r.append((R.min() / nrm, R.max() / nrm))
    r = np.array(r)
    say(f"  gap D={D:8.1e}:  min_v delta^2lam/<v,v> = {r[:,0].min():+.6e}   (-1/(2D) = {-1/(2*D):+.3e})"
        f"   max_v = {r[:,1].max():+.6e}   sum_i delta^2lam = 0 (by tr)")

say("  (2-level check, N=2, pure off-diagonal v with <v,v>=1: delta^2lam_+ = +1/(2d), delta^2lam_- = -1/(2d),")
say("   g'(d)*delta^2lam_- = -g'(d)/(2d) -> -1/12 as d->0 : FINITE ; the pair ledger 4g'(d)|W|^2/d -> 2/3)")
say("  => delta^2 lam has NO gap-independent lower bound;  the composite (g'' term + g' term) is bounded below.")

say("=" * 104)
say("[5] small-X expansion: Phi'' - N/12 = O(a0^2).  pointwise sharp coefficient c_p = inf (Phi''-N/12)/a0^2")
say("    (N=2 exact, radial formula o09 S1: Phi''-1/6 = rho^2[u^2/120+(1-u^2)/360], a0=rho u")
say("     => Phi''-1/6 = rho^2/360 + a0^2/180 >= a0^2/120  [T] exact, rho^2>=a0^2)")
for N in (2, 3, 4):
    Gb = su_gens(N)
    best = (1e9, None, None)
    for eps in (0.02, 0.05, 0.1, 0.2):
        for _ in range(120):
            x = rng.normal(size=N); x = x - x.mean()
            x = x / np.linalg.norm(x) * eps * np.sqrt(N / (N - 1.0))
            lam = np.sort(x)[::-1]
            dirs = list(Gb)
            for _t in range(8):                                   # random traceless Hermitian directions
                A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
                vv = A + A.conj().T
                vv = vv - np.trace(vv) / N * np.eye(N)
                dirs.append(vv / np.sqrt(2 * np.real(np.trace(vv @ vv))))
            al = np.diag(lam).astype(complex)
            dirs.append(al / np.sqrt(2 * np.sum(lam ** 2)))        # aligned with X (a0 maximal)
            for B in dirs:
                a0 = abs(2.0 * np.sum(lam * np.real(np.diag(B))))     # a0 = <v,X>, <A,B> := 2tr AB
                if a0 < 1e-8:
                    continue
                q = (phi2(lam, B) - N / 12.0) / a0 ** 2
                if q < best[0]:
                    best = (q, eps, np.sort(lam)[::-1])
    say(f"  N={N}: min (Phi''-N/12)/a0^2 over small-X battery = {best[0]:.8f}  (at eps={best[1]:g}, "
        f"spectrum {np.array2string(best[2], precision=4)})")
say(f"  analytic aligned-direction value (1/120)*sum d^4/|X|^4 (N=2 -> exactly 1/120 = {1/120:.8f})")

say("=" * 104)
say("[6] conditioned slice curve E[Phi''|a0] : SU(2) single plaquette, H1A-EXACT model (as o09 [4])")
RMAX = TWO_PI - 1e-3
NS = 2000
HH = 1e-4
_GL = {}


def gl(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = (0.5 * (x + 1.0), 0.5 * w)
    return _GL[n]


def slice_a(a, beta, theta):
    if abs(a) >= RMAX - 1e-5:
        return None
    smax = np.sqrt(max(RMAX ** 2 - a ** 2, 0.0))
    nd, wt = gl(NS)
    s = smax * nd
    wgt = smax * wt
    c, sth = np.cos(theta / 2), np.sin(theta / 2)

    def fields(x3):
        rho = np.maximum(np.sqrt(s ** 2 + x3 ** 2), 1e-12)
        u = x3 / rho
        sig = 1.0 - np.cos(rho / 2) * c + np.sin(rho / 2) * sth * u
        Phi = -2.0 * np.log(np.sin(rho / 2) / (rho / 2))
        return sig, Phi

    sig, Phi = fields(a)
    sp, Pp = fields(a + HH)
    sm, Pm = fields(a - HH)
    s2p, P2p = fields(a + 2 * HH)
    s2m, P2m = fields(a - 2 * HH)
    dPhi = (Pp - Pm) / (2 * HH)
    ddPhi = (-P2p + 16 * Pp - 30 * Phi + 16 * Pm - P2m) / (12 * HH ** 2)
    w = np.exp(-(beta * sig + Phi)) * s * wgt
    Z = w.sum()
    m = lambda q: float((w * q).sum() / Z)
    return dict(Phi2=m(ddPhi), varPhi=m(dPhi ** 2) - m(dPhi) ** 2, varF=m(dPhi ** 2) - m(dPhi) ** 2)


a_grid = np.linspace(-2.0, 2.0, 17)
for beta in (2.5, 4.30661):
    for th, tag in ((np.pi, "theta=pi"), (0.0, "theta=0 ")):
        vals = np.array([slice_a(float(a), beta, th)["Phi2"] for a in a_grid])
        i0 = int(np.argmin(np.abs(a_grid)))
        A0 = vals[i0]
        msk = np.abs(a_grid) <= 1.5
        c_fit = np.sum((a_grid[msk] ** 2) * (vals[msk] - A0)) / np.sum((a_grid[msk] ** 2) ** 2)
        j = int(np.argmin(vals))
        say(f"  beta={beta:.5f} {tag}: E[Phi''|0] = {A0:.6f}   min over a-grid = {vals.min():.6f} at a = {a_grid[j]:+.3f}"
            f"   argmin |a| = {a_grid[i0]:+.2f}")
        say(f"      fitted  E[Phi''] = {A0:.5f} + c_eff a^2 ,  c_eff = {c_fit:+.6f}   (N=2: N/12 = {1/6:.6f})")
        s = "      curve a: " + " ".join(f"{a:+.2f}:{v:.5f}" for a, v in zip(a_grid, vals))
        say(s)

say("=" * 104)
say("[7] comparison with the o06 multi-link gap  -beta E[C_geo]  (o06 log [4], beta=2.5, sum_l ||v_l||^2 = 2)")
say("  min H_vv = -1.893109  => per unit ||v||^2 : -0.9466     (o06 [4])")
say(f"  Phi'' floor per unit ||v||^2 : +N/12 = +0.1667 ; Phi'' GROWTH at |a0|<=1.5 : <= +0.08 (see [6])")
say("  => Phi'' growth is O(a0^2) with c_eff ~ 3.5e-2 and saturates at the flat-chart edge;  beta E[C_geo] is O(beta)")
say("     and unbounded  =>  the growth CANNOT cancel it (ratio ~ 10x at beta=2.5, worsens with beta).")
say("=" * 104)
LOG.close()
