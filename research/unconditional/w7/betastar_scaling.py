"""
β*(N) for the single-plaquette H1 margin.  numpy only, deterministic where possible.

Model (canonical SU(N) lift of H1A-EXACT.md §2.1, N=2 is *identical* to H1A):
  X = sum_{a off-Cartan} xi^a T^a   (Cartan-transverse slice == conditioning <v,xi>=0 for all Cartan v)
  U = exp(iX),   S = beta (1 - ReTr(U B)/N),   B = centre element (det=1 phase fixed); N=2 -> diag(-i,i) ~ H1A's i*sigma_3
  measure e^{-S} J dX,  J = prod_{i<j} [sin(mu_ij/2)/(mu_ij/2)]^2 ,  mu = eig(X)   (de Wit - 't Hooft)
  domain |mu_i - mu_j| <= 2 pi   (injectivity domain of exp; N=2 -> |xi| <= 2 pi, H1A's cut)
  v = T_3 = block sigma_3/2
  delta_max = <H>_w - Var(G)_w + <Phi''>_w - Var(Phi')_w - 2 Cov(G,Phi')_w     [H1A §1.3 identity]
Everything except the beta factor is beta-independent  ->  build one pool of samples, scan beta for free.
"""
import numpy as np
import time

np.set_printoptions(precision=10)


def su_gens(N):
    G = []
    for i in range(N):
        for j in range(i + 1, N):
            A = np.zeros((N, N), complex); A[i, j] = 1; A[j, i] = 1; G.append(A / 2)
            Bm = np.zeros((N, N), complex); Bm[i, j] = -1j; Bm[j, i] = 1j; G.append(Bm / 2)
    for k in range(1, N):
        D = np.zeros((N, N), complex)
        for i in range(k):
            D[i, i] = 1.0
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        G.append(D)
    return G


def centre_B(N, conj=False):
    w = np.exp(2j * np.pi / N)
    d = np.array([w ** k for k in range(N)], complex) * np.exp(-1j * np.pi * (N - 1) / N)
    B = np.diag(d)
    return B.conj() if conj else B


def pool(N, nwant, seed=1, chunk=200000, conj=False, h=1e-3, verbose=True):
    """uniform-in-box + rejection on the domain; per-sample beta-independent factors."""
    rng = np.random.default_rng(seed)
    G = su_gens(N)
    off = G[:N * (N - 1)]
    Tv = G[N * (N - 1)]
    B = centre_B(N, conj=conj)
    A = 2 * np.pi
    iu = np.triu_indices(N, 1)
    store = {k: [] for k in ('J', 'Phi', 'dPhi', 'ddPhi', 'f', 'Ag', 'Ch')}
    got = 0
    tried = 0
    while got < nwant:
        xi = rng.uniform(-A, A, size=(chunk, N * (N - 1)))
        tried += chunk
        X = np.zeros((chunk, N, N), complex)
        for a, T in enumerate(off):
            X += xi[:, a, None, None] * T[None]
        mu, V = np.linalg.eigh(X)
        d = (mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]]
        keep = np.all(np.abs(d) <= 2 * np.pi, axis=1)
        if not keep.any():
            continue
        X, mu, V, d = X[keep], mu[keep], V[keep], d[keep]
        half = d / 2
        J = np.prod((np.sin(half) / half) ** 2, axis=1)
        Phi = -np.log(J)
        e = np.exp(1j * mu)
        U = (V * e[:, None, :]) @ V.conj().transpose(0, 2, 1)
        f = np.real(np.trace(U @ B, axis1=1, axis2=2)) / N
        # NB: d/dxi exp(iX) != i T exp(iX) unless [T,X]=0 -> all v-derivatives by finite differences
        Pp = np.empty_like(Phi); Pm = np.empty_like(Phi)
        fp = np.empty_like(Phi); fm = np.empty_like(Phi)
        for sgn, dst, fdst in ((+1, Pp, fp), (-1, Pm, fm)):
            Xs = X + sgn * h * Tv
            mu2, V2 = np.linalg.eigh(Xs)
            d2 = (mu2[:, :, None] - mu2[:, None, :])[:, iu[0], iu[1]] / 2
            dst[:] = -np.log(np.prod((np.sin(d2) / d2) ** 2, axis=1))
            U2 = (V2 * np.exp(1j * mu2)[:, None, :]) @ V2.conj().transpose(0, 2, 1)
            fdst[:] = np.real(np.trace(U2 @ B, axis1=1, axis2=2)) / N
        dPhi = (Pp - Pm) / (2 * h)
        ddPhi = (Pp - 2 * Phi + Pm) / h ** 2
        Ag = (fp - fm) / (2 * h)                       # d_v f   ->  G = -beta*Ag
        Ch = (fp - 2 * f + fm) / h ** 2                # d_v^2 f ->  H = -beta*Ch
        m = min(nwant - got, len(J))
        store['J'].append(J[:m]); store['Phi'].append(Phi[:m]); store['dPhi'].append(dPhi[:m])
        store['ddPhi'].append(ddPhi[:m]); store['f'].append(f[:m]); store['Ag'].append(Ag[:m])
        store['Ch'].append(Ch[:m])
        got += m
    P = {k: np.concatenate(v) for k, v in store.items()}
    P['acc'] = got / tried
    if verbose:
        print(f"  N={N}: {got} valid / {tried} draws  (acc {P['acc']:.2e})", flush=True)
    return P


def delta(P, beta, nblk=25):
    f, J, Ag, Ch = P['f'], P['J'], P['Ag'], P['Ch']
    logw = -beta * (1 - f) + np.log(J)
    logw -= logw.max()
    n = len(logw)
    B_ = nblk
    idx = np.array_split(np.arange(n), B_)
    vals = dict(H=-beta * Ch, G=-beta * Ag, P2=P['ddPhi'], D=P['dPhi'])

    def avg(w, v):
        return (w * v).sum() / w.sum()

    w = np.exp(logw)
    res = {}
    subs = {k: np.array([avg(w[ix], vals[k][ix]) for ix in idx]) for k in vals}
    for k in vals:
        mu = avg(w, vals[k])
        err = np.sqrt((B_ - 1) / B_ * ((subs[k] - subs[k].mean()) ** 2).sum())
        res[k] = (mu, err)
    H, eH = res['H']; Gm, eG = res['G']
    varG, eVG = res['G2'] if 'G2' in res else (0, 0), 0
    # second moments
    subs2 = {k: np.array([avg(w[ix], vals[k][ix] ** 2) for ix in idx]) for k in ('G', 'D')}
    mG2 = avg(w, vals['G'] ** 2); eG2 = np.sqrt((B_ - 1) / B_ * ((subs2['G'] - subs2['G'].mean()) ** 2).sum())
    mD2 = avg(w, vals['D'] ** 2); eD2 = np.sqrt((B_ - 1) / B_ * ((subs2['D'] - subs2['D'].mean()) ** 2).sum())
    mGD = avg(w, vals['G'] * vals['D'])
    subGD = np.array([avg(w[ix], (vals['G'] * vals['D'])[ix]) for ix in idx])
    eGD = np.sqrt((B_ - 1) / B_ * ((subGD - subGD.mean()) ** 2).sum())
    P2, eP2 = res['P2']; Dm, eD = res['D']
    varG = mG2 - Gm ** 2; varD = mD2 - Dm ** 2; cov = mGD - Gm * Dm
    Wgeo = P2 - varD
    Xc = -2 * cov
    dl = H - varG + Wgeo + Xc
    e = eH + eG2 + 2 * abs(Gm) * eG + eP2 + eD2 + 2 * abs(Dm) * eD + 2 * eGD
    return dict(H=H, varG=varG, Phi2=P2, varD=varD, Wgeo=Wgeo, X=Xc, delta=dl, err=e,
                varPhiPrime=varD, cov=cov, G=Gm, D=Dm)


def n2_exact(n=400001):
    r = np.linspace(1e-9, 2 * np.pi, n); dr = r[1] - r[0]
    w = np.sin(r / 2) ** 2 / r; w = w / (w.sum() * dr)
    ddPhi = (2 / r - 1 / np.tan(r / 2)) / r
    W = (w * ddPhi).sum() * dr
    fun = np.sin(r / 2) / r
    m1 = (w * fun).sum() * dr; m2 = (w * fun ** 2).sum() * dr
    Cin = np.trapezoid((1 - np.cos(r)) / r, r)
    return W, m2 - m1 ** 2, Cin


def bstar_from(betas, deltas):
    b = np.asarray(betas, float); d = np.asarray(deltas, float)
    idx = np.where(np.diff(np.sign(d)) != 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    return b[i] - d[i] * (b[i + 1] - b[i]) / (d[i + 1] - d[i])


if __name__ == "__main__":
    t0 = time.time()
    W, V0, Cin = n2_exact()
    print("=" * 90)
    print("[0] N=2 exact (1D deterministic quadrature over the slice)")
    print(f"  Cin(2pi)=int_0^2pi (1-cos t)/t dt = {Cin:.12f}")
    print(f"  W_geo = <Phi''>_Haar          = {W:.10f}")
    print(f"  1/(2 Cin(2pi))                = {1/(2*Cin):.10f}   (closed form, rel.diff {abs(W-1/(2*Cin))/W:.1e})")
    print(f"  V0 = Var(sin(r/2)/r)          = {V0:.10f}")
    print(f"  beta*(2)                      = {np.sqrt(W/V0):.6f}")
    print(f"  doc anchors: W=0.2051153  V0=0.0110593  beta*=4.30661")
    print("=" * 90)

    pools = {}
    for N, nw in ((2, 400000), (3, 500000), (4, 200000)):
        pools[N] = pool(N, nw, seed=10 + N)
    print(f"  (pool build {time.time()-t0:.1f}s)")

    print("=" * 90)
    print("[1] delta(beta) from the pools (weight e^-S J exact, beta-free pool)")
    for N in (2, 3, 4):
        print(f"--- N={N} ---")
        for beta in (0.0, 2.0, 4.0, 4.30661, 6.0, 8.0, 10.0, 12.0, 14.0):
            r = delta(pools[N], beta)
            print(f"  b={beta:7.4f}  delta={r['delta']:+9.5f} +- {r['err']:.5f}"
                  f" | <H>={r['H']:+8.4f} -V={-r['varG']:+8.5f} <Phi''>={r['Phi2']:+.4f}"
                  f" -V(Phi')={-r['varD']:+.5f} X={r['X']:+.5f}", flush=True)
    print("=" * 90)
    print(f"elapsed {time.time()-t0:.1f}s")
