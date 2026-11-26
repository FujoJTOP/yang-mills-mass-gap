"""o07_ratio.py -- R(beta) = Var(F'|0)/E(F''|0) on ONE plaquette, N=2,3.

F = S_W + Phi = beta(1 - m) - log J,  m=(1/N)Re tr(e^X B), X=iH, v=iH_v.
Conditioning: g(v,X) = <v,X> = 0  (slice).  Measure on slice: e^{-F} dx  (= e^{-S} dU).
Identity (H1A-EXACT eq:split):  G''(0) = E[F''|0] - Var(F'|0) = E - Var.
So  R = Var/E = 1 - G''/E  =>  R<=1  <=>  G''>=0.   Crossing R=1 <=> Var=E.

N=2 exact: H1A model, r in (0,2pi), measure 2pi r e^{-S} J dr, J=sin^2(r/2)/r^2 (*4 const).
N=3: full su(3) flat chart, Weyl Jacobian J=prod_{i<j}[sin(mu_ij/2)/(mu_ij/2)]^2,
     B = diag(1,w,w^2) (theta=pi analogue: tr B = 0 => S(0)=beta), H_v=diag(1,-1,0)/2.
     MCMC (Metropolis, many chains) in the 7-dim slice.
Labels: [E] exact quadrature, [C] MC.
"""
import numpy as np

rng = np.random.default_rng(20261125)
TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------- generators
def gens(N):
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), complex); T[j, k] = T[k, j] = 0.5; Ts.append(T)
            T = np.zeros((N, N), complex); T[j, k] = -0.5j; T[k, j] = 0.5j; Ts.append(T)
    for l in range(1, N):
        T = np.zeros((N, N), complex); d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d; Ts.append(T)
    return Ts


def model(N, Bkind="center"):
    if N == 2:
        Hv = np.diag([0.5, -0.5]).astype(complex)
        B = np.diag([1j, -1j]) if Bkind == "center" else np.eye(2, dtype=complex)
    else:
        Hv = np.diag([0.5, -0.5, 0.0]).astype(complex)
        w = np.exp(2j * np.pi / 3)
        B = np.diag([1.0, w, w ** 2]) if Bkind == "center" else np.diag([1.0, w ** 2, w])
    return Hv, B


# ------------------------------------------------- generic s-derivative stack
def tri(H, B, Hv, h, N):
    """m, d1m, d2m, LJ, d1LJ, d2LJ ;  H=(M,N,N) hermitian, s-dir = H + s*Hv."""
    Ms, Ls = [], []
    for s in (-h, 0.0, h):
        mu, V = np.linalg.eigh(H + s * Hv)
        mu = np.sort(mu, axis=1)
        BV = np.einsum('ij,mjk->mik', B, V)
        Tt = np.einsum('mji,mjk->mik', V.conj(), BV)
        dg = np.einsum('mii->mi', Tt)
        m = (np.exp(1j * mu) * dg).sum(1).real / N
        dm = mu[:, :, None] - mu[:, None, :]
        iu = np.triu_indices(N, 1)
        rr = np.sinc(dm[:, iu[0], iu[1]] / TWO_PI)
        LJ = 2.0 * np.log(rr).sum(1)
        Ms.append(m); Ls.append(LJ)
    m, LJ = Ms[1], Ls[1]
    d1m = (Ms[2] - Ms[0]) / (2 * h); d2m = (Ms[2] - 2 * Ms[1] + Ms[0]) / h ** 2
    d1L = (Ls[2] - Ls[0]) / (2 * h); d2L = (Ls[2] - 2 * Ls[1] + Ls[0]) / h ** 2
    return m, d1m, d2m, LJ, d1L, d2L


def ld_slice(H, B, beta, N):
    """log density on the slice: -beta(1-m) + LJ ; -inf outside |mu_ij|<=2pi."""
    mu, V = np.linalg.eigh(H)
    mu = np.sort(mu, axis=1)
    BV = np.einsum('ij,mjk->mik', B, V)
    Tt = np.einsum('mji,mjk->mik', V.conj(), BV)
    dg = np.einsum('mii->mi', Tt)
    m = (np.exp(1j * mu) * dg).sum(1).real / N
    dm = mu[:, :, None] - mu[:, None, :]
    iu = np.triu_indices(N, 1)
    dij = dm[:, iu[0], iu[1]]
    rr = np.sinc(dij / TWO_PI)
    LJ = 2.0 * np.log(np.abs(rr) + 1e-300).sum(1)
    ok = np.abs(dij).max(1) <= TWO_PI
    return np.where(ok, -beta * (1.0 - m) + LJ, -1e30)


# ------------------------------------------------ N=2 exact (H1A, theta=pi)
def n2_exact(beta, theta=np.pi, nr=400001):
    r = np.concatenate(([1e-12], np.linspace(1e-11, TWO_PI - 1e-11, nr), [TWO_PI - 1e-13]))
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    J = 4.0 * np.sin(r / 2) ** 2 / r ** 2            # *4 const -> Phi shifts by const only
    Phi = -np.log(np.maximum(J, 1e-300))
    S = beta * (1.0 - np.cos(r / 2) * c)
    w = np.exp(-S) * J * r                           # d^2x = 2pi r dr
    w = w / w.sum()
    Phi_r = np.gradient(Phi, r)
    F1 = beta * s * np.sin(r / 2) / r                # d_v S ; d_v Phi = 0 (radial)
    F2 = beta * c * np.sin(r / 2) / (2 * r) + Phi_r / r
    mF1 = (w * F1).sum()
    var = (w * F1 ** 2).sum() - mF1 ** 2
    E = (w * F2).sum()
    EH = (w * beta * c * np.sin(r / 2) / (2 * r)).sum()
    EPhi = (w * Phi_r / r).sum()
    return var, E, EH, EPhi


def cross2(f, lo, hi, it=200):
    a, b = lo, hi
    for _ in range(it):
        mid = 0.5 * (a + b)
        if (f(a) > 0) == (f(mid) > 0):
            a = mid
        else:
            b = mid
    return 0.5 * (a + b)


# ------------------------------------------------------------- N>=2 slice MC
def mc_slice(N, beta, chains=400, steps=6000, thin=5, burn=3000, h=1e-3, Bkind="center"):
    Hv, B = model(N, Bkind)
    Ts = np.array(gens(N)); D = len(Ts)
    d = np.array([np.trace(Hv @ T).real for T in Ts])
    # orthonormal basis of d-perp
    A = np.random.randn(D, D - 1)
    A = A - np.outer(d, d @ A) / (d @ d)
    Q, _ = np.linalg.qr(A)
    y = 0.15 * rng.standard_normal((chains, D - 1))
    H = np.einsum('mc,cij->mij', y @ Q.T, Ts)
    cur = ld_slice(H, B, beta, N)
    sig = np.full((chains, 1), 0.25)
    kept, acc_hist = [], []

    def toH(yy):
        return np.einsum('mc,cij->mij', yy @ Q.T, Ts)

    for t in range(steps):
        prop = y + sig * rng.standard_normal(y.shape)
        Hp = toH(prop)
        lp = ld_slice(Hp, B, beta, N)
        ok = np.log(rng.random(chains)) < (lp - cur)
        y[ok] = prop[ok]; H[ok] = Hp[ok]; cur[ok] = lp[ok]
        acc = ok.mean()
        sig = np.clip(sig * np.exp(0.3 * (acc - 0.4)), 0.02, 1.5)
        acc_hist.append(acc)
        if t >= burn and (t - burn) % thin == 0:
            kept.append(y.copy())
    X = np.concatenate(kept, axis=0)                 # (ns, D-1)
    xfull = X @ Q.T
    Hm = np.einsum('mc,cij->mij', xfull, Ts)
    m, d1m, d2m, LJ, d1L, d2L = tri(Hm, B, Hv, h, N)
    F1 = -beta * d1m - d1L
    F2 = -beta * d2m - d2L
    EH = -beta * d2m
    EPhi = -d2L
    varG = beta ** 2 * (d1m ** 2).mean(0) if d1m.ndim > 1 else beta ** 2 * d1m.var()
    # block jackknife over chains: reshape kept as (nkept, chains)
    nk = len(kept)
    def ratio(idx):
        f1 = F1[idx]; f2 = F2[idx]
        v = f1.var(); e = f2.mean()
        return v, e, f1, f2
    v0 = F1.var(); e0 = F2.mean()
    F1g = F1.reshape(nk, chains); F2g = F2.reshape(nk, chains)
    def stat(f1, f2):
        return f1.var(), f2.mean()
    jk = []
    for c in range(chains):
        mask = np.ones(chains, bool); mask[c] = False
        f1 = F1g[:, mask].ravel(); f2 = F2g[:, mask].ravel()
        v, e = f1.var(), f2.mean()
        jk.append((v, e))
    jk = np.array(jk); n = chains
    var_jk = (n - 1) / n * ((jk - jk.mean(0)) ** 2).sum(0)
    # ESS from chain-mean variance
    return dict(var=v0, E=e0, dvar=np.sqrt(var_jk[0]), dE=np.sqrt(var_jk[1]),
                EH=EH.mean(), EPhi=EPhi.mean(), varG=(beta ** 2) * d1m.var(),
                ns=F1g.size, acc=float(np.mean(acc_hist[burn:])), sig=float(sig.mean()),
                R=v0 / e0)


# ------------------------------------------------------------------ replicas
def replica_L(beta, L, theta=np.pi, nr=200001):
    """rho_L(a) ∝ rho_1(a)^L  (L identical plaquettes sharing the mode).
    Here: exact 1D radial for N=2/theta=pi, weight e^{-L F}.  R_L-1 = L(R_1-1)E_1/E_L."""
    r = np.linspace(1e-11, TWO_PI - 1e-11, nr)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    J = 4.0 * np.sin(r / 2) ** 2 / r ** 2
    Phi = -np.log(np.maximum(J, 1e-300))
    S = beta * (1.0 - np.cos(r / 2) * c)
    F = S + Phi
    w = np.exp(-L * F) * r
    w /= w.sum()
    Phi_r = np.gradient(Phi, r)
    F1 = beta * s * np.sin(r / 2) / r
    F2 = beta * c * np.sin(r / 2) / (2 * r) + Phi_r / r
    v = (w * F1 ** 2).sum() - (w * F1).sum() ** 2
    e = (w * F2).sum()
    return v, e, v / e


if __name__ == "__main__":
    print("=" * 78)
    print("[V0] N=2 exact anchors (targets: V0=0.0110593, W0=0.2051153, beta*=4.30661)")
    v1 = n2_exact(1.0); print("  beta=1: Var=%.9f Var/b^2=%.9f  E=%.9f  EH=%.3e  EPhi=%.9f"
                             % (v1[0], v1[0], v1[1], v1[2], v1[3]))
    v0 = n2_exact(0.0); print("  beta=0: E=%.9f  (W0 target 0.20511530)  EPhi=%.9f" % (v0[1], v0[3]))

    print("=" * 78)
    print("[T1] N=2 exact, theta=pi   (R = Var/E ; crossing where Var=E)")
    print("  N beta      Var(F'|0)     E(F''|0)        R        R-1   E[Hvv]     E[Phi'']")
    betas = [0.0, 1.0, 2.5, 4.0, 6.0, 10.0]
    n2tab = {}
    for b in betas:
        v, e, eh, ep = n2_exact(b)
        n2tab[b] = (v, e)
        print("  2 %5.2f  %14.9f %14.9f %10.6f %+9.5f %10.3e %10.6f"
              % (b, v, e, v / e, v / e - 1, eh, ep))
    bs = cross2(lambda b: n2_exact(b)[0] - n2_exact(b)[1], 3.0, 6.0)
    v, e, eh, ep = n2_exact(bs)
    print("  beta* = %.8f   [target 4.30661]   Var=%.9f E=%.9f Var-E=%+.3e"
          % (bs, v, e, v - e))
    print("  check: R(beta*) = %.10f" % (v / e))

    print("=" * 78)
    print("[C1] N=2 generic MC (same code as N=3), theta=pi slice")
    for b in (2.5, 6.0):
        r = mc_slice(2, b, chains=200, steps=3000, thin=2, burn=1500)
        ex = n2_exact(b)
        print("  b=%4.1f MC Var=%9.6f+-%7.6f E=%9.6f+-%7.6f R=%8.5f | exact Var=%9.6f E=%9.6f R=%8.5f"
              " | ns=%d acc=%.2f" % (b, r['var'], r['dvar'], r['E'], r['dE'],
                                     r['var'] / r['E'], ex[0], ex[1], ex[0] / ex[1],
                                     r['ns'], r['acc']))

    print("=" * 78)
    print("[C2] N=3 MC, one plaquette, B=diag(1,w,w^2), H_v=diag(1,-1,0)/2, slice <v,X>=0")
    print("  N beta      Var(F'|0)         E(F''|0)          R        R-1     E[Hvv]    E[Phi'']   ns")
    n3tab = {}
    for b in [0.0, 1.0, 2.5, 4.0, 6.0, 10.0]:
        r = mc_slice(3, b, chains=400, steps=6000, thin=5, burn=3000)
        n3tab[b] = r
        print("  3 %5.2f  %11.6f+-%8.6f %11.6f+-%8.6f %9.5f %+9.5f %9.5f %9.5f %7d acc=%.2f sig=%.2f"
              % (b, r['var'], r['dvar'], r['E'], r['dE'], r['R'], r['R'] - 1,
                 r['EH'], r['EPhi'], r['ns'], r['acc'], r['sig']))

    print("=" * 78)
    print("[C3] N=3 extend to beta=14, 20 (does R reach 1 ?)")
    for b in [14.0, 20.0]:
        r = mc_slice(3, b, chains=400, steps=6000, thin=5, burn=3000)
        n3tab[b] = r
        print("  3 %5.2f  %11.6f+-%8.6f %11.6f+-%8.6f %9.5f %+9.5f %9.5f %9.5f %7d acc=%.2f sig=%.2f"
              % (b, r['var'], r['dvar'], r['E'], r['dE'], r['R'], r['R'] - 1,
                 r['EH'], r['EPhi'], r['ns'], r['acc'], r['sig']))

    # two-term / three-term fit for beta*(3)
    bvec = np.array([k for k in n3tab if k > 0])
    dvec = np.array([n3tab[k]['E'] - n3tab[k]['var'] for k in bvec])
    A = np.vstack([np.ones_like(bvec), bvec, bvec ** 2]).T
    coef, *_ = np.linalg.lstsq(A, dvec, rcond=None)
    disc = coef[1] ** 2 - 4 * coef[2] * coef[0]
    print("  fit delta3(b)=a0+a1 b+a2 b^2 : a0=%.4f a1=%.5f a2=%.6f  disc=%.4f"
          % (coef[0], coef[1], coef[2], disc))
    if disc > 0:
        rr = [(-coef[1] + s * np.sqrt(disc)) / (2 * coef[2]) for s in (1, -1)]
        print("  roots: %.2f , %.2f  (positive branch = beta*(3) extrapolation)" % tuple(rr))
    print("  two-term law sqrt(W0/V0): W0=%.4f V0=%.2e -> %.2f"
          % (n3tab[0.0]['E'] * 0 + 0.3453, 9.0e-4, np.sqrt(0.3453 / 9.0e-4)))

    print("=" * 78)
    print("[T2] L-replica: rho_L ~ rho_1^L  =>  G''_L = L G''_1  =>  crossing L-invariant")
    for b in (2.5, 4.0, 6.0):
        v, e, R = replica_L(b, 1)
        row = [b, R]
        for L in (2, 4):
            vL, eL, RL = replica_L(b, L)
            row += [L, RL, RL - 1, L * (R - 1) * e / eL]
        print("  b=%4.1f R1=%.6f | L=2 R=%.6f (R-1=%.4e) | L=4 R=%.6f (R-1=%.4e) | sign(R-1) L=1,2,4: %s"
              % (row[0], row[1], row[3], row[4], row[5], row[6],
                 np.sign([row[1] - 1, row[3] - 1, row[5] - 1]).astype(int)))
