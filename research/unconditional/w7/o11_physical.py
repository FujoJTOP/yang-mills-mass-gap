"""
O11-PHYS: Delta := E[kappa|0] - Var(F|0) at beta = 2.5 (physical point).

  kappa = d_v^2 V,  F = d_v V,  V = S_W + Phi,  Phi = -log J   (flat-graph Haar Jacobian)
  Delta = <kappa|0> - Var(F|0) = G''(0)   [GPP-STRUCT.md eq:split first line, identity]

(1) single plaquette SU(N), N = 2,3 at beta=2.5:
      N=2: exact 1D quadrature over the conditioning slice (deterministic) + MC pools
      N=3: MC pools (two independent groups), weights e^{-beta(1-f)} J exact
(2) L=2 4D SU(2) Wilson lattice, beta=2.5: Metropolis on the exact conditioning
    hyperplane {<v,X>=0}, flat coordinates X_link = log U_link.
    two coarse modes: (A) k=0 uniform link direction  (B) single-link direction.

numpy only.  Run:  py -3.12 o11_physical.py
"""
import numpy as np
import time
import betastar_scaling as bs

BETA = 2.5
T0 = time.time()


def P(*a):
    print(*a, flush=True)


# =====================================================================
# (1) single plaquette
# =====================================================================
def jk(P_, beta, nblk=25):
    """Delta and its block-jackknife sigma, plus the pieces E[H], Var(G), <Phi''>, Var(Phi'), X."""
    f, J, Ag, Ch, Dd, Dp = P_['f'], P_['J'], P_['Ag'], P_['Ch'], P_['ddPhi'], P_['dPhi']
    lw = -beta * (1 - f) + np.log(np.maximum(J, 1e-300))
    lw -= lw.max()
    w = np.exp(lw)
    kap = -beta * Ch + Dd          # kappa = H + Phi''
    F = -beta * Ag + Dp            # F     = G + Phi'
    n = len(w)
    idx = np.array_split(np.arange(n), nblk)
    keep = [np.concatenate([idx[j] for j in range(nblk) if j != i]) for i in range(nblk)]

    def mom(sel):
        ww = w[sel]
        Z = ww.sum()
        mF = (ww * F[sel]).sum() / Z
        return dict(Z=Z,
                    k=(ww * kap[sel]).sum() / Z,
                    F=mF, F2=(ww * F[sel] ** 2).sum() / Z,
                    H=(ww * (-beta * Ch[sel])).sum() / Z,
                    vG=(ww * (-beta * Ag[sel]) ** 2).sum() / Z - ((ww * (-beta * Ag[sel])).sum() / Z) ** 2,
                    Dd=(ww * Dd[sel]).sum() / Z,
                    vD=(ww * Dp[sel] ** 2).sum() / Z - ((ww * Dp[sel]).sum() / Z) ** 2,
                    nD=(ww * Dp[sel]).sum() / Z,
                    mGD=(ww * (-beta * Ag[sel]) * Dp[sel]).sum() / Z)

    full = mom(slice(None))
    dl = lambda m: m['k'] - (m['F2'] - m['F'] ** 2)
    d0 = dl(full)
    ds = np.array([dl(mom(s)) for s in keep])
    sig = np.sqrt((nblk - 1) / nblk * ((ds - ds.mean()) ** 2).sum())
    return dict(delta=d0, sig=sig, ess=full['Z'] ** 2 / (w ** 2).sum(),
                H=full['H'], vG=full['vG'], Dd=full['Dd'], vD=full['vD'],
                kap=full['k'], varF=full['F2'] - full['F'] ** 2, F=full['F'])


def n2_exact_beta(beta):
    """deterministic N=2 quadrature: E[H]=0 exactly (centre B), Delta = W0 - beta^2 V0."""
    W, V0, Cin = bs.n2_exact()
    return W, V0, Cin, W - beta ** 2 * V0


# =====================================================================
# (2) L=2 4D SU(2) lattice
# =====================================================================
def lat_tables(L):
    V = L ** 4
    tau = [(L ** k) for k in range(4)]

    def site(c):
        return c[0] * tau[0] + c[1] * tau[1] + c[2] * tau[2] + c[3] * tau[3]

    def link(mu, s):
        return mu * V + s

    plaq = []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for s in range(V):
                c = [(s // tau[k]) % L for k in range(4)]
                em = list(c); em[mu] = (em[mu] + 1) % L
                en = list(c); en[nu] = (en[nu] + 1) % L
                plaq.append((link(mu, s), link(nu, site(em)), link(mu, site(en)), link(nu, s)))
    tab = np.array(plaq)
    near = [np.where((tab == l).any(1))[0] for l in range(4 * V)]
    return tab, near


def expi(x):
    r = np.linalg.norm(x, axis=1)
    c = np.cos(r / 2)
    s = np.where(r > 1e-10, np.sin(r / 2) / np.where(r > 1e-10, r, 1.0), 0.5)
    U = np.zeros((x.shape[0], 2, 2), complex)
    U[:, 0, 0] = c + 1j * s * x[:, 2]
    U[:, 0, 1] = 1j * s * (x[:, 0] - 1j * x[:, 1])
    U[:, 1, 0] = 1j * s * (x[:, 0] + 1j * x[:, 1])
    U[:, 1, 1] = c - 1j * s * x[:, 2]
    return U


def phi_of(r):
    q = np.where(r > 1e-10, np.sin(r / 2) / np.where(r > 1e-10, r / 2, 1.0), 1.0)
    return -2 * np.log(np.abs(q))


def tr_plaq(U, tab, rows):
    A = U[tab[rows, 0]]
    Bm = U[tab[rows, 1]]
    C = U[tab[rows, 2]].conj().transpose(0, 2, 1)
    D = U[tab[rows, 3]].conj().transpose(0, 2, 1)
    Qt = A @ Bm @ C @ D
    return 0.5 * np.real(Qt[:, 0, 0] + Qt[:, 1, 1])


def Vtot(x, tab, beta):
    U = expi(x)
    tr = tr_plaq(U, tab, np.arange(len(tab)))
    return beta * np.sum(1.0 - tr) + phi_of(np.linalg.norm(x, axis=1)).sum()


def chain(x0, tab, near, beta, v, nsweep, eps, rng, h=1e-3, burn=400):
    x = x0 - (x0 * v).sum() * v
    U = expi(x)
    tr = tr_plaq(U, tab, np.arange(len(tab)))
    ph = phi_of(np.linalg.norm(x, axis=1)).sum()
    nlink = x.shape[0]
    Fs, Ks = [], []
    nacc = 0
    for sw in range(nsweep + burn):
        for _ in range(nlink):
            l = rng.integers(nlink)
            xt = x.copy()
            xt[l] += rng.uniform(-eps, eps, 3)
            t = (xt * v).sum()
            xt -= t * v                                  # exact projection back onto <v,X>=0
            rows = near[l]
            Ut = U.copy()
            Ut[l] = expi(xt[l:l + 1])[0]
            trn = tr_plaq(Ut, tab, rows)
            phn = phi_of(np.linalg.norm(xt, axis=1)).sum()
            dS = beta * np.sum(tr[rows] - trn) + (phn - ph)
            if dS <= 0 or rng.random() < np.exp(-dS):
                x, U, ph = xt, Ut, phn
                tr[rows] = trn
                nacc += 1
        if sw < burn:
            continue
        Vp = Vtot(x + h * v, tab, beta)
        V0 = Vtot(x, tab, beta)
        Vm = Vtot(x - h * v, tab, beta)
        Fs.append((Vp - Vm) / (2 * h))
        Ks.append((Vp - 2 * V0 + Vm) / h ** 2)
    Fs = np.array(Fs); Ks = np.array(Ks)
    return Fs, Ks, nacc / ((nsweep + burn) * nlink)


def jk_lat(K, F, nblk=25):
    n = len(K)
    idx = np.array_split(np.arange(n), nblk)

    def dl(sel):
        return K[sel].mean() - F[sel].var()

    d0 = dl(slice(None))
    ds = np.array([dl(np.concatenate([idx[j] for j in range(nblk) if j != i])) for i in range(nblk)])
    sig = np.sqrt((nblk - 1) / nblk * ((ds - ds.mean()) ** 2).sum())
    return d0, sig


def ac(x, maxlag=200):
    x = x - x.mean()
    v = (x * x).sum()
    return np.array([1.0] + [abs((x[:-k] * x[k:]).sum() / v) for k in range(1, maxlag + 1)])


# =====================================================================
if __name__ == "__main__":
    P("=" * 96)
    P(f"[O11-PHYS] Delta = E[kappa|0] - Var(F|0)   beta = {BETA}   ({time.strftime('%Y-%m-%d %H:%M:%S')})")
    P("=" * 96)

    # ---------------- (1) single plaquette ----------------
    W, V0, Cin, dex = n2_exact_beta(BETA)
    P("\n[1a] single plaquette SU(2): EXACT (deterministic quadrature, 4e5 pts)")
    P(f"   <Phi''> = W0 = {W:.10f}   Var(G)/beta^2 = V0 = {V0:.10f}   E[H] = 0 (exact, centre B)")
    P(f"   Delta = W0 - beta^2 V0 = {dex:+.10f}     (quadrature err < 1e-9)")

    pool_tab = {}
    for N, nw, seeds in ((2, 500000, (101, 102)), (3, 500000, (201, 202))):
        P(f"\n[1b] single plaquette SU({N}): MC pools, {nw} valid samples per group, 2 groups")
        for sd in seeds:
            t = time.time()
            Pp = bs.pool(N, nw, seed=sd, verbose=True)
            r = jk(Pp, BETA)
            P(f"   N={N} seed={sd}: acc={Pp['acc']:.3e}  ESS={r['ess']:.3e}  "
              f"E[H]={r['H']:+.5f}  Var(G)={r['vG']:.6f}  <Phi''>={r['Dd']:.5f}  Var(Phi')={r['vD']:.2e}")
            P(f"            Var(F)={r['varF']:.6f}  E[kappa]={r['kap']:.6f}  "
              f"Delta={r['delta']:+.6f} +- {r['sig']:.6f}  ({r['delta']/max(r['sig'],1e-12):.1f} sigma)"
              f"   [{time.time()-t:.1f}s]")
            pool_tab[(N, sd)] = r
        a, b = pool_tab[(N, seeds[0])], pool_tab[(N, seeds[1])]
        ad = 0.5 * (a['delta'] + b['delta']); sd_ = 0.5 * np.hypot(a['sig'], b['sig'])
        P(f"   >>> N={N}: 2-group mean Delta = {ad:+.6f} +- {sd_:.6f}"
          f"   (group spread {abs(a['delta']-b['delta']):.6f})  = {ad/max(sd_,1e-12):.1f} sigma")

    # ---------------- (2) L=2 lattice ----------------
    L = 2
    tab, near = lat_tables(L)
    nlink = 4 * L ** 4
    P(f"\n[2] L={L} 4D SU(2): {L**4} sites, {nlink} links, {len(tab)} plaquettes, beta={BETA}")
    rng = np.random.default_rng(7)

    modes = {}
    vA = np.zeros((nlink, 3)); vA[:, 2] = 1.0 / np.sqrt(nlink)     # k=0 uniform, unit norm
    vB = np.zeros((nlink, 3)); vB[0, 2] = 1.0                     # one link, unit norm
    modes['A_k0'] = vA
    modes['B_local'] = vB

    for name, v in modes.items():
        P(f"\n--- mode {name}: |v|=1, <v,X>=0 imposed exactly ---")
        for sd in ((11, 12) if name == 'A_k0' else (21,)):
            x0 = rng.normal(0, 0.5, size=(nlink, 3))
            t = time.time()
            Fs, Ks, acc = chain(x0, tab, near, BETA, v, 1200, 0.35, np.random.default_rng(sd))
            d0, sg = jk_lat(Ks, Fs)
            tau = 1 + 2 * ac(Ks, 100).cumsum()[-1] if False else None
            acfK = ac(Ks, 100); tauK = 1 + 2 * acfK[1:].sum()
            acfF = ac(Fs, 100); tauF = 1 + 2 * acfF[1:].sum()
            neff = len(Ks) / max(tauK, tauF, 1.0)
            P(f"   seed={sd}: acc={acc:.3f}  Nsamp={len(Ks)}  tau_K={tauK:.1f} tau_F={tauF:.1f}  Neff~{neff:.0f}")
            P(f"            E[kappa]={Ks.mean():.4f} (sem {Ks.std()/np.sqrt(neff):.4f})  "
              f"Var(F)={Fs.var():.4f}  E[F]={Fs.mean():+.4f}")
            P(f"            Delta={d0:+.4f} +- {sg:.4f} (block-jk, 25 blocks)  "
              f"= {d0/max(sg,1e-12):.1f} sigma   [{time.time()-t:.1f}s]")
            modes.setdefault(name, []).append((d0, sg, neff, Ks.mean(), Fs.var()))

    P("\n" + "=" * 96)
    P(f"total elapsed {time.time()-T0:.1f}s")
