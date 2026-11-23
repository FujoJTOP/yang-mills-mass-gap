"""
Analyse SU(2) production data: Creutz ratios chi(R,T) with binned-jackknife errors,
static potential V(R), string tension sigma*a^2, and 0++ glueball mass.
"""
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def creutz_from(w):
    """w: dict (R,T)->value. Return dict (R,T)->chi for valid entries."""
    out = {}
    for (R, T) in list(w.keys()):
        a = w.get((R, T)); b = w.get((R-1, T-1)); c = w.get((R, T-1)); d = w.get((R-1, T))
        if None in (a, b, c, d):
            continue
        if min(a, b, c, d) <= 0:
            out[(R, T)] = np.nan
        else:
            out[(R, T)] = -np.log(a * b / (c * d))
    return out


def binned_jackknife(Warr, nbin=20):
    """Warr (nR,nT,ncfg). Returns means dict, chi dict, err dict."""
    nR, nT, n = Warr.shape
    means = {(R+1, T+1): Warr[R, T].mean() for R in range(nR) for T in range(nT)}
    chi0 = creutz_from(means)
    nbin = min(nbin, n)
    edges = np.linspace(0, n, nbin+1).astype(int)
    chis = []
    for b in range(nbin):
        idx = np.r_[0:edges[b], edges[b+1]:n]
        if len(idx) < 2:
            continue
        m = {(R+1, T+1): Warr[R, T][idx].mean() for R in range(nR) for T in range(nT)}
        chis.append(creutz_from(m))
    err = {}
    for key in chi0:
        vals = np.array([c.get(key, np.nan) for c in chis])
        vals = vals[np.isfinite(vals)]
        if len(vals) < 2:
            err[key] = np.nan
        else:
            mm = vals.mean()
            err[key] = np.sqrt((len(vals)-1)/len(vals) * np.sum((vals-mm)**2))
    return means, chi0, err


def static_potential(Warr):
    """V(R) from -ln[W(R,T+1)/W(R,T)] at each T; returns dict R->(T,V)."""
    nR, nT, n = Warr.shape
    Vt = {}
    for R in range(nR):
        for T in range(nT-1):
            a = Warr[R, T].mean(); b = Warr[R, T+1].mean()
            if a > 0 and b > 0:
                Vt[(R+1, T+1)] = -np.log(b/a)
    return Vt


def glueball(O):
    """O shape (ncfg, L). Return C(t) and jackknife err (t=0..L-1)."""
    ncfg, L = O.shape
    Obar = O.mean()

    def corr(Oo):
        C = np.zeros(L)
        ob = Oo.mean()
        for t in range(L):
            s = np.mean([np.mean(Oo[:, t0]*Oo[:, (t0+t) % L]) for t0 in range(L)])
            C[t] = s - ob*ob
        return C

    C0 = corr(O)
    nbin = min(20, ncfg)
    edges = np.linspace(0, ncfg, nbin+1).astype(int)
    Cs = []
    for b in range(nbin):
        idx = np.r_[0:edges[b], edges[b+1]:ncfg]
        if len(idx) < 2:
            continue
        Cs.append(corr(O[idx]))
    Cs = np.array(Cs)
    err = np.sqrt((len(Cs)-1)/len(Cs) * np.sum((Cs-Cs.mean(0))**2, axis=0))
    return C0, err


def cosh_fit_global(C, Lt, tmax=None):
    """Fit C(t)=A[exp(-m t)+exp(-m(Lt-t))]+b, t=0..tmax. Returns m."""
    if tmax is None:
        tmax = Lt//2
    ts = np.arange(0, tmax+1)
    y = C[:tmax+1]
    mask = np.isfinite(y)
    ts, y = ts[mask], y[mask]
    if len(ts) < 3:
        return np.nan, np.nan
    best = (np.inf, np.nan, None)
    for m in np.linspace(0.05, 3.0, 590):
        f = np.exp(-m*ts) + np.exp(-m*(Lt-ts))
        A = np.vstack([f, np.ones_like(f)]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = np.sum((A@coef - y)**2)
        if r < best[0]:
            best = (r, m, coef)
    return best[1], np.sqrt(best[0]/max(len(ts)-2, 1))


def local_mass(C, t):
    if t-1 < 0 or t+1 >= len(C):
        return np.nan
    a = C[t-1]+C[t+1]
    b = 2*C[t]
    if b <= 0 or a <= 0:
        return np.nan
    return np.arccosh(a/b)


def analyse_file(path, nbin=20):
    d = np.load(path)
    beta = float(d['beta']); L = int(d['L']); ncfg = int(d['ncfg'])
    Warr = d['W']
    means, chi, err = binned_jackknife(Warr, nbin=nbin)
    smears = [int(x) for x in d['smears']]
    return dict(beta=beta, L=L, ncfg=ncfg, Warr=Warr, means=means, chi=chi,
                err=err, smears=smears, d=d)


def fmt_chi_table(res, Rmax_show=None):
    chi, err = res['chi'], res['err']
    Rs = sorted(set(R for (R, T) in chi))
    Ts = sorted(set(T for (R, T) in chi))
    lines = []
    for R in Rs:
        cells = []
        for T in Ts:
            if (R, T) not in chi:
                continue
            v = chi[(R, T)]; e = err[(R, T)]
            if np.isfinite(v) and np.isfinite(e):
                cells.append(f"{v:+.4f}±{e:.4f}")
            else:
                cells.append("  --  ")
        lines.append(f"R={R}: " + " | ".join(cells))
    return "\n".join(lines)


def sigma_from_creutz(res, rel_err_max=0.4):
    """Plateau value of sigma*a^2 from diagonal chi(R,R) and best off-diagonal.
    Returns dict with chosen (R,T), value, err, and plateau list."""
    chi, err = res['chi'], res['err']
    cand = []
    for (R, T), v in chi.items():
        if R < 2 or T < 2:
            continue
        e = err[(R, T)]
        if np.isfinite(v) and np.isfinite(e) and v > 0 and e < rel_err_max * v:
            cand.append((R, T, v, e))
    if not cand:
        return None
    # prefer largest min(R,T), then smallest error
    cand.sort(key=lambda c: (-min(c[0], c[1]), c[3]))
    best = cand[0]
    return dict(R=best[0], T=best[1], sigma=best[2], err=best[3],
                plate=[c for c in cand if min(c[0], c[1]) >= min(best[0], best[1]) - 1])


def correlated_area_fit(Warr, nboot=200, seed=0):
    """Bootstrap: fit ln W = c - sigma*R*T - p*(R+T) over R,T>=2, W>0.
    Returns sigma mean, std."""
    nR, nT, n = Warr.shape
    rng = np.random.default_rng(seed)
    sig = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        Wb = Warr[:, :, idx].mean(axis=-1)
        xs, ar, pe, y = [], [], [], []
        for R in range(1, nR+1):
            for T in range(1, nT+1):
                v = Wb[R-1, T-1]
                if v > 0 and R >= 2 and T >= 2:
                    ar.append(R*T); pe.append(R+T); y.append(np.log(v))
        if len(y) < 4:
            continue
        A = np.vstack([np.ones(len(y)), np.array(pe), np.array(ar)]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        sig.append(-coef[2])
    sig = np.array(sig)
    return sig.mean(), sig.std()


if __name__ == '__main__':
    import glob
    files = sorted(glob.glob(os.path.join(HERE, 'prod_L8_b*.npz')))
    allres = []
    for f in files:
        r = analyse_file(f)
        allres.append(r)
        print("="*90)
        print(f"beta={r['beta']}  L={r['L']}  ncfg={r['ncfg']}")
        print(fmt_chi_table(r))
        sc = sigma_from_creutz(r)
        if sc:
            print(f"  => plateau Creutz: chi({sc['R']},{sc['T']}) = "
                  f"{sc['sigma']:.4f} ± {sc['err']:.4f}  (sigma*a^2)")
            pl = [(c[0], c[1], round(c[2], 4), round(c[3], 4)) for c in sc['plate']]
            print(f"     plateau set (R,T,chi,err): {pl}")
        sf, se = correlated_area_fit(r['Warr'])
        print(f"  => bootstrap area-law fit (R,T>=2): sigma*a^2 = {sf:.4f} ± {se:.4f}")
        # potential differences at largest usable T
        Vt = static_potential(r['Warr'])
        print("V(R) = -ln[W(R,T+1)/W(R,T)]  at T=2,3,4,5:")
        for R in range(1, r['Warr'].shape[0]+1):
            vs = [Vt.get((R, T), np.nan) for T in (2, 3, 4, 5)]
            print(f"  R={R}: " + " ".join(f"{v:.3f}" if np.isfinite(v) else " --  " for v in vs))
        # glueball per smearing
        print("glueball 0++ (spatial plaquette op):")
        for n in r['smears']:
            O = r['d'][f'O_{n}']
            C, Ce = glueball(O)
            Ct = C/C[0] if C[0] != 0 else C
            m4, _ = cosh_fit_global(C, r['L'], tmax=4)
            ms = [local_mass(C, t) for t in (1, 2, 3)]
            print(f"  smear={n:2d}  C/C0(t=0..4)=" +
                  " ".join(f"{Ct[t]:+.3f}" for t in range(5)) +
                  f"  m*a(cosh,t<=4)={m4:.3f}  meff(1,2,3)=" +
                  " ".join(f"{x:.2f}" if np.isfinite(x) else "--" for x in ms))
    json.dump([{k: v for k, v in r.items() if k not in ('Warr', 'd', 'means', 'chi', 'err',
               'smears')} for r in allres], open(os.path.join(HERE, 'sigma_summary.json'), 'w'))
