"""Analysis: autocorrelation-aware error bars (binning + jackknife), all
observables, plots, results.csv."""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'runs')
PLOTS = os.path.join(HERE, 'plots')


# ----------------------------------------------------------------------
# statistics helpers
# ----------------------------------------------------------------------
def tau_int(x, c=6.0):
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return 0.5
    s = 0.0
    tau = 0.5
    W = max(1, min(n // 10, 2000))
    for t in range(1, W):
        rho = np.dot(x[:-t], x[t:]) / (n - t) / v
        s += rho
        tau = 0.5 + s
        if t >= c * tau:
            break
    return max(tau, 0.5)


def bins_of(n, b):
    b = max(b, int(np.ceil(n / 400.0)))   # cap the number of jackknife bins
    nb = n // b
    return [np.arange(i * b, (i + 1) * b) for i in range(nb)]


def jackknife_est(arrays, stat, b):
    """arrays: dict name -> (n,) or (n,...) per-config data. stat(dict)->float."""
    n = len(next(iter(arrays.values())))
    bs = bins_of(n, b)
    full = stat(arrays)
    vals = []
    for k in range(len(bs)):
        idx = np.concatenate([bs[j] for j in range(len(bs)) if j != k])
        vals.append(stat({key: a[idx] for key, a in arrays.items()}))
    vals = np.array(vals)
    nb = len(bs)
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


# ----------------------------------------------------------------------
# load
# ----------------------------------------------------------------------
def load_runs():
    out = []
    for f in sorted(glob.glob(os.path.join(RUNS, '*.npz'))):
        if os.path.basename(f).startswith('block_'):
            continue
        d = dict(np.load(f, allow_pickle=True))
        d['file'] = os.path.basename(f)
        out.append(d)
    return out


def group_runs(runs):
    g = {}
    for d in runs:
        key = (int(d['L']), round(float(d['beta']), 4))
        g.setdefault(key, []).append(d)
    return g


def merge(ds):
    """Concatenate per-config measurements of independent chains."""
    keys = ['plaq', 'mag_slice', 'poly_corr', 'poly_mean', 'defect',
            'plaq_min', 'plaq_max']
    out = {k: np.concatenate([d[k] for d in ds]) for k in keys}
    blocks = set()
    for d in ds:
        for k in d.files if hasattr(d, 'files') else d.keys():
            if k.startswith('blk') and k[3:].isdigit():
                blocks.add(int(k[3:]))
    for b in sorted(blocks):
        out[f'blk{b}'] = np.concatenate([d[f'blk{b}'] for d in ds])
    out['loops'] = np.concatenate([d['loops'] for d in ds])
    out['L'] = int(ds[0]['L'])
    out['beta'] = float(ds[0]['beta'])
    out['rmax'] = int(ds[0]['rmax'])
    out['nchains'] = len(ds)
    return out


# ----------------------------------------------------------------------
# observables
# ----------------------------------------------------------------------
def plaquette_obs(d):
    return float(np.mean(d['plaq']))


def loop_obs(arrays):
    return np.mean(arrays['loops'], axis=0)


def corr_connected(d):
    """C(t) = <S(t0)S(t0+t)> - <S>^2 for the magnetic density slice sums."""
    S = d['mag_slice']
    n, L = S.shape
    m = S.mean()
    C = np.zeros(L)
    for t in range(L):
        C[t] = np.mean(S * np.roll(S, -t, axis=1))
    C -= m * m
    return C


def fit_xi(C, L, tmin=1, tmax=None):
    """Fit C(t) = A (e^{-m t} + e^{-m (L-t)}) on t in [tmin,tmax]."""
    if tmax is None:
        tmax = L // 2
    ts = np.arange(tmin, tmax + 1)
    y = C[ts]
    ms = np.linspace(0.05, 3.0, 3000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    return float(ms[k]), float(A[k]), ts, y


def m_eff(C, L):
    out = []
    for t in range(1, L // 2 + 1):
        cm = C[(t - 1) % L]
        c0 = C[t % L]
        cp = C[(t + 1) % L]
        if cm > 1e-16 and c0 > 1e-16 and cp > 1e-16:
            out.append((t, np.arccosh((cm + cp) / (2 * c0))))
        else:
            out.append((t, np.nan))
    return out


def r4_from_sums(sums):
    """sums: (..., 7) = [S1,S2,S3,S4,count,min,max]; returns (mean,kappa2,kappa4,R4)."""
    s1 = np.sum(sums[..., 0])
    s2 = np.sum(sums[..., 1])
    s3 = np.sum(sums[..., 2])
    s4 = np.sum(sums[..., 3])
    n = np.sum(sums[..., 4])
    m1 = s1 / n
    m2 = s2 / n
    m3 = s3 / n
    m4 = s4 / n
    c2 = m2 - m1 * m1
    mu4 = m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 - 3 * m1 ** 4
    c4 = mu4 - 3 * c2 * c2
    r4 = c4 / c2 ** 2 if c2 > 0 else np.nan
    return m1, c2, c4, r4


def kappa_results():
    f = os.path.join(HERE, 'kappa_results.npz')
    if not os.path.exists(f):
        return {}
    d = dict(np.load(f, allow_pickle=True))
    out = {}
    for k in d:
        if k in ('keys',):
            continue
        out.setdefault(tuple(d[k][:3].astype(int)) if d[k].ndim else k, None)
    return d


# ----------------------------------------------------------------------
def main():
    os.makedirs(PLOTS, exist_ok=True)
    runs = load_runs()
    g = group_runs(runs)
    rows = []

    summary = {}
    for (L, beta), ds in sorted(g.items()):
        d = merge(ds)
        n = len(d['plaq'])
        tauP = tau_int(d['plaq'])
        b = max(1, int(np.ceil(2 * tauP)))
        P, eP = jackknife_est({'p': d['plaq']}, lambda a: float(np.mean(a['p'])), b)
        Lw, eL = jackknife_est({'loops': d['loops']}, loop_obs, b)
        C = corr_connected(d)

        def xi_stat(a, L=L):
            Cc = corr_connected({'mag_slice': a['s']})
            mm = fit_xi(Cc, L, 1, L // 2)[0]
            return 1.0 / mm

        xi_val, xi_err = jackknife_est({'s': d['mag_slice']}, xi_stat, b)

        def m_stat(a, L=L):
            Cc = corr_connected({'mag_slice': a['s']})
            return fit_xi(Cc, L, 1, L // 2)[0]

        m, m_err = jackknife_est({'s': d['mag_slice']}, m_stat, b)
        A = fit_xi(C, L)[1]
        rows += [
            dict(obs='plaquette', L=L, beta=beta, value=P, err=eP, n=n, tau=tauP, binsize=b),
            dict(obs='xi_over_a', L=L, beta=beta, value=xi_val, err=xi_err, n=n, tau=tauP, binsize=b),
            dict(obs='mass_a', L=L, beta=beta, value=m, err=m_err, n=n, tau=tauP, binsize=b),
        ]
        for k, thr in enumerate([0.0, 0.1, 0.2, 1.0 / 3.0]):
            rho, erho = jackknife_est({'x': d['defect'][:, k]},
                                      lambda a: float(np.mean(a['x'])), b)
            rows.append(dict(obs=f'defect_frac_below_{thr:.3f}', L=L, beta=beta,
                             value=rho, err=erho, n=n, tau=tauP, binsize=b))
        for R in range(1, Lw.shape[0] + 1):
            for T in range(1, Lw.shape[1] + 1):
                rows.append(dict(obs=f'W_{R}_{T}', L=L, beta=beta,
                                 value=Lw[R - 1, T - 1], err=eL[R - 1, T - 1],
                                 n=n, tau=tauP, binsize=b))
        r4 = {}
        for k in d:
            if k.startswith('blk') and k[3:].isdigit():
                bb = int(k[3:])
                def st(a, bb=bb):
                    return r4_from_sums(a['x'])[3]
                val, ev = jackknife_est({'x': d[k]}, st, b)
                m1, c2, c4, _ = r4_from_sums(d[k])
                r4[bb] = dict(mean=m1, kappa2=c2, kappa4=c4, R4=val, err=ev)
                rows.append(dict(obs=f'R4_block{bb}', L=L, beta=beta, value=val,
                                 err=ev, n=n, tau=tauP, binsize=b))
        # Polyakov cross-check
        def poly_mass(a):
            Pc = np.mean(a['pc'], axis=0)
            Pm = a['pm']
            C2 = Pc.copy()
            C2[1:] -= np.abs(np.mean(Pm[:, 0] + 1j * Pm[:, 1])) ** 2
            rr = np.arange(len(C2))
            ok = (rr >= 1) & (rr <= len(C2) // 2) & (C2 > 0)
            if ok.sum() < 2:
                return np.nan
            return -np.polyfit(rr[ok], np.log(C2[ok]), 1)[0]
        try:
            pm, epm = jackknife_est({'pc': d['poly_corr'], 'pm': d['poly_mean']},
                                    poly_mass, b)
        except Exception:
            pm, epm = np.nan, np.nan
        rows.append(dict(obs='polyakov_mass_a', L=L, beta=beta, value=pm, err=epm,
                         n=n, tau=tauP, binsize=b))
        summary[(L, beta)] = dict(d=d, P=P, eP=eP, tauP=tauP, b=b, C=C, m=m, em=m_err,
                                  xi=xi_val, exi=xi_err,
                                  loops=Lw, eL=eL, r4=r4, n=n, poly_mass=pm)

    # ---------------- results.csv ----------------
    with open(os.path.join(HERE, 'results.csv'), 'w') as f:
        f.write('observable,L,beta,value,error,n_configs,tau_int_plaq,bin_size\n')
        for r in rows:
            f.write(f"{r['obs']},{r['L']},{r['beta']:.4f},{r['value']:.8g},"
                    f"{r['err']:.4g},{r['n']},{r['tau']:.3f},{r['binsize']}\n")
    print("wrote results.csv")

    # ---------------- plots ----------------
    # 1 plaquette vs beta
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for L in sorted({k[0] for k in summary}):
        xs, ys, es = [], [], []
        for (LL, beta), s in sorted(summary.items()):
            if LL == L:
                xs.append(beta); ys.append(s['P']); es.append(s['eP'])
        ax.errorbar(xs, ys, yerr=es, marker='o', ls='-', label=f'L={L}', capsize=3)
    bb = np.linspace(0.5, 3, 100)
    ax.plot(bb, bb / 18.0, 'k--', lw=1, label=r'strong coupling $\beta/18$')
    ax.set_xlabel(r'$\beta$'); ax.set_ylabel(r'$\langle \frac{1}{3} \mathrm{Re\,tr}\,U_p\rangle$')
    ax.set_title('Average plaquette'); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, 'plaquette_vs_beta.png'), dpi=140)
    plt.close(fig)

    # 2 Wilson loops and static potential
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for (L, beta), s in sorted(summary.items()):
        if L != 6:
            continue
        lw = np.maximum(s['loops'], 1e-12)
        for R in range(lw.shape[0]):
            axs[0].plot(np.arange(1, lw.shape[1] + 1), -np.log(lw[R]), 'o-',
                        label=f'R={R+1}, b={beta}')
    axs[0].set_xlabel('T'); axs[0].set_ylabel(r'$-\log W(R,T)$')
    axs[0].set_title('Wilson loops (L=6)'); axs[0].legend(fontsize=7); axs[0].grid(alpha=.3)
    for (L, beta), s in sorted(summary.items()):
        if L != 6:
            continue
        lw = np.maximum(s['loops'], 1e-14)
        V = []
        for R in range(lw.shape[0]):
            Ts = np.arange(1, lw.shape[1] + 1)
            # V(R) from ratio W(R,T+1)/W(R,T) = exp(-V) (large T)
            if len(Ts) >= 2:
                V.append(-np.log(lw[R, -1] / lw[R, -2]))
            else:
                V.append(np.nan)
        axs[1].plot(np.arange(1, len(V) + 1), V, 'o-', label=f'b={beta}')
    axs[1].set_xlabel('R'); axs[1].set_ylabel(r'$aV(R)$')
    axs[1].set_title('Static potential (crude, L=6)'); axs[1].legend(fontsize=7)
    axs[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, 'wilson_loops.png'), dpi=140)
    plt.close(fig)

    # 3 correlation functions
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for (L, beta), s in sorted(summary.items()):
        C = s['C']; ts = np.arange(1, L // 2 + 1)
        axs[0].plot(ts, C[ts], 'o-', label=f'L={L}, b={beta}')
    axs[0].set_yscale('log'); axs[0].set_xlabel('t'); axs[0].set_ylabel('C(t)')
    axs[0].set_title('connected plaquette correlator'); axs[0].legend(fontsize=7)
    axs[0].grid(alpha=.3)
    for (L, beta), s in sorted(summary.items()):
        me = m_eff(s['C'], L)
        if len(me):
            axs[1].plot([a for a, b in me], [b for a, b in me], 'o-',
                        label=f'L={L}, b={beta}')
    axs[1].set_xlabel('t'); axs[1].set_ylabel(r'$a\,m_{\rm eff}$')
    axs[1].set_title('effective mass'); axs[1].legend(fontsize=7); axs[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, 'correlators.png'), dpi=140)
    plt.close(fig)

    # 4 R4 vs blocking
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for (L, beta), s in sorted(summary.items()):
        bs = sorted(s['r4'])
        ax.errorbar(bs, [s['r4'][b]['R4'] for b in bs],
                    yerr=[s['r4'][b]['err'] for b in bs],
                    marker='o', ls='-', capsize=3, label=f'L={L}, b={beta}')
    ax.axhline(0, color='k', lw=.8)
    ax.set_xlabel('block size b (in lattice units)')
    ax.set_ylabel(r'$R_4=\kappa_4/\kappa_2^2$')
    ax.set_title('connected 4-point ratio of blocked plaquette')
    ax.legend(fontsize=7); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, 'R4_blocking.png'), dpi=140)
    plt.close(fig)

    # 5 defect density
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for (L, beta), s in sorted(summary.items()):
        d = s['d']
        ax.plot([0, .1, .2, 1 / 3],
                [np.mean(d['defect'][:, k]) for k in range(4)], 'o-',
                label=f'L={L}, b={beta}')
    ax.set_yscale('log'); ax.set_xlabel('threshold on $\\frac{1}{3}$Re tr $U_p$')
    ax.set_ylabel(r'defect fraction'); ax.set_title('defect density')
    ax.legend(fontsize=7); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, 'defects.png'), dpi=140)
    plt.close(fig)

    # ---------------- console summary ----------------
    with open(os.path.join(HERE, 'summary.txt'), 'w') as f:
        for (L, beta), s in sorted(summary.items()):
            if L == 6 or True:
                f.write(f"L={L} beta={beta:.3f}: <P>={s['P']:.5f}({s['eP']:.5f}) "
                        f"tau_P={s['tauP']:.2f} bin={s['b']} N={s['n']} "
                        f"a*m={s['m']:.4f} xi/a={1/s['m']:.3f}\n")
                for bb in sorted(s['r4']):
                    r = s['r4'][bb]
                    f.write(f"    block {bb}: R4={r['R4']:+.4f}({r['err']:.4f}) "
                            f"kappa2={r['kappa2']:.3e} kappa4={r['kappa4']:.3e}\n")
                f.write(f"    defects(0,.1,.2,1/3)="
                        f"{np.round(np.mean(s['d']['defect'],axis=0),5)}\n")
                f.write(f"    polyakov mass (a): {s['poly_mass']:.4f}\n")
    print(open(os.path.join(HERE, 'summary.txt')).read())


if __name__ == '__main__':
    main()
