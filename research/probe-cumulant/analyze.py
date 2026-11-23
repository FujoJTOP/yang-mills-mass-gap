"""Analysis of probe.py output.

The CPB-W tilt is   e^{-xi V} mu_G,   mu_G(alpha) ∝ e^{-alpha Q},
    V(alpha)_x = S_x - alpha Q_x,   S_x = beta sum_{p∋x} (1 - Re tr U_p/3),
                                    Q_x = beta sum_{p∋x} tr(A_p^2)/6,
so that at xi=1, for EVERY alpha, e^{-V(alpha)} mu_G(alpha) = e^{-S_W}.
alpha is therefore the free choice of Gaussian reference; alpha=1 is the
"canonical" quadratic part of the Wilson action.

Measured per (L, beta), at b=1 and b=2 block means:
  * k2, k3, k4, R3, R4 for a grid of alpha
  * the zero-separation lower bound on C_cum:  C_cum >= |k3(V_x,V_x,V_x)|
  * three-point function k3(d1,d2), |d|_inf<=1, and the local density
        c3_loc = sum_{d1,d2} |k3(d1,d2)|
    with the MC noise bias of the absolute-value sum estimated and subtracted
  * collinear decay k2(r), k3(0,r,s)
  * worst-configuration / tail analysis
  * control: the blocked plaquette FIELD (expected Gaussian, k3 ~ 0)
Errors: jackknife over 100 config blocks.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'runs')

ALPHAS = (0.0, 0.25, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5)
NBLOCK = 100


def block_mean(A, b):
    """A shape (...,L,L,L,L) -> (...,(L/b)^4) unweighted block means."""
    pre = A.shape[:-4]
    L1 = A.shape[-4]
    n = L1 // b
    R = A.reshape(*pre, n, b, n, b, n, b, n, b)
    nd = len(pre)
    return R.mean(axis=(nd + 1, nd + 3, nd + 5, nd + 7)).reshape(*pre, n ** 4)


def disp_set(L, r):
    out = []
    for a in range(-r, r + 1):
        for b_ in range(-r, r + 1):
            for c in range(-r, r + 1):
                for d in range(-r, r + 1):
                    out.append((a % L, b_ % L, c % L, d % L))
    seen, uniq = set(), []
    for v in out:
        if v not in seen:
            seen.add(v)
            uniq.append(v)
    return uniq


def shifted(A, d):
    for ax, s in zip(range(A.ndim - 4, A.ndim), d):
        A = np.roll(A, -s, axis=ax)
    return A


def raw_moments(x):
    x = x.astype(np.float64)
    return np.array([x.sum(), (x * x).sum(), (x ** 3).sum(), (x ** 4).sum()])


def k_from_raw(s, n):
    m1, m2, m3, m4 = s / n
    k2 = m2 - m1 * m1
    k3 = m3 - 3 * m1 * m2 + 2 * m1 ** 3
    c4 = m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 - 3 * m1 ** 4   # central 4th moment
    k4 = c4 - 3 * k2 * k2                                     # 4th cumulant
    return m1, k2, k3, k4


def jack_err(rep):
    """rep shape (nb, ...): jackknife error of the mean."""
    nb = rep.shape[0]
    return np.sqrt((nb - 1) / nb * np.sum((rep - rep.mean(axis=0)) ** 2, axis=0))


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(RUNS, 'L*_b*.npz')))
    files = [f for f in files if 'test' not in os.path.basename(f)]
    results = []
    for f in files:
        d = np.load(f)
        L = int(d['L']); beta = float(d['beta']); N = d['Sx'].shape[0]
        Sx = d['Sx'].astype(np.float64)
        Qx = d['Qx'].astype(np.float64)
        plaq = d['plaq']
        nb = min(NBLOCK, N)
        blocks = np.array_split(np.arange(N), nb)
        rec = {'file': os.path.basename(f), 'L': L, 'beta': beta, 'N': N,
               'plaq': float(plaq.mean())}
        print(f"\n===== {rec['file']}  L={L} beta={beta}  N={N}  <P>={plaq.mean():.5f}")

        for b in (1, 2):
            if L % b:
                continue
            Sb = block_mean(Sx, b); Qb = block_mean(Qx, b)
            n = L // b
            nv = n ** 4
            ntot = N * nv
            ns = np.array([len(ix) * nv for ix in blocks], float)
            tot = ns.sum()

            # ---------- alpha scan ----------
            sums = {a: np.zeros((nb, 4)) for a in ALPHAS}
            for j, idx in enumerate(blocks):
                for a in ALPHAS:
                    sums[a][j] = raw_moments(Sb[idx] - a * Qb[idx])
            Stot = {a: sums[a].sum(axis=0) for a in ALPHAS}
            print(f"  b={b}:")
            tbl = []
            for a in ALPHAS:
                m1, k2, k3, k4 = k_from_raw(Stot[a], ntot)
                rep = np.array([k_from_raw(Stot[a] - sums[a][j], ntot - ns[j])[2]
                                for j in range(nb)])
                ek3 = float(jack_err(rep))
                R3 = k3 / k2 ** 1.5 if k2 > 0 else np.nan
                R4 = k4 / k2 ** 2 if k2 > 0 else np.nan
                rec[f'b{b}_a{a}_k2'] = float(k2)
                rec[f'b{b}_a{a}_k3'] = float(k3)
                rec[f'b{b}_a{a}_ek3'] = float(ek3)
                rec[f'b{b}_a{a}_R4'] = float(R4)
                tbl.append((a, k2, k3, ek3, R3))
            for a, k2, k3, ek3, R3 in tbl:
                print(f"    alpha={a:4.2f}: var={k2:9.4f}  k3={k3:10.4f}+-{ek3:7.4f}"
                      f"  R3={R3:7.3f}")
            amin = min(tbl, key=lambda t: abs(t[2]))
            aa = np.array([t[0] for t in tbl]); kk = np.array([t[2] for t in tbl])
            ee = np.array([max(t[3], 1e-12) for t in tbl]); vv = np.array([t[1] for t in tbl])
            try:
                c3 = np.polyfit(aa, kk, 3, w=1.0 / ee)
                c2 = np.polyfit(aa, vv, 2)
                ag = np.linspace(0.0, 1.5, 1501)
                i = int(np.argmin(np.abs(np.polyval(c3, ag))))
                astar, kstar = float(ag[i]), float(np.polyval(c3, ag[i]))
                vstar = float(np.polyval(c2, astar))
            except Exception:
                astar, kstar, vstar = amin[0], amin[2], amin[1]
            print(f"    -> min |k3| at alpha={amin[0]:.2f}: k3={amin[2]:.4f}"
                  f"+-{amin[3]:.4f}; cubic fit alpha*={astar:.3f}, k3={kstar:+.4f}, "
                  f"var={vstar:.4f}")
            rec[f'b{b}_alpha_min'] = float(amin[0])
            rec[f'b{b}_k3_min'] = float(amin[2])
            rec[f'b{b}_alpha_star'] = astar
            rec[f'b{b}_k3_star'] = kstar
            rec[f'b{b}_var_star'] = vstar

            # ---------- field control P = 1 - S/(6 beta) ----------
            Pb = 1.0 - Sb / (6.0 * beta)
            fs = np.zeros((nb, 4))
            for j, idx in enumerate(blocks):
                fs[j] = raw_moments(Pb[idx])
            Ftot = fs.sum(axis=0)
            m1P, k2P, k3P, k4P = k_from_raw(Ftot, ntot)
            repP = np.array([k_from_raw(Ftot - fs[j], ntot - ns[j])[2] for j in range(nb)])
            ekP = float(jack_err(repP))
            print(f"    FIELD P: var={k2P:.6f}  k3={k3P:+.2e}+-{ekP:.2e}"
                  f"  R3={k3P/k2P**1.5:+.4f}   [control: ~0]")
            rec[f'b{b}_field_k3'] = float(k3P)
            rec[f'b{b}_field_ek3'] = ekP
            rec[f'b{b}_field_var'] = float(k2P)

            # ---------- three-point function at alpha=1 ----------
            V = Sb - Qb
            mu = V.mean()
            X = (V - mu).reshape(N, n, n, n, n)
            D = disp_set(n, 1)
            nD = len(D)
            rr = min(3, n // 2)
            Dcol = []
            for ax in range(4):
                for r in range(rr + 1):
                    v = [0, 0, 0, 0]; v[ax] = r
                    Dcol.append(tuple(v))
            B3 = np.zeros((nb, nD, nD)); B3c = np.zeros((nb, len(Dcol), len(Dcol)))
            B2c = np.zeros((nb, len(Dcol)))
            k3c = np.zeros(N); maxabs = np.zeros(N)
            for j, idx in enumerate(blocks):
                for c in idx:
                    x = X[c].ravel()
                    k3c[c] = np.mean(x ** 3)
                    maxabs[c] = np.abs(V[c]).max()
                    Sh = np.stack([shifted(X[c], t).ravel() for t in D])
                    B3[j] += Sh @ (Sh * x[None, :]).T
                    Shc = np.stack([shifted(X[c], t).ravel() for t in Dcol])
                    B2c[j] += Shc @ x
                    B3c[j] += Shc @ (Shc * x[None, :]).T
            tot3 = B3.sum(axis=0)
            rep = (tot3[None] - B3) / (tot - ns)[:, None, None]
            sig = jack_err(rep)                       # per-pair jackknife sigma
            B3t = tot3 / tot
            full = np.abs(B3t).sum()
            bias = np.sqrt(2.0 / np.pi) * sig.sum()   # E|noise| in the abs sum
            ec3 = float(jack_err(np.array([np.abs(r).sum() for r in rep])))
            strong = np.abs(B3t) > 2.0 * sig
            c3sig = float(np.abs(B3t[strong]).sum())
            i0 = D.index((0, 0, 0, 0))
            print(f"    THREE-POINT (|d|<=1, {nD}^2 pairs): k3(0,0)={B3t[i0,i0]:9.4f}")
            print(f"      c3_abs={full:9.3f}+-{ec3:.3f}   noise bias={bias:8.3f}"
                  f"   c3_corr={full-bias:9.3f}   c3(>2sig)={c3sig:9.3f}"
                  f"   [{int(strong.sum())}/{nD*nD} pairs]")
            rec[f'b{b}_k3_000'] = float(B3t[i0, i0])
            rec[f'b{b}_c3_abs'] = float(full)
            rec[f'b{b}_c3_corr'] = float(full - bias)
            rec[f'b{b}_c3_2sig'] = c3sig
            rec[f'b{b}_ec3'] = ec3

            # ---------- decay ----------
            B2t = B2c.sum(axis=0) / tot
            B3ct = B3c.sum(axis=0) / tot
            nax = 4
            C2 = [float(np.mean([B2t[ax * (rr + 1) + r] for ax in range(nax)]))
                  for r in range(rr + 1)]
            C3c = {}
            for r in range(rr + 1):
                for s in range(rr + 1):
                    C3c[(r, s)] = float(np.mean(
                        [B3ct[ax * (rr + 1) + r, ax * (rr + 1) + s] for ax in range(nax)]))
            rec[f'b{b}_C2'] = C2
            rec[f'b{b}_C3col'] = {f'{k[0]},{k[1]}': v for k, v in C3c.items()}
            print(f"    C2(r=0..{rr}) = {np.round(C2,5)}")
            print(f"    C3 collinear: " + "  ".join(
                f"({r},{s})={C3c[(r,s)]:+.3e}" for r in range(rr + 1)
                for s in range(r, rr + 1)))
            if C2[0] > 0 and len(C2) > 1 and abs(C2[1]) > 0:
                print(f"    C2(1)/C2(0)={C2[1]/C2[0]:+.5f} -> xi_V <= "
                      f"{1.0/max(1e-12,-np.log(abs(C2[1]/C2[0]))):.2f} sites")

            # ---------- worst configuration / tail ----------
            order = np.argsort(np.abs(k3c))
            n1 = max(1, N // 100)
            frac_top1 = k3c[order[-n1:]].sum() / k3c.sum()
            frac_top5 = k3c[order[-max(1, N // 20):]].sum() / k3c.sum()
            wc = int(np.argmax(maxabs))
            k3site = np.mean(X.reshape(N, -1) ** 3, axis=0)
            print(f"    WORST: max|V|={maxabs.max():.2f} (config {wc}); k3_c mean="
                  f"{k3c.mean():+.4f} top1%={frac_top1:.3f} top5%={frac_top5:.3f}")
            print(f"    per-site k3: sup={k3site.max():+.3f} std={k3site.std():.3f}"
                  f" (site-noise ~{np.std(X**3)/np.sqrt(N):.3f})")
            rec[f'b{b}_maxabsV'] = float(maxabs.max())
            rec[f'b{b}_k3_top1_share'] = float(frac_top1)
            rec[f'b{b}_k3_top5_share'] = float(frac_top5)
            rec[f'b{b}_k3site_sup'] = float(k3site.max())
            rec[f'b{b}_k3site_std'] = float(k3site.std())
            rec[f'b{b}_Vmin'] = float(V.min()); rec[f'b{b}_Vmax'] = float(V.max())

        # ---------- per-plaquette (activity-normalisation) ----------
        if 'pmom' in d.files:
            AL2 = list(d['alphas']); pm = d['pmom']; fmo = d['fmom']
            npl = 6 * L ** 4
            ns_p = np.array([len(ix) * npl for ix in blocks], float)
            tot_p = ns_p.sum()
            print("  PER-PLAQUETTE (6L^4 values/config):")
            for k, a in enumerate(AL2):
                sums = np.array([pm[ix, k, :].sum(axis=0) for ix in blocks])
                Sp = sums.sum(axis=0)
                m1, k2, k3, k4 = k_from_raw(Sp, tot_p)
                rep = np.array([k_from_raw(Sp - sums[j], tot_p - ns_p[j])[2]
                                for j in range(nb)])
                ek3 = float(jack_err(rep))
                if a in (0.0, 0.5, 0.6, 1.0, 1.5):
                    print(f"    alpha={a:4.2f}: var={k2:8.4f} k3={k3:9.4f}+-{ek3:6.4f}"
                          f" R3={k3/k2**1.5:7.3f}")
                rec[f'p_a{a}_k2'] = float(k2); rec[f'p_a{a}_k3'] = float(k3)
                rec[f'p_a{a}_ek3'] = ek3
            # field control on the plaquette field p_p = 1-Re tr/3
            fsum = np.array([fmo[ix].sum(axis=0) for ix in blocks])
            Ft = fsum.sum(axis=0)
            m1, k2, k3, k4 = k_from_raw(Ft, tot_p)
            rep = np.array([k_from_raw(Ft - fsum[j], tot_p - ns_p[j])[2]
                            for j in range(nb)])
            print(f"    FIELD p_p: var={k2:.6f} k3={k3:+.2e}+-{jack_err(rep):.2e}")
            rec['p_field_k3'] = float(k3)
            rec['p_maxabsV'] = float(d['pmax'].max())

        # white-noise check: k3(block mean) vs k3(site)/b^8
        if 'b1_a1.0_k3' in rec and 'b2_a1.0_k3' in rec:
            pred = rec['b1_a1.0_k3'] / 2.0 ** 8
            print(f"  white-noise check: k3(b=2)={rec['b2_a1.0_k3']:+.5f} vs "
                  f"k3(b=1)/256={pred:+.5f}  ratio={rec['b2_a1.0_k3']/pred:.3f}")
            rec['white_ratio'] = float(rec['b2_a1.0_k3'] / pred)
        results.append(rec)

    with open(os.path.join(HERE, 'results.json'), 'w') as fh:
        json.dump(results, fh, indent=1)
    print('\nwrote results.json')


if __name__ == '__main__':
    main()
