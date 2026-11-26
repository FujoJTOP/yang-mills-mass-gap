"""GH3 analysis -- reads gh3_<tag>.npz (no MC). Tables, GEVP, p-ratio, cost model.

Run: py -3.12 gh3_analyze.py L6
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
OBSN = ['W11_p0', 'W11_px', 'W11_pxxx', 'W12_p0', 'W22_p0']
LNAMES = ['raw', 'a10n1', 'a20n1', 'a30n1', 'a50n1', 'a10n3']


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def corr_pair(Si, Sj, Lt):
    a = Si - Si.mean(axis=0)
    b = Sj - Sj.mean(axis=0)
    C = np.zeros(Lt)
    for t in range(Lt):
        C[t] = np.mean(np.sum(a * np.roll(b, -t, axis=1), axis=1)) / Lt
    return C


def jack_stats(S, Lt, nb=20):
    """returns C(t) mean and JK error on C(t) for a single slice-sum array."""
    N = S.shape[0]
    C = corr_pair(S, S, Lt)
    idx = np.array_split(np.arange(N), nb)
    vals = np.zeros((nb, Lt))
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        vals[b] = corr_pair(S[keep], S[keep], Lt)
    m = vals.mean(axis=0)
    e = np.sqrt((nb - 1) / nb * np.sum((vals - m) ** 2, axis=0))
    return C, e, vals


def gewp_lam(Ct, t0):
    """C is (Lt, n, n). Return lambda_n(t) sorted desc, plus rank/eig."""
    G = 0.5 * (Ct[t0] + Ct[t0].T)
    w, V = np.linalg.eigh(G)
    keep = w > 1e-2 * w.max()
    W = V[:, keep] / np.sqrt(w[keep])[None, :]
    out = {}
    for t in range(t0 + 1, Ct.shape[0]):
        M = 0.5 * (Ct[t] + Ct[t].T)
        M = W.T @ M @ W
        out[t] = np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]
    return out, int(keep.sum()), w


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'L6'
    d = np.load(os.path.join(HERE, f'gh3_{tag}.npz'))
    S = d['S']                    # (N, nlev, NOBS, Lt)
    Lt = int(d['Lt']); Ls = int(d['Ls']); beta = float(d['beta'])
    N = int(d['n']); gap = int(d['gap']); nch = int(d['chains'])
    nlev = S.shape[1]
    al = d['lev_alpha']; ns = d['lev_n']
    P("=" * 104)
    P(f"GH3 analysis  tag={tag}  beta={beta}  Ls={Ls} Lt={Lt}  N={N} ({nch}x{N//nch}) gap={gap}")
    P("=" * 104)

    # ---------------- per (level, obs) correlation table ----------------
    tab = {}
    for lev in range(nlev):
        P("")
        P(f"### {LNAMES[lev]:>6}  alpha={al[lev]:.2f} n={ns[lev]}")
        P(f"{'obs':>9} {'C(0)':>9} | " + " | ".join(f"Phi{t}" for t in range(1, 6)) +
          " |  " + " ".join(f"SN{t:<4}" for t in range(1, 6)) +
          " | " + " ".join(f"m{t:<2}" for t in range(1, 5)))
        for k in range(5):
            Sm = S[:, lev, k, :]
            C, e, vals = jack_stats(Sm, Lt)
            if C[0] <= 0:
                continue
            ph = [C[t] / C[0] for t in range(1, 6)]
            err = [e[t] / C[0] for t in range(1, 6)]
            sn = [C[t] / e[t] if e[t] > 0 else float('nan') for t in range(1, 6)]
            me = []
            eme = []
            for t in range(1, 5):
                if C[t] > 0 and C[t + 1] > 0:
                    me.append(np.log(C[t] / C[t + 1]))
                    eme.append(np.hypot(e[t] / C[t], e[t + 1] / C[t + 1]))
                else:
                    me.append(float('nan')); eme.append(float('nan'))
            P(f"{OBSN[k]:>9} {C[0]:9.4f} | " + " | ".join(f"{v:6.4f}" for v in ph) +
              " | " + " ".join(f"{v:5.1f}" for v in sn) +
              " | " + " ".join(f"{v:5.2f}" if np.isfinite(v) else "  nan" for v in me))
            P(f"{'  +-':>9} {'':>9} | " + " | ".join(f"{v:6.4f}" for v in err) +
              " | " + " " * 29 + " | " + " ".join(f"{v:5.2f}" if np.isfinite(v) else "  nan" for v in eme))
            tab[(lev, k)] = dict(C0=float(C[0]), phi=ph, err=err, snr=sn, meff=me, meff_err=eme)

    # ---------------- smearing response (is the response sign-flipped vs beta=2.2?) ------
    P("")
    P("### smearing response of the 1x1 operator (Phi1 and C(0) vs alpha)")
    P(f"{'level':>7} {'C(0)':>9} {'Phi1':>8} {'+-':>7} {'C(1)=C0*Phi1':>13} {'C0/C0_raw':>10} {'C1/C1_raw':>10}")
    base = tab[(0, 0)]
    c1r = base['C0'] * base['phi'][0]
    for lev in range(nlev):
        r = tab[(lev, 0)]
        c1 = r['C0'] * r['phi'][0]
        P(f"{LNAMES[lev]:>7} {r['C0']:9.4f} {r['phi'][0]:8.5f} {r['err'][0]:7.5f} {c1:13.4f} "
          f"{r['C0']/base['C0']:10.3f} {c1/c1r:10.3f}")

    # ---------------- p-ratio: is the one-step correlation spatially extended? ----------
    P("")
    P("### momentum ratio  rho_p = Phi1(p) / Phi1(p=0)   [site-local contact -> rho=1 ;")
    P("    propagating state of mass m -> rho ~ e^{-(sqrt(m^2+p^2)-m)}]   p=pi: 0.043 (m<<pi), 0.22 (m=2.6)")
    P(f"{'level':>7} {'rho_px':>8} {'+-':>7} {'rho_pxxx':>9} {'+-':>7}  {'m_eff(1) p=0':>12} {'pi_pred':>8} {'ppp_pred':>9}")
    for lev in range(nlev):
        a0 = tab[(lev, 0)]['phi'][0]; e0 = tab[(lev, 0)]['err'][0]
        a1 = tab[(lev, 1)]['phi'][0]; e1 = tab[(lev, 1)]['err'][0]
        a2 = tab[(lev, 2)]['phi'][0]; e2 = tab[(lev, 2)]['err'][0]
        r1 = a1 / a0; er1 = r1 * np.hypot(e1 / a1, e0 / a0)
        r2 = a2 / a0; er2 = r2 * np.hypot(e2 / a2, e0 / a0)
        m = tab[(lev, 0)]['meff'][0]
        pr1 = np.exp(-(np.sqrt(m * m + np.pi ** 2) - m)) if np.isfinite(m) else float('nan')
        pr2 = np.exp(-(np.sqrt(m * m + 3 * np.pi ** 2) - m)) if np.isfinite(m) else float('nan')
        P(f"{LNAMES[lev]:>7} {r1:8.4f} {er1:7.4f} {r2:9.4f} {er2:7.4f}  {m:12.3f} {pr1:8.4f} {pr2:9.4f}")

    # ---------------- GEVP over the 6 W11_p0 levels, metric t0 = 0,1,2 ------------------
    P("")
    P("### GEVP over 6 W11_p0 levels (raw + 5 smearing levels): lambda_n of C(t) vs C(t0)")
    Ct = np.zeros((Lt, nlev, nlev))
    for i0 in range(nlev):
        for j0 in range(nlev):
            Ct[:, i0, j0] = corr_pair(S[:, i0, 0, :], S[:, j0, 0, :], Lt)
    for t0 in (0, 1, 2):
        lam, rank, w = gewp_lam(Ct, t0)
        P(f"  t0={t0}  rank={rank}  Gram eig=" + " ".join(f"{v:.3g}" for v in np.sort(w)[::-1]))
        for t in sorted(lam):
            l = lam[t]
            if l[0] <= 0:
                P(f"    t={t}  lam_max={l[0]:.3e} (<=0, no plateau)")
                continue
            P(f"    t={t}  lam_max={l[0]:.4e}  m_eff={-np.log(l[0])/(t-t0):7.4f}"
              + (f"  lam2={l[1]:.3e} m2={(-np.log(l[1])/(t-t0) if l[1]>0 else float('nan')):7.3f}"
                 if len(l) > 1 else ""))
    # jackknife error on the t0=1 GEVP m_eff
    nb = 20
    idx = np.array_split(np.arange(N), nb)
    for t0 in (0, 1):
        acc = {t: [] for t in range(t0 + 1, 7)}
        for b in range(nb):
            keep = np.concatenate([idx[j] for j in range(nb) if j != b])
            Cb = np.zeros((Lt, nlev, nlev))
            for i0 in range(nlev):
                for j0 in range(nlev):
                    Cb[:, i0, j0] = corr_pair(S[keep, i0, 0, :], S[keep, j0, 0, :], Lt)
            lb, _, _ = gewp_lam(Cb, t0)
            for t in acc:
                acc[t].append(-np.log(lb[t][0]) / (t - t0) if lb[t][0] > 0 else np.nan)
        P(f"  [JK] t0={t0} GEVP m_eff(t) with 20-block delete-one error:")
        for t in acc:
            v = np.array(acc[t], float)
            m = np.nanmean(v); e = np.sqrt((nb - 1) / nb * np.nansum((v - m) ** 2))
            P(f"       t={t}  m_eff={m:7.4f} +- {e:6.4f}")
    json.dump({f"{k[0]}_{k[1]}": v for k, v in tab.items()},
              open(os.path.join(HERE, f'gh3_analyze_{tag}.json'), 'w'), indent=1)
    open(os.path.join(HERE, f'gh3_analyze_{tag}.log'), 'w').write("\n".join(OUT) + "\n")
    P("")
    P(f"wrote gh3_analyze_{tag}.log")


if __name__ == '__main__':
    main()
