"""GH3 -- GEVP subset-robustness scan, off the npz (no MC).

Question: the paper's basis was {raw, APE alpha=0.5 n=1} (2 operators, window t=1,2).
Does a plateau need the extra levels, and is it stable under level subsets?

Run: py -3.12 gh3_gevp_scan.py L6 L8
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
LNAMES = ['raw', 'a10n1', 'a20n1', 'a30n1', 'a50n1', 'a10n3']
SUBSETS = [[0, 4], [0, 4, 5], [0, 1, 2, 3, 4, 5], [0, 1, 3, 4], [0, 3], [1, 3, 5], [0, 3, 5]]


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def corr_pair(Si, Sj, Lt):
    a = Si - Si.mean(axis=0)
    b = Sj - Sj.mean(axis=0)
    return np.array([np.mean(np.sum(a * np.roll(b, -t, axis=1), axis=1)) / Lt for t in range(Lt)])


def lam_max(Ct, t0):
    G = 0.5 * (Ct[t0] + Ct[t0].T)
    w, V = np.linalg.eigh(G)
    keep = w > 1e-2 * w.max()
    W = V[:, keep] / np.sqrt(w[keep])[None, :]
    r = {}
    for t in range(t0 + 1, Ct.shape[0]):
        M = W.T @ (0.5 * (Ct[t] + Ct[t].T)) @ W
        r[t] = (np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1], int(keep.sum()))
    return r


def main():
    for tag in sys.argv[1:]:
        f = os.path.join(HERE, f'gh3_{tag}.npz')
        if not os.path.exists(f):
            P(f"missing {f}")
            continue
        d = np.load(f)
        S = d['S']; Lt = int(d['Lt']); N = int(d['n']); Ls = int(d['Ls'])
        nlev = S.shape[1]
        P("=" * 108)
        P(f"GEVP subset scan   tag={tag}  beta={d['beta']}  Ls={Ls} Lt={Lt}  N={N}")
        P("=" * 108)
        nb = 20
        idx = np.array_split(np.arange(N), nb)
        for sub in SUBSETS:
            if max(sub) >= nlev:
                continue
            ns_ = len(sub)
            def build(rows):
                Ct = np.zeros((Lt, ns_, ns_))
                for a, i0 in enumerate(sub):
                    for b, j0 in enumerate(sub):
                        Ct[:, a, b] = corr_pair(S[rows, i0, 0, :], S[rows, j0, 0, :], Lt)
                return Ct
            Ct = build(np.arange(N))
            r0 = lam_max(Ct, 0)
            rank = r0[1][1] if 1 in r0 else 0
            P("")
            P(f"subset {{{', '.join(LNAMES[i] for i in sub)}}}  ({ns_} ops, Gram rank={rank})")
            P(f"{'t0':>3} {'m_eff(t) t=1..6  [20-block JK error]':>0}")
            for t0 in (0,):
                r = r0
                ms, es = [], []
                for t in range(1, 7):
                    l = r[t][0]
                    ms.append(-np.log(l[0]) / (t - t0) if l[0] > 0 else np.nan)
                    # jackknife
                    acc = []
                    for b in range(nb):
                        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
                        Cb = build(keep)
                        lb = lam_max(Cb, t0)[t][0]
                        acc.append(-np.log(lb[0]) / (t - t0) if lb[0] > 0 else np.nan)
                    v = np.array(acc, float)
                    m = np.nanmean(v)
                    es.append(np.sqrt((nb - 1) / nb * np.nansum((v - m) ** 2)))
                P("    t=1..6 : " + " ".join(f"{m:6.3f}" for m in ms))
                P("     JK    : " + " ".join(f"{e:6.3f}" for e in es))
                # lambda2 gap
                l2 = [r[t][0][1] / r[t][0][0] if len(r[t][0]) > 1 and r[t][0][0] > 0 else np.nan
                      for t in range(1, 7)]
                P("    lam2/lam1: " + " ".join(f"{v:6.3f}" for v in l2))
        open(os.path.join(HERE, f'gh3_gevp_scan_{tag}.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
