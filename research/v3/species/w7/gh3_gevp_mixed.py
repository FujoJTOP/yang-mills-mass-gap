"""GH3 -- mixed-geometry GEVP, off the npz (no MC).

The smearing-level family is degenerate (adjacent alpha overlap ~0.8, Gram rank 3/6), so a
GEVP built on it cannot separate states.  Here we build GEVP bases that mix *geometry*
(1x1, 1x2, 2x2 spatial loops, i.e. obs 0/3/4) with smearing level, which adds genuinely
independent directions.

Run: py -3.12 gh3_gevp_mixed.py L6 L8
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
LNAMES = ['raw', 'a10n1', 'a20n1', 'a30n1', 'a50n1', 'a10n3']
OBSN = ['W11_p0', 'W11_px', 'W11_pxxx', 'W12_p0', 'W22_p0']

BASES = {
    'B1 paper: 11(raw)+11(a50n1)':        [(0, 0), (4, 0)],
    'B2 geometry: 11+12+22 (raw)':        [(0, 0), (0, 3), (0, 4)],
    'B3 geom+smear: 11r,11a50,12r,22r':   [(0, 0), (4, 0), (0, 3), (0, 4)],
    'B4 geom x2: +22(a30n1),12(a30n1)':   [(0, 0), (0, 3), (0, 4), (3, 3), (3, 4)],
    'B5 smearing only: 6 alpha levels':   [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
}


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def corr_pair(Si, Sj, Lt):
    a = Si - Si.mean(axis=0)
    b = Sj - Sj.mean(axis=0)
    return np.array([np.mean(np.sum(a * np.roll(b, -t, axis=1), axis=1)) / Lt for t in range(Lt)])


def build(S, rows, pairs, Lt):
    n = len(pairs)
    Ct = np.zeros((Lt, n, n))
    for a, (la, oa) in enumerate(pairs):
        for b, (lb, ob) in enumerate(pairs):
            Ct[:, a, b] = corr_pair(S[rows, la, oa, :], S[rows, lb, ob, :], Lt)
    return Ct


def gevp(Ct, t0):
    G = 0.5 * (Ct[t0] + Ct[t0].T)
    w, V = np.linalg.eigh(G)
    keep = w > 1e-2 * w.max()
    W = V[:, keep] / np.sqrt(w[keep])[None, :]
    rank = int(keep.sum())
    lam = {}
    for t in range(t0 + 1, Ct.shape[0]):
        M = W.T @ (0.5 * (Ct[t] + Ct[t].T)) @ W
        lam[t] = np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]
    return lam, rank, w


def main():
    for tag in sys.argv[1:]:
        f = os.path.join(HERE, f'gh3_{tag}.npz')
        if not os.path.exists(f):
            P(f"missing {f}")
            continue
        d = np.load(f)
        S = d['S']; Lt = int(d['Lt']); N = int(d['n']); Ls = int(d['Ls'])
        P("=" * 110)
        P(f"mixed-geometry GEVP   tag={tag}  beta={d['beta']}  Ls={Ls} Lt={Lt}  N={N}")
        P("=" * 110)
        nb = 20
        idx = np.array_split(np.arange(N), nb)
        for name, pairs in BASES.items():
            if any(o >= S.shape[2] or l >= S.shape[1] for l, o in pairs):
                continue
            Ct = build(S, np.arange(N), pairs, Lt)
            lam, rank, w = gevp(Ct, 0)
            P("")
            P(f"{name}   ({len(pairs)} ops, Gram rank={rank}/"
              f"{len(pairs)}, kappa={np.sort(w)[0]/np.sort(w)[-1]:.2e})")
            P("    ops: " + ", ".join(f"{LNAMES[l]}/{OBSN[o]}" for l, o in pairs))
            ms, es, mono = [], [], []
            for t in range(1, 7):
                l = lam[t]
                ms.append(-np.log(l[0]) / t if l[0] > 0 else np.nan)
                acc = []
                for b in range(nb):
                    keep = np.concatenate([idx[j] for j in range(nb) if j != b])
                    lb, _, _ = gevp(build(S, keep, pairs, Lt), 0)
                    acc.append(-np.log(lb[t][0]) / t if lb[t][0] > 0 else np.nan)
                v = np.array(acc, float)
                m = np.nanmean(v)
                es.append(np.sqrt((nb - 1) / nb * np.nansum((v - m) ** 2)))
                mono.append(l[0])
            P("    m_eff t=1..6: " + " ".join(f"{m:6.3f}" for m in ms))
            P("    JK         : " + " ".join(f"{e:6.3f}" for e in es))
            P("    lam_max(t) : " + " ".join(f"{v:8.2e}" for v in mono) +
              ("   [MONOTONE: single-exponential OK]" if all(mono[i] > mono[i + 1] for i in range(5))
               else "   [NON-monotone -> GEVP lambda is not a clean exponential]"))
            if len(pairs) > 1:
                P("    lam2/lam1  : " + " ".join(
                    f"{(lam[t][1]/lam[t][0] if lam[t][0] > 0 else np.nan):6.3f}" for t in range(1, 7)))
        open(os.path.join(HERE, f'gh3_gevp_mixed_{tag}.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
