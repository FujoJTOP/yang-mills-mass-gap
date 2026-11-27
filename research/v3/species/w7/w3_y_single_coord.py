"""W3 -- y = m2^2 - m0^2 as a SINGLE coordinate, off the npz (zero new MC).

Read-only analysis of the archived gh3 correlator matrices:
  S[n, level=6, obs=5, Lt]  (level alpha = lev_alpha, obs = W11_p0/W11_px/W11_pxxx/W12_p0/W22_p0)

Estimator is byte-identical to gh3_gevp_mixed.py (whiten at t0=0, keep w>1e-2*wmax,
lam sorted desc, r(t) = lam2/lam1, Delta = -ln r, y = 2*m1*Delta + Delta^2, m1=1.546).

Three readings:
  (A) the 5 archived bases  -> validation against gh3_gevp_mixed_L{6,8}.log
  (B) the full operator pool (30 ops = 6 levels x 5 obs), one number per volume, no basis choice
  (C) the ensemble of ALL size-3 subsets of the 8 logged ops -> spread of y at fixed L
"""
import os, itertools, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M1 = 1.546
OBSN = ['W11_p0', 'W11_px', 'W11_pxxx', 'W12_p0', 'W22_p0']
LNAMES = ['raw', 'a10n1', 'a20n1', 'a30n1', 'a50n1', 'a10n3']
BASES = {
    'B1': [(0, 0), (4, 0)],
    'B2': [(0, 0), (0, 3), (0, 4)],
    'B3': [(0, 0), (4, 0), (0, 3), (0, 4)],
    'B4': [(0, 0), (0, 3), (0, 4), (3, 3), (3, 4)],
    'B5': [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
}
NB = 20


def corr_blocks(S, idx, Lt):
    """Ct for one config block, full 30x30, shape (Lt,30,30). S: [n,6,5,Lt]."""
    X = S[idx].reshape(len(idx), 30, Lt).astype(np.float64)
    X = X - X.mean(axis=0, keepdims=True)
    Ct = np.zeros((Lt, 30, 30))
    for t in range(Lt):
        Z = np.roll(X, -t, axis=2)
        Ct[t] = np.einsum('cai,cbi->ab', X, Z, optimize=True) / (len(idx) * Lt)
    return Ct


def lam_of(Ct, t0=0):
    G = 0.5 * (Ct[t0] + Ct[t0].T)
    w, V = np.linalg.eigh(G)
    keep = w > 1e-2 * w.max()
    W = V[:, keep] / np.sqrt(w[keep])[None, :]
    out = {}
    for t in range(t0 + 1, Ct.shape[0]):
        M = W.T @ (0.5 * (Ct[t] + Ct[t].T)) @ W
        out[t] = np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]
    return out, int(keep.sum()), w


def r_of(lam, t):
    l = lam[t]
    if len(l) < 2 or l[0] <= 0:
        return np.nan
    return l[1] / l[0]


def y_of(r):
    if not np.isfinite(r) or r <= 0:
        return np.nan
    D = -np.log(r)
    return 2 * M1 * D + D * D


def jk_r(Ctb, sel, t):
    """delete-one-block JK over the block-averaged Ct restricted to operator indices sel."""
    vals = []
    for b in range(len(Ctb)):
        Ct = np.mean([Ctb[j] for j in range(len(Ctb)) if j != b], axis=0)
        lam, _, _ = lam_of(Ct[:, sel][:, :, sel])
        vals.append(r_of(lam, t))
    v = np.array(vals, float)
    m = np.nanmean(v)
    e = np.sqrt((len(v) - 1) / len(v) * np.nansum((v - m) ** 2))
    return m, e


def main():
    res = {}
    for tag in ['L6', 'L8']:
        d = np.load(os.path.join(HERE, f'gh3_{tag}.npz'))
        S = d['S']
        Lt = int(d['Lt'])
        N = S.shape[0]
        idx = np.array_split(np.arange(N), NB)
        Ctb = [corr_blocks(S, ix, Lt) for ix in idx]
        Ct = np.mean(Ctb, axis=0)
        print(f"# {tag}: N={N} Lt={Lt} Ls={int(d['Ls'])} beta={d['beta']}")

        # ---- (A) archived 5 bases: byte-identical reproduction
        for name, pairs in BASES.items():
            sel = [l * 5 + o for l, o in pairs]
            lam, rank, w = lam_of(Ct[:, sel][:, :, sel])
            rr = [r_of(lam, t) for t in range(1, 7)]
            print(f"  {name} rank={rank}/{len(sel)} kappa={np.sort(w)[0]/np.sort(w)[-1]:.2e} "
                  f"r(t=1..6)=" + " ".join(f"{v:6.3f}" for v in rr))
        # ---- (B) full pool
        lam, rank, w = lam_of(Ct)
        print(f"  FULL 30 ops: kept rank={rank}  kappa={np.sort(w)[0]/np.sort(w)[-1]:.2e}")
        for t in (1, 2, 3):
            rr = r_of(lam, t)
            print(f"    t={t}: r={rr:.4f}  Delta={-np.log(rr):.4f}  y={y_of(rr):.4f}")
        # JK for full pool
        for t in (1, 2, 3):
            m, e = jk_r(Ctb, list(range(30)), t)
            print(f"    t={t}: JK r={m:.4f}+-{e:.4f}  sigma_y={(2*M1+2*(-np.log(m)))*e/m:.3f}")
        res[tag] = {}

        # ---- (C) all size-3 subsets of the 8 logged ops
        ops8 = sorted({l * 5 + o for pairs in BASES.values() for l, o in pairs})
        print(f"  #8 logged ops -> {len(ops8)}  all triples = {len(list(itertools.combinations(ops8,3)))}")
        for t in (1, 2, 3):
            ys, rs = [], []
            for tri in itertools.combinations(ops8, 3):
                sel = list(tri)
                lam3, _, _ = lam_of(Ct[:, sel][:, :, sel])
                rr = r_of(lam3, t)
                rs.append(rr)
                ys.append(y_of(rr))
            rs = np.array(rs); ys = np.array(ys)
            good = np.isfinite(ys)
            print(f"    t={t}: n={good.sum()}/{len(ys)} r>0 | r med={np.nanmedian(rs):.4f} "
                  f"[{np.nanmin(rs):.4f},{np.nanmax(rs):.4f}] | "
                  f"y med={np.nanmedian(ys):.3f} q16={np.nanpercentile(ys[good],16):.3f} "
                  f"q84={np.nanpercentile(ys[good],84):.3f} [{np.nanmin(ys):.3f},{np.nanmax(ys):.3f}]")
            res[tag][t] = dict(y_med=float(np.nanmedian(ys)), y_min=float(np.nanmin(ys)),
                               y_max=float(np.nanmax(ys)), r_med=float(np.nanmedian(rs)),
                               n=int(good.sum()), ntot=int(len(ys)))
    print(json.dumps(res))


if __name__ == '__main__':
    main()
