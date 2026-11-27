"""W3b -- the basis-PAIRED ratio y(8)/y(6): contamination cancels, basis bias drops out.

For every size-3 subset of the 8 logged operators, compute y = 2*m1*D + D^2 at L=6 and L=8
(the SAME subset at both volumes), with 20-block delete-one JK, and form ratio = y8/y6.
Collapse template 6pi^2/L^2 predicts 6^2/8^2 = 0.5625 ; constant branch predicts 1.0.

Also reports the 5 archived bases and the full 30-operator pool.
"""
import os, itertools, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M1 = 1.546
NB = 20
BASES = {
    'B1': [(0, 0), (4, 0)],
    'B2': [(0, 0), (0, 3), (0, 4)],
    'B3': [(0, 0), (4, 0), (0, 3), (0, 4)],
    'B4': [(0, 0), (0, 3), (0, 4), (3, 3), (3, 4)],
    'B5': [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
}


def corr_blocks(S, idx, Lt):
    X = S[idx].reshape(len(idx), 30, Lt).astype(np.float64)
    X = X - X.mean(axis=0, keepdims=True)
    Ct = np.zeros((Lt, 30, 30))
    for t in range(Lt):
        Z = np.roll(X, -t, axis=2)
        Ct[t] = np.einsum('cai,cbi->ab', X, Z, optimize=True) / (len(idx) * Lt)
    return Ct


def lam_of(Ct):
    G = 0.5 * (Ct[0] + Ct[0].T)
    w, V = np.linalg.eigh(G)
    keep = w > 1e-2 * w.max()
    W = V[:, keep] / np.sqrt(w[keep])[None, :]
    M = W.T @ (0.5 * (Ct[1] + Ct[1].T)) @ W
    return np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]


def y_from_sel(Ctb, sel):
    vals = []
    for b in range(len(Ctb)):
        Ct = np.mean([Ctb[j] for j in range(len(Ctb)) if j != b], axis=0)
        l = lam_of(Ct[:, sel][:, :, sel])
        vals.append(l[1] / l[0] if (len(l) > 1 and l[0] > 0 and l[1] > 0) else np.nan)
    v = np.array(vals, float)
    m = np.nanmean(v)
    e = np.sqrt((len(v) - 1) / len(v) * np.nansum((v - m) ** 2))
    D = -np.log(m) if m > 0 else np.nan
    y = 2 * M1 * D + D * D if np.isfinite(D) else np.nan
    sy = (2 * M1 + 2 * D) * e / m if (m > 0 and np.isfinite(D)) else np.nan
    return y, sy, m


def load(tag):
    d = np.load(os.path.join(HERE, f'gh3_{tag}.npz'))
    S = d['S']; Lt = int(d['Lt']); N = S.shape[0]
    idx = np.array_split(np.arange(N), NB)
    return [corr_blocks(S, ix, Lt) for ix in idx]


def main():
    C6, C8 = load('L6'), load('L8')
    out = {}
    print("=== paired ratio on the 5 archived bases ===")
    for name, pairs in BASES.items():
        sel = [l * 5 + o for l, o in pairs]
        y6, s6, r6 = y_from_sel(C6, sel)
        y8, s8, r8 = y_from_sel(C8, sel)
        rat = y8 / y6 if (y6 == y6 and y6 != 0) else np.nan
        srat = abs(rat) * np.hypot(s8 / y8, s6 / y6) if np.isfinite(rat) else np.nan
        print(f"{name}: r6={r6:.3f} r8={r8:.3f} y6={y6:7.3f}+-{s6:5.3f} y8={y8:7.3f}+-{s8:5.3f} "
              f"ratio={rat:6.3f}+-{srat:5.3f}")
        out[name] = dict(y6=y6, s6=s6, y8=y8, s8=s8, ratio=rat, srat=srat)

    print("\n=== full 30-op pool ===")
    for tag, C in (('L6', C6), ('L8', C8)):
        y, s, r = y_from_sel(C, list(range(30)))
        print(f"{tag}: r={r:.4f} y={y:.3f}+-{s:.3f}")
        out['FULL_' + tag] = dict(y=y, s=s, r=r)

    print("\n=== all size-3 subsets of the 8 logged ops (paired) ===")
    ops8 = sorted({l * 5 + o for pairs in BASES.values() for l, o in pairs})
    rats, ys6, ys8, ok = [], [], [], 0
    for tri in itertools.combinations(ops8, 3):
        sel = list(tri)
        y6, s6, _ = y_from_sel(C6, sel)
        y8, s8, _ = y_from_sel(C8, sel)
        if np.isfinite(y6) and np.isfinite(y8) and y6 > 0 and y8 > 0:
            rat = y8 / y6
            srat = rat * np.hypot(s8 / y8, s6 / y6)
            rats.append((rat, srat)); ys6.append(y6); ys8.append(y8); ok += 1
    rats = np.array(rats)
    r = rats[:, 0]
    print(f"n={ok}/120 usable | ratio: med={np.median(r):.3f} mean={r.mean():.3f} "
          f"q16={np.percentile(r,16):.3f} q84={np.percentile(r,84):.3f} "
          f"min={r.min():.3f} max={r.max():.3f}")
    for lo, hi, lab in ((0, 0.5625, '<0.5625 (超塌)'), (0.5625, 0.65, '塌缩支 0.5625~8%'),
                        (0.65, 1.0, '中间'), (1.0, 9, '>=1 (常/升)')):
        print(f"   {lab}: {int(((r >= lo) & (r < hi)).sum())}/120")
    print(f"   within 2sig of 0.5625: {int((np.abs(r-0.5625) <= 2*rats[:,1]).sum())}/120")
    print(f"   within 2sig of 1.0000: {int((np.abs(r-1.0) <= 2*rats[:,1]).sum())}/120")
    # weighted mean ratio
    w = 1 / rats[:, 1] ** 2
    rm = float((r * w).sum() / w.sum()); srm = float(1 / np.sqrt(w.sum()))
    print(f"   inverse-variance mean ratio = {rm:.3f} +- {srm:.3f}  "
          f"(chi2={float((((r-rm)**2)*w).sum()):.1f}/{len(r)-1} dof)")
    print(f"   y6 med={np.median(ys6):.3f} y8 med={np.median(ys8):.3f}")
    out['ENS'] = dict(n=ok, med=float(np.median(r)), q16=float(np.percentile(r, 16)),
                      q84=float(np.percentile(r, 84)), wmean=rm, swmean=srm,
                      chi2=float((((r - rm) ** 2) * w).sum()))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
