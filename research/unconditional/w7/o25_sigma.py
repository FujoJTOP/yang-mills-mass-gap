"""O25-L3 auxiliary: tau-corrected sigma for Delta = E[kappa]-Var(F), from the saved chains.

usage: py -3.12 o25_sigma.py
"""
import glob
import numpy as np


def ac(x, maxlag=300):
    x = x - x.mean()
    v = (x * x).sum()
    m = min(maxlag, len(x) // 4)
    return np.array([1.0] + [abs((x[:-k] * x[k:]).sum() / v) for k in range(1, m + 1)])


def tau_of(x):
    a = ac(x)
    return 1.0 + 2.0 * a[1:].sum()


rows = []
for f in sorted(glob.glob("o25_L*_*.npy")):
    d = np.load(f)
    K, F = d[0], d[1]
    tK, tF = tau_of(K), tau_of(F)
    n = len(K)
    nK, nF = n / max(tK, 1.0), n / max(tF, 1.0)
    # delta-method: Delta = E[K] - Var(F); d Var(F) = (2/N) sum (F-<F>)(F^2-<F^2>) ... use the
    # block-jackknife value as the primary; here the tau-corrected independent-sample estimate:
    g = (F - F.mean()) ** 2
    sg2 = np.sqrt(K.var() / nK + g.var() / nF)
    s0 = sg2 / np.sqrt(2.0)          # Var(F) and E[K] are ~independently driven (conservative /2^0.5)
    rows.append((f, len(K), K.mean(), F.var(), K.mean() - F.var(), tK, tF,
                 nK, nF, s0, K.mean() / (16 if "_L2_" in f else 81),
                 F.var() / (16 if "_L2_" in f else 81)))
    print(f"{f:26s} n={len(K):5d} EK={K.mean():+.4f} VF={F.var():.4f} "
          f"D={K.mean()-F.var():+.4f} tauK={tK:6.1f} tauF={tF:6.1f} "
          f"Neff={min(nK,nF):6.0f} sig_tau={s0:.4f}", flush=True)

for lab, pref in (("L=2", "o25_L2_"), ("L=3", "o25_L3_")):
    sel = [r for r in rows if pref in r[0]]
    if not sel:
        continue
    m = len(sel)
    D = np.mean([r[4] for r in sel])
    s = np.sqrt(np.mean([r[9] ** 2 for r in sel])) / np.sqrt(m)
    EK = np.mean([r[2] for r in sel]); VF = np.mean([r[3] for r in sel])
    print(f"\n{lab}: {m} groups, EK={EK:+.4f} VF={VF:.4f} Delta={D:+.4f} +- {s:.4f} "
          f"({D/max(s,1e-12):+.2f} sigma, tau-corrected)  EK/L4={EK/(16 if lab=='L=2' else 81):+.5f} "
          f"VF/L4={VF/(16 if lab=='L=2' else 81):+.5f}", flush=True)
