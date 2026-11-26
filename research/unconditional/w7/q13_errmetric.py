"""Q13 errmetric — why does sigma_tau (Geyer/delta method) sit ~2.5x above the block/scatter errors?

For each chain:
    s_kappa = sqrt(Var(kappa)(1+2 tau_kappa)/n)
    s_F_old = sqrt(2 Var(F)^2 (1+2 tau_F)/n)          <-- p19_combine.py / q13_combine.py form
    s_F_new = sqrt(Var(z)(1+2 tau_z)/n),  z=(F-mean F)^2   <-- correct decorrelation of the
                                                             sample-variance FLUCTUATION
Then the combined errors, plus a batch-means check with batches of ns/4 sweeps (>> tau_K).
Also the 8-chain independent scatter (tau-free) and its leave-one-chain-out jackknife.

usage: py -3.12 q13_errmetric.py [tag]
"""
import os
import sys

import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else 'b6'
d = np.load(os.path.join(os.environ.get('TEMP', '.'), f'q13_{TAG}.npz'))
C, ns = d['cols'], int(d['ns'])
G = C.shape[1]


def acf(x, ml=2000):
    y = x - x.mean()
    v = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / v for k in range(1, min(ml, len(y) - 2))])


def gey(x, ml=2000):
    r = np.concatenate(([1.0], acf(x, ml)))
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    tot = 0.0
    for si in s:
        if si <= 0:
            break
        tot += si
    return max(2 * tot - 1.0, 1.0)


def gpp(c):
    return c[:, 0].mean() - c[:, 1].var()


print(f"=== Q13 errmetric  beta={float(d['beta']):g}  ns={ns}  G={G} ===")
print(f"{'g':>2} {'tauK':>7} {'tauF':>6} {'tau_z':>7} {'VarK':>7} {'VarF':>6} {'Var(z)':>7} "
      f"{'s_k':>7} {'s_F_old':>8} {'s_F_new':>8}")
sk, so, sn = [], [], []
for g in range(G):
    k, f = C[:, g, 0], C[:, g, 1]
    z = (f - f.mean()) ** 2
    tk, tf, tz = gey(k), gey(f), gey(z)
    vk, vf, vz = k.var(ddof=1), f.var(ddof=1), z.var(ddof=1)
    a = np.sqrt(vk * (1 + 2 * tk) / ns)
    b = np.sqrt(2 * vf ** 2 * (1 + 2 * tf) / ns)
    c = np.sqrt(vz * (1 + 2 * tz) / ns)
    sk.append(a), so.append(b), sn.append(c)
    print(f"{g:>2} {tk:7.1f} {tf:6.1f} {tz:7.1f} {vk:7.4f} {vf:6.4f} {vz:7.2f} "
          f"{a:7.4f} {b:8.4f} {c:8.4f}")
sk, so, sn = np.array(sk), np.array(so), np.array(sn)
tot_old = np.sqrt((sk ** 2 + so ** 2).sum()) / G
tot_new = np.sqrt((sk ** 2 + sn ** 2).sum()) / G
print(f"\nsigma_tau OLD (2Var(F)^2, tau_F) = {tot_old:.4f}   <- the 1.3-sigma-of-p19 form")
print(f"sigma_tau NEW (Var(z), tau_z)    = {tot_new:.4f}")
print(f"   kappa piece = {np.sqrt((sk ** 2).sum()) / G:.4f}   "
      f"Var(F) piece old = {np.sqrt((so ** 2).sum()) / G:.4f} -> new = {np.sqrt((sn ** 2).sum()) / G:.4f}")

n4 = ns // 4
seg = np.array([[gpp(C[i * n4:(i + 1) * n4, g, :]) for i in range(4)] for g in range(G)])
print(f"\nbatch means k=4/chain ({n4} sweeps ~ {n4 / 200:.0f} tauK): "
      f"mean={seg.mean():+.4f} se={seg.std(ddof=1) / np.sqrt(seg.size):.4f}")

v = np.array([gpp(C[:, g, :]) for g in range(G)])
s = v.std(ddof=1)
print(f"8 independent chains: mean={v.mean():+.4f} s={s:.4f} SEM={s / np.sqrt(G):.4f} "
      f"-> {v.mean() / (s / np.sqrt(G)):.1f} sigma")
print(f"leave-one-chain-out jackknife SE = {s / np.sqrt(G):.4f} (same thing, exact for iid chains)")
print(f"all 8 chains positive: {bool((v > 0).all())}   min={v.min():+.4f} max={v.max():+.4f}   "
      f"sign-test p={2.0 ** -G:.2e}")
# 95% upper edge of the SEM itself (chi2 with G-1 dof)
import scipy.stats as st
lo, hi = np.sqrt((G - 1) * s ** 2 / st.chi2.ppf(0.975, G - 1)), np.sqrt((G - 1) * s ** 2 / st.chi2.ppf(0.025, G - 1))
print(f"SEM 95% interval from the chain scatter: [{lo / np.sqrt(G):.4f}, {hi / np.sqrt(G):.4f}] "
      f"-> sigma in [{v.mean() / (hi / np.sqrt(G)):.1f}, {v.mean() / (lo / np.sqrt(G)):.1f}]")
