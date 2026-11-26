"""P19 combine — G''(0), Delta_dyn, Var(F) at fixed L=2 with two independent errors:
sigma_jk (delete-one-block jackknife, 25 blocks/chain) and sigma_tau (Geyer integrated
autocorrelation time), then the 4-group combination.

usage: py -3.12 p19_combine.py <tag>
reads %TEMP%/p19_<tag>.npz written by p19_sweep.py
"""
import os
import sys

import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else 'b2.5'
d = np.load(os.path.join(os.environ.get('TEMP', '.'), f'p19_{TAG}.npz'))
C = d['cols']
beta, ns, G = float(d['beta']), int(d['ns']), C.shape[1]
NB = 25


def acf(x, ml=600):
    y = x - x.mean()
    v = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / v for k in range(1, min(ml, len(y) - 2))])


def geyer(x, ml=600):
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


def f_gpp(c):
    return c[:, 0].mean() - c[:, 1].var()


def f_dyn(c):
    return c[:, 2].mean() - c[:, 3].var()


def f_var(c):
    return c[:, 1].var()


def block_jk(c, fn, nb=NB):
    n = len(c)
    idx = np.array_split(np.arange(n), nb)
    ds = np.array([fn(np.delete(c, b, axis=0)) for b in idx])
    return ds.mean(), np.sqrt((nb - 1) / nb * ((ds - ds.mean()) ** 2).sum())


def tau_err(c, fn):
    """delta method on the two ingredients, each with its own Geyer factor."""
    n = len(c)
    k, f = c[:, 0], c[:, 1]
    tk, tf = geyer(k), geyer(f)
    if fn is f_gpp or fn is f_dyn:
        a, b = (c[:, 0], c[:, 1]) if fn is f_gpp else (c[:, 2], c[:, 3])
        sk = np.sqrt(a.var(ddof=1) * (1 + 2 * geyer(a)) / n)
        sf = np.sqrt(2 * b.var(ddof=1) ** 2 * (1 + 2 * geyer(b)) / n)
        return np.sqrt(sk ** 2 + sf ** 2)
    v = f.var(ddof=1)
    return np.sqrt(2 * v ** 2 * (1 + 2 * tf) / n)


print(f"=== P19  beta={beta:g}  tag={TAG}  L=2  SU(2)  Nsweep/chain={ns}  chains={G}  "
      f"seeds={d['seeds']}  eps={float(d['eps']):g}  domain-reject={float(d['domrej']):.4f}  "
      f"wall={float(d['el']):.0f}s ({1e3*float(d['el'])/(ns+800):.2f} ms/sweep, {G} chains) ===")
print(f"acc_link={np.round(d['acc'],3)}")

rows = {}
for nm, fn in (('Gpp(0)', f_gpp), ('Delta_dyn', f_dyn), ('Var(F)', f_var)):
    rows[nm] = []
    for g in range(G):
        c = C[:, g, :]
        val, sjk = block_jk(c, fn)
        stau = tau_err(c, fn)
        rows[nm].append((val, sjk, stau))
    v = np.array([r[0] for r in rows[nm]])
    sjk = np.array([r[1] for r in rows[nm]])
    stau = np.array([r[2] for r in rows[nm]])
    m = v.mean()
    S_jk = np.sqrt((sjk ** 2).sum()) / G
    S_tau = np.sqrt((stau ** 2).sum()) / G
    dg = G * m - (G - 1) * v
    S_grp = np.sqrt((G - 1) / G * ((dg - dg.mean()) ** 2).sum())
    print(f"\n--- {nm} ---")
    print(f"   per chain: " + "  ".join(f"{r[0]:+.4f}" for r in rows[nm]))
    print(f"   {nm} = {m:+.4f} +- {S_jk:.4f}(jk) / {S_tau:.4f}(tau) / {S_grp:.4f}(grp)   "
          f"spread={v.std(ddof=1):.4f}  chi2/dof={np.sum((v-m)**2/sjk**2)/(G-1):.2f}")

tk = np.array([geyer(C[:, g, 0]) for g in range(G)])
tf = np.array([geyer(C[:, g, 1]) for g in range(G)])
print(f"\ntauK={np.round(tk,1)}  tauF={np.round(tf,1)}  Neff={np.round(ns/np.maximum(tk,tf)).astype(int)}")
kp = C[:, :, 0].ravel()
print(f"kappa pooled: mean={kp.mean():.4f} median={np.median(kp):.4f} max={kp.max():.3e} "
      f"min={kp.min():.3e} std={kp.std():.4f}")
gpp = [r[0] for r in rows['Gpp(0)']]
print(f"\nVERDICT beta={beta:g}: Gpp(0)={np.mean(gpp):+.4f} "
      f"(jk {np.sqrt(sum(r[1]**2 for r in rows['Gpp(0)']))/G:.4f}) "
      f"| Delta_dyn={np.mean([r[0] for r in rows['Delta_dyn']]):+.4f} "
      f"| Var(F)={np.mean([r[0] for r in rows['Var(F)']]):+.4f}")
