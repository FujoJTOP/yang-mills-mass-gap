"""Q13 combine — G''(0) at L=2, SU(2), beta=6 with 8 independent chains.

Two independent errors per chain (same conventions as p19_combine.py):
  sigma_jk  = delete-one-block jackknife (25 and 40 blocks; block 800 / 500 sweeps > tau_K)
  sigma_tau = delta method with Geyer integrated autocorrelation time
then the G-group combination (sqrt(sum s^2)/G, plus a group jackknife over chains),
the sigma level, and the sweep count needed for 5 sigma.

usage: py -3.12 q13_combine.py [tag]
reads %TEMP%/q13_<tag>.npz written by q13_beta6.py
"""
import os
import sys

import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else 'b6'
d = np.load(os.path.join(os.environ.get('TEMP', '.'), f'q13_{TAG}.npz'))
C = d['cols']
beta, ns, G = float(d['beta']), int(d['ns']), C.shape[1]
ML = 2000


def acf(x, ml=ML):
    y = x - x.mean()
    v = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / v for k in range(1, min(ml, len(y) - 2))])


def geyer(x, ml=ML):
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


def block_jk(c, fn, nb):
    n = len(c)
    idx = np.array_split(np.arange(n), nb)
    ds = np.array([fn(np.delete(c, b, axis=0)) for b in idx])
    return ds.mean(), np.sqrt((nb - 1) / nb * ((ds - ds.mean()) ** 2).sum())


def tau_err(c, fn, ml=ML):
    n = len(c)
    if fn is f_gpp or fn is f_dyn:
        a, b = (c[:, 0], c[:, 1]) if fn is f_gpp else (c[:, 2], c[:, 3])
        sk = np.sqrt(a.var(ddof=1) * (1 + 2 * geyer(a, ml)) / n)
        sf = np.sqrt(2 * b.var(ddof=1) ** 2 * (1 + 2 * geyer(b, ml)) / n)
        return np.sqrt(sk ** 2 + sf ** 2)
    v = c[:, 1].var(ddof=1)
    return np.sqrt(2 * v ** 2 * (1 + 2 * geyer(c[:, 1], ml)) / n)


print(f"=== Q13  beta={beta:g}  tag={TAG}  L=2  SU(2)  Nsweep/chain={ns}  chains={G}  "
      f"seeds={d['seeds']}  eps={float(d['eps']):g}  domain-reject={float(d['domrej']):.4f}  "
      f"wall={float(d['el']):.0f}s ({1e3 * float(d['el']) / (ns + int(d['burn'])):.2f} ms/sweep, "
      f"{G} chains) ===")
print(f"acc_link={np.round(d['acc'], 3)}")

res = {}
for nm, fn in (('Gpp(0)', f_gpp), ('Delta_dyn', f_dyn), ('Var(F)', f_var)):
    rows = []
    for g in range(G):
        c = C[:, g, :]
        v25, s25 = block_jk(c, fn, 25)
        v40, s40 = block_jk(c, fn, 40)
        rows.append((v25, s25, s40, tau_err(c, fn)))
    v = np.array([r[0] for r in rows])
    s25 = np.array([r[1] for r in rows])
    s40 = np.array([r[2] for r in rows])
    st = np.array([r[3] for r in rows])
    m = v.mean()
    S25 = np.sqrt((s25 ** 2).sum()) / G
    S40 = np.sqrt((s40 ** 2).sum()) / G
    St = np.sqrt((st ** 2).sum()) / G
    dg = G * m - (G - 1) * v
    Sgr = np.sqrt((G - 1) / G * ((dg - dg.mean()) ** 2).sum())
    res[nm] = (m, S25, S40, St, Sgr, v, s25)
    print(f"\n--- {nm} ---")
    print("   per chain: " + "  ".join(f"{r[0]:+.4f}" for r in rows))
    print(f"   {nm} = {m:+.4f} +- {S25:.4f}(jk25) / {S40:.4f}(jk40) / {St:.4f}(tau) / "
          f"{Sgr:.4f}(grp)   spread={v.std(ddof=1):.4f}  "
          f"chi2/dof={np.sum((v - m) ** 2 / s25 ** 2) / (G - 1):.2f}")

m, S25, S40, St, Sgr, v, s25 = res['Gpp(0)']

# ---- drift / equilibration check: first half vs second half of every chain ----
h1 = C[:ns // 2, :, :]
h2 = C[ns // 2:, :, :]
a1 = np.array([f_gpp(h1[:, g, :]) for g in range(G)])
a2 = np.array([f_gpp(h2[:, g, :]) for g in range(G)])
sa1 = np.array([block_jk(h1[:, g, :], f_gpp, 12)[1] for g in range(G)])
sa2 = np.array([block_jk(h2[:, g, :], f_gpp, 12)[1] for g in range(G)])
dh = a2 - a1
sd = np.sqrt((sa1 ** 2 + sa2 ** 2).sum()) / G
print(f"half-split: 1st={a1.mean():+.4f} 2nd={a2.mean():+.4f} diff={dh.mean():+.4f} +- {sd:.4f} "
      f"({dh.mean() / sd:+.1f} sigma)")

# ---- batch means: k batches of ns//k sweeps per chain (valid because batch >> tau_K) ----
KBM = 32
nbm = ns // KBM


def batch_means(fn, k=KBM):
    seg = np.array([[fn(C[i * nbm:(i + 1) * nbm, g, :]) for i in range(k)] for g in range(G)])
    flat = seg.ravel()
    return flat.mean(), flat.std(ddof=1) / np.sqrt(flat.size), seg


print(f"\n--- batch means ({KBM} batches x {nbm} sweeps/chain, {G * KBM} batch values) ---")
for nm, fn in (('Gpp(0)', f_gpp), ('Delta_dyn', f_dyn), ('Var(F)', f_var)):
    mb, sb, seg = batch_means(fn)
    perchn = seg.std(axis=1, ddof=1) / np.sqrt(KBM)
    print(f"   {nm} = {mb:+.4f} +- {sb:.4f}(bm)   per-chain sem="
          + " ".join(f"{q:.3f}" for q in perchn))
    if nm == 'Gpp(0)':
        BM = (mb, sb)

print(f"\nsigma level (Gpp(0)/error):  {m / S25:.1f} (jk25)   {m / S40:.1f} (jk40)   "
      f"{m / St:.1f} (tau)   {BM[0] / BM[1]:.1f} (batch-means)")

# ---- tau / Neff table ----
tk = np.array([geyer(C[:, g, 0]) for g in range(G)])
tf = np.array([geyer(C[:, g, 1]) for g in range(G)])
tk6 = np.array([geyer(C[:, g, 0], 600) for g in range(G)])
tf6 = np.array([geyer(C[:, g, 1], 600) for g in range(G)])
neff = ns / np.maximum(tk, tf)
print(f"tauK={np.round(tk, 1)}  (ml=600: {np.round(tk6, 1)})")
print(f"tauF={np.round(tf, 1)}  (ml=600: {np.round(tf6, 1)})")
print(f"Neff per chain={np.round(neff).astype(int)}  total Neff={neff.sum():.0f}")
kp = C[:, :, 0].ravel()
print(f"kappa pooled: mean={kp.mean():.4f} median={np.median(kp):.4f} max={kp.max():.3e} "
      f"min={kp.min():.3e} std={kp.std():.4f}  frac(kappa<0)={np.mean(kp < 0):.3f}")

# ---- sweeps needed for 5 sigma (errors scale as 1/sqrt(N)) ----
print("\n--- needs for 5 sigma (same G=%d chains, errors ~ ns^-1/2) ---" % G)
for nm in ('jk25', 'jk40', 'tau', 'bm'):
    s = {'jk25': S25, 'jk40': S40, 'tau': St, 'bm': BM[1]}[nm]
    k = (5.0 * s / m) ** 2
    print(f"   from sigma_{nm}={s:.4f}: ns_needed = {ns * k:.3g} sweeps/chain "
          f"(x{k:.2f} of today, wall ~ {float(d['el']) * k / 60:.0f} min)")

# ---- finite-difference-in-a alternative: cost estimate (analytic, uses measured Var(F), tau_F) ----
nn = ns
varF = C[:, :, 1].var(ddof=1)
for aa in (0.1, 0.2, 0.3):
    rel = np.sqrt((np.exp(aa ** 2 * varF) - 1) * (1 + 2 * tf.mean()) / nn)
    print(f"   FD-in-a at a={aa}: reweight rel.err={rel:.4f} -> sigma_Gpp="
          f"{2 * rel / aa ** 2:.2f}  (vs jk {S25:.3f})")

print(f"\nVERDICT beta={beta:g}: Gpp(0)={m:+.4f} (jk25 {S25:.4f} / tau {St:.4f}) "
      f"| Delta_dyn={res['Delta_dyn'][0]:+.4f} | Var(F)={res['Var(F)'][0]:+.4f} "
      f"| sigma={m / S25:.1f}(jk) {m / St:.1f}(tau)")
