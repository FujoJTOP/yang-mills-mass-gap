#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""p14_cov_an.py -- statistics-only pass over the saved p14 records.

Fixes the one real bug of the online pass: pooling the 4 chains with reshape() INTERLEAVES
them, so a Geyer IPS on the pooled series sees acf(1)~0 and returns the tau=1 floor.
Here tau is estimated PER CHAIN, and four independent error bars are reported:
  jk     : jackknife over 20 blocks/chain (=200-sweep contiguous blocks, jk on 100 blocks)
  tau    : sqrt(Var(G Phi') tau_int / n) per chain, combined over chains
  batch  : scatter of 20 batch means per chain, combined
  group  : 4 fully independent chains (seed to seed)
usage: py -3.12 p14_cov_an.py [npz ...]
"""
import sys, glob
import numpy as np

P = lambda *a: print(*a, flush=True)


def acf(y, ml):
    y = y - y.mean()
    v = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / v for k in range(1, ml + 1)])


def tau_ips(y, ml=3000):
    r = np.concatenate(([1.0], acf(y, min(ml, len(y) // 2))))
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    s = s[s > 0]
    if s.size == 0:
        return 1.0
    return max(2.0 * min(s.sum(), 10 * s[0]), 1.0)      # Sokal-style cap at 10*S(0)


def per_chain(y):
    """y=(n,4): G, Phi', Phi'', H -> dict of Cov and friends (exact zero means by parity)"""
    G, Pp = y[:, 0], y[:, 1]
    z = (G - G.mean()) * (Pp - Pp.mean())
    n = y.shape[0]
    nb = 20
    k = n // nb
    bm = z[:k * nb].reshape(nb, k).mean(1)
    tG, tP, tZ = tau_ips(G), tau_ips(Pp), tau_ips(z)
    mG, mP = G.mean(), Pp.mean()
    cov = z.mean()
    vG, vP = G.var(), Pp.var()
    # jackknife over the nb batch means
    jk = (z.sum() - bm * k) / (n - k)
    s_jk = np.sqrt((nb - 1) / nb * ((jk - jk.mean()) ** 2).sum())
    return dict(cov=cov, vG=vG, vP=vP, vF=(G + Pp).var(), mG=mG, mP=mP,
                Phi2=y[:, 2].mean(), H=y[:, 3].mean(),
                tG=tG, tP=tP, tZ=tZ, s_tau=np.sqrt(z.var() * tZ / n),
                s_jk=s_jk, s_bm=bm.std(ddof=1) / np.sqrt(nb))


if __name__ == "__main__":
    files = sys.argv[1:] or ["p14_main.npz"]
    st = []
    for f in files:
        d = np.load(f)
        r = d['rec']
        P(f"[{f}] ns={int(d['ns'])} chains={r.shape[1]} eps={float(d['eps'])} "
          f"acc={np.round(d['acc'],3)} seeds={d['seeds']}")
        for g in range(r.shape[1]):
            st.append(per_chain(r[:, g, :4]))
    G = len(st)
    ns = None
    keys = ["cov", "vG", "vP", "vF", "mG", "mP", "Phi2", "H", "tG", "tP", "tZ", "s_tau", "s_jk", "s_bm"]
    P(f"\n{G} independent chains ({files})")
    P("  chain  " + "  ".join(f"{k:>9s}" for k in keys))
    for g, s in enumerate(st):
        P(f"  {g:5d}  " + "  ".join(f"{s[k]:+9.4f}" for k in keys))
    sw = np.array([s['cov'] for s in st])
    mean = sw.mean()
    s_grp = sw.std(ddof=1) / np.sqrt(G)
    s_jk = np.sqrt(np.mean([s['s_jk'] ** 2 for s in st]) / G)
    s_tau = np.sqrt(np.mean([s['s_tau'] ** 2 for s in st]) / G)
    s_bm = np.sqrt(np.mean([s['s_bm'] ** 2 for s in st]) / G)
    P(f"\n  Cov(G_v,Phi') = {mean:+.6f}  +- {s_jk:.6f}(jk)  +- {s_tau:.6f}(tau)  "
      f"+- {s_bm:.6f}(batch)  +- {s_grp:.6f}(group, {G} chains)")
    P(f"  -> /3sigma: jk {mean/s_jk:+.2f}  tau {mean/s_tau:+.2f}  batch {mean/s_bm:+.2f}  "
      f"group {mean/s_grp:+.2f}")
    V = np.array([[s['vG'], s['vP'], s['vF'], s['Phi2'], s['H'], s['mG'], s['mP']] for s in st])
    nm = ["VarG", "VarPhi'", "Var(G+Phi')", "E[Phi'']", "E[H_vv]", "E[G]", "E[Phi']"]
    P("\n  " + "   ".join(f"{n}={v.mean():+.4f}(+-{v.std(ddof=1)/np.sqrt(G):.4f})"
                          for n, v in zip(nm, V.T)))
    P(f"  corr = {mean/np.sqrt(V[:,0].mean()*V[:,1].mean()):+.4f}   "
      f"sigma_G sigma_Phi' = {np.sqrt(V[:,0].mean()*V[:,1].mean()):.4f}  (C-S: |Cov| <= this)")
    P(f"  W_full = E[Phi'']+VarG-Var(G+Phi') = {V[:,3].mean()+V[:,0].mean()-V[:,2].mean():+.4f}")
    P(f"  cross part -VarPhi'-2Cov = {-(V[:,1].mean()+2*mean):+.4f}  "
      f"[ -VarPhi'={-V[:,1].mean():+.4f} , -2Cov={-2*mean:+.4f} ]")
    P(f"  2|Cov|/Var(G+Phi') = {2*abs(mean)/V[:,2].mean()*100:.2f}%   "
      f"2|Cov|/E[Phi''] = {2*abs(mean)/V[:,3].mean()*100:.2f}%")
    P(f"  bound: Cov >= {mean-3*s_tau:+.5f} (3 sigma, tau)   Cov <= {mean+3*s_tau:+.5f}")
