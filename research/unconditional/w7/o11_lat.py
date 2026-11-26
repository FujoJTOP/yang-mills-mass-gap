"""O11 lattice worker: one MH chain on the exact conditioning hyperplane {<v,X>=0}.

usage: py -3.12 o11_lat.py <mode:A|B> <seed> <nsweep> <eps>
prints a single result line (also written by the caller's redirect).
"""
import sys
import numpy as np
import time
import o11_physical as o

mode = sys.argv[1]
seed = int(sys.argv[2])
ns = int(sys.argv[3])
eps = float(sys.argv[4])

tab, near = o.lat_tables(2)
nlink = 64
v = np.zeros((nlink, 3))
if mode == 'A':
    v[:, 2] = 1.0 / np.sqrt(nlink)
else:
    v[0, 2] = 1.0

t = time.time()
x0 = np.random.default_rng(1000 + seed).normal(0, 0.6, (nlink, 3))
Fs, Ks, acc = o.chain(x0, tab, near, 2.5, v, ns, eps, np.random.default_rng(seed), burn=500)
d0, sg = o.jk_lat(Ks, Fs, nblk=20)
aK = o.ac(Ks, 300); aF = o.ac(Fs, 300)
tauK = 1 + 2 * aK[1:].sum(); tauF = 1 + 2 * aF[1:].sum()
neff = len(Ks) / max(tauK, tauF, 1.0)
sem_simple = np.sqrt(Ks.var() / neff)
print(f"mode={mode} seed={seed} eps={eps} nsweep={ns} acc={acc:.3f} "
      f"E[kappa]={Ks.mean():+.4f} varF={Fs.var():.4f} E[F]={Fs.mean():+.4f} "
      f"tauK={tauK:.1f} tauF={tauF:.1f} Neff={neff:.0f} "
      f"Delta={d0:+.4f} sigma_jk={sg:.4f} ratio={d0/max(sg,1e-12):+.1f} "
      f"t={time.time()-t:.1f}s", flush=True)
np.save(f"o11_lat_{mode}_{seed}.npy", np.array([Ks, Fs]))
