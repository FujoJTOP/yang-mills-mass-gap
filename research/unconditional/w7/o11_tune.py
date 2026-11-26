import numpy as np, time
import o11_physical as o
tab, near = o.lat_tables(2)
v = np.zeros((64,3)); v[:,2] = 1/8
for eps in (0.5, 0.9, 1.5, 2.5):
    Fs,Ks,acc = o.chain(np.random.default_rng(5).normal(0,0.6,(64,3)), tab, near, 2.5, v, 400, eps, np.random.default_rng(9), burn=200)
    aK=o.ac(Ks,200); aF=o.ac(Fs,200)
    tK=1+2*aK[1:].sum(); tF=1+2*aF[1:].sum()
    print(f"eps={eps} acc={acc:.2f} K={Ks.mean():+.3f} varF={Fs.var():.3f} E[F]={Fs.mean():+.3f} tauK={tK:.1f} tauF={tF:.1f} Neff={len(Ks)/max(tK,tF):.0f} d={Ks.mean()-Fs.var():+.3f}")
