import sys
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import numpy as np
import gpp_core as C
from w7_final import width_scan

rng = np.random.default_rng(0)
for tag, fn in [('SU2', 'partA_su2_N4000.npz'), ('SU3', 'partA_su3_N1500.npz')]:
    d = np.load(fn)
    names = [k[:-2] for k in d.files if k.endswith('_a')]
    kaps, ses, EHs, VGs = [], [], [], []
    Ws = (0.10, 0.15, 0.20, 0.30, 0.50)
    acc_ev = []
    for n in names:
        a = d[n + '_a']; H = d[n + '_H']; G = d[n + '_G']
        kap, s, lam, _, inv = C.veff_curvature(a, rng=rng, nboot=200)
        kaps.append(kap); ses.append(s); EHs.append(H.mean()); VGs.append(G.var())
        sc = width_scan(a, H, G, ws=Ws)
        acc_ev.append([eh - vg for (w, eh, vg, c) in sc])
    kaps = np.array(kaps); ses = np.array(ses)
    w = 1.0 / ses ** 2
    mean = (kaps * w).sum() / w.sum(); err = 1.0 / np.sqrt(w.sum())
    print(f'{tag}: {len(names)} modes  kappa range [{kaps.min():+.3f},{kaps.max():+.3f}]  n_pos={int((kaps>0).sum())}/{len(names)}')
    print(f'   inv-var mean Gpp = {mean:+.4f} +- {err:.4f}  (mean SE {ses.mean():.3f}, E[H]={np.mean(EHs):.3f}, VarG={np.mean(VGs):.3f})')
    ev = np.nanmean(np.array(acc_ev, float), axis=0)
    print('   mode-avg Ev(w)=E[H|w]-Var(G|w):', ' '.join(f'{ww:.2f}:{e:+.3f}' for ww, e in zip(Ws, ev)))
