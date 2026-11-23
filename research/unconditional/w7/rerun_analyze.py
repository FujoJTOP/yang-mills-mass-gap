"""Aggregate the rerun_*.npz into the final W7-RERUN numbers."""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import rerun_core as RC
import rerun_run as R


def summarize(fn, group, label):
    d = np.load(fn)
    A, H, G, P1, P2 = d['a'], d['H'], d['G'], d['P1'], d['P2']
    Nm = A.shape[1]
    rng = np.random.default_rng(7000)
    kk, kkse, kh, inv, br, W, EP2, EH, VG, tot = ([] for _ in range(10))
    for j in range(Nm):
        r = R.analyse_mode(A[:, j], H[:, j], G[:, j], P1[:, j], P2[:, j], rng, bw=0.30)
        c = r['cond']
        kk.append(r['kkde']); kkse.append(r['skde']); kh.append(r['kap'])
        inv.append(r['inv']); br.append(c['bracket']); W.append(c['W'])
        EP2.append(c['EP2']); EH.append(c['EH']); VG.append(c['VG']); tot.append(c['tot'])
    kk = np.array(kk); kkse = np.array(kkse)
    # pooled Gpp = mean over modes; SE = combination of mode SE and mode scatter
    mean = kk.mean(); se_tot = np.sqrt(kkse.mean()**2 + kk.var(ddof=1)/len(kk))
    print(f"{label:14s} N={A.shape[0]:5d} modes={Nm:2d} | "
          f"Gpp(KDE)={mean:+.3f}+-{se_tot:.3f}  [mode range {kk.min():+.3f}..{kk.max():+.3f}]  "
          f"1/Var={np.mean(inv):.3f}")
    print(f"{'':14s}   kappa_hist={np.mean(kh):+.3f}  "
          f"E[H]-Var={np.mean(br):+.3f}  E[Phi2]={np.mean(EP2):+.3f}  "
          f"W={np.mean(W):+.3f}  bracket+W={np.mean(tot):+.3f}  "
          f"E[H]={np.mean(EH):.3f} VarG={np.mean(VG):.3f}")
    return dict(label=label, N=A.shape[0], gpp=mean, se=se_tot, inv=np.mean(inv),
                kh=np.mean(kh), br=np.mean(br), EP2=np.mean(EP2), W=np.mean(W),
                tot=np.mean(tot), EH=np.mean(EH), VG=np.mean(VG))


if __name__ == '__main__':
    import os
    os.chdir('D:/Dev/math/research/unconditional/w7')
    res = []
    print("=" * 100)
    for fn, grp, lab in [('rerun_su2_L4.npz', 'su2', 'SU(2) L=4'),
                         ('rerun_su2_L6.npz', 'su2', 'SU(2) L=6'),
                         ('rerun_su3_L4.npz', 'su3', 'SU(3) L=4'),
                         ('rerun_su3_L6.npz', 'su3', 'SU(3) L=6')]:
        if os.path.exists(fn):
            res.append(summarize(fn, grp, lab))
        else:
            print(f"{lab}: MISSING {fn}")
    print("=" * 100)
