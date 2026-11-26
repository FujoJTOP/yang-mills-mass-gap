import numpy as np
import p11_volt as p

R = {}
for L in (2, 3):
    d = np.load(f"p11_L{L}_L{L}.npz")
    K, F = d["K"], d["F"]
    G = K.shape[1]
    S = [p.group_stats(K[:, g], F[:, g]) for g in range(G)]
    D = np.mean([r["d"] for r in S])
    DS = np.array([r["d"] for r in S])
    Sjk = np.sqrt(np.sum([r["s_jk"] ** 2 for r in S])) / G
    Stau = np.sqrt(np.sum([r["s_tau"] ** 2 for r in S])) / G
    Sgrp = np.sqrt((G - 1) / G * ((G * D - (G - 1) * DS - D) ** 2).sum())
    R[L] = dict(D=D, Sjk=Sjk, Stau=Stau, Sgrp=Sgrp, DS=DS, Np=6 * L ** 4, ns=K.shape[0],
                sK=np.sqrt(np.mean([r["s_jk"] ** 2 for r in S])),
                mK=np.mean([r["mK"] for r in S]), vF=np.mean([r["vF"] for r in S]),
                tK=np.mean([r["tK"] for r in S]), tF=np.mean([r["tF"] for r in S]),
                neff=np.mean([r["neff"] for r in S]), acc=float(np.mean(d["acc"])),
                domrej=float(d["domrej"]), wall=float(d["el"]), eps=float(d["eps"]))
    print(f"L={L} Np={6*L**4} NL={4*L**4} ns={K.shape[0]}/chain eps={float(d['eps'])} "
          f"acc={float(np.mean(d['acc'])):.3f} domrej={float(d['domrej']):.5f} wall={float(d['el']):.0f}s")
    print(f"   E[kappa]={R[L]['mK']:+.4f}+-{R[L]['sK']:.4f}  Var(F)={R[L]['vF']:.4f}  "
          f"tauK={R[L]['tK']:.1f} tauF={R[L]['tF']:.1f} Neff={R[L]['neff']:.0f}")
    print(f"   Delta={D:+.4f} +- {Sjk:.4f}(jk) +- {Stau:.4f}(tau) +- {Sgrp:.4f}(grp)  "
          f"Delta/Np={D/6/L**4:+.3e}  L^2*Delta={L**2*D:+.3f}+-{L**2*Sjk:.3f}")
    print(f"   per-chain Delta={np.round(DS,4)}   sign +: {int((DS>0).sum())}/{G}")

D2, D3 = R[2]["D"], R[3]["D"]
lr = np.log(1.5)
print(f"\nsign: L=2 {np.sign(D2):+.0f}  L=3 {np.sign(D3):+.0f}  same={np.sign(D2)==np.sign(D3)}")
print(f"ratio D3/D2 = {D3/D2:.3f}")
for lab, s2, s3 in (("jk", R[2]["Sjk"], R[3]["Sjk"]), ("tau", R[2]["Stau"], R[3]["Stau"]),
                    ("grp", R[2]["Sgrp"], R[3]["Sgrp"])):
    k = np.log(D3 / D2) / lr
    sk = np.sqrt((s2 / D2) ** 2 + (s3 / D3) ** 2) / lr
    kn = np.log(9 * D3 / (4 * D2)) / lr
    print(f"k(raw) from sigma_{lab}: {k:+.2f} +- {sk:.2f}   k(L^2-normalised): {kn:+.2f} +- {sk:.2f}")
print(f"extensivity: E[kappa] {R[2]['mK']:+.4f} -> {R[3]['mK']:+.4f} (ratio {R[3]['mK']/R[2]['mK']:.3f}), "
      f"Var(F) {R[2]['vF']:.4f} -> {R[3]['vF']:.4f} (ratio {R[3]['vF']/R[2]['vF']:.3f})  [L^0 would be 1.00]")
print(f"Delta/Np ratio {R[3]['D']/486/(R[2]['D']/96):.3f}  (L^-4-extensive would be 96/486={96/486:.3f})")
