#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""q09_cov_cs.py -- Cauchy-Schwarz lower bound on Cov(G,Phi') and whether it suffices.

Three panels:
  [1] re-analysis of w7/p14_main.npz  (L=2, SU(2), beta=2.5, o24 normalization)
  [2] window bookkeeping: the three candidate Var(Phi') numbers and where each comes from
  [3] exact algebra: what CS gives, what it costs, and the thresholds Var(Phi') must beat
Also prints the two structural checks ([T]): Phi''(0)=N/12 chain, and the Lambda-spread
killed-route estimate for a pathwise bound on |Phi'|.

usage: py -3.12 q09_cov_cs.py [npz]
"""
import sys
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
    return max(2.0 * min(s.sum(), 10 * s[0]), 1.0)


# ---------------------------------------------------------------- [1] p14 data
def panel1(path):
    d = np.load(path)
    rec = d["rec"]                       # (nsweeps, nchains, 4): G, Phi', Phi'', H
    ns, nc, _ = rec.shape
    P("=" * 78)
    P("[1] p14_main.npz  (nsweeps=%d, nchains=%d)  cols: 0=G  1=Phi'  2=Phi''  3=H" % (ns, nc))
    rows = []
    for c in range(nc):
        G, Pp = rec[:, c, 0], rec[:, c, 1]
        vG, vP = G.var(), Pp.var()
        cov = ((G - G.mean()) * (Pp - Pp.mean())).mean()
        rows.append((vG, vP, cov, cov / np.sqrt(vG * vP)))
    rows = np.array(rows)
    vG_a, vP_a, cv_a, cr_a = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]
    for c in range(nc):
        P("    chain %d  VarG=%9.5f  VarPhi'=%9.6f  Cov=%+9.6f  corr=%+9.5f  CS=%8.5f"
          % (c, vG_a[c], vP_a[c], cv_a[c], cr_a[c], np.sqrt(vG_a[c] * vP_a[c])))
    # pooled (means subtracted per chain -> genuine zero means by parity)
    G = rec[:, :, 0] - rec[:, :, 0].mean(0)
    Pp = rec[:, :, 1] - rec[:, :, 1].mean(0)
    vG, vP = G.var(), Pp.var()
    cov = (G * Pp).mean()
    cs = np.sqrt(vG * vP)
    P("    pooled   VarG=%9.5f  VarPhi'=%9.6f  Cov=%+9.6f  corr=%+9.5f  CS=%8.5f"
      % (vG, vP, cov, cov / cs, cs))
    P("    ratio VarPhi'/VarG = %.5f  (%.3f%%)" % (vP / vG, 100 * vP / vG))
    # tau-based sigma on Cov, per chain
    z = G * Pp
    # 4 independent chains, nc*n sweeps total: sigma^2 = mean_c(var_c tau_c / n) / nc
    v_tau = np.array([z[:, c].var() * tau_ips(z[:, c]) / ns for c in range(nc)])
    sg = np.sqrt(v_tau.mean() / nc)
    P("    per-chain sigma_tau(Cov) = %s" % np.array2string(np.sqrt(v_tau), precision=5))
    P("    sigma_tau(Cov, 4 chains pooled) = %.6f  =>  Cov in [%+.5f, %+.5f] (3 sigma)"
      % (sg, cov - 3 * sg, cov + 3 * sg))
    P("    CS slack vs 3-sigma measured: %.2fx  (1/|corr| = %.1f)"
      % (3 * sg / cs, 1 / abs(cov / cs) if cov else float("inf")))
    return dict(vG=vG, vP=vP, cov=cov, cs=cs, sig=sg, corr=cov / cs)


# ------------------------------------------------- [2] the three Var(Phi') values
# sources are cited in Q09-COVLOW.md; all L=2 / SU(2) / beta=2.5 unless noted
WINDOWS = {
    "A_prompt": dict(vG=5.4477, vP=3.6e-3, EH=6.3315,
                     src="VarG: P02-DYNDYN.md:48 ; E[H]: P02:46 ; VarPhi'=3.6e-3: "
                         "O09-PHICURV.md:80 (single plaquette, a=0.75)  <-- WINDOW MISMATCH"),
    "B_lattice_P02": dict(vG=5.4477, vP=0.162409, EH=6.3315, EHp=None,
                          src="VarPhi'=Var(G+Phi')-Var(G)-2Cov algebra, P02-DYNDYN.md:54"),
    "C_lattice_o27": dict(vG=5.7307, vP=0.14140, EH=6.319132, EHp=None,
                          src="o27: VarG=5.7307, Var(G+Phi')=5.8291, Cov=-0.0215 "
                              "(O27-ADJUDICATE.md:76-78)"),
}


def panel2():
    P("=" * 78)
    P("[2] window bookkeeping.  Var(G+Phi') = VarG + VarPhi' + 2Cov  (identity, [T])")
    for k, w in WINDOWS.items():
        vG, vP = w["vG"], w["vP"]
        sG, sP = np.sqrt(vG), np.sqrt(vP)
        cs = sG * sP
        P("  %-14s VarG=%8.4f  VarPhi'=%.6f  sigmaG=%.4f sigmaPhi=%.4f  CS=%.5f"
          % (k, vG, vP, sG, sP, cs))
        P("                 ratio VarPhi'/VarG = %.5f   1/(ratio) = %.1f" % (vP / vG, vG / vP))
        P("                 src: %s" % w["src"])
    # reconciliations
    B, C = WINDOWS["B_lattice_P02"], WINDOWS["C_lattice_o27"]
    P("  cross-check B: Var(G+Phi') = %.6f + %.6f + 2*(%.6f) = %.6f  (P02:48 gives 5.580327)"
      % (B["vG"], B["vP"], -0.014886, B["vG"] + B["vP"] - 2 * 0.014886))
    P("  cross-check C: Var(G+Phi') = %.6f + %.6f + 2*(%.6f) = %.6f  (o27 gives 5.8291)"
      % (C["vG"], C["vP"], -0.0215, C["vG"] + C["vP"] - 2 * 0.0215))
    P("  A vs B: VarPhi' differs by %.1fx  ->  3.6e-3 is NOT a lattice number"
      % (B["vP"] / WINDOWS["A_prompt"]["vP"]))


# ------------------------------------------------------------- [3] exact algebra
def panel3():
    P("=" * 78)
    P("[3] algebra.  Cov>=-VarG/2 via CS needs VarPhi'<=VarG/4 ; H1 needs the other side.")
    P("    Var(F)=VarG+VarPhi'+2Cov,  F=G+Phi',  H1: W_full=E[H]+VarG-Var(F)>=0")
    for k, w in WINDOWS.items():
        vG, vP, EH = w["vG"], w["vP"], w["EH"]
        cs = np.sqrt(vG * vP)
        d = EH - vG
        # (i) the requested inequality  Cov >= -VarG/2
        thr_A = vG / 4.0
        P("  --- %s ---" % k)
        P("    (i) Cov >= -VarG/2 = %+.4f : CS gives Cov >= %+.5f  [T]  slack %.2fx"
          % (-vG / 2, -cs, (vG / 2) / cs))
        P("        equals VarPhi' <= VarG/4 = %.5f ; have %.6f  => %s (%.1fx margin)"
          % (thr_A, vP, "PASS" if vP <= thr_A else "FAIL", thr_A / vP))
        # (ii) H1 needs an UPPER bound on Cov
        cov_max_needed = (EH - vG - vP) / 2.0
        P("    (ii) H1 one-sided: need Cov <= (E[H]-VarG-VarPhi')/2 = %.5f ; CS gives Cov <= %.5f  => %s"
          % (cov_max_needed, cs, "SUFFICES" if cs <= cov_max_needed else "TOO LOOSE by %.2fx" % (cs / cov_max_needed)))
        # threshold on VarPhi' for H1
        sG = np.sqrt(vG)
        x = -sG + np.sqrt(sG * sG + d)          # sigma_Phi' <= -sigmaG+sqrt(sigmaG^2+Delta)
        P("        threshold: VarPhi' <= %.5f gives (sigmaG+sigmaPhi)^2 <= E[H]-?  ->  VarPhi' <= %.5f"
          % (x * x, x * x))
        P("        Delta=E[H]-VarG=%.5f ; (sigmaG+sigmaPhi)^2=%.5f <= E[H]=%.4f ? %s"
          % (d, (sG + np.sqrt(vP)) ** 2, EH, "YES" if (sG + np.sqrt(vP)) ** 2 <= EH else "NO"))
        # what Cov must be for H1 to be at risk: Var(F)>E[H]
        P("        Var(F) at Cov=+CS: %.5f  vs E[H]=%.4f  (excess %+.5f)"
          % (vG + vP + 2 * cs, EH, vG + vP + 2 * cs - EH))
    P("  NOTE: (i) is satisfied by ANY window, but it only forces Var(F) >= VarPhi'+2Cov+0 >= 0;")
    P("        it carries no information about H1. The binding side for H1 is (ii).")


# --------------------------------------------------------- [T] structural checks
def panel4():
    P("=" * 78)
    P("[4] structural checks ([T] unless noted)")
    # Phi''(0) = 1/6 per link; N=64 links at L=2 -> N/12 = 5.333 needs sum c^2 >= N/2
    N = 4 * 2 ** 4
    phi2_0 = 1.0 / 6.0
    P("    N(links,L=2)=%d ; Phi''(0)=%.10f per link (Taylor of -cot(r/2)+2/r)"
      % (N, phi2_0))
    P("    Phi''>= (1/6)*sum_l c_l^2 ; with sum c_l^2 >= N/2 => Phi'' >= N/12 = %.4f  (E[H]=6.3315 OK)"
      % (N / 12.0))
    rr = np.linspace(0.30, 5.657, 400001)
    lam = (-1.0 / np.tan(rr / 2) + 2.0 / rr) / rr      # Lambda(r)=phi'(r)/r
    P("    Lambda(r)=phi'(r)/r on r in [0.30,5.657]: min=%.5f max=%.5f spread=%.5f"
      % (lam.min(), lam.max(), lam.max() - lam.min()))
    P("    pathwise route |Phi'|<=sum r_l|Lambda_l-mu| <= spread*sum r_l : spread=%.4f, sum r ~ %.1f"
      % (lam.max() - lam.min(), N * 3.06))
    bound = (lam.max() - lam.min()) * N * 3.06
    P("    => Var(Phi') <= %.1f^2 = %.0f  (need <= %.4f) : PATHWISE ROUTE DEAD by %.0fx"
      % (bound, bound ** 2, WINDOWS["B_lattice_P02"]["vG"] / 4, bound ** 2 / WINDOWS["B_lattice_P02"]["vG"] * 4))
    # exactly-aligned limit
    P("    exactly aligned (phi'=mu*r for all l) => Phi'=mu*sum r_l c_l=0 exactly => Var(Phi')=0 [T]")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else r"p14_main.npz"
    panel1(path)
    panel2()
    panel3()
    panel4()
    P("=" * 78)
    P("done.")
