"""Q06 - does G''(beta)=S-D satisfy an integrable Riccati-type inequality in beta?

S(beta) := E[H+Phi'']  (= E[kappa])      D(beta) := Var(G+Phi') (= Var(F))
                            G''(beta) = S-D      (H1 <=> S>D)

Measure: dmu_beta = exp(-beta*h + Phi),  h = sum_P (1-ReTr U_P /2)  (Wilson),  Phi = -log J
(beta-free).  Flow identity for any observable O:

    d/dbeta E[O] = E[d_beta O] - Cov(O,h)                                     [T]

with the pointwise-affine structure  V = beta*h + Phi  =>  F = beta*g + F_Phi,
kappa = beta*h2 + kappa_Phi,  g := d_t h,  h2 := d_t^2 h  (both beta-free random variables):

    S' = E[h2] - Cov(kappa,h)                       (1)
    D' = 2 Cov(F,g) - Cov(F^2,h) + 2 E[F] Cov(F,h)  (2)

Both are 3rd-cumulant objects -> not closed in (S,D).  Samples from p19 (%TEMP%/p19_b*.npz,
L=2, mod A, 4 chains x 5000 sweeps) give g = F_W/beta and h2 = kappa_W/beta directly, so the
response terms E[h2] and 2Cov(F,g) are measurable; the Cov(.,h) terms are not (h is not stored)
and are inferred here as residuals against finite-difference S',D'.

usage: py -3.12 q06_riccati.py            (reads %TEMP%/p19_b{1,2.5,4,6}.npz, prints table)
"""
import os

import numpy as np

TEMP = os.environ.get('TEMP', '.')
BETAS = [1.0, 2.5, 4.0, 6.0]
TAGS = ['b1', 'b2.5', 'b4', 'b6']
NB = 25
# p19 published numbers (research/v3/species/P19-BETA-AT-L2.md sec.0) for the cross-check
PUB = dict(G=[0.2835, 0.4116, 0.5108, 0.5913], S=[0.5723, 1.9348, 3.3529, 5.3242],
           D=[0.2888, 1.5232, 2.8421, 4.7329],
           sG=[0.0049, 0.0421, 0.0939, 0.1732], sD=[0.0058, 0.0426, 0.0948, 0.1716])


def jk(c, fn, nb=NB):
    """delete-one-block jackknife on one chain's (ns,6) array -> (value, sigma_jk)."""
    idx = np.array_split(np.arange(len(c)), nb)
    ds = np.array([fn(np.delete(c, b, axis=0)) for b in idx])
    return ds.mean(), np.sqrt((nb - 1) / nb * ((ds - ds.mean()) ** 2).sum())


def comb(rows):
    """4-group combination: mean, sigma_jk = sqrt(sum s_jk^2)/G, group error."""
    v = np.array([r[0] for r in rows])
    s = np.array([r[1] for r in rows])
    G = len(v)
    m = v.mean()
    dg = G * m - (G - 1) * v
    sgrp = np.sqrt((G - 1) / G * ((dg - dg.mean()) ** 2).sum())
    return m, np.sqrt((s ** 2).sum()) / G, sgrp


def d1(xs, ys, i, st):
    """2nd-order Lagrange derivative at node xs[i] from stencil st (3 indices) -> (value, w).

    Fornberg first-derivative weights; exact for quadratics, including the endpoint case i in st.
    """
    n = 3
    w = np.zeros(n)
    for j in range(n):
        if st[j] == i:
            continue
        t = 1.0 / (xs[st[j]] - xs[i])
        for q in range(n):
            if q != j and st[q] != i:
                t *= (xs[i] - xs[st[q]]) / (xs[st[j]] - xs[st[q]])
        w[j] = t
    for j in range(n):
        if st[j] == i:
            w[j] = sum(1.0 / (xs[i] - xs[st[q]]) for q in range(n) if st[q] != i)
    return float(w @ np.array([ys[j] for j in st])), w


def _selfcheck():
    """derivative stencils must be exact on quadratics (they are 3-point Lagrange)."""
    xs = np.array([1.0, 2.5, 4.0, 6.0])
    for i in range(4):
        st = (0, 1, 2) if i < 2 else (1, 2, 3)
        v, _ = d1(xs, xs ** 2, i, st)
        assert abs(v - 2 * xs[i]) < 1e-9, (i, v)


FNS = {
    'S': lambda z, b: z[:, 0].mean(),
    'D': lambda z, b: z[:, 1].var(),
    'G': lambda z, b: z[:, 0].mean() - z[:, 1].var(),
    'Rsp': lambda z, b: 2.0 * float(((z[:, 1] - z[:, 1].mean())
                                     * (z[:, 3] / b - (z[:, 3] / b).mean())).mean()),
    'h2': lambda z, b: z[:, 2].mean() / b,
    'kPhi': lambda z, b: z[:, 4].mean(),
    'A': lambda z, b: (z[:, 3] / b).var(),
    'B': lambda z, b: float(((z[:, 3] / b - (z[:, 3] / b).mean())
                             * (z[:, 5] - z[:, 5].mean())).mean()),
    'C': lambda z, b: z[:, 5].var(),
}


def main():
    _selfcheck()
    print(__doc__)
    per = {k: [] for k in FNS}                 # per[k] = [ [4 chains at beta_0], [.. beta_1], ...]
    zmax = 0.0
    for beta, tag in zip(BETAS, TAGS):
        C = np.load(os.path.join(TEMP, f'p19_{tag}.npz'))['cols']      # (ns,G,6)
        G = C.shape[1]
        rows = {k: [] for k in FNS}
        for g in range(G):
            c = C[:, g, :]
            F, kap, kW, FW, FP = c[:, 1], c[:, 0], c[:, 2], c[:, 3], c[:, 5]
            zmax = max(zmax, float(np.abs(F - (FW + FP)).max()),
                       float(np.abs(kap - (kW + c[:, 4])).max()))
            for k, f in FNS.items():
                rows[k].append(jk(c, lambda z, f=f, b=beta: f(z, b)))
        for k in FNS:
            per[k].append(rows[k])
        print(f"[Q06] beta={beta:g}  chains={G}  ns={C.shape[0]}  "
              f"affine check max|F-(F_W+F_Phi)|,|kappa-(kappa_W+kappa_Phi)| = {zmax:.1e}")
    M = {k: [comb(x) for x in v] for k, v in per.items()}
    b = np.array(BETAS)
    S = np.array([m[0] for m in M['S']]); sS = np.array([m[1] for m in M['S']])
    D = np.array([m[0] for m in M['D']]); sD = np.array([m[1] for m in M['D']])
    Gp = np.array([m[0] for m in M['G']])
    Rsp = np.array([m[0] for m in M['Rsp']]); sRsp = np.array([m[1] for m in M['Rsp']])
    h2 = np.array([m[0] for m in M['h2']])
    kPhi = np.array([m[0] for m in M['kPhi']])
    A = np.array([m[0] for m in M['A']]); B = np.array([m[0] for m in M['B']])
    Cc = np.array([m[0] for m in M['C']])
    R = D / S

    print("\n--- reproduction of p19 (P19-BETA-AT-L2.md sec.0) ---")
    print(f"{'beta':>5} {'S_jk':>8} {'S_pub':>8} {'D_jk':>8} {'D_pub':>8} {'G_jk':>8} {'G_pub':>8} {'maxrel':>9}")
    for i in range(4):
        rel = max(abs(S[i] / PUB['S'][i] - 1), abs(D[i] / PUB['D'][i] - 1), abs(Gp[i] / PUB['G'][i] - 1))
        print(f"{b[i]:5.1f} {S[i]:8.4f} {PUB['S'][i]:8.4f} {D[i]:8.4f} {PUB['D'][i]:8.4f} "
              f"{Gp[i]:+8.4f} {PUB['G'][i]:+8.4f} {rel:9.2e}")

    print("\n--- affine structure: which beta-power carries D and S ---")
    print(f"{'beta':>5} {'b^2Var(g)/D':>12} {'2bCov(g,Phi)/D':>16} {'Var(Phi)/D':>11} "
          f"{'Var_b(g)':>9} {'Var_b(F_Phi)':>12} {'E[h2]':>7} {'E[kap_Phi]':>10} {'S-afit':>8}")
    for i in range(4):
        print(f"{b[i]:5.1f} {b[i]**2 * A[i] / D[i]:12.3f} {2 * b[i] * B[i] / D[i]:16.3f} "
              f"{Cc[i] / D[i]:11.3f} {A[i]:9.4f} {Cc[i]:12.4f} {h2[i]:7.4f} {kPhi[i]:10.4f} "
              f"{S[i] - (b[i] * h2[i] + kPhi[i]):+8.4f}")

    st = [(0, 1, 2), (0, 1, 2), (1, 2, 3), (1, 2, 3)]
    sSall = np.sqrt(sS ** 2 + sD ** 2)                 # S = G + D, conservative (independent)
    Sp, sSp, Dp, sDp = [], [], [], []
    for i in range(4):
        v, w = d1(b, S, i, st[i]); Sp.append(v)
        sSp.append(float(np.sqrt(((w * sSall[list(st[i])]) ** 2).sum())))
        v, w = d1(b, D, i, st[i]); Dp.append(v)
        sDp.append(float(np.sqrt(((w * sD[list(st[i])]) ** 2).sum())))
    Sp, Dp, sSp, sDp = map(np.array, (Sp, Dp, sSp, sDp))
    gap = Dp - Sp
    # gap must equal -G' identically (S-D = G); propagate p19's own sigma(G'') through the stencil
    Gpd, sGpd = [], []
    for i in range(4):
        v, w = d1(b, Gp, i, st[i]); Gpd.append(v)
        sGpd.append(float(np.sqrt(((w * np.array(PUB['sG'])[list(st[i])]) ** 2).sum())))
    Gpd, sGpd = np.array(Gpd), np.array(sGpd)
    print(f"\n  identity check  gap=-G': max|gap+G'| = {np.abs(gap + Gpd).max():.2e}"
          f" ; G'={np.round(Gpd, 4)} +- {np.round(sGpd, 4)} (p19 sigma)")
    q = b * Dp / D
    mm = b * Sp / S
    print("\n--- pointwise derivatives (non-uniform 3-pt stencils; sn = sqrt(sS^2+sD^2)) ---")
    print(f"{'beta':>5} {'S':>8} {'Sprime':>8} {'+-jk':>6} {'Dprime':>8} {'+-jk':>6} "
          f"{'Dprime-Sprime':>13} {'q=bDp/D':>8} {'m=bSp/S':>8} {'m-q*R':>7} {'Dp<=Sp':>7}")
    for i in range(4):
        print(f"{b[i]:5.1f} {S[i]:8.4f} {Sp[i]:8.4f} {sSp[i]:6.4f} {Dp[i]:8.4f} {sDp[i]:6.4f} "
              f"{gap[i]:+13.4f} {q[i]:8.3f} {mm[i]:8.3f} {mm[i] - q[i] * R[i]:+7.4f} "
              f"{('yes' if gap[i] <= 0 else 'NO'):>7}  gap/sig(G')={gap[i] / sGpd[i]:+.2f}")
    print("  interval averages:")
    for i in range(3):
        db = b[i + 1] - b[i]
        dd, ds = (D[i + 1] - D[i]) / db, (S[i + 1] - S[i]) / db
        print(f"    [{b[i]:4.1f},{b[i+1]:4.1f}] mid={0.5 * (b[i] + b[i + 1]):5.2f} dD={dd:+.4f} "
              f"dS={ds:+.4f} diff={dd - ds:+.4f}  {'D/S up' if dd * S[i] > D[i] * ds else 'D/S down'}")

    print("\n--- 3rd cumulants inferred from (1),(2) ---")
    print(f"{'beta':>5} {'E[h2]':>8} {'2Cov(F,g)':>10} {'+-jk':>6} {'Cov(kap,h)':>11} {'Cov(F^2,h)':>11}")
    for i in range(4):
        print(f"{b[i]:5.1f} {h2[i]:8.4f} {Rsp[i]:10.4f} {sRsp[i]:6.4f} {h2[i] - Sp[i]:+11.4f} "
              f"{Rsp[i] - Dp[i]:+11.4f}")

    print("\n--- beta-dependence of D/S: competing closures ---")
    c0, p0 = np.polyfit(np.log(b), np.log(1 - R), 1)
    c, p = float(np.exp(c0)), float(-p0)
    print(f"  [P1] 1-D/S=c*beta^-p   c={c:.4f} p={p:+.4f} -> D/S->1 (slow); "
          f"R(10)={1 - c * 10**-p:.3f} R(20)={1 - c * 20**-p:.3f} R(100)={1 - c * 100**-p:.3f}")
    sR = float(np.polyfit(np.log(b), np.log(R), 1)[0])
    print(f"  [P2] D/S~beta^s        s={sR:+.4f} -> unbounded, crossing exists")
    Dq = np.polyfit(b, D, 2); Sl = np.polyfit(b, S, 1)
    zq = np.roots([Dq[0], Dq[1] - Sl[0], Dq[2] - Sl[1]])
    zq = sorted(z.real for z in zq if abs(z.imag) < 1e-9 and z.real > 0)
    print(f"  [P3] structural D=Ab^2+Bb+C, S=lb+e: D={np.round(Dq, 5)} res={abs(D - np.polyval(Dq, b)).max():.1e}"
          f" | S={np.round(Sl, 5)} res={abs(S - np.polyval(Sl, b)).max():.1e}")
    print(f"       D=S at beta={np.round(zq, 2)};  D'=S' at beta={(Sl[0] - Dq[1]) / (2 * Dq[0]):.2f}"
          f"   D(6)pred={np.polyval(Dq, 6):.4f} (data {D[3]:.4f})")
    e = float((q - mm)[-1])
    print(f"  [P4] R=Rs*(b/6)^(q-m), q-m(b=6)={e:+.3f} -> D=S at beta="
          f"{(6 * (1 / R[-1])**(1 / e) if e > 0 else float('nan')):.2f}")
    zG = sorted(z.real for z in np.roots(np.polyfit(b, Gp, 2)) if abs(z.imag) < 1e-9 and z.real > 1)
    print(f"  [P5] G''=S-D quadratic in beta: zeros>1 at {np.round(zG, 2)} (p19_fit quotes 17.31)")
    gmid = 0.5 * (b[1:] + b[:-1])
    gint = np.array([(D[i + 1] - D[i]) / (b[i + 1] - b[i]) - (S[i + 1] - S[i]) / (b[i + 1] - b[i])
                     for i in range(3)])
    gl = np.polyfit(gmid, gint, 1)
    print(f"  [P6] interval gap {np.round(gint, 4)} vs midpoint -> zero at beta={-gl[1]/gl[0]:.2f}")

    print("\n--- verdict ---")
    print(f"  (a) pointwise D'<=S' on [1,6]: {'yes' if (gap <= 0).all() else 'NO at beta=' + str(b[gap > 0])}"
          f" ; interval-averaged: {'yes' if (gint < 0).all() else 'NO'}"
          f" ; provable from (S,D): NO (3rd cumulants)")
    print(f"  (b) D/S={np.round(R, 3)} (monotone up, <1 on [1,6]); 1-D/S ~ beta^-{p:.2f} (p<1 => no")
    print(f"      crossing in window); structural closure puts D=S at beta={np.round(zq[1] if len(zq) > 1 else zq[0], 1)}")
    print("  (c) Cov(kap,h), Cov(F^2,h) required -> not closed in (S,D); flow identity is single-beta,")
    print("      so one run storing h closes them exactly (no hierarchy).")
    print("  acid test: free field S=D and every term above vanishes -> comparison is boundary-tight.")


if __name__ == '__main__':
    main()
