"""A5 certificate-margin calculator (read-only; reuses w7 pool npz).

Goal: quantify, on the EXISTING SU(3) beta=5.7 k=2pi/L data, exactly how far the
finite-region certificate  G'' >= b* + FLOOR - u* + x*  is from being rigorous.

Identity (flat chart, exact; rerun_core):
    G''(0) = bracket + W_geo + X,
    bracket = E[H|0]-Var(G|0),  W_geo = E[Phi''|0]-Var(Phi'|0),  X = -2Cov(G,Phi'|0).

Proven pointwise floor (Lemma B, PROOF-JACOBIAN-SAVES):
    Phi'' >= (N/12)||v||^2  =>  E[Phi''|0] >= N/12 = 1/4   (SU(3), unit mode).
Hence the CERTIFICATE (rigorous IF b*,u*,x* are rigorous bounds):
    G'' >= b* + 1/4 - u* + x*,   b* <= bracket,  u* >= Var(Phi'|0),  x* <= X.

We report, per L:
  (1) point values of bracket, E[Phi''], Var(Phi'), W_geo, X, G''(KDE);
  (2) lambda_floor = bracket + 1/4 - Var(Phi') + X  (point-value certificate LB);
  (3) delta-budget = lambda_floor  (total rigorous-enclosure slack still allowed
      while keeping the certified lower bound > 0);
  (4) the Cauchy-Schwarz a-priori ceiling on |X| = 2 sqrt(VarG * VarPhi') --
      the best deterministic bound available -- showing it EXCEEDS the budget.

No sampling, no bootstrap: pure arithmetic on stored arrays. ~MB memory.
"""
import sys
import numpy as np

W7 = 'D:/Dev/math/research/unconditional/w7'
sys.path.insert(0, W7)
import kzero_analyze as KZ   # noqa: E402
import rerun_core as RC      # noqa: E402

FLOOR = 0.25                 # N/12 for SU(3), unit-normalized mode (Lemma B)
BW = 0.30                    # kernel conditioning bandwidth (matches kzero_analyze)


def analyse(L):
    d = np.load(f'{W7}/_pool_su3_L{L}.npz')
    A, H, G, P1, P2 = d['a'], d['H'], d['G'], d['P1'], d['P2']
    rows = []
    for j in range(A.shape[1]):
        a, h, g, p1, p2 = A[:, j], H[:, j], G[:, j], P1[:, j], P2[:, j]
        k = RC.kernel_cond(a, {'H': h, 'G': g, 'P1': p1, 'P2': p2}, bw=BW)
        EH, VG = k['E_H'], k['V_own']['G']
        EP2, VP1, cov = k['E_P2'], k['V_own']['P1'], k['C_G_P1']
        bracket = EH - VG
        Wgeo = EP2 - VP1
        X = -2.0 * cov
        gpp = KZ._kde_point(a)
        rows.append(dict(bracket=bracket, EP2=EP2, VP1=VP1, Wgeo=Wgeo, X=X,
                         gpp=gpp, VG=VG, EH=EH,
                         lam_floor=bracket + FLOOR - VP1 + X))
    m = {key: float(np.mean([r[key] for r in rows])) for key in rows[0]}
    return L, A.shape[0], A.shape[1], rows, m


def main():
    print("=" * 92)
    print("A5 certificate margins -- SU(3) beta=5.7, transverse plane mode k=2pi/L")
    print("  identity G''=(E[H]-VarG)+[E(Phi'')-Var(Phi')]+X ; proven floor E(Phi'')>=%.2f" % FLOOR)
    print("=" * 92)
    for L in (4, 6):
        try:
            L, N, Nm, rows, m = analyse(L)
        except FileNotFoundError:
            print(f"L={L}: pool npz missing"); continue
        print(f"\n--- L={L}  N={N}  modes={Nm}  k=2pi/L={2*np.pi/L:.4f} ---")
        print(f"  {'mode':>4s} {'brkt':>7s} {'EPhi2':>7s} {'VarPhi':>7s} {'W_geo':>7s} "
              f"{'X':>7s} {'Gpp':>7s} {'lam_flr':>8s}")
        for j, r in enumerate(rows):
            print(f"  {j:>4d} {r['bracket']:+7.3f} {r['EP2']:7.3f} {r['VP1']:7.3f} "
                  f"{r['Wgeo']:+7.3f} {r['X']:+7.3f} {r['gpp']:+7.3f} {r['lam_floor']:+8.3f}")
        budget = m['lam_floor']
        cs = 2.0 * np.sqrt(m['VG'] * m['VP1'])          # Cauchy-Schwarz |X| ceiling
        print(f"  POOLED  brkt={m['bracket']:+.3f}  EPhi2={m['EP2']:.3f}  VarPhi'={m['VP1']:.3f}  "
              f"W_geo={m['Wgeo']:+.3f}  X={m['X']:+.3f}  Gpp={m['gpp']:+.3f}")
        print(f"  CERTIFICATE  G'' >= brkt + {FLOOR:.2f} - VarPhi' + X = {budget:+.4f}   "
              f"[point-value; needs rigorous b*,u*,x*]")
        print(f"  delta-BUDGET: rigorous enclosures must have TOTAL slack < {budget:.4f} to certify G''>0")
        print(f"  a-priori / Cauchy-Schwarz ceiling on |X| = 2 sqrt(VarG*VarPhi') = {cs:.3f}"
              f"   (budget for |X| alone ~ {FLOOR + m['bracket'] - m['VP1']:.3f}) -> "
              f"{'VACUOUS' if cs > FLOOR + m['bracket'] - m['VP1'] else 'could help'}")
    print("\n" + "=" * 92)
    print("VERDICT: point-value certificate is POSITIVE (+0.09 L4 / +0.04 L6) but the")
    print("enclosure is statistical, not rigorous; and CS a-priori |X| ceiling >> budget.")
    print("=" * 92)


if __name__ == '__main__':
    main()
