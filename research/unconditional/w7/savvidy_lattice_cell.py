"""
savvidy_lattice_cell.py -- SAVVIDY-LATTICE-CELL witness  (2026-11-24)

Empirical test of the "softened chi-floor" prediction in the negative-mode-active
region, using EXISTING data only (no MC):

  * research/unconditional/w7/k0_su2_b2.4_L8_n2.npz  -- per-config, per-mode
    (a=|<v,x>|, H=d^2 S/ds^2, G=dS/ds, P1=d(-logJ)/ds, P2=d^2(-logJ)/ds^2)
    at SU(2) beta=2.4, L=8, k=2*pi/L and 4*pi/L, N=400 configs.
  * research/lattice-cpb/runs/L6_b2.600_s51.npz  -- plaquette "defect" time series
    (negative-plaquette fraction) at beta=2.6, L=6.

Exact identity used (k0_continuum_run.py docstring):
    G''(0) = (E[H] - Var(G)) + W ,  W = E[P2] - Var(P1) - 2*Cov(G,P1)
    => S = E[H],  D = Var(G) - W,  G'' = S - D,  D/S = defect fraction.

Question: is D/S bounded (<1) and positive, consistent with a softened floor?
Label: [V-num] (measured, not proved).  No QED.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))  # D:/Dev/math
K0 = os.path.join(ROOT, "research", "unconditional", "w7", "k0_su2_b2.4_L8_n2.npz")
L6 = os.path.join(ROOT, "research", "lattice-cpb", "runs", "L6_b2.600_s51.npz")


def mode_stats(H, G, P1, P2, j):
    h, g, p1, p2 = H[:, j], G[:, j], P1[:, j], P2[:, j]
    EH = float(h.mean()); Vg = float(g.var())
    Ep2 = float(p2.mean()); Vp1 = float(p1.var())
    Cov = float((g * p1).mean() - g.mean() * p1.mean())
    W = Ep2 - Vp1 - 2.0 * Cov
    return EH, Vg, W, EH - Vg + W, Vg - W         # S, VarG, W, Gpp, D


def block(mask, H, G, P1, P2):
    h, g = H[:, mask].ravel(), G[:, mask].ravel()
    p1, p2 = P1[:, mask].ravel(), P2[:, mask].ravel()
    EH = float(h.mean()); Vg = float(g.var())
    Ep2 = float(p2.mean()); Vp1 = float(p1.var())
    Cov = float((g * p1).mean() - g.mean() * p1.mean())
    W = Ep2 - Vp1 - 2.0 * Cov
    S, D = EH, Vg - W
    return S, Vg, W, S - D, D, D / S


def main():
    d = np.load(K0, allow_pickle=True)
    H, G, P1, P2 = d["H"], d["G"], d["P1"], d["P2"]
    kv = np.asarray(d["kvals"], float)
    names = [str(x) for x in d["names"]]
    print("=" * 74)
    print("SAVVIDY-LATTICE-CELL  (SU(2) beta=%.2f L=%d N=%d)" %
          (float(d["beta"]), int(d["L"]), H.shape[0]))
    print("=" * 74)
    print("per mode:  S=E[H]  VarG  W  G''=S-D  D  D/S")
    for j in range(H.shape[1]):
        S, Vg, W, Gpp, D = mode_stats(H, G, P1, P2, j)
        print("  %-15s k=%.3f  S=%.3f VarG=%.3f W=%+.3f G''=%+.3f D=%.3f D/S=%.3f"
              % (names[j], kv[j], S, Vg, W, Gpp, D, D / S))
    print()
    print("aggregate by k:")
    for lab, mask in [("k=2pi/L (n=1)", kv < 1.0),
                      ("k=4pi/L (n=2)", kv >= 1.0),
                      ("ALL", np.ones_like(kv, bool))]:
        S, Vg, W, Gpp, D, ds = block(mask, H, G, P1, P2)
        print("  %-14s S=%.4f VarG=%.4f W=%+.4f  G''=%+.4f  D=%.4f  D/S=%.4f"
              % (lab, S, Vg, W, Gpp, D, ds))
    Gpp_all = np.array([mode_stats(H, G, P1, P2, j)[3] for j in range(H.shape[1])])
    DS_all = np.array([mode_stats(H, G, P1, P2, j)[4] / mode_stats(H, G, P1, P2, j)[0]
                       for j in range(H.shape[1])])
    print()
    print("  per-mode spread:  G'' in [%+.3f, %+.3f]   D/S in [%.3f, %.3f]"
          % (Gpp_all.min(), Gpp_all.max(), DS_all.min(), DS_all.max()))
    print("  n(mode) with G''<0 (D/S>1): %d / %d"
          % (int((Gpp_all < 0).sum()), Gpp_all.size))

    if os.path.exists(L6):
        e = np.load(L6, allow_pickle=True)
        de, mn = e["defect"], e["plaq_min"]
        print()
        print("=" * 74)
        print("negative-mode ubiquity  (SU(2) beta=%.2f L=%d)"
              % (float(e["beta"]), int(e["L"])))
        print("  frac(configs with a negative plaquette) = %.3f" % float((mn < 0).mean()))
        print("  plaq_min mean = %.3f   defect(thr) mean = %s"
              % (float(mn.mean()), np.round(de.mean(0), 3)))
    print()
    print("[V-num] D/S is BOUNDED (<=~1.1, no divergence) but O(1) (~0.94-1.00),")
    print("        i.e. G''=S-D is a near-cancellation (~0), not robustly positive.")
    print("        Softest measured mode is k=2pi/L; k->0 (n=0) is NOT in this data.")
    print("        No QED. Clay 0%.")


if __name__ == "__main__":
    main()
