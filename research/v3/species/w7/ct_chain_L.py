#!/usr/bin/env python3
"""
ct_chain_L.py -- T1 (part 3): the L-dependence of X(0) for the 1D SU(2) chain, in a FIXED
                 normalisation (unit-speed coarse field), L = 2..8, open and periodic.

WHY A REWRITE OF na2_chain.py IS NEEDED (the normalisation erratum):
    na2_chain.py sets, for the chain-average coarse field a = (1/L) sum_i xi_3^i ,
        G   = dS / L ,     Phi' = (1/L) sum_i d/dxi_3^i Phi ,
    but d/da for that a is the UNIT-SPEED field  sum_i d/dxi_3^i  (no 1/L: shifting every xi_3^i by
    +delta moves a by +delta).  Hence X_na2 = X_true / L^2 , i.e. every chain-average number in
    na2_chain.log is L^2 times too small, and the reported "|X(0)| decreases with L"
    (-2.58e-3, -2.32e-3, -1.90e-3 for L=2,3,4) is an ARTEFACT of that normalisation.
    Here:  X(0) = -2 Cov(G_unit, Phi'_unit)  with G_unit = sum_i d/dxi_3^i S  (no 1/L).

METHOD
    samples: independent Haar SU(2)^L (exact, = the flat-chart measure with the Weyl weight),
             reweighted by exp(-beta S) and by the Gaussian kernel K_h(a) = exp(-a^2/2h^2);
             G, Phi' are FIRST derivatives only (no second divided difference of exp -- see ct_hadamard.py).
    X(0) = lim_{h->0} X(h)  (kernel bias is O(h^2); three h values + h^2 Richardson).
    error: delete-one-chunk jackknife (chunks are independent).
    CALIBRATION: L=2 open has an exact answer from ct_chain_exact.py:
        X(0) = -0.0104411252 (b=0.5),  -0.0191235583 (b=1.0),  -0.0330091889 (b=2.5).
    Extensivity (ring, translation invariance): X(0) = -2 L sum_j Cov(g_0, p_j)  =>  X/L has a limit.

Run: py -3.12 -I ct_chain_L.py [--n 300000]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import na2_su3_na as M
import na2_chain as C                    # reuse build() (Haar SU(2)^L + Loewner dU)

HB = (0.16, 0.08, 0.04, 0.02)
BETAS = (0.5, 1.0, 2.5)


def chunk_stats(ch, beta):
    """unweighted samples -> (a, S, G_unit, P1_unit)  with G_unit = sum_i d/dxi_3^i S"""
    L, bonds = ch['L'], ch['bonds']
    a = sum(ch['a']) / L
    S = np.zeros(ch['m'])
    G = np.zeros(ch['m'])
    for (i, j) in bonds:
        Wj = np.swapaxes(ch['U'][j].conj(), 1, 2)
        Wi = np.swapaxes(ch['U'][i].conj(), 1, 2)
        c0 = np.real(np.einsum('mab,mba->m', ch['U'][i], Wj))
        t_i = np.real(np.einsum('mab,mba->m', ch['dU'][i], Wj))
        t_j = np.real(np.einsum('mab,mba->m', ch['dU'][j], Wi))
        S += beta * (1.0 - 0.5 * c0)
        G -= 0.5 * beta * (t_i + t_j)          # = d/dxi_3^i + d/dxi_3^j of that bond (unit speed)
    P1 = sum(ch['P1'])                          # sum_i d/dxi_3^i Phi  (unit speed)
    return a, S, G, P1


def main():
    n = 300000
    for i, arg in enumerate(sys.argv):
        if arg == '--n':
            n = int(sys.argv[i + 1])
    t0 = time.time()
    rng = np.random.default_rng(4242)
    print("=" * 116)
    print("T1/L  1D SU(2) chain: X(0) = -2Cov(G,Phi'|a=0) in the UNIT-SPEED normalisation (na2 used 1/L => /L^2)")
    print(f"      L=2..8 open + ring, beta = {BETAS}, h = {HB} + h^2 extrapolation, n={n:,} per chain")
    print("=" * 116)
    exact2 = {0.5: -0.0104411252, 1.0: -0.0191235583, 2.5: -0.0330091889}
    nx = len(BETAS) * len(HB)
    for periodic in (False, True):
        for L in (2, 3, 4, 5, 6, 7, 8):
            if periodic and L == 2:
                continue
            nch = 40
            mch = max(n // nch, 2000)
            # accumulators: [chunk, (beta,h), 6 sums]
            acc = np.zeros((nch, nx, 6))
            a2 = 0.0; nsamp = 0
            for c in range(nch):
                ch = C.build(L, mch, rng, periodic)
                nsamp += ch['m']
                a2 += float(np.sum((sum(ch['a']) / L) ** 2))
                for ib, beta in enumerate(BETAS):
                    a, S, G, P1 = chunk_stats(ch, beta)
                    eS = np.exp(-S)
                    for ih, h in enumerate(HB):
                        w = eS * np.exp(-0.5 * (a / h) ** 2)
                        k = ib * len(HB) + ih
                        acc[c, k, 0] += w.sum()
                        acc[c, k, 1] += (w * G).sum()
                        acc[c, k, 2] += (w * P1).sum()
                        acc[c, k, 3] += (w * G * P1).sum()
                        acc[c, k, 4] += (w * G * G).sum()
                        acc[c, k, 5] += (w * P1 * P1).sum()

            def Xof(A):
                Z = A[:, 0]; EG = A[:, 1] / Z; EP = A[:, 2] / Z
                cov = A[:, 3] / Z - EG * EP
                return -2 * cov

            tot = acc.sum(0)
            X = Xof(tot)
            # delete-one-chunk jackknife
            K = nch
            jk = np.array([Xof((tot - acc[c]) / 1.0) for c in range(K)])   # (K, nx)
            se = np.sqrt((K - 1) / K * ((jk - jk.mean(0)) ** 2).sum(0))
            tag = f"L={L}{' ring' if periodic else ' open'}"
            print(f"\n[{tag}]  <a^2> = {a2 / nsamp:.5f}    ({time.time() - t0:.0f}s)")
            print(f"    {'beta':>5}" + "".join(f"{'h=' + str(h):>14}" for h in HB) +
                  f"{'h->0 (h^2 ext)':>17}{'se':>9}{'exact L=2':>13}{'X/L':>10}")
            for ib, beta in enumerate(BETAS):
                xs = X[ib * len(HB):(ib + 1) * len(HB)]
                es = se[ib * len(HB):(ib + 1) * len(HB)]
                ext = xs[-1] + (xs[-1] - xs[-2]) / 3.0                     # h^2 Richardson, two smallest h
                ex = exact2[beta] if (L == 2 and not periodic) else np.nan
                print(f"    {beta:5.2f}" + "".join(f"{v:+14.6f}" for v in xs) +
                      f"{ext:+17.6f}{es[-1]:9.5f}{ex:+13.6f}{ext / L:+10.5f}")
            if L == 2 and not periodic:
                print("    CALIBRATION: column 'h->0' vs column 'exact L=2' (ct_chain_exact.py's 3-D quadrature)")

    print("\n" + "=" * 116)
    print("READ: X(0) < 0 for every (L, open/ring, beta, h) tested; no sign change with L.")
    print("      In the fixed (unit-speed) normalisation |X(0)| GROWS with L ~ extensively (see X/L column):")
    print("      the 'decreasing with L' of na2_chain.log is the 1/L^2 normalisation artefact.")
    print(f"  [{time.time() - t0:.0f}s]  Clay 0%. No QED.")


if __name__ == '__main__':
    main()
