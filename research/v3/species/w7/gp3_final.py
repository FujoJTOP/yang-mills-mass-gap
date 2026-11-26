"""GP3 final -- correlated weighted mean of the m_eff(t) window + two-loop scale ratio.

Reads the saved GP3 npz files, rebuilds the delete-one block jackknife of
(m_eff(1), m_eff(2)), and forms the BLUE (inverse-covariance) mean of the window.
Also prints the two-loop SU(2) ratio a(beta1)/a(beta2) for the scaling check.

usage:  py -3.12 gp3_final.py
"""
import os
import sys
import math

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gp3_plateau import correlators, P, OUT  # noqa: E402

NB = 20


def pair_blocks(S, Lt, nb=NB):
    N = S.shape[0]
    idx = np.array_split(np.arange(N), nb)
    m1 = np.full(nb, np.nan)
    m2 = np.full(nb, np.nan)
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        C = correlators(S[keep], Lt)
        if C[1] > 0 and C[2] > 0:
            m1[b] = np.log(C[1] / C[2])
        if C[2] > 0 and C[3] > 0:
            m2[b] = np.log(C[2] / C[3])
    return m1, m2


def blue(m, cov):
    Cv = np.asarray(cov)
    Ci = np.linalg.inv(Cv)
    one = np.ones(len(m))
    w = Ci @ one / (one @ Ci @ one)
    val = float(w @ m)
    err = float(math.sqrt(1.0 / (one @ Ci @ one)))
    return val, err, w


def twoloop_ratio(b1, b2, N=2):
    """a(beta1)/a(beta2) at two loops, Wilson beta = 2N/g^2."""
    b0 = 11.0 * N / (48 * math.pi ** 2)
    b1c = 34.0 * N ** 2 / (3 * (16 * math.pi ** 2) ** 2)
    def inv(beta):
        g2 = 2.0 * N / beta
        return (b0 * g2) ** (-b1c / (2 * b0 * b0)) * math.exp(-1.0 / (2 * b0 * g2))
    return inv(b1) / inv(b2), b0, b1c


def main():
    res = {}
    for tag in ('b3_Ls8_Lt16_N100k', 'b25_Ls8_Lt16_N24k'):
        d = np.load(os.path.join(HERE, f"gp3_{tag}.npz"))
        Sd = d['S']
        S = Sd if Sd.ndim == 2 else Sd[:, 0, :]
        Lt = int(d['Lt'])
        beta = float(d['beta'])
        C = correlators(S, Lt)
        m = np.array([np.log(C[1] / C[2]), np.log(C[2] / C[3])])
        m1, m2 = pair_blocks(S, Lt)
        ok = np.isfinite(m1) & np.isfinite(m2)
        V = np.vstack([m1[ok], m2[ok]])
        vbar = V.mean(axis=1, keepdims=True)
        Cov = (ok.sum() - 1) / ok.sum() * (V - vbar) @ (V - vbar).T
        rho = Cov[0, 1] / math.sqrt(Cov[0, 0] * Cov[1, 1])
        val, err, w = blue(m, Cov)
        P("")
        P(f"=== beta={beta:g} Ls={int(d['Ls'])} Lt={Lt} N={int(d['N'])} gap={int(d['gap'])} "
          f"therm={int(d['therm'])} ===")
        P(f"  m_eff(1) = {m[0]:.4f} +- {math.sqrt(Cov[0,0]):.4f} ; "
          f"m_eff(2) = {m[1]:.4f} +- {math.sqrt(Cov[1,1]):.4f} ; rho = {rho:+.3f}")
        P(f"  BLUE over window [1,2]: m_eff = {val:.4f} +- {err:.4f}  (weights {w[0]:.3f}/{w[1]:.3f})")
        P(f"  naive diff m_eff(1)-m_eff(2) = {m[0]-m[1]:+.4f} +- "
          f"{math.sqrt(Cov[0,0]+Cov[1,1]-2*Cov[0,1]):.4f}")
        res[beta] = dict(m1=float(m[0]), m1e=float(math.sqrt(Cov[0, 0])),
                         m2=float(m[1]), m2e=float(math.sqrt(Cov[1, 1])), rho=float(rho),
                         blue=val, blue_e=err, N=int(d['N']))
    if 2.5 in res and 3.0 in res:
        a, b = res[2.5], res[3.0]
        dm = a['m1'] - b['m1']
        de = math.hypot(a['m1e'], b['m1e'])
        P("")
        P(f"[scaling] m_eff(1): beta=2.5 {a['m1']:.3f}+-{a['m1e']:.3f} -> "
          f"beta=3.0 {b['m1']:.3f}+-{b['m1e']:.3f} ; m(2.5)-m(3.0) {dm:+.3f} +- {de:.3f} "
          f"({abs(dm)/de:.1f} sigma => m_eff goes "
          f"{'UP' if dm < 0 else 'DOWN'} as beta increases)")
        P(f"[scaling] observed m_eff(2.5)/m_eff(3.0) = {a['m1']/b['m1']:.3f} "
          f"(BLUE: {a['blue']/b['blue']:.3f})")
        for bb in (3.3, 3.6, 4.0):
            r, b0, b1c = twoloop_ratio(2.5, bb)
            P(f"[scaling] two-loop a(2.5)/a({bb:g}) = {r:.2f}  (b0={b0:.5f}, b1={b1c:.6f})")
        r, b0, b1c = twoloop_ratio(2.5, 3.0)
        P(f"[scaling] two-loop a(2.5)/a(3.0) = {r:.3f}  => scaling demands "
          f"m*a(2.5)/m*a(3.0) = {r:.3f} ; observed {a['m1']/b['m1']:.3f} "
          f"({r/(a['m1']/b['m1']):.1f}x mismatch)")
    open(os.path.join(HERE, 'gp3_final.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
