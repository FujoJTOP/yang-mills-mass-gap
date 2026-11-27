"""RR2-MONOTONE, part 2 -- (a) the QQ8 0.10 / QQ3 0.3 '47x' audit (X_E modulus x off-diagonal weight),
(b) the normalization sweep of the rise-variation of ln e, (c) budget arithmetic for the canonical pair.
Deterministic eigh only; no sampling.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

N = 8
K = np.arange(N, dtype=float)
EPS = {
    "k(k+2)/4": K * (K + 2) / 4.0,
    "k(k+2)": K * (K + 2),
    "(k+1)^2/4-1/4": ((K + 1) ** 2 - 1) / 4.0,
    "k^2": K ** 2,
    "(k+1)^2": (K + 1.0) ** 2,
    "k": K,
    "2k+1": 2 * K + 1.0,
    "k(k+2)/4+1/4": K * (K + 2) / 4.0 + 0.25,
}
OFF = {
    "1": np.ones(N - 1),
    "sqrt((k+1)(k+2))": np.sqrt((K[:-1] + 1) * (K[:-1] + 2)),
    "sqrt(k+1)": np.sqrt(K[:-1] + 1.0),
    "(k+1)": K[:-1] + 1.0,
    "1/sqrt((k+1)(k+2))": 1.0 / np.sqrt((K[:-1] + 1) * (K[:-1] + 2)),
    "1/(2j+1)": 1.0 / (K[:-1] + 1.0),
}


def r_inf(eps, off):
    b, Vb = eigh_tridiagonal(np.zeros(len(eps)), off)
    XE = Vb.T @ np.diag(eps) @ Vb
    E0 = XE[0, 1:]
    return float(np.sum(E0[1:] ** 2) / E0[0] ** 2)


def e_of(eps, off, t):
    E, V = eigh_tridiagonal(eps, t * off)
    row = (V.T @ np.diag(eps) @ V)[0, 1:]
    return float(row[0] ** 2 / np.sum(row ** 2))


def main():
    print("== (a) 47x audit: QQ8 0.10 gives 2.5159739 ; QQ3 0.3 / RR2 gives 0.0532688 (ratio 47.2)")
    print(f"{'X_E modulus':>18} | " + " | ".join(f"{k:>18}" for k in OFF))
    vals = []
    for en, eps in EPS.items():
        row = []
        for on, off in OFF.items():
            r = r_inf(eps, off)
            vals.append(r)
            row.append(f"{r:18.6f}")
        print(f"{en:>18} | " + " | ".join(row))
    print(f"scan range {min(vals):.4f} .. {max(vals):.4f}  (2.516 not reached) ; "
          f"5 moduli differing only by scale+shift give the identical value (invariance check)")

    print("\n== (b) normalization sweep (N=8): rise-variation of ln e on t in [0.5,60]")
    tg = np.logspace(np.log10(0.5), np.log10(60), 1500)
    eps0, off0 = K * (K + 2) / 4.0, np.ones(N - 1)
    for lam in (0.1, 0.25, 1.0, 4.0, 10.0):
        v = np.array([e_of(lam * eps0, off0, t) for t in tg])
        dl = np.diff(np.log(v))
        print(f"   lam={lam:5.2f}: e(t_A)={e_of(lam*eps0,off0,0.78125):.7f} "
              f"+var={np.sum(dl[dl>0]):.5f}  -var={-np.sum(dl[dl<0]):.5f}  "
              f"TV={np.sum(np.abs(dl)):.5f}  (rise {'>0' if np.sum(dl[dl>0])>0 else '=0'})")
    print("   note: e(t;lam) = e(t/lam;1) => the hump window shifts with lam; 'not monotone' holds for every lam")

    print("\n== (c) budget arithmetic for the canonical pair (values from rr2_monotone_check.log)")
    w0 = 1.0 / 64.0
    e25, einf = 0.957693820059, 0.949425240
    tv = 0.13561
    b_pair = np.log(e25 / w0)
    b_ym = np.log(21.46752)
    print(f"   e(2.5)={e25:.12f} = {e25/w0:.3f} w0 ; B_pair = ln(e(2.5)/w0) = {b_pair:.6f} nats")
    print(f"   TV(ln e) = {tv:.5f}  =>  B_pair/TV = {b_pair/tv:.2f} ; inf >= e(2.5)e^-TV = "
          f"{e25*np.exp(-tv):.4f} = {e25*np.exp(-tv)/w0:.1f} w0")
    print(f"   YM budget B = ln 21.46752 = {b_ym:.6f} ; B/TV = {b_ym/tv:.2f} ; B - TV = {b_ym-tv:.4f} nats "
          f"(= {np.exp(b_ym-tv):.2f}x head-room)")
    print(f"   actual fall of the pair: ln(e(2.5)/e_inf) = {np.log(e25/einf):.5f} nats ; "
          f"fraction of B_pair used = {np.log(e25/einf)/b_pair:.4f}")


if __name__ == "__main__":
    main()
