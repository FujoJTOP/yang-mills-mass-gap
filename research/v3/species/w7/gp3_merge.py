"""GP3 merge -- pool independent GP3 points with identical parameters and re-run the verdict.

Only same (beta, Ls, Lt, gap, therm, operator) runs may be pooled; each run is an
independent set of chains, so concatenating the config axis keeps the chain-major
layout and the delete-one block jackknife stays valid (blocks are 5e3-config
contiguous stretches, >> tau_int).

usage:  py -3.12 gp3_merge.py --files gp3_b3_Ls8_Lt16_N60k.npz,gp3_b3_Ls8_Lt16_N40k.npz \
                             --tag b3_Ls8_Lt16_N100k
"""
import os
import sys
import json
import argparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gp3_plateau import P, OUT, jack_corr, table, plateau_scan  # noqa: E402


def tau_int(x):
    x = np.asarray(x, float)
    n = len(x)
    y = x - x.mean()
    v = float(np.dot(y, y)) / n
    if v <= 0:
        return 0.5
    t = 0.5
    for lag in range(1, min(n // 4, 300)):
        c = float(np.dot(y[:-lag], y[lag:])) / (n - lag) / v
        if c <= 0:
            break
        t += c
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--files', required=True)
    ap.add_argument('--tag', default='merged')
    ap.add_argument('--tmax', type=int, default=6)
    a = ap.parse_args()

    Ss, plaqs, walls = [], [], []
    meta = None
    for f in a.files.split(','):
        d = np.load(os.path.join(HERE, f.strip()))
        m = (float(d['beta']), int(d['Ls']), int(d['Lt']), int(d['gap']), int(d['therm']))
        if meta is not None and m != meta:
            raise SystemExit(f"parameter mismatch: {m} vs {meta}")
        meta = m
        Sf = d['S']
        Ss.append(Sf[:, 0, :] if Sf.ndim == 3 else Sf)
        plaqs.append(d['plaq'])
        walls.append(float(d['wall']))
        P(f"[in] {f.strip()}  N={Ss[-1].shape[0]}  <plaq>={plaqs[-1].mean():.5f}")
    S = np.concatenate(Ss, axis=0)
    plaq = np.concatenate(plaqs, axis=0)
    beta, Ls, Lt, gap, therm = meta
    N = S.shape[0]
    P("")
    P("=" * 96)
    P(f"GP3-MERGED  beta={beta:g}  Ls={Ls}^3 x Lt={Lt}  N={N}  gap={gap}  therm={therm}  "
      f"operator=run3/b01 lev0 raw (/3 over 3 spatial planes)")
    P("=" * 96)
    P(f"[E] <1/2 Re tr U_p> = {plaq.mean():.8f} +- {plaq.std(ddof=1)/np.sqrt(N):.8f}")
    tau = float(np.mean([tau_int(p) for p in plaqs]))
    P(f"[E] tau_int(plaq) at gap={gap}: {tau:.2f} cfg  => N_eff = {N/(2*tau):.0f}")
    C, sig, meff, meff_e = jack_corr(S, Lt)
    s = table(f"MERGED beta={beta:g} Ls={Ls} Lt={Lt} N={N}", C, sig, meff, meff_e, Lt, a.tmax)
    tmax = min(a.tmax, Lt - 2)
    wins = plateau_scan(meff, meff_e, s, 0, tmax)
    wins_nocut = plateau_scan(meff, meff_e, np.full(Lt, 1e9), 0, tmax)
    P("")
    hi = max([t for t in range(min(7, Lt)) if s[t] >= 5.0], default=-1)
    P(f"[S/N>=5] highest t with S/N>=5: {hi}")
    P(f"[plateau A] windows, S/N>=5 required: {wins if wins else 'NONE'}")
    P(f"[plateau B] windows, S/N cut removed (geometry only): "
      f"{[(w[0], w[1], round(float(w[2]), 2)) for w in wins_nocut] if wins_nocut else 'NONE'}")
    # drift over the only 2-step candidate [1,2]
    if np.isfinite(meff[1]) and np.isfinite(meff[2]):
        d = meff[1] - meff[2]
        sd = np.hypot(meff_e[1], meff_e[2]) if np.isfinite(meff_e[2]) else np.nan
        P(f"[drift] m_eff(1)-m_eff(2) = {d:+.3f} +- {sd:.3f}  ({abs(d)/sd:.1f} sigma)"
          if np.isfinite(sd) else
          f"[drift] m_eff(1)-m_eff(2) = {d:+.3f}  (sigma(m_eff(2)) undefined)")
    P(f"[cost] to reach S/N(2)=5: N >= {N*(5.0/s[2])**2:.0f} configs ; "
      f"to reach 3-sigma drift at the observed gap: N ~ {N*(3*sd/abs(d))**2:.0f}")
    np.savez(os.path.join(HERE, f"gp3_{a.tag}.npz"), S=S, plaq=plaq,
             beta=np.float64(beta), Ls=np.int64(Ls), Lt=np.int64(Lt), N=np.int64(N),
             gap=np.int64(gap), therm=np.int64(therm), wall=np.float64(sum(walls)),
             C=C, sigC=sig, meff=meff, meff_e=meff_e, tau_plaq=np.float64(tau))
    json.dump(dict(beta=beta, Ls=Ls, Lt=Lt, N=N, gap=gap, therm=therm, wall=sum(walls),
                   plaq=float(plaq.mean()),
                   plaq_err=float(plaq.std(ddof=1)/np.sqrt(N)), tau_plaq=tau,
                   C=C.tolist(), sigmaC=sig.tolist(), meff=meff.tolist(),
                   meff_e=meff_e.tolist(), snr=s.tolist(),
                   plateau_A=wins, plateau_B=[[int(w[0]), int(w[1]), float(w[2])]
                                              for w in wins_nocut]),
              open(os.path.join(HERE, f"gp3_{a.tag}.json"), 'w'), indent=1)
    open(os.path.join(HERE, f"gp3_{a.tag}.log"), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
