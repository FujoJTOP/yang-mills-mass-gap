"""GH3 step 0 -- FREE diagnosis on archived S(t) arrays (no new MC).

Reads:
  run1_b25_l6.npz            beta=2.5, Ls=6,  Lt=12, N=6000, levels raw + APE n=1..4  (alpha=0.5)
  run3_L8_b2.5_N12k.npz      beta=2.5, Ls=8,  Lt=12, N=12000, levels raw + TAPE n=1,2 (alpha=0.5)
  f1sq_b2.2_L6.npz           beta=2.2, Ls=6,  Lt=6,  N=4000, levels s=0,2,4,8          (alpha=0.5)

Prints, per level:
  C(0), Phi(t)=C(t)/C(0) for t=1..Lt/2, JK error, S/N(t), m_eff(t)=ln C(t)/C(t+1),
  and the noise-model check sigma_model(t>=1) = C(0)*sqrt(2/(N*Lt)).

No file outside this directory is written; nothing is modified.
Run: py -3.12 gh3_diag.py   ->  gh3_diag.log
"""
import os, json
import numpy as np
from run3_apebase import correlators, jack

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def report(tag, S, Lt, names, beta, Ls, N):
    P("=" * 100)
    P(f"{tag}   beta={beta}  Ls={Ls}  Lt={Lt}  N={N}   levels={len(names)}")
    P("=" * 100)
    rows = {}
    for lev, nm in enumerate(names):
        Sm = S[:, lev, :]
        C = correlators(Sm, Lt)
        if C[0] <= 0:
            P(f"{nm}: C(0)={C[0]:.4g} non-positive")
            continue
        # jackknife errors on every C(t)
        err = np.zeros(Lt)
        for t in range(Lt):
            _, e, _ = jack(Sm, Lt, 20, lambda Cc, t=t: Cc[t])
            err[t] = e
        sig_model = C[0] * np.sqrt(2.0 / (N * Lt))
        P(f"\n-- {nm}:  C(0)={C[0]:.4f}   sigma_model(t>=1)={sig_model:.3f}")
        hdr = "   t    C(t)/C(0)        +-        S/N     m_eff      +-"
        P(hdr)
        for t in range(1, Lt // 2 + 1):
            phi = C[t] / C[0]
            e = err[t] / C[0]
            snr = C[t] / err[t] if err[t] > 0 else float('nan')
            # m_eff from C(t)->C(t+1), error by linear propagation
            if t + 1 <= Lt - 1 and C[t + 1] > 0 and C[t] > 0:
                me = np.log(C[t] / C[t + 1])
                eme = np.hypot(err[t] / C[t], err[t + 1] / C[t + 1])
            else:
                me, eme = float('nan'), float('nan')
            P(f"{t:4d} {phi:12.5f} {e:9.5f} {snr:8.2f} {me:9.4f} {eme:8.4f}")
        rows[nm] = dict(C0=float(C[0]), phi=[float(C[t] / C[0]) for t in range(Lt)],
                        err=[float(err[t] / C[0]) for t in range(Lt)],
                        snr=[float(C[t] / err[t]) if err[t] > 0 else None for t in range(Lt)],
                        sigma_model=float(sig_model))
    return rows


def main():
    res = {}

    # ---------------- dataset A: beta=2.5, L=6, 5 levels
    f = os.path.join(HERE, 'run1_b25_l6.npz')
    if os.path.exists(f):
        d = np.load(f)
        S = d['S']
        res['b25_L6'] = report('A  RUN1  beta=2.5 L=6', S, int(d['Lt']),
                               ['raw', 'APE1', 'APE2', 'APE3', 'APE4'],
                               float(d['beta']), int(d['Ls']), int(d['n']))

    # ---------------- dataset B: beta=2.5, L=8, 3 levels
    f = os.path.join(HERE, 'run3_L8_b2.5_N12k.npz')
    if os.path.exists(f):
        d = np.load(f)
        S = d['S']
        res['b25_L8'] = report('B  RUN3  beta=2.5 L=8', S, int(d['Lt']),
                               ['raw', 'TAPE1', 'TAPE2'],
                               float(d['beta']), int(d['Ls']), int(d['n']))

    # ---------------- dataset C: beta=2.2, L=6, s=0,2,4,8 (T=6)
    f = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'unconditional', 'f1sq_b2.2_L6.npz'))
    if os.path.exists(f):
        d = np.load(f)
        S = np.stack([d['S0'], d['S2'], d['S4'], d['S8']], axis=1)
        res['b22_L6'] = report('C  F1SQ  beta=2.2 L=6', S, int(S.shape[2]),
                               ['s0', 's2', 's4', 's8'],
                               float(d['beta']), int(d['L']), int(d['N']))
    else:
        P(f"TOOL-NOTE: {f} missing")

    json.dump(res, open(os.path.join(HERE, 'gh3_diag.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'gh3_diag.log'), 'w').write("\n".join(OUT) + "\n")
    P("\nwrote gh3_diag.log / gh3_diag.json")


if __name__ == '__main__':
    main()
