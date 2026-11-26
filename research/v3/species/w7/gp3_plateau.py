"""GP3 -- does the transfer-matrix machine (tm_gap family) have an effective-mass plateau?

Single question: at beta = 3.0 (finer than the 2.5 batch), does m_eff(t) = ln[C(t)/C(t+1)]
show a t-window where it is flat within 1 sigma AND S/N >= 5 ?

Operator: run3/b01 `lev0 = raw` -- zero-momentum per-slice spatial plaquette sum,
    S(t) = (1/3) sum_{space} sum_{3 spatial planes} (1/2) Re tr U_p
The normalisation is /3 (average over the 3 spatial planes), NOT /4 and NOT /6.
m_eff(t) is normalisation-free.

Core is IMPORTED UNCHANGED from run3_apebase (same heat bath / same _chains / same
correlators); nlev = 1 so no APE level is ever built.  Nothing in the existing files
is modified.

usage:
  py -3.12 gp3_plateau.py run --beta 3.0 --Ls 8 --Lt 16 --n 48000 --gap 5 --therm 900
  py -3.12 gp3_plateau.py analyze            # re-read every b01_b*.npz at zero MC cost
"""
import os
import sys
import json
import time
import argparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run3_apebase import _chains, index_tables, hot, correlators  # noqa: E402

OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def jack_corr(Sm, Lt, nb=20):
    """C(t), sigma_C(t) (delete-one block jackknife), m_eff(t) = ln C(t)/C(t+1) +- 1 sigma."""
    N = Sm.shape[0]
    C = correlators(Sm, Lt)
    idx = np.array_split(np.arange(N), nb)
    Cs = np.empty((nb, Lt))
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        Cs[b] = correlators(Sm[keep], Lt)
    sig = np.sqrt((nb - 1.0) / nb * np.nansum((Cs - C) ** 2, axis=0))
    meff = np.full(Lt - 1, np.nan)
    meff_e = np.full(Lt - 1, np.nan)
    for t in range(Lt - 1):
        if C[t] > 0 and C[t + 1] > 0:
            meff[t] = np.log(C[t] / C[t + 1])
        v = np.array([np.log(c[t] / c[t + 1]) if (c[t] > 0 and c[t + 1] > 0) else np.nan
                      for c in Cs])
        ok = np.isfinite(v)
        if ok.sum() >= nb - 2:
            mm = v[ok].mean()
            meff_e[t] = np.sqrt((nb - 1.0) / nb * np.nansum((v[ok] - mm) ** 2))
    return C, sig, meff, meff_e


def snr(C, sig):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(sig > 0, np.abs(C) / sig, np.inf)


def table(tag, C, sig, meff, meff_e, Lt, tmax=6):
    s = snr(C, sig)
    P("")
    P(f"--- {tag}  (Lt={Lt}) ---")
    P(f"{'t':>3} {'C(t)/C(0)':>13} {'+-':>10} {'S/N(t)':>8} {'m_eff(t)':>9} {'+-':>8} "
      f"{'[M] cosh m':>11}")
    for t in range(0, min(tmax, Lt - 1) + 1):
        c = C[t] / C[0]
        ce = sig[t] / abs(C[0])
        me = meff[t] if t < Lt - 1 else np.nan
        mee = meff_e[t] if t < Lt - 1 else np.nan
        # cosh effective mass (periodic time) as cross-check, only if defined
        cm, c0, cp = C[(t - 1) % Lt], C[t % Lt], C[(t + 1) % Lt]
        mc = np.nan
        if t >= 1 and min(c0, cm + cp) > 1e-14 and (cm + cp) / (2 * c0) >= 1.0:
            mc = float(np.arccosh((cm + cp) / (2 * c0)))
        P(f"{t:>3} {c:>13.5e} {ce:>10.2e} {s[t]:>8.2f} {me:>9.3f} {mee:>8.3f} {mc:>11.3f}")
    return s


def plateau_scan(meff, meff_e, s, tlo=0, thi=None):
    """Task definition: window [t1,t2] with >=2 consecutive t, all m_eff(t) mutually
    consistent within 1 sigma, and S/N>=5 across the window.  Returns list of windows."""
    n = len(meff)
    thi = n - 1 if thi is None else thi
    wins = []
    for t1 in range(tlo, thi + 1):
        for t2 in range(t1 + 1, thi + 1):
            ts = list(range(t1, t2 + 1))
            if any(not np.isfinite(meff[t]) or not np.isfinite(meff_e[t]) for t in ts):
                continue
            if any(s[t] < 5.0 for t in ts):
                continue
            worst = 0.0
            ok = True
            for a in range(len(ts)):
                for b in range(a + 1, len(ts)):
                    i, j = ts[a], ts[b]
                    d = abs(meff[i] - meff[j])
                    sd = np.hypot(meff_e[i], meff_e[j])
                    if d > sd:
                        ok = False
                        break
                    worst = max(worst, d / sd if sd > 0 else np.inf)
                if not ok:
                    break
            if ok:
                wins.append((t1, t2, worst))
    return wins


def point(beta, Ls, Lt, n, gap, therm, chains, threads, seed, alpha=0.5, nlev=1):
    import numba
    numba.set_num_threads(threads)
    per = max(1, n // chains)
    N = per * chains
    plus, minus = index_tables(Ls, Lt)
    U = hot(chains, Ls, Lt, seed)
    rng = np.empty((chains, 2), np.uint64)
    ss = np.uint64(seed)
    for c in range(chains):
        ss = ss * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
        rng[c, 0] = ss | np.uint64(1)
        s1 = ss ^ np.uint64(0x9E3779B97F4A7C15)
        rng[c, 1] = s1 if s1 != np.uint64(0) else np.uint64(0x2545F4914F6CDD1D)
    S = np.zeros((N, nlev, Lt))
    plaq = np.zeros(N)

    P(f"[T] sites={Ls**3*Lt} links={4*Ls**3*Lt} U={U.nbytes*chains/1e6:.1f} MB  "
      f"threads={threads} chains={chains} per={per}")
    t0 = time.time()
    _chains(U[:1], rng[:1], plus, minus, beta / 2.0, 2, gap, 2, alpha, nlev, Lt,
            np.zeros((2, nlev, Lt)), np.zeros(2))
    t_comp = time.time() - t0
    t1 = time.time()
    _chains(U, rng, plus, minus, beta / 2.0, per, gap, therm, alpha, nlev, Lt, S, plaq)
    wall = time.time() - t1
    tot_sweeps = N * gap + therm
    P(f"[E] {N} cfg in {wall:.1f}s ({N/wall:.1f} cfg/s, {tot_sweeps/wall:.0f} sweeps/s)  "
      f"compile {t_comp:.1f}s")
    P(f"[E] <1/2 Re tr U_p> = {plaq.mean():.5f} +- {plaq.std(ddof=1)/np.sqrt(N):.5f}  "
      f"N={N} gap={gap} therm={therm}")
    return S[:, 0, :], plaq, wall, t_comp, N


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    r = sub.add_parser('run')
    r.add_argument('--beta', type=float, default=3.0)
    r.add_argument('--Ls', type=int, default=8)
    r.add_argument('--Lt', type=int, default=16)
    r.add_argument('--n', type=int, default=48000)
    r.add_argument('--gap', type=int, default=5)
    r.add_argument('--therm', type=int, default=900)
    r.add_argument('--chains', type=int, default=8)
    r.add_argument('--threads', type=int, default=8)
    r.add_argument('--seed', type=int, default=20261126)
    r.add_argument('--tag', type=str, default='')
    a = sub.add_parser('analyze')
    a.add_argument('--tmax', type=int, default=6)
    a.add_argument('--files', type=str, default='b01_b3.npz,b01_b4.npz,b01_b6.npz,b01_b10.npz')
    args = ap.parse_args()

    if args.mode == 'analyze':
        for f in args.files.split(','):
            f = f.strip()
            if not f:
                continue
            p = os.path.join(HERE, f)
            if not os.path.exists(p):
                P(f"[MISS] {f}")
                continue
            d = np.load(p)
            S = d['S'][:, 0, :]
            Lt = int(d['Lt'])
            C, sig, meff, meff_e = jack_corr(S, Lt)
            tag = (f"{f}  beta={float(d['beta']):g} N={int(d['N'])} Ls=8 "
                   f"<plaq>={float(np.mean(d['plaq'])):.5f}")
            s = table(tag, C, sig, meff, meff_e, Lt, args.tmax)
            wins = plateau_scan(meff, meff_e, s, 0, min(args.tmax, Lt - 2))
            P(f"[verdict] S/N>=5 up to t={int(np.max(np.nonzero(s[:min(args.tmax, Lt-1)+1] >= 5.0)[0])) if np.any(s[:min(args.tmax, Lt-1)+1] >= 5.0) else -1}"
              f" ; plateau windows(1sig & S/N>=5): {wins if wins else 'NONE'}")
        open(os.path.join(HERE, 'gp3_analyze_preexisting.log'), 'w').write("\n".join(OUT) + "\n")
        return

    tag = args.tag or f"b{args.beta:g}_Ls{args.Ls}_Lt{args.Lt}"
    P("=" * 96)
    P(f"GP3  SU(2) Wilson  beta={args.beta:g}  Ls={args.Ls}^3 x Lt={args.Lt}  "
      f"operator=run3/b01 lev0 raw (normalised /3 over the 3 spatial planes)")
    P(f"N={args.n} ({args.chains} chains x {args.n//args.chains})  gap={args.gap}  "
      f"therm={args.therm}  seed={args.seed}")
    P("=" * 96)
    T0 = time.time()
    S, plaq, wall, t_comp, N = point(args.beta, args.Ls, args.Lt, args.n, args.gap,
                                     args.therm, args.chains, args.threads, args.seed)
    C, sig, meff, meff_e = jack_corr(S, args.Lt)
    s = table(f"beta={args.beta:g} Ls={args.Ls} Lt={args.Lt} N={N}", C, sig, meff, meff_e,
              args.Lt, 6)
    tmax = min(6, args.Lt - 2)
    wins = plateau_scan(meff, meff_e, s, 0, tmax)
    P("")
    P(f"[S/N>=5] highest t with S/N>=5: "
      f"{max([t for t in range(min(7, args.Lt)) if s[t] >= 5.0], default=-1)}")
    P(f"[plateau] windows (>=2 consecutive t, 1 sigma, S/N>=5): {wins if wins else 'NONE'}")
    P(f"[wall] total {time.time()-T0:.1f}s (MC {wall:.1f}s)")
    np.savez(os.path.join(HERE, f"gp3_{tag}.npz"), S=S, plaq=plaq,
             beta=np.float64(args.beta), Ls=np.int64(args.Ls), Lt=np.int64(args.Lt),
             N=np.int64(N), gap=np.int64(args.gap), therm=np.int64(args.therm),
             chains=np.int64(args.chains), seed=np.int64(args.seed),
             wall=np.float64(wall), C=C, sigC=sig, meff=meff, meff_e=meff_e)
    json.dump(dict(beta=args.beta, Ls=args.Ls, Lt=args.Lt, N=N, gap=args.gap,
                   therm=args.therm, chains=args.chains, seed=args.seed, wall=wall,
                   plaq=float(plaq.mean()),
                   plaq_err=float(plaq.std(ddof=1) / np.sqrt(N)),
                   C=C.tolist(), sigmaC=sig.tolist(), meff=meff.tolist(),
                   meff_e=meff_e.tolist(), snr=s.tolist(), plateau=wins),
              open(os.path.join(HERE, f"gp3_{tag}.json"), 'w'), indent=1)
    open(os.path.join(HERE, f"gp3_{tag}.log"), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
