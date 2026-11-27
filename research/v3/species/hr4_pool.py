"""HR4-PTEST-POOL -- pool every beta=2.5, Ls=6, Lt=20 file (identical lattice+operator)
into ONE correlator, then run the two-point test.  ZERO MC, one file in memory at a time.

Pool used: gh3_probe.npz (N=160000, seed 20261128) + gh3_x2_{1..8}.npz (N=100000 each,
seeds 20261327..20261334)  ->  N = 960000.  These nine files share Ls=6, Lt=20, beta=2.5,
the same 6 APE levels x 5 observables.  Per-file per-slice means are removed before the
correlation (that is what run3_apebase.correlators does), blocks are pooled across files.

Also reports the noise-floor diagnostic: the t at which |C(t)|/C(0) becomes smaller than
its own error, because that is what decides whether the test is resolvable at all.

Run: py -3.12 hr4_pool.py            (~90 s, 4.6 GB read, no MC)
"""
import os
import sys
import json
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.join(HERE, 'w7')
FILES = ['gh3_probe.npz'] + [f'gh3_x2_{i}.npz' for i in range(1, 9)]
NBF = 10          # jackknife blocks per file  -> 90 blocks
CHUNK = 5000

OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def main():
    t00 = time.time()
    acc = None
    cnt = []
    Lt = None
    K = None
    shape_mid = None
    Ns = []
    nb_off = 0
    nbtot = NBF * len(FILES)
    P("=" * 116)
    P("HR4 PTEST-POOL:  beta=2.5, Ls=6, Lt=20, pooled over 9 files (N = 960000)")
    P("=" * 116)
    for jf, fn in enumerate(FILES):
        f = os.path.join(W7, fn)
        if not os.path.exists(f):
            P(f"[MISS] {fn}")
            continue
        d = np.load(f, mmap_mode='r')
        with d:
            S = d['S']
            Lt = int(d['Lt'])
            Ls = int(d['Ls'])
            beta = float(d['beta'])
            N = S.shape[0]
            shape_mid = S.shape[1:-1]
            K = int(np.prod(shape_mid)) if len(shape_mid) else 1
            mu = np.zeros((K, Lt))
            for i0 in range(0, N, CHUNK):
                mu += np.asarray(S[i0:i0 + CHUNK]).reshape(-1, K, Lt).sum(axis=0)
            mu /= N
            a_j = np.zeros((NBF, K, Lt))
            c_j = np.zeros(NBF, dtype=np.int64)
            for i0 in range(0, N, CHUNK):
                X = np.asarray(S[i0:i0 + CHUNK]).reshape(-1, K, Lt) - mu
                F = np.fft.rfft(X, axis=2)
                ac = np.fft.irfft(F.real ** 2 + F.imag ** 2, n=Lt, axis=2)
                blk = (np.arange(i0, i0 + X.shape[0]) * NBF) // N
                for b in np.unique(blk):
                    m = blk == b
                    a_j[b] += ac[m].sum(axis=0)
                    c_j[b] += int(m.sum())
            if acc is None:
                acc = np.zeros((nbtot, K, Lt))
                cnt = np.zeros(nbtot, dtype=np.int64)
            acc[nb_off:nb_off + NBF] = a_j
            cnt[nb_off:nb_off + NBF] = c_j
            nb_off += NBF
            Ns.append(N)
            P(f"  [{jf+1}/9] {fn:<18} N={N:>7}  beta={beta:g} Ls={Ls} Lt={Lt}  "
              f"beta*={time.time()-t00:.1f}s")
        del d
    Ntot = int(sum(Ns))
    P(f"  pooled N = {Ntot}   blocks = {nbtot}   streams = {K}   shape_mid = {shape_mid}")

    tot = acc.sum(axis=0)
    C = tot / (Ntot * Lt)
    nb = nbtot

    def jk(fn):
        v = fn(C)
        vs = np.array([fn((tot - acc[b]) / ((Ntot - cnt[b]) * Lt)) for b in range(nb)])
        ok = np.isfinite(vs)
        if ok.sum() < nb - 2 or not np.isfinite(v):
            return v, float('nan')
        m = vs[ok].mean()
        return v, float(np.sqrt((nb - 1.0) / nb * np.sum((vs[ok] - m) ** 2)))

    L = 6
    t1, t2 = L // 2, L
    P("")
    P("-" * 116)
    P(f"  L={L}  tau=3 and 6 ;  inequality  C(0)C(6) >= C(3)^2")
    P(f"  {'k':>3} {'lev':>4} {'obs':>4} {'C(0)':>11} {'C(3)':>12} {'+-':>9}"
      f" {'C(6)':>12} {'+-':>9} {'C(3)^2/C(0)':>13} {'D':>11} {'+-sD':>9}"
      f" {'D/sD':>7} {'K':>10} {'verdict':>7}")
    rows = []
    for k in range(K):
        lev, obs = (divmod(k, shape_mid[1]) if len(shape_mid) >= 2 else (0, k))
        v0, e0 = jk(lambda Cc, k=k: Cc[k, 0])
        v1, e1 = jk(lambda Cc, k=k: Cc[k, t1])
        v2, e2 = jk(lambda Cc, k=k: Cc[k, t2])
        vD, eD = jk(lambda Cc, k=k: Cc[k, 0] * Cc[k, t2] - Cc[k, t1] ** 2)
        sig = vD / eD if eD > 0 else np.nan
        Kk = vD / (v1 * v1) if v1 != 0 else np.nan
        rows.append(dict(k=k, lev=lev, obs=obs, C0=v0, Ct1=v1, Ct2=v2, s1=e1, s2=e2,
                         rhs=v1 * v1 / v0, D=vD, sD=eD, sig=float(sig), K=float(Kk),
                         ok=bool(vD >= 0)))
        P(f"  {k:>3} {lev:>4} {obs:>4} {v0:11.4f} {v1:12.4e} {e1:9.1e}"
          f" {v2:12.4e} {e2:9.1e} {v1*v1/v0:13.4e} {vD:11.4e} {eD:9.2e}"
          f" {sig:7.2f} {Kk:10.2f} {'PASS' if vD >= 0 else 'FAIL':>7}")
    npass = sum(1 for r in rows if r['ok'])
    P(f"  streams: PASS {npass}/{K}")
    sigs = [r['sig'] for r in rows if np.isfinite(r['sig'])]
    P(f"  D/sD over {len(sigs)} streams: min {min(sigs):+.2f}  median "
      f"{np.median(sigs):+.2f}  max {max(sigs):+.2f}   (chi2-ish: mean {np.mean(sigs):+.2f})")

    # ---------------- noise floor: where does C(t) stop being resolvable?
    P("")
    P("-" * 116)
    P("  NOISE FLOOR: |C(t)/C(0)| vs its own error, canonical stream k=0 (lev0 raw W11_p0)")
    P(f"  {'t':>3} {'C(t)/C(0)':>13} {'+-':>10} {'S/N':>8}")
    snr = []
    for t in range(Lt):
        vt, et = jk(lambda Cc, t=t: Cc[0, t])
        v0, e0 = jk(lambda Cc: Cc[0, 0])
        snr.append(abs(vt) / et if et > 0 else np.inf)
        P(f"  {t:>3} {vt/v0:13.5e} {et/v0:10.2e} {snr[-1]:8.2f}")
    tlast = max([t for t in range(Lt) if snr[t] >= 3.0], default=-1)
    P(f"  S/N >= 3 up to t = {tlast} ;  t = L = {L} has S/N = {snr[L]:.2f}")
    P(f"  => the two-point test at tau = L needs S/N at t = L; that is the binding constraint.")

    # ---------------- K(t) = Var_nu(e^{-st})/C(t)^2 on the RESOLVED times
    P("")
    P("-" * 116)
    P("  K(t) := (C(0)C(2t) - C(t)^2)/C(t)^2 = Var_nu(e^{-st})/C(t)^2  (>=0 always)")
    P("  measured only where C(2t) is above its own error -- that is the only usable window")
    P(f"  {'t':>3} {'C(t)/C(0)':>12} {'C(2t)/C(0)':>12} {'S/N(2t)':>8} {'K(t)':>10} {'+-':>9}")
    Krows = []
    for t in range(1, min(5, Lt // 2)):
        vt, et = jk(lambda Cc, t=t: Cc[0, t])
        v2, e2 = jk(lambda Cc, t=t: Cc[0, 2 * t])
        v0, e0 = jk(lambda Cc: Cc[0, 0])
        vK, eK = jk(lambda Cc, t=t: (Cc[0, 0] * Cc[0, 2 * t] - Cc[0, t] ** 2) / Cc[0, t] ** 2)
        ok = abs(v2) / e2 >= 3.0 if e2 > 0 else False
        P(f"  {t:>3} {vt/v0:12.5e} {v2/v0:12.5e} {abs(v2)/e2:8.2f} {vK:10.4f} {eK:9.4f}"
          f"   {'<-resolved' if ok else ''}")
        Krows.append(dict(t=t, Ct=float(vt / v0), C2t=float(v2 / v0), snr2t=float(abs(v2) / e2),
                          K=float(vK), sK=float(eK), resolved=bool(ok)))
    P("")
    P("  COST for the test itself: the margin it compares against is C(L/2)^2/C(0);")
    s6 = rows[0]['s2'] / rows[0]['C0']
    marg = rows[0]['rhs'] / rows[0]['C0']
    P(f"    measured sigma(C(6)/C(0)) = {s6:.2e} at N = {Ntot}")
    P(f"    the margin C(3)^2/C(0)^2 = {marg:.2e}  ->  need sigma below that:")
    fac = (s6 / marg) ** 2
    P(f"    N needed ~ {Ntot*fac:.2e}  (x{fac:.0f})  => the test is not a data route at L=6.")
    json.dump(dict(N=Ntot, blocks=nb, Lt=Lt, Ls=6, beta=2.5, files=FILES,
                   rows=rows, snr=snr, tlast=int(tlast), K=Krows),
              open(os.path.join(HERE, 'hr4_pool.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'hr4_pool.log'), 'w').write("\n".join(OUT) + "\n")
    P(f"\n[ok] wall {time.time()-t00:.1f}s ; wrote hr4_pool.json / .log")


if __name__ == '__main__':
    main()
