"""HR4-PTEST -- two-point complete-monotonicity test, run on the existing npz (ZERO MC).

Test (HZ1-PLATEAU.md:19,143,295 [10]):
      C_L(L)  >=  C_L(L/2)^2 / C_L(0),        L = spatial box size,  tau = L.
Equivalently  D := C(0)C(L) - C(L/2)^2 >= 0 ,  K := D/C(L/2)^2 = R - 1 >= 0.

Estimator is copied EXACTLY from run3_apebase.correlators (per-time-slice mean removed,
periodic time wrap, /Lt).  Errors: delete-one-block jackknife (nb=20) on the *combination*
D (correlated), not quadrature.

Memory discipline: one file at a time, S read in chunks through mmap, the full array is
never materialised.  No MC, no fits, no new random numbers.

Tier of faithfulness (reflection term in the wrapped correlator e^{-m(Lt-t)}):
    A  Lt >= 2L   the cylinder value is recovered to O(e^{-m(Lt-L)} / e^{-mL}) <= 1
    B  L <= Lt-1 < 2L   t=L exists but is reflection-dominated -> the *sign* is still
                        guaranteed by positivity, the *magnitude* is not a P readout
    C  Lt <= L-1  no data at t=L  -> NOT TESTABLE (nothing is computed)

Run:  py -3.12 hr4_ptest.py --set core
      py -3.12 hr4_ptest.py --set probe      (gh3_probe + one gh3_x2 seed, 1.25 GB read)
      py -3.12 hr4_ptest.py --set x2         (all eight gh3_x2 seeds, 3.8 GB read)
"""
import os
import sys
import json
import time
import argparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.join(HERE, 'w7')
sys.path.insert(0, HERE)
rd = None

OUT = []
ROWS = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def wrapped_correlator(path, chunk=4000, nb=20):
    """C[k,t] (and per-block jackknife) exactly like run3_apebase.correlators."""
    d = np.load(path, mmap_mode='r')
    with d:
        S = d['S']
        Lt = int(d['Lt'])
        N = S.shape[0]
        K = int(np.prod(S.shape[1:-1])) if S.ndim > 2 else 1
        shape_mid = S.shape[1:-1]
        # ---- pass 1: global per-slice mean (exactly what correlators() subtracts)
        mu = np.zeros((K, Lt))
        for i0 in range(0, N, chunk):
            X = np.asarray(S[i0:i0 + chunk]).reshape(-1, K, Lt)
            mu += X.sum(axis=0)
        mu /= N
        # ---- pass 2: circular autocorrelation via FFT, block-accumulated
        acc = np.zeros((nb, K, Lt))
        cnt = np.zeros(nb, dtype=np.int64)
        for i0 in range(0, N, chunk):
            X = np.asarray(S[i0:i0 + chunk]).reshape(-1, K, Lt)
            n = X.shape[0]
            X = X - mu
            F = np.fft.rfft(X, axis=2)
            ac = np.fft.irfft(F.real ** 2 + F.imag ** 2, n=Lt, axis=2)  # sum_s x_s x_{s+t}
            blk = (np.arange(i0, i0 + n) * nb) // N
            for b in np.unique(blk):
                m = blk == b
                acc[b] += ac[m].sum(axis=0)
                cnt[b] += int(m.sum())
        return acc, cnt, N, K, Lt, shape_mid, d


def jack(acc, cnt, N, Lt, fn, nb=20):
    """value and delete-one-block jackknife error of fn(C[k,t] vector)."""
    tot = acc.sum(axis=0)
    C = tot / (N * Lt)
    v = fn(C)
    vs = np.empty(nb)
    for b in range(nb):
        Cb = (tot - acc[b]) / ((N - cnt[b]) * Lt)
        vs[b] = fn(Cb)
    ok = np.isfinite(vs)
    if ok.sum() < nb - 2 or not np.isfinite(v):
        return v, float('nan'), C
    m = vs[ok].mean()
    e = np.sqrt((nb - 1.0) / nb * np.sum((vs[ok] - m) ** 2))
    return v, e, C


def tier(L, Lt):
    if L > Lt - 1:
        return 'C'
    if Lt >= 2 * L:
        return 'A'
    return 'B'


def do_file(path, tag, lev=None, obs=None, Ls_override=None, label=None):
    f = os.path.join(W7, path)
    if not os.path.exists(f):
        P(f"[MISS] {path}")
        return None
    t0 = time.time()
    acc, cnt, N, K, Lt, shape_mid, d = wrapped_correlator(f)
    d = np.load(f, mmap_mode='r')
    with d:
        Ls = int(d['Ls']) if 'Ls' in d.files else Ls_override
        beta = float(d['beta']) if 'beta' in d.files else float('nan')
        seed = int(d['seed']) if 'seed' in d.files else -1
        gap = int(d['gap']) if 'gap' in d.files else -1
    if Ls is None:
        P(f"[SKIP] {path}: no Ls stored and none supplied")
        return None
    L = Ls
    T = tier(L, Lt)
    t1, t2 = L // 2, L
    hdr = (f"{label or path}  Ls={Ls} Lt={Lt} N={N} beta={beta:g} seed={seed} gap={gap}  "
           f"streams={K} shape_mid={shape_mid}  tier={T}  wall={time.time()-t0:.1f}s")
    P("")
    P("-" * 118)
    P(hdr)
    P("-" * 118)
    if T == 'C':
        P(f"  NOT TESTABLE: t = L = {L} exceeds Lt-1 = {Lt-1}; the sequence stops before the")
        P(f"  matching time.  Nothing is computed (no extrapolation, no substitution).")
        ROWS.append(dict(file=path, Ls=Ls, Lt=Lt, N=N, beta=beta, seed=seed, K=K,
                         tier=T, testable=False))
        return None

    def C_sel(C, k, t):
        return C[k, t]

    rows = []
    for k in range(K):
        v0, e0, C = jack(acc, cnt, N, Lt, lambda Cc, k=k: Cc[k, 0])
        v1, e1, _ = jack(acc, cnt, N, Lt, lambda Cc, k=k: Cc[k, t1])
        v2, e2, _ = jack(acc, cnt, N, Lt, lambda Cc, k=k: Cc[k, t2])
        vD, eD, _ = jack(acc, cnt, N, Lt,
                         lambda Cc, k=k: Cc[k, 0] * Cc[k, t2] - Cc[k, t1] ** 2)
        vR, eR, _ = jack(acc, cnt, N, Lt,
                         lambda Cc, k=k: Cc[k, 0] * Cc[k, t2] / Cc[k, t1] ** 2)
        rhs = v1 * v1 / v0 if v0 != 0 else np.nan
        K_ = vD / (v1 * v1) if v1 != 0 else np.nan
        sig = vD / eD if (eD > 0 and np.isfinite(eD)) else np.nan
        rows.append(dict(k=k, C0=v0, Ct1=v1, Ct2=v2, rhs=rhs, D=vD, sD=eD,
                         K=K_, R=vR, sR=eR, sig=sig,
                         ok=bool(vD >= 0), viol_sig=float(-sig) if vD < 0 else 0.0))
    # which stream is the canonical target operator?
    kc = 0
    if lev is not None and obs is not None and len(shape_mid) >= 2:
        kc = int(lev) * shape_mid[1] + int(obs)
    P(f"  canonical stream k={kc} (lev,obs index from shape_mid={shape_mid})")
    P(f"  {'k':>3} {'C(0)':>12} {'C(L/2)':>12} {'C(L)':>12} {'C(L/2)^2/C(0)':>15}"
      f" {'D':>12} {'+-sD':>10} {'K':>10} {'R':>12} {'verdict':>8}")
    for r in rows:
        mark = '*' if r['k'] == kc else ' '
        P(f" {mark}{r['k']:>2} {r['C0']:12.5e} {r['Ct1']:12.5e} {r['Ct2']:12.5e}"
          f" {r['rhs']:15.5e} {r['D']:12.4e} {r['sD']:10.2e} {r['K']:10.4f}"
          f" {r['R']:12.4f} {'PASS' if r['ok'] else 'FAIL':>8}")
    npass = sum(1 for r in rows if r['ok'])
    worst = min(rows, key=lambda r: r['D'] / r['sD'] if r['sD'] > 0 else np.inf)
    P(f"  streams: PASS {npass}/{K} ; worst D/sD = {worst['D']/worst['sD']:+.2f} (k={worst['k']})")
    if kc < len(rows):
        r = rows[kc]
        P(f"  * canonical: C(0)={r['C0']:.5e}  C(L/2)={r['Ct1']:.5e}  C(L)={r['Ct2']:.5e}")
        P(f"    rhs=C(L/2)^2/C(0)={r['rhs']:.5e}   D={r['D']:+.4e} +- {r['sD']:.2e}"
          f"   R={r['R']:.4f} +- {r['sR']:.4f}   D/sD={r['sig']:+.2f}")
        P(f"    P-readout C(L)/C(0) = {r['Ct2']/r['C0']:.4e}")
    ROWS.append(dict(file=path, Ls=Ls, Lt=Lt, N=N, beta=beta, seed=seed, K=K, tier=T,
                     testable=True, npass=npass,
                     canon=rows[kc] if kc < len(rows) else None,
                     worst=dict(D=worst['D'], sD=worst['sD'], k=worst['k']),
                     allrows=[{kk: (None if isinstance(vv, float) and not np.isfinite(vv) else vv)
                               for kk, vv in r.items()} for r in rows]))
    return rows


def selfcheck():
    """the FFT/block machinery must reproduce run3_apebase.correlators bit-for-bit."""
    P("=" * 118)
    P("[0] self-check: chunked FFT estimator vs run3_apebase.correlators (same file)")
    P("=" * 118)
    sys.path.insert(0, W7)
    from run3_apebase import correlators
    f = os.path.join(W7, 'gh2_b2.5_Ls6_Lt12.npz')
    d = np.load(f)
    S = d['S'][:, 0, :]
    Lt = int(d['Lt'])
    Cref = correlators(S, Lt)
    del S
    acc, cnt, N, K, Lt2, _, dd = wrapped_correlator(f)
    with dd:
        pass
    Cmine = acc.sum(axis=0)[0] / (N * Lt2)
    err = np.max(np.abs(Cmine - Cref) / np.maximum(np.abs(Cref), 1e-300))
    P(f"  max rel diff over t=0..{Lt-1}: {err:.3e}   (must be < 1e-12)")
    assert err < 1e-10, err
    P("  OK")


CORE = [
    ('gh2_b2.5_Ls6_Lt20.npz', None, None, None),
    ('gh2_b2.5_Ls8_Lt20.npz', None, None, None),
    ('gh2_b2.5_Ls6_Lt12.npz', None, None, None),
    ('gh2_b2.5_Ls8_Lt12.npz', None, None, None),
    ('gh3_L6.npz', 0, 0, None),
    ('gh3_L8.npz', 0, 0, None),
    ('run1_b25_l6.npz', None, None, None),
    ('run1_b2.5_Ls8Lt12.npz', None, None, None),
    ('run3_L8_b2.5_N12k.npz', None, None, None),
    ('run3_bench.npz', None, None, None),
    ('gp3_b25_Ls8_Lt16_N24k.npz', None, None, None),
    ('gp3_b3_Ls8_Lt16_N100k.npz', None, None, None),
    ('gw7_b2.59_Ls10.npz', None, None, None),
    ('gw7_b2.66_Ls12.npz', None, None, None),
    ('b01_b2.2.npz', None, None, 8),
    ('b01_b2.75.npz', None, None, 8),
    ('b01_b3.npz', None, None, 8),
    ('b01_b4.npz', None, None, 8),
    ('b01_b6.npz', None, None, 8),
    ('b01_b10.npz', None, None, 8),
]
PROBE = [('gh3_probe.npz', 0, 0, None)]
X2 = [(f'gh3_x2_{i}.npz', 0, 0, None) for i in range(1, 9)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', default='core', choices=['core', 'probe', 'x2', 'one'])
    ap.add_argument('--file', default=None)
    a = ap.parse_args()
    P("HR4-PTEST: two-point complete-monotonicity test on existing npz  (zero MC)")
    if a.set == 'core':
        selfcheck()
        seq = CORE
    elif a.set == 'probe':
        seq = PROBE
    elif a.set == 'x2':
        seq = X2
    else:
        seq = [(a.file, 0, 0, None)]
    for (path, lev, obs, lso) in seq:
        do_file(path, None, lev, obs, Ls_override=lso)
    tag = a.set
    json.dump(ROWS, open(os.path.join(HERE, f'hr4_ptest_{tag}.json'), 'w'), indent=1, default=str)
    logp = os.path.join(HERE, f'hr4_ptest_{tag}.log')
    mode = 'a' if (a.set != 'core' and os.path.exists(logp)) else 'w'
    open(logp, mode).write("\n".join(OUT) + "\n")
    P(f"\n[ok] wrote hr4_ptest_{tag}.json / .log")


if __name__ == '__main__':
    main()
