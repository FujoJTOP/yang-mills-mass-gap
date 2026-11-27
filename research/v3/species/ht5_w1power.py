"""HT5-W1POWER -- the L-power of w1^{YM} = f1^2/M0 on the 0++ carrier.

ZERO new MC.  Reads archived S(t) arrays / json records only.
Estimator (unchanged, RUN1/RUN3/GS1/HQ7 convention):
    C(t) = correlated slice operator, zero momentum
    w1   = A_at(C, Lt, m=1.546) / C(0)
    errors = 20-block delete-one jackknife over configurations.

Route A (primary): beta=2.5, Ls in {6,8}, Lt in {12,20}, N=48000  (gh2 A npz + gh2 B json)
Route B (independent): beta=2.5, Ls in {6,8}, Lt=12, APE levels 0..5 (gh3_L6/L8 npz)

EXCLUSIONS (see HT5-W1POWER.md sec.1): every stored w1 field with w1>1 or
w1_err<1e-12 (the registered correlators() axis bug) is screened out before use;
b01_* w1 columns are not used at all.
"""
import os
import sys
import json
import math

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.join(HERE, 'w7')
sys.path.insert(0, W7)
from run3_apebase import correlators, A_at, fit_A  # noqa: E402

M = 1.546
NB = 20
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def jk_from_C(Cs):
    def jk(fn):
        v = np.array([fn(Cs[b]) for b in range(len(Cs))], float)
        m = float(np.nanmean(v))
        e = float(np.sqrt((NB - 1) / NB * np.nansum((v - m) ** 2)))
        return m, e
    return jk


def cells_from_npz(path):
    d = np.load(path)
    Lt = int(d['Lt'])
    S = d['S']
    Sm = S[:, 0, :] if S.ndim == 3 else S
    C = correlators(Sm, Lt)
    N = Sm.shape[0]
    idx = np.array_split(np.arange(N), NB)
    Cs = [correlators(np.concatenate([Sm[i] for j, i in enumerate(idx) if j != b]), Lt)
          for b in range(NB)]
    jk = jk_from_C(Cs)
    w1, w1e = jk(lambda c: A_at(c, Lt, M) / c[0])
    wf, wfe = jk(lambda c: fit_A(c, Lt, ts=np.arange(1, min(7, Lt)))[1] / c[0])
    C0e = float(np.sqrt((NB - 1) / NB * np.nansum(
        (np.array([c[0] for c in Cs]) - C[0]) ** 2)))
    return dict(src=os.path.basename(path), beta=float(d['beta']), Ls=int(d['Ls']), Lt=Lt,
                N=N, C0=float(C[0]), C0e=C0e, w1=w1, w1e=w1e, w1f=wf, w1fe=wfe)


def cells_gh3(path, lev):
    d = np.load(path, mmap_mode='r')
    Lt = int(d['Lt'])
    S = np.asarray(d['S'][:, lev, 0, :])     # obs index 0 = p=0 carrier (gh3_analyze.py:78)
    C = correlators(S, Lt)
    N = S.shape[0]
    idx = np.array_split(np.arange(N), NB)
    Cs = [correlators(np.concatenate([S[i] for j, i in enumerate(idx) if j != b]), Lt)
          for b in range(NB)]
    jk = jk_from_C(Cs)
    w1, w1e = jk(lambda c: A_at(c, Lt, M) / c[0])
    wf, wfe = jk(lambda c: fit_A(c, Lt, ts=np.arange(1, min(7, Lt)))[1] / c[0])
    return dict(src=os.path.basename(path) + f':lev{lev}', beta=float(d['beta']),
                Ls=int(d['Ls']), Lt=Lt, N=N, C0=float(C[0]),
                w1=w1, w1e=w1e, w1f=wf, w1fe=wfe)


def comb(a, b):
    wa, wb = 1 / a['w1e'] ** 2, 1 / b['w1e'] ** 2
    m = (wa * a['w1'] + wb * b['w1']) / (wa + wb)
    return m, 1 / math.sqrt(wa + wb)


def q_of(w6, e6, w8, e8):
    r = w8 / w6
    ln = math.log(r) / math.log(8 / 6)
    se = math.hypot(e6 / w6, e8 / w8) / math.log(8 / 6)
    return r, math.hypot(e6 / w6, e8 / w8), ln, se


def main():
    P("=" * 100)
    P("HT5  w1^{YM} = f1^2/M0  on the 0++ carrier : its L-power  (zero new MC)")
    P(f"estimator A_at(C,Lt,m={M})/C(0), t=1..3 window, {NB}-block delete-one jackknife")
    P("=" * 100)

    # ---------------- screening / exclusion ledger ----------------
    P("")
    P("[0] SCREEN of every stored w1 field in research/v3/species/w7/*.json")
    bad = []
    for f in sorted(os.listdir(W7)):
        if not f.endswith('.json'):
            continue
        try:
            d = json.load(open(os.path.join(W7, f)))
        except Exception:
            continue
        if not isinstance(d, list):
            continue
        for r in d:
            if not isinstance(r, dict) or 'w1' not in r:
                continue
            w = r.get('w1')
            e = r.get('w1_err', r.get('w1e'))
            tag = None
            if w is not None and abs(w) > 1.0:
                tag = 'w1>1'
            if e is not None and abs(e) < 1e-12:
                tag = (tag + '+') if tag else ''
                tag += 'w1_err<1e-12'
            if tag:
                bad.append((f, r.get('beta'), r.get('Ls'), r.get('Lt'), w, e, tag))
    for b in bad:
        P(f"    DROP {b[0]:22s} beta={b[1]} Ls={b[2]} Lt={b[3]} w1={b[4]:.5f} "
          f"w1_err={b[5]:.2e}  [{b[6]}]")
    P(f"    => {len(bad)} record(s) dropped; every b01_* w1 column is excluded wholesale")
    P("       (HQ7-TWOPILLAR.md:22 / GP3-PLATEAU.md:154 / GM3-NUMERICS.md:119).")

    # ---------------- route A ----------------
    P("")
    P("[1] ROUTE A  beta=2.5, N=48000/cell  (A = archived npz, B = gh2_ltB.json)")
    A = {}
    for Ls in (6, 8):
        for Lt in (12, 20):
            c = cells_from_npz(os.path.join(W7, f'gh2_b2.5_Ls{Ls}_Lt{Lt}.npz'))
            A[(Ls, Lt)] = c
            P(f"    A Ls={Ls} Lt={Lt:2d} N={c['N']:6d} C(0)={c['C0']:8.4f}+-{c['C0e']:.4f} "
              f"w1={c['w1']:.5f}+-{c['w1e']:.5f}  w1_fit={c['w1f']:.5f}+-{c['w1fe']:.5f}")
    P("    NOTE: the four archived npz ARE the gh2_ltB replica (C(0) digit-identical to")
    P("          gh2_ltB.json) -> the second replica must be taken from gh2_ltA.json, not ltB.")
    B = {(r['Ls'], r['Lt']): r for r in json.load(open(os.path.join(W7, 'gh2_ltA.json')))}
    for k in sorted(A):
        P(f"    A-vs-npz check Ls={k[0]} Lt={k[1]:2d}: npz w1={A[k]['w1']:.5f} "
          f"ltB json w1={[r['w1'] for r in json.load(open(os.path.join(W7,'gh2_ltB.json'))) if (r['Ls'],r['Lt'])==k][0]:.5f}"
          f"  (identical => npz == ltB)")
    comb_c = {}
    for k in sorted(A):
        w, e = comb(A[k], dict(w1e=B[k]['w1_err'], w1=B[k]['w1']))
        comb_c[k] = dict(w1=w, w1e=e, C0=(A[k]['C0'] + B[k]['C0']) / 2)
        P(f"    AB Ls={k[0]} Lt={k[1]:2d} N=48000 : w1={w:.5f}+-{e:.5f}")
    P("")
    P(f"    cross-check vs GH2-LTART.md:197-200 (0.3354/0.3491/0.3558/0.3568):")
    P("      " + "  ".join(f"{comb_c[k]['w1']:.4f}" for k in sorted(comb_c)
                           if k[1] == 12) + "   (Lt=12) ;  "
      + "  ".join(f"{comb_c[k]['w1']:.4f}" for k in sorted(comb_c) if k[1] == 20) + "  (Lt=20)")

    # ---------------- route B ----------------
    P("")
    P("[2] ROUTE B  independent runs gh3 (beta=2.5, Lt=12), APE levels 0..5")
    G = {}
    for lev in range(6):
        a = cells_gh3(os.path.join(W7, 'gh3_L6.npz'), lev)
        b = cells_gh3(os.path.join(W7, 'gh3_L8.npz'), lev)
        G[lev] = (a, b)
        flag = '' if (abs(a['w1']) <= 1 and abs(b['w1']) <= 1) else '  [OUT OF RANGE w1>1 -> dropped]'
        P(f"    lev{lev}: Ls=6 N={a['N']:6d} C0={a['C0']:7.4f} w1={a['w1']:.5f}+-{a['w1e']:.5f} | "
          f"Ls=8 N={b['N']:6d} C0={b['C0']:7.4f} w1={b['w1']:.5f}+-{b['w1e']:.5f}{flag}")

    # ---------------- power ----------------
    P("")
    P("[3] L-POWER  q = dln(w1)/dln(L)   over Ls = 6 -> 8   (ln(8/6) = 0.287682)")
    rows = []
    for Lt in (12, 20):
        w6, e6 = comb_c[(6, Lt)]['w1'], comb_c[(6, Lt)]['w1e']
        w8, e8 = comb_c[(8, Lt)]['w1'], comb_c[(8, Lt)]['w1e']
        r, dr, q, dq = q_of(w6, e6, w8, e8)
        rows.append((f'A Lt={Lt}', q, dq, r))
        P(f"    A Lt={Lt:2d} : w1(8)/w1(6) = {r:.4f}+-{dr:.4f}  =>  q = {q:+.4f}+-{dq:.4f}")
    wa = 1 / rows[0][2] ** 2
    wb = 1 / rows[1][2] ** 2
    qc = (wa * rows[0][1] + wb * rows[1][1]) / (wa + wb)
    eq = 1 / math.sqrt(wa + wb)
    chi2 = wa * (rows[0][1] - qc) ** 2 + wb * (rows[1][1] - qc) ** 2
    P(f"    A combined (Lt=12,20 independent) : q = {qc:+.4f}+-{eq:.4f}   "
      f"chi2/1dof = {chi2:.2f}")
    P(f"      2 sigma: [{qc - 2 * eq:+.4f}, {qc + 2 * eq:+.4f}]")
    for tgt, name in ((-3.0, 'L^-3 (vacuous)'), (0.0, 'L^0 (content)'),
                      (-1.0, 'L^-1 (free-field extensive)'), (3.0, 'L^+3')):
        P(f"      vs {name:28s}: {(qc - tgt) / eq:+.1f} sigma")
    # weighted-across-Lt version (independent samples)
    w6 = (comb_c[(6, 12)]['w1'] / comb_c[(6, 12)]['w1e'] ** 2
          + comb_c[(6, 20)]['w1'] / comb_c[(6, 20)]['w1e'] ** 2)
    d6 = (1 / comb_c[(6, 12)]['w1e'] ** 2 + 1 / comb_c[(6, 20)]['w1e'] ** 2)
    w8 = (comb_c[(8, 12)]['w1'] / comb_c[(8, 12)]['w1e'] ** 2
          + comb_c[(8, 20)]['w1'] / comb_c[(8, 20)]['w1e'] ** 2)
    d8 = (1 / comb_c[(8, 12)]['w1e'] ** 2 + 1 / comb_c[(8, 20)]['w1e'] ** 2)
    r, dr, q2, eq2 = q_of(w6 / d6, 1 / math.sqrt(d6), w8 / d8, 1 / math.sqrt(d8))
    P(f"    A all-Lt pooled : w1(6)={w6/d6:.5f}+-{1/math.sqrt(d6):.5f} "
      f"w1(8)={w8/d8:.5f}+-{1/math.sqrt(d8):.5f} => q = {q2:+.4f}+-{eq2:.4f}")
    P(f"      2 sigma: [{q2 - 2 * eq2:+.4f}, {q2 + 2 * eq2:+.4f}]")
    P(f"      vs L^-3: {(q2 + 3) / eq2:+.1f} sigma ; vs L^0: {(q2 - 0) / eq2:+.1f} sigma ;"
      f" vs L^-1: {(q2 + 1) / eq2:+.1f} sigma")

    P("")
    P("[4] ROUTE B power, per APE level (each level = one 0++ carrier)")
    qs = []
    for lev in range(6):
        a, b = G[lev]
        r, dr, q, dq = q_of(a['w1'], a['w1e'], b['w1'], b['w1e'])
        if abs(a['w1']) > 1 or abs(b['w1']) > 1:
            P(f"    lev{lev}: DROPPED (w1 out of [0,1])")
            continue
        qs.append((q, dq))
        P(f"    lev{lev}: ratio {r:.4f}+-{dr:.4f}  q = {q:+.4f}+-{dq:.4f}")
    wsum = sum(1 / d ** 2 for _, d in qs)
    qm = sum(q / d ** 2 for q, d in qs) / wsum
    P(f"    B pooled over 6 levels : q = {qm:+.4f}+-{1/math.sqrt(wsum):.4f}")
    raw = qs[0]
    P(f"    B lev0 only (raw op)   : q = {raw[0]:+.4f}+-{raw[1]:.4f}")

    P("")
    P("[5] C(0) power on the same cells (denominator M0):")
    for Lt in (12, 20):
        c6, c8 = comb_c[(6, Lt)]['C0'], comb_c[(8, Lt)]['C0']
        P(f"    Lt={Lt:2d}: C(0)(6)={c6:.4f} C(0)(8)={c8:.4f}  ratio={c8/c6:.4f} "
          f"(L^3={2.3704:.4f})  => p_C0={math.log(c8/c6)/math.log(8/6):+.3f}")
    P("    => M0 ~ L^3 (volume) while w1 ~ L^0  =>  f1^2 = w1*M0 ~ L^{3+q} ~ V")

    P("")
    P("[6] SYSTEMATIC: Lt truncation shifts q")
    dd = 0.0137 - 0.0010   # GH2-LTART.md:242-244  delta(6) - delta(8)
    wbar = 0.345
    P(f"    delta(6)-delta(8) = +0.0127+-0.0099 (GH2-LTART.md:244) => dq = "
      f"{dd / (wbar * math.log(8 / 6)):.3f} at 1 sigma of that difference")
    P(f"    q(Lt->inf) = q(20) = {rows[1][1]:+.4f}+-{rows[1][2]:.4f}  (no correction needed)")
    P("=" * 100)
    open(os.path.join(HERE, 'ht5_w1power.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    os.chdir(HERE)
    main()
