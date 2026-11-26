"""GH2 pooling + verdict arithmetic.  Reads the per-job records that gh2_ltscan.py wrote
(gh2_lt*.json -- one entry per (beta, Ls, Lt, seed) job), pools independent replicate jobs by
inverse variance, and prints the verdict arithmetic used in GH2-LTART.md.

Why JSON and not the npz: jobs A and B share the same output npz name (gh2_b<beta>_Ls<Ls>_Lt<Lt>.npz),
so job B overwrote job A's arrays.  The per-job JSON records survive and carry every number quoted
here (they are the same numbers gh2_ltscan.py printed).  The npz of the last job are kept as raw
arrays.  Recorded as a bookkeeping defect in GH2-LTART.md S4.

Run:  py -3.12 gh2_pool.py
"""
import os
import sys
import json
import glob

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M = 1.546
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def pool(rows):
    """Inverse-variance pooling of independent replicate jobs at the same (beta, Ls, Lt)."""
    v = np.array([r['w1'] for r in rows])
    w = np.array([1 / r['w1_err'] ** 2 for r in rows])
    mw = float(np.sum(w * v) / np.sum(w))
    ew = float(1 / np.sqrt(np.sum(w)))
    chi2 = float(np.sum(w * (v - mw) ** 2))
    vR = np.array([r['Rbar'] for r in rows])
    wR = np.array([1 / r['Rbar_err'] ** 2 for r in rows])
    mR = float(np.sum(wR * vR) / np.sum(wR))
    eR = float(1 / np.sqrt(np.sum(wR)))
    chi2R = float(np.sum(wR * (vR - mR) ** 2))
    vk = np.array([r['kappa'] for r in rows])
    wk = np.array([1 / r['kappa_err'] ** 2 for r in rows])
    mk = float(np.sum(wk * vk) / np.sum(wk))
    ek = float(1 / np.sqrt(np.sum(wk)))
    return dict(w1=mw, w1e=ew, chi2=chi2, Rbar=mR, Rbar_e=eR, chi2R=chi2R,
                kappa=mk, kappa_e=ek, nj=len(rows), N=int(sum(r['N'] for r in rows)),
                C0=float(np.mean([r['C0'] for r in rows])),
                plaq=float(np.mean([r['plaq'] for r in rows])),
                C1o0=float(np.mean([r['C1o0'] for r in rows])))


def main():
    files = sorted(glob.glob(os.path.join(HERE, 'gh2_lt*.json')))
    rows = []
    for f in files:
        for r in json.load(open(f)):
            r['job'] = os.path.basename(f)
            rows.append(r)
    P("=" * 112)
    P("GH2 POOL  --  L_t artefact test, beta = 2.5, SU(2) Wilson, raw operator, m = 1.546 pinned")
    P("estimator = run3_apebase.correlators + A_at, window t=1,2,3, 20-block delete-one jackknife")
    P(f"job files = {[os.path.basename(f) for f in files]}  ;  {len(rows)} point-jobs")
    P("=" * 112)

    P("")
    P("--- per-job ---")
    P(f"{'job':>9} {'Ls':>3} {'Lt':>3} {'N':>6} {'<plaq>':>9} {'C(0)':>9} {'C1/C0':>8} {'w1':>8} "
      f"{'+-':>7} {'Rbar':>7} {'+-':>6} {'kap':>6} {'S/N2':>5} {'s':>6}")
    for r in sorted(rows, key=lambda r: (r['Ls'], r['Lt'], r['seed'])):
        P(f"{r['job'].replace('gh2_','').replace('.json',''):>9} {r['Ls']:>3d} {r['Lt']:>3d} "
          f"{r['N']:>6d} {r['plaq']:>9.5f} {r['C0']:>9.4f} {r['C1o0']:>8.5f} {r['w1']:>8.4f} "
          f"{r['w1_err']:>7.4f} {r['Rbar']:>7.4f} {r['Rbar_err']:>6.4f} {r['kappa']:>6.3f} "
          f"{r['snr2']:>5.1f} {r['sec']:>6.1f}")

    P("")
    P("--- C(t)/C(0) with jackknife errors, t = 0..4 ---")
    for r in sorted(rows, key=lambda r: (r['Ls'], r['Lt'], r['seed'])):
        P(f"  {r['job']} Ls={r['Ls']} Lt={r['Lt']}: " + "  ".join(
            "1.00000" if t == 0 else f"{r['r'][str(t)]:.5f}+-{r['r_e'][str(t)]:.5f}"
            for t in range(0, 5)))
        P(f"      m_eff(t) = ln C(t)/C(t+1), t=0..4 : " + "  ".join(
            f"{r['meff'][str(t)]:+.3f}" for t in range(0, 5)))
        P(f"      m_fit(t=1..6) = {r['m_fit']:.3f} ; w1_nowrap = {r['w1_nowrap']:.5f} "
          f"(wrap term's weight in the fit = {100*(r['w1_nowrap']/r['w1']-1):+.3f}%)")

    keys = sorted(set((r['beta'], r['Ls'], r['Lt']) for r in rows))
    pooled = {}
    P("")
    P("--- pooled by (Ls, Lt), inverse variance ---")
    P(f"{'beta':>5} {'Ls':>3} {'Lt':>3} {'jobs':>4} {'N_tot':>6} {'<plaq>':>9} {'C(0)':>9} "
      f"{'C1/C0':>8} {'w1':>8} {'+-':>7} {'Rbar':>7} {'+-':>6} {'kap':>6} {'chi2_w':>7}")
    for k in keys:
        sub = [r for r in rows if (r['beta'], r['Ls'], r['Lt']) == k]
        pp = pool(sub)
        pooled[k] = pp
        P(f"{k[0]:>5.3f} {k[1]:>3d} {k[2]:>3d} {pp['nj']:>4d} {pp['N']:>6d} {pp['plaq']:>9.5f} "
          f"{pp['C0']:>9.4f} {pp['C1o0']:>8.5f} {pp['w1']:>8.4f} {pp['w1e']:>7.4f} "
          f"{pp['Rbar']:>7.4f} {pp['Rbar_e']:>6.4f} {pp['kappa']:>6.3f} {pp['chi2']:>7.2f}")

    def g(k):
        return pooled[k] if k in pooled else None

    P("")
    P("=" * 112)
    P("VERDICT ARITHMETIC")
    P("=" * 112)
    P("")
    P("(1) L_s 6->8 difference at each FIXED L_t")
    for Lt in sorted(set(k[2] for k in keys)):
        A, B = g((2.5, 6, Lt)), g((2.5, 8, Lt))
        if A and B:
            d, e = A['w1'] - B['w1'], np.hypot(A['w1e'], B['w1e'])
            dR, eR = B['Rbar'] - A['Rbar'], np.hypot(A['Rbar_e'], B['Rbar_e'])
            P(f"  Lt={Lt:2d}: w1(6)-w1(8) = {d:+.4f} +- {e:.4f}  ({100*d/A['w1']:+.2f}% of w1(6), "
              f"{abs(d)/e:.2f} sigma)   [N={A['N']}+{B['N']}]")
            P(f"          Rbar(8)-Rbar(6) = {dR:+.4f} +- {eR:.4f}  ({100*dR/A['Rbar']:+.2f}% of Rbar, "
              f"{abs(dR)/eR:.2f} sigma)")
    P("")
    P("(2) L_t drift at each FIXED L_s   <-- this IS the truncation the paper is accused of")
    for Ls in sorted(set(k[1] for k in keys)):
        ks = sorted([k for k in keys if k[1] == Ls], key=lambda k: k[2])
        for i in range(len(ks) - 1):
            A, B = pooled[ks[i]], pooled[ks[i + 1]]
            d, e = B['w1'] - A['w1'], np.hypot(A['w1e'], B['w1e'])
            dR, eR = B['Rbar'] - A['Rbar'], np.hypot(A['Rbar_e'], B['Rbar_e'])
            P(f"  Ls={Ls}: Lt {ks[i][2]:2d}->{ks[i+1][2]:2d}  dw1 = {d:+.4f} +- {e:.4f} "
              f"({100*d/A['w1']:+.2f}%, {abs(d)/e:.2f} sigma)")
            P(f"           dRbar = {dR:+.4f} +- {eR:.4f} ({100*dR/A['Rbar']:+.2f}%, "
              f"{abs(dR)/eR:.2f} sigma)")
    P("")
    P("(3) differential truncation, the only thing that can fake an L_s dependence")
    d6 = g((2.5, 6, 20))['w1'] - g((2.5, 6, 12))['w1']
    e6 = np.hypot(g((2.5, 6, 20))['w1e'], g((2.5, 6, 12))['w1e'])
    d8 = g((2.5, 8, 20))['w1'] - g((2.5, 8, 12))['w1']
    e8 = np.hypot(g((2.5, 8, 20))['w1e'], g((2.5, 8, 12))['w1e'])
    P(f"  delta(6) = w1(6,20)-w1(6,12) = {d6:+.4f} +- {e6:.4f}")
    P(f"  delta(8) = w1(8,20)-w1(8,12) = {d8:+.4f} +- {e8:.4f}")
    P(f"  delta(6)-delta(8) = {d6-d8:+.4f} +- {np.hypot(e6,e8):.4f}  "
      f"= {100*(d6-d8)/g((2.5,6,12))['w1']:+.2f}% of w1(6,12) "
      f"({abs(d6-d8)/np.hypot(e6,e8):.2f} sigma)")
    P("")
    P("(4) archived anchors (paper rows) vs this work")
    arch = {(6, 12): (0.37292, 0.01753, 2.5164, 0.0460, 6000, 'RUN1.md:47'),
            (8, 12): (0.34273, 0.01189, 2.6075, 0.0355, 12000, 'RUN3.md:51')}
    for (Ls, Lt), (w, we, R, Re, N, src) in arch.items():
        A = g((2.5, Ls, Lt))
        ds = A['w1'] - w
        es = np.hypot(we, A['w1e'])
        P(f"  Ls={Ls} Lt={Lt} ({src}): archive w1={w:.4f}+-{we:.4f} (N={N}) vs "
          f"{A['w1']:.4f}+-{A['w1e']:.4f} (N={A['N']}) -> shift {ds:+.4f} +- {es:.4f} "
          f"({abs(ds)/es:.2f} sigma)")
        dRs = A['Rbar'] - R
        eRs = np.hypot(Re, A['Rbar_e'])
        P(f"          archive Rbar={R:.4f}+-{Re:.4f} vs {A['Rbar']:.4f}+-{A['Rbar_e']:.4f} -> "
          f"shift {dRs:+.4f} +- {eRs:.4f} ({abs(dRs)/eRs:.2f} sigma)")
    P("")
    P("(5) the paper's number vs this work, in the paper's own units")
    paper_d = 0.37292 - 0.34273
    paper_e = np.hypot(0.01753, 0.01189)
    for Lt in sorted(set(k[2] for k in keys)):
        A, B = g((2.5, 6, Lt)), g((2.5, 8, Lt))
        d, e = A['w1'] - B['w1'], np.hypot(A['w1e'], B['w1e'])
        dd, ee = paper_d - d, np.hypot(paper_e, e)
        P(f"  Lt={Lt:2d}: paper {paper_d:+.4f}+-{paper_e:.4f} vs this work {d:+.4f}+-{e:.4f} "
          f"-> difference {dd:+.4f} +- {ee:.4f} ({abs(dd)/ee:.2f} sigma)")
    P("")
    P("(6) L_t-independence of C(0)  (a light mode would move the denominator)")
    for Ls in sorted(set(k[1] for k in keys)):
        vals = [(k[2], pooled[k]['C0']) for k in keys if k[1] == Ls]
        c0s = [v[1] for v in vals]
        P(f"  Ls={Ls}: C(0) at Lt={[v[0] for v in vals]} -> {['%.4f' % c for c in c0s]} "
          f"spread {100*(max(c0s)-min(c0s))/np.mean(c0s):.2f}%")
    P("")
    P(f"MC wall clock total (both jobs, from their own [wall] lines) = 241 + 241 = 482 s "
      f"(cap 600 s)")

    json.dump(dict(per_job=rows, pooled={str(k): v for k, v in pooled.items()}),
              open(os.path.join(HERE, 'gh2_pool.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'gh2_pool.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
