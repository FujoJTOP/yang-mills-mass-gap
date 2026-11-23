"""Consolidate tm_gap_gevp.py extractions into a tidy table + CSV."""
import glob
import numpy as np
import tm_gap_gevp as G


def pick(ts, full, errs):
    """return list of (t,val,err) finite."""
    return [(t, float(full[i]), float(errs[i])) for i, t in enumerate(ts)
            if np.isfinite(full[i]) and errs[i] > 0]


def main():
    rows = []
    print(f"{'point':<12}{'L':>3}{'beta':>6} | "
          f"{'GEVP[4,8] t=L/2':>22} | {'GEVP[4,8] best-S/N<=L/2':>26} | {'single[8] t=L/2':>18}")
    print("-" * 100)
    for f in sorted(glob.glob('tm_gap_data/L*_b*.npz')):
        tag = f.split('\\')[-1].split('/')[-1].replace('.npz', '')
        r48 = G.analyse(f, t0=1, basis=[4, 8], verbose=False)
        r08 = G.analyse(f, t0=1, basis=[0, 4, 8], verbose=False)
        r8 = G.analyse(f, t0=1, basis=[8], verbose=False)
        L = r48['L']; beta = r48['beta']
        # t=L/2 value from GEVP [4,8]
        tL = L // 2
        iL = r48['ts'].index(tL)
        vL, eL = r48['full'][iL], r48['errs'][iL]
        # best S/N over t<=L/2 for [4,8]
        cand = pick(r48['ts'], r48['full'], r48['errs'])
        best = min(cand, key=lambda x: x[2] / max(abs(x[1]), 1e-9)) if cand else (None, np.nan, np.nan)
        # single [8] at t=L/2
        iL8 = r8['ts'].index(tL)
        v8, e8 = r8['full'][iL8], r8['errs'][iL8]
        print(f"{tag:<12}{L:>3}{beta:>6.1f} | "
              f"{vL:>9.3f} +- {eL:<8.3f} | "
              f"t={best[0]} {best[1]:>6.3f} +- {best[2]:<8.3f} | "
              f"{v8:>9.3f} +- {e8:<6.3f}")
        rows.append(dict(tag=tag, L=L, beta=beta,
                         gevp48_tL_t=tL, gevp48_tL=vL, gevp48_tL_err=eL,
                         gevp48_best_t=best[0], gevp48_best=best[1], gevp48_best_err=best[2],
                         single8_tL=v8, single8_tL_err=e8))
    # write CSV
    with open('tm_gap_results.csv', 'w') as fh:
        fh.write("tag,L,beta,gevp48_tL_t,gevp48_tL,gevp48_tL_err,"
                 "gevp48_best_t,gevp48_best,gevp48_best_err,single8_tL,single8_tL_err\n")
        for r in rows:
            fh.write(f"{r['tag']},{r['L']},{r['beta']},{r['gevp48_tL_t']},"
                     f"{r['gevp48_tL']:.4f},{r['gevp48_tL_err']:.4f},"
                     f"{r['gevp48_best_t']},{r['gevp48_best']:.4f},{r['gevp48_best_err']:.4f},"
                     f"{r['single8_tL']:.4f},{r['single8_tL_err']:.4f}\n")
    print("\nwrote tm_gap_results.csv")


if __name__ == '__main__':
    main()
