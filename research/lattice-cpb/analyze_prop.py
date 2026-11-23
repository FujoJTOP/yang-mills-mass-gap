"""Final analysis of gluon_prop_results.json -> numbers for the report."""
import json
import numpy as np
from scipy.optimize import least_squares


def linfit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return a, b   # 1/D = a k2 + b ; m2=b/a, Z=1/a, D0=1/b


def fit_GS(x, y):
    def res(p):
        Z, m1s, m2s, m3s = p
        return Z * (x + m1s) / (x * x + m2s * x + m3s ** 2) - y
    best = None
    for m1g in (0.05, 0.5, 2.0):
        try:
            r = least_squares(res, [max(y[0] * (x[0] + m1g), 1e-3), m1g, 1.0, 0.5],
                              bounds=(1e-9, np.inf), max_nfev=20000)
            if best is None or r.cost < best.cost:
                best = r
        except Exception:
            pass
    Z, m1s, m2s, m3s = best.x
    D0 = Z * m1s / (m3s ** 2)
    return D0, np.sqrt(m1s)


data = json.load(open('D:/Dev/math/research/lattice-cpb/gluon_prop_results.json'))
print("=" * 78)
print("FREE CONTROLS (input: massless 1/k^2)")
for r in data:
    if r['kind'] != 'free':
        continue
    sh = np.array(r['ushells']); D = np.array(r['meanD'])
    print(f"  L={r['L']}:  " + "  ".join(f"k2={s:g}:D={d:.3f}" for s, d in zip(sh[:4], D[:4])))
    print(f"          (D*k^2 = {np.round(D[:4]*sh[:4],3)}  ~ const => massless)")

print("=" * 78)
print("MC RESULTS   (fit 1/D = a k^2 + b over RELIABLE shells k^2>=2)")
hdr = f"{'beta':>5} {'L':>2} {'ncfg':>4} {'plaq':>5} | "
for r in data:
    if r['kind'] != 'mc':
        continue
    sh = np.array(r['ushells']); D = np.array(r['meanD']); Dk = np.array(r['meanDk'])
    # reliable: k2 >= 2, and D>0
    sel = (sh >= 2 - 1e-9) & (D > 1e-6)
    x = sh[sel]; y = 1.0 / D[sel]
    a, b = linfit(x, y)
    m2 = b / a; Z = 1 / a
    D0 = 1 / b if b > 1e-6 else np.inf
    try:
        D0g, m1g = fit_GS(x, y * 0)  # placeholder
    except Exception:
        D0g, m1g = np.nan, np.nan
    try:
        D0g, m1g = fit_GS(x, D[sel])
    except Exception:
        D0g, m1g = np.nan, np.nan
    # D*k^2 diagnostic
    dk2_lo = D[sel][0] * x[0]           # at smallest reliable k2 (=2)
    dk2_hi = D[sel][-1] * x[-1]         # at largest shell
    print(f"\n beta={r['beta']:.1f} L={r['L']} ncfg={r['ncfg']} plaq={r['plaq']:.3f} "
          f"div/rms={r['div_rms_ratio']:.1e}")
    print("   k2      D(k)        D*k2      1/D")
    for s, d in zip(sh, D):
        flag = "  <-- k2=1 unreliable (excluded)" if s < 1.5 else ""
        print(f"  {s:5.2f}  {d: .4e}  {d*s: .4f}   {1/d: .4f}{flag}")
    print(f"   massive fit (k2>=2): m2={m2:+.4f}  m={np.sqrt(max(m2,0)):.4f}  "
          f"Z={Z:.3f}  D0={D0:.3f}  intercept b={b:+.4f}")
    print(f"   D*k2  at k2=2 -> {dk2_lo:.2f},  at k2={x[-1]:.0f} -> {dk2_hi:.2f}   "
          f"({'DECREASING => steeper than 1/k^2' if dk2_hi < dk2_lo else 'increasing => massive pull'})")
    print(f"   Gribov-Stingl: D0={D0g:.2f}  m1={m1g:.3f}")
    # khat-projected for comparison
    selK = (sh >= 2 - 1e-9) & (Dk > 1e-6)
    aK, bK = linfit(sh[selK], 1.0 / Dk[selK])
    print(f"   [khat-proj] m2={bK/aK:+.4f}  D0={(1/bK if bK>1e-6 else float('inf')):.3f}")
print("=" * 78)

# volume trend of D at comparable k2 across L (beta 2.5)
print("VOLUME CHECK at beta=2.5 (D(k) at shared momenta, L=4 vs 6 vs 8):")
for ktarget in (2.0, 4.0, 6.0):
    row = []
    for r in data:
        if r['kind'] != 'mc' or abs(r['beta'] - 2.5) > 1e-6:
            continue
        sh = np.array(r['ushells']); D = np.array(r['meanD'])
        i = np.argmin(np.abs(sh - ktarget))
        if abs(sh[i] - ktarget) < 1e-6:
            row.append(f"L{r['L']}:{D[i]:.3f}")
    print(f"  k2={ktarget}:  " + "   ".join(row))
