"""
SM6 soft-mode lattice scan witness  (SM6: softest-mode conditioned curvature at beta=2.5).

Reads ONLY local campaign per-mode data (no downloads, no raw links exist for SU(2) @ beta~2.5):
  partA_su2_N4000.npz   SU(2) beta=2.5 L=4  48 fixed transverse plane-wave modes (a,H,G)
  partA_su2_N800.npz    SU(2) beta=2.5 L=6  48 fixed transverse plane-wave modes (a,H,G)
  k0_su2_b2.4_L8_n2.npz SU(2) beta=2.4 L=8  12 modes (6 @ k=pi/4, 6 @ k=pi/2), full (a,H,G,P1,P2)

These are FIXED probe directions v (plane waves in the UNFIXED flat chart), NOT a Hessian
eigenbasis, and there are NO stored cross-couplings <v_i, Hess v_j>.  Therefore:

  * We cannot diagonalise the block Hessian to get its true lowest eigenvector.
  * We CAN evaluate, on every stored direction v, the gauge-dependent conditioned curvature
        G''(v;a0*) = -d^2/da^2 log P_v(a)|_{a0*},   a=<v,x>,
    and take the MINIMUM over the sampled set.  By the Rayleigh-quotient inequality that
    sampled minimum is an UPPER bound on the true softest block-Hessian eigenvalue
        lambda_min(block Hess) = min_{|v|=1} v^T (Hess V) v  <=  v^T (Hess V) v  for every sampled v,
    so  min_sampled G''  >=  lambda_min.  A positive sampled minimum therefore does NOT certify a
    positive true softest margin; a NEGATIVE sampled value certifies lambda_min < 0.

Estimators per mode (identical conventions to COMPUTE-GPP-K0-CONTINUUM / h5_head_reanalysis.py):
  a0*    = KDE mode of P_v(a)                (coarse minimiser; symmetric measure -> ~0)
  kappa  = curvature of -log P_v at a0*      (KDE gold standard)
  bracket= E[H|a0*] - Var(G|a0*)             (Wilson action bracket; equipartition ~0)
  Wgeo   = E[Phi''|a0*] - Var(Phi'|a0*)       (only k0; proven >= 0 pointwise)
  X      = -2 Cov(G,Phi'|a0*)                 (only k0; sign-indefinite cross term)
  Gid    = bracket + Wgeo + X                 (exact identity; only k0)
Errors: config bootstrap (B=80) of the per-mode mean and of the per-mode MIN.
NOTE: the min is non-smooth; its bootstrap SE is a spread indicator, not a Gaussian error.

Run:  py -3.12 -I sm6_soft_scan.py
"""
import os, json
import numpy as np
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "unconditional", "w7")
B = 80
RNG = np.random.default_rng(20261124)


def _mode_kde(a, bw=0.20, ngrid=161, span=4.0):
    a = np.asarray(a, float); sd = a.std()
    if sd <= 0:
        return 0.0
    kde = gaussian_kde(a, bw_method=bw)
    xs = np.linspace(-span * sd, span * sd, ngrid)
    return float(xs[int(np.argmax(kde(xs)))])


def _kappa_kde(a, a0, bw=0.20, rfit=0.8, ngrid=41):
    a = np.asarray(a, float); sd = a.std()
    if sd <= 0:
        return 0.0
    kde = gaussian_kde(a, bw_method=bw)
    xs = np.linspace(a0 - rfit * sd, a0 + rfit * sd, ngrid)
    p = np.maximum(kde(xs), 1e-300)
    y = -np.log(p)
    A = np.vstack([np.ones_like(xs), (xs - a0), (xs - a0) ** 2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(2.0 * coef[2])


def _kernel_cond(a, cols, a0, bw=0.30):
    a = np.asarray(a, float); sd = a.std()
    h = bw * sd if sd > 0 else 1.0
    w = np.exp(-0.5 * ((a - a0) / h) ** 2); ws = w.sum()
    if ws <= 0:
        ws = 1.0
    means = {k: float(np.sum(w * np.asarray(x, float)) / ws) for k, x in cols.items()}
    out = {"E_" + k: v for k, v in means.items()}
    for k, x in cols.items():
        for l, y in cols.items():
            if k <= l:
                xm = np.asarray(x, float) - means[k]
                ym = np.asarray(y, float) - means[l]
                out["C_%s_%s" % (k, l)] = float(np.sum(w * xm * ym) / ws)
    return out


def _per_mode(A, H, G, P1=None, P2=None, idx=None, kdebw=0.20, bwc=0.30):
    """Return arrays over modes: a0, kappa, bracket, EH, [Wgeo, X, Gid]."""
    if idx is None:
        idx = np.arange(len(A))
    Nm = A.shape[1]
    a0 = np.empty(Nm); kap = np.empty(Nm); br = np.empty(Nm); EH = np.empty(Nm)
    Wg = np.full(Nm, np.nan); Xc = np.full(Nm, np.nan); Gid = np.full(Nm, np.nan)
    for j in range(Nm):
        a = A[idx, j]
        m0 = _mode_kde(a, bw=kdebw)
        kap[j] = _kappa_kde(a, m0, bw=kdebw)
        cols = {"H": H[idx, j], "G": G[idx, j]}
        if P1 is not None:
            cols["P1"] = P1[idx, j]; cols["P2"] = P2[idx, j]
        c = _kernel_cond(a, cols, m0, bw=bwc)
        EH[j] = c["E_H"]
        br[j] = c["E_H"] - c["C_G_G"]
        a0[j] = m0
        if P1 is not None:
            Wg[j] = c["E_P2"] - c["C_P1_P1"]
            Xc[j] = -2.0 * c["C_G_P1"]
            Gid[j] = br[j] + Wg[j] + Xc[j]
    return dict(a0=a0, kappa=kap, bracket=br, EH=EH, Wgeo=Wg, X=Xc, Gid=Gid)


def _boot_minmean(A, H, G, P1=None, P2=None, key="kappa"):
    N = len(A)
    mins = np.empty(B); means = np.empty(B)
    for b in range(B):
        idx = RNG.integers(0, N, N)
        pm = _per_mode(A, H, G, P1, P2, idx)
        vals = pm[key]
        mins[b] = np.nanmin(vals); means[b] = np.nanmean(vals)
    return float(np.std(mins, ddof=1)), float(np.std(means, ddof=1))


def do_partA(path, label):
    d = np.load(path, allow_pickle=False)
    modes = sorted({k[:-2] for k in d.files if k.endswith("_H")})
    A = np.column_stack([d[m + "_a"] for m in modes])
    H = np.column_stack([d[m + "_H"] for m in modes])
    G = np.column_stack([d[m + "_G"] for m in modes])
    pm = _per_mode(A, H, G)
    N, Nm = A.shape
    res = dict(label=label, group="su2", beta=2.5, N=N, Nm=Nm, has_phi=False, modes=modes)
    for key in ("kappa", "bracket", "EH"):
        v = pm[key]
        res[key] = dict(mean=float(v.mean()), min=float(v.min()), max=float(v.max()),
                        median=float(np.median(v)), argmin=int(np.argmin(v)),
                        argmin_name=modes[int(np.argmin(v))])
    # bootstrap SE of mean/min of kappa
    se_min, se_mean = _boot_minmean(A, H, G, key="kappa")
    res["kappa"]["se_boot_min"] = se_min
    res["kappa"]["se_boot_mean"] = se_mean
    res["a0_mean"] = float(pm["a0"].mean())
    res["inv_var_mean"] = float(np.mean([1.0 / A[:, j].var() for j in range(Nm)]))
    # sorted per-mode kappa (softest first) for the record
    res["kappa_sorted"] = [float(x) for x in np.sort(pm["kappa"])]
    return res, pm


def do_k0(path, label):
    d = np.load(path, allow_pickle=False)
    A, H, G, P1, P2 = d["a"], d["H"], d["G"], d["P1"], d["P2"]
    beta = float(d["beta"]); L = int(d["L"]); kv = np.asarray(d["kvals"])
    names = [str(x) for x in d["names"]]
    pm = _per_mode(A, H, G, P1, P2)
    N, Nm = A.shape
    res = dict(label=label, group="su2", beta=beta, L=L, N=N, Nm=Nm, has_phi=True,
               names=names, kvals=kv.tolist())
    for key in ("kappa", "bracket", "Wgeo", "X", "Gid", "EH"):
        v = pm[key]
        res[key] = dict(mean=float(v.mean()), min=float(v.min()), max=float(v.max()),
                        median=float(np.median(v)), argmin=int(np.argmin(v)),
                        argmin_name=names[int(np.argmin(v))],
                        argmin_k=float(kv[int(np.argmin(v))]))
    for key in ("kappa", "Gid"):
        se_min, se_mean = _boot_minmean(A, H, G, P1, P2, key=key)
        res[key]["se_boot_min"] = se_min
        res[key]["se_boot_mean"] = se_mean
    res["a0_mean"] = float(pm["a0"].mean())
    res["inv_var_mean"] = float(np.mean([1.0 / A[:, j].var() for j in range(Nm)]))
    res["kappa_sorted"] = [float(x) for x in np.sort(pm["kappa"])]
    res["Gid_sorted"] = [float(x) for x in np.sort(pm["Gid"])]
    return res, pm


def show(r):
    print("\n=== %s ===" % r["label"])
    print("  N=%d modes=%d  a0*(mean)=%+.3f  1/Var=%+.3f  %s"
          % (r["N"], r["Nm"], r["a0_mean"], r["inv_var_mean"],
             "FULL DECOMP" if r["has_phi"] else "a,H,G only (no W/X split)"))
    for key in ("EH", "bracket", "Wgeo", "X", "kappa", "Gid"):
        if key not in r or r[key] is None or np.isnan(r[key]["mean"]):
            continue
        x = r[key]
        line = ("  %-8s mean=%+.3f  min=%+.3f (%s)  max=%+.3f  median=%+.3f"
                % (key, x["mean"], x["min"], x["argmin_name"], x["max"], x["median"]))
        if "se_boot_min" in x:
            line += "  | bootSE(min)=%.3f bootSE(mean)=%.3f" % (x["se_boot_min"], x["se_boot_mean"])
        print(line)
    print("  kappa sorted (softest->stiffest):",
          ", ".join("%+.2f" % v for v in r["kappa_sorted"][:12]),
          "..." if len(r["kappa_sorted"]) > 12 else "")


if __name__ == "__main__":
    out = {}
    specs = [
        ("partA_su2_N4000.npz", "SM6  SU(2) beta=2.5 L=4  (48 fixed transverse plane waves, k=pi/2)", "partA"),
        ("partA_su2_N800.npz",  "SM6  SU(2) beta=2.5 L=6  (48 fixed transverse plane waves, k=pi/3)", "partA"),
        ("k0_su2_b2.4_L8_n2.npz", "SM6  SU(2) beta=2.4 L=8  (12 modes, full decomposition)", "k0"),
    ]
    for fn, lab, kind in specs:
        p = os.path.join(DATA, fn)
        if not os.path.exists(p):
            print("MISSING", p); continue
        r, _ = do_partA(p, lab) if kind == "partA" else do_k0(p, lab)
        show(r); out[lab] = r
    with open(os.path.join(HERE, "sm6_soft_scan_results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote sm6_soft_scan_results.json")
