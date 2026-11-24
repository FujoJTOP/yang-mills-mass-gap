"""
H5-HEAD reanalysis witness (H5: lattice test of G''(v;a0*) >= delta_W at beta=2.5).

Reads ONLY local campaign data (no downloads):
  partA_su2_N4000.npz  SU(2) beta=2.5 L=4  48 modes  (a,H,G)
  partA_su2_N800.npz   SU(2) beta=2.5 L=6  48 modes  (a,H,G)
  k0_su2_b2.4_L8_n2.npz SU(2) beta=2.4 L=8 12 modes  (a,H,G,P1,P2)  FULL decomposition
  k0_su3_b5.7_L8_n2.npz SU(3) beta=5.7 L=8 12 modes  (a,H,G,P1,P2)  FULL decomposition

For every mode:
  a0*      = KDE mode of P(a)             (the COARSE MINIMISER, V3/(*_loc))
  V''(a0*) = curvature of -log P(a) at a0* (KDE gold standard; full marginal curvature)
  bracket  = E[H|a0*] - Var(G|a0*)        (Wilson-only when only a,H,G are stored)
  W_geo    = E[Phi''|a0*] - Var(Phi'|a0*)  (needs P1,P2; proven >= 0 pointwise)
  X        = -2 Cov(G,Phi'|a0*)            (needs P1,P2)
  G''_id   = bracket + W_geo + X           (exact identity, needs P1,P2)
One-sided bounds (k0 datasets only): keep ONLY the proven W_geo>=0, treat bracket+X
statistically:  LB_rig = (bracket+X) - 1.645 SE.  If < 0, 0 is not excluded.

Object caveat (honest): a = <v, log U> with v a transverse plane wave is the
GAUGE-DEPENDENT gluon amplitude in a flat chart. It is NOT the gauge-invariant
conditioned coarse curvature the head inequality needs. partA stores no Phi
derivatives, so no rigorous-style bound is formable there at all.
"""
import os, json
import numpy as np
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "unconditional", "w7")


def _kappa_kde(a, a0, bw=0.20, rfit=0.8, ngrid=41):
    """Curvature of -log KDE at a0 (quadratic fit of -log p over |a-a0|<rfit*sd)."""
    a = np.asarray(a, float)
    sd = a.std()
    kde = gaussian_kde(a, bw_method=bw)
    xs = np.linspace(a0 - rfit * sd, a0 + rfit * sd, ngrid)
    p = np.maximum(kde(xs), 1e-300)
    y = -np.log(p)
    A = np.vstack([np.ones_like(xs), (xs - a0), (xs - a0) ** 2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return 2.0 * coef[2]


def _mode_kde(a, bw=0.20, ngrid=161, span=4.0):
    """KDE mode (argmax of P(a)) on a grid over +-span*sd."""
    a = np.asarray(a, float)
    sd = a.std()
    kde = gaussian_kde(a, bw_method=bw)
    xs = np.linspace(-span * sd, span * sd, ngrid)
    return float(xs[int(np.argmax(kde(xs)))])


def _kernel_cond(a, cols, a0, bw=0.30):
    """Nadaraya-Watson conditional moments at a=a0 (Gaussian kernel)."""
    a = np.asarray(a, float)
    sd = a.std()
    h = bw * sd
    w = np.exp(-0.5 * ((a - a0) / h) ** 2)
    ws = w.sum()
    means = {k: float(np.sum(w * np.asarray(x, float)) / ws) for k, x in cols.items()}
    out = {"E_" + k: v for k, v in means.items()}
    for k, x in cols.items():
        for l, y in cols.items():
            if k <= l:
                xm = np.asarray(x, float) - means[k]
                ym = np.asarray(y, float) - means[l]
                out["C_%s_%s" % (k, l)] = float(np.sum(w * xm * ym) / ws)
    return out


# ---------------------------------------------------------------------------
def reanalyze_partA(path, label, bw=0.30, kdebw=0.20, nblk=20):
    """partA: only a,H,G.  ->  a0*, V''(a0*) KDE, bracket, 1/Var.  No W_geo/X split."""
    d = np.load(path, allow_pickle=True)
    modes = sorted({k[:-2] for k in d.files if k.endswith("_H")})
    A = np.column_stack([d[m + "_a"] for m in modes])
    H = np.column_stack([d[m + "_H"] for m in modes])
    G = np.column_stack([d[m + "_G"] for m in modes])
    N, Nm = A.shape
    edges = np.linspace(0, N, nblk + 1).astype(int)

    def per_mode(sl):
        rows = []
        for j in range(Nm):
            a = A[sl, j]
            a0 = _mode_kde(a, bw=kdebw)
            kappa = _kappa_kde(a, a0, bw=kdebw)
            cond = _kernel_cond(a, {"H": H[sl, j], "G": G[sl, j]}, a0, bw=bw)
            bracket = cond["E_H"] - cond["C_G_G"]
            rows.append(dict(a0=a0, kappa=kappa, bracket=bracket,
                             inv=1.0 / a.var(), Wfull=kappa - bracket))
        return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}

    full = per_mode(slice(None))
    leave = {k: [] for k in full}
    for b in range(nblk):
        keep = np.ones(N, bool); keep[edges[b]:edges[b + 1]] = False
        s = per_mode(keep)
        for k in full:
            leave[k].append(s[k])
    se = {k: float(np.sqrt((nblk - 1) / nblk * np.sum((np.array(v) - np.mean(v)) ** 2)))
          for k, v in leave.items()}
    return dict(label=label, object="gauge-dependent gluon (a,H,G only)", N=N, Nm=Nm,
                full=full, se=se, nblk=nblk, has_phi=False)


def reanalyze_k0(path, label, bw=0.30, kdebw=0.20, nblk=20):
    """k0: full a,H,G,P1,P2 -> exact identity + rigorous-style one-sided bound at a0*."""
    d = np.load(path)
    A, H, G, P1, P2 = d["a"], d["H"], d["G"], d["P1"], d["P2"]
    N, Nm = A.shape
    beta = float(d["beta"]); L = int(d["L"]); group = str(d["group"])
    edges = np.linspace(0, N, nblk + 1).astype(int)

    def stats(sl):
        rows = []
        for j in range(Nm):
            a = A[sl, j]
            a0 = _mode_kde(a, bw=kdebw)
            c = _kernel_cond(a, {"H": H[sl, j], "G": G[sl, j],
                                 "P1": P1[sl, j], "P2": P2[sl, j]}, a0, bw=bw)
            bracket = c["E_H"] - c["C_G_G"]
            Wgeo = c["E_P2"] - c["C_P1_P1"]
            X = -2.0 * c["C_G_P1"]
            rows.append(dict(a0=a0, bracket=bracket, Wgeo=Wgeo, X=X,
                             Gid=bracket + Wgeo + X, brX=bracket + X,
                             kappa=_kappa_kde(a, a0, bw=kdebw), inv=1.0 / a.var()))
        return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}

    full = stats(slice(None))
    leave = {k: [] for k in full}
    for b in range(nblk):
        keep = np.ones(N, bool); keep[edges[b]:edges[b + 1]] = False
        s = stats(keep)
        for k in full:
            leave[k].append(s[k])
    se = {k: float(np.sqrt((nblk - 1) / nblk * np.sum((np.array(v) - np.mean(v)) ** 2)))
          for k, v in leave.items()}
    lb_rig = full["brX"] - 1.645 * se["brX"]
    lb_stat = full["kappa"] - 1.645 * se["kappa"]
    return dict(label=label, object="gauge-dependent gluon (full decomposition)",
                beta=beta, L=L, group=group, N=N, Nm=Nm, full=full, se=se,
                LB_rig=lb_rig, LB_stat=lb_stat, nblk=nblk, has_phi=True)


def show(r):
    f, s = r["full"], r["se"]
    print("\n=== %s ===" % r["label"])
    print("  object:", r["object"], "| N=%d modes=%d nblk=%d" % (r["N"], r["Nm"], r["nblk"]))
    print("  a0* = %+.4f   (coarse minimiser; symmetric => ~0)" % f["a0"])
    if r["has_phi"]:
        print("  bracket=%+.3f+-%.3f  Wgeo=%+.3f+-%.3f  X=%+.3f+-%.3f  Gid=%+.3f+-%.3f"
              % (f["bracket"], s["bracket"], f["Wgeo"], s["Wgeo"], f["X"], s["X"],
                 f["Gid"], s["Gid"]))
        print("  KDE G''(a0*)=%+.3f+-%.3f   1/Var=%+.3f" % (f["kappa"], s["kappa"], f["inv"]))
        print("  RIGOROUS-style 95%% 1-sided LB (Wgeo>=0 only, brX stat) = %+.3f" % r["LB_rig"])
        print("  STATISTICAL  95%% 1-sided LB (trusts full G'')          = %+.3f" % r["LB_stat"])
    else:
        print("  bracket=%+.3f+-%.3f  KDE G''(a0*)=%+.3f+-%.3f  1/Var=%+.3f  Wfull(G''-br)=%+.3f"
              % (f["bracket"], s["bracket"], f["kappa"], s["kappa"], f["inv"], f["Wfull"]))
        print("  (no Phi-derivative columns stored -> W_geo/X NOT separable -> no rigorous-style LB)")


if __name__ == "__main__":
    out = {}
    specs = [
        ("partA_su2_N4000.npz", "PART A  SU(2) beta=2.5 L=4 (k=pi/2)  EXACT HEAD", "partA"),
        ("partA_su2_N800.npz",  "PART A  SU(2) beta=2.5 L=6 (k=pi/3)  EXACT HEAD", "partA"),
        ("k0_su2_b2.4_L8_n2.npz", "K0      SU(2) beta=2.4 L=8 (k=pi/4,pi/2) FULL DECOMP", "k0"),
        ("k0_su3_b5.7_L8_n2.npz", "K0      SU(3) beta=5.7 L=8 (k=pi/4,pi/2) FULL DECOMP", "k0"),
    ]
    for fn, lab, kind in specs:
        p = os.path.join(DATA, fn)
        if not os.path.exists(p):
            print("MISSING", p); continue
        r = reanalyze_partA(p, lab) if kind == "partA" else reanalyze_k0(p, lab)
        show(r); out[lab] = r
    with open(os.path.join(HERE, "h5_head_reanalysis_results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote h5_head_reanalysis_results.json")
