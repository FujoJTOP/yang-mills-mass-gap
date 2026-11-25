#!/usr/bin/env python3
"""
wg_fss_fit.py -- finite-size scaling fit to the gauge-invariant chi_O data.

Data (identical estimator / measure / APE protocol, beta=2.5, SU(2)):
  L   = 4,   6,   8,   12,  16
  chi = .327650, .377436, .410459, .468797, .586532
  se  = .022273, .021339, .057637, .027880, .054732

Models compared (weighted least squares, w = 1/se^2):
  M(p)  chi(L) = chi_inf + c / L^p        p fixed at 2, 1, 0;  and free p
  M(log) chi(L) = A + B * ln(L)           (no finite limit; A is the L=1 intercept)
Optional saturating-log: chi(L) = chi_inf + c / ln(L)   (bounded, quoted as a check)

p=0 in the ansatz degenerates to a constant: chi = chi_inf + c  ->  chi = <chi>
(the pure plateau / weighted mean).

No downloads, numpy only. Run: py -3.12 -I wg_fss_fit.py
"""
import numpy as np

# ----------------------------------------------------------------------------- data
L   = np.array([4.0, 6.0, 8.0, 12.0, 16.0])
chi = np.array([0.32765014859553365, 0.37743611373983316, 0.4104591281817326,
                0.46879746110233433, 0.5865319300536103])
se  = np.array([0.0222733807527846, 0.02133943765327901, 0.057637361481588344,
                0.027879999106099587, 0.05473246685404735])
w   = 1.0 / se**2
n   = len(L)


def chi2_and_cov(resid, J, w):
    """chi^2, dof, and parameter covariance = inv(J^T W J)."""
    chi2 = float(np.sum(w * resid**2))
    dof  = n - J.shape[1]
    A    = J.T @ (w[:, None] * J)          # normal matrix
    cov  = np.linalg.inv(A)
    return chi2, dof, cov


def fit_fixed_p(p):
    """chi = chi_inf + c / L^p,  p fixed. Linear in (chi_inf, c)."""
    X = np.column_stack([np.ones(n), L**(-p)])
    WX = w[:, None] * X
    A  = X.T @ WX
    b  = X.T @ (w * chi)
    theta = np.linalg.solve(A, b)
    resid = chi - X @ theta
    chi2, dof, cov = chi2_and_cov(resid, X, w)
    perr = np.sqrt(np.diag(cov))
    return dict(p=p, chi_inf=theta[0], chi_inf_err=perr[0],
                c=theta[1], c_err=perr[1], chi2=chi2, dof=dof,
                chi2_dof=chi2/dof)


def fit_free_p(pmin=0.05):
    """chi = chi_inf + c / L^p with p free. Grid p; linear solve (chi_inf,c) each."""
    ps = np.linspace(pmin, 8.0, 79501)
    best = None
    for p in ps:
        X = np.column_stack([np.ones(n), L**(-p)])
        A = X.T @ (w[:, None] * X)
        b = X.T @ (w * chi)
        theta = np.linalg.solve(A, b)
        r = chi - X @ theta
        c2 = float(np.sum(w * r**2))
        # favour approach-from-below (c<0) is automatic if data rises; don't constrain
        if best is None or c2 < best[0]:
            best = (c2, p, theta)
    c2, p, theta = best
    X = np.column_stack([np.ones(n), L**(-p)])
    # curvature-based p error: d2chi2/dp2 at min -> sigma_p = sqrt(2/(d2chi2))
    dp = 1e-3
    def chi2_at(pp):
        Xp = np.column_stack([np.ones(n), L**(-pp)])
        Ap = Xp.T @ (w[:, None] * Xp)
        bp = Xp.T @ (w * chi)
        th = np.linalg.solve(Ap, bp)
        rr = chi - Xp @ th
        return float(np.sum(w * rr**2))
    d2 = (chi2_at(p+dp) - 2*c2 + chi2_at(p-dp)) / dp**2
    p_err = np.sqrt(2.0/d2) if d2 > 0 else np.nan
    chi2, dof, cov = chi2_and_cov(chi - X @ theta, X, w)
    perr = np.sqrt(np.diag(cov))
    return dict(p=p, p_err=p_err, chi_inf=theta[0], chi_inf_err=perr[0],
                c=theta[1], c_err=perr[1], chi2=chi2, dof=dof,
                chi2_dof=chi2/dof)


def fit_log():
    """chi = A + B ln L  (unbounded; no finite thermodynamic limit)."""
    X = np.column_stack([np.ones(n), np.log(L)])
    WX = w[:, None] * X
    A  = X.T @ WX
    b  = X.T @ (w * chi)
    theta = np.linalg.solve(A, b)
    resid = chi - X @ theta
    chi2, dof, cov = chi2_and_cov(resid, X, w)
    perr = np.sqrt(np.diag(cov))
    return dict(A=theta[0], A_err=perr[0], B=theta[1], B_err=perr[1],
                chi2=chi2, dof=dof, chi2_dof=chi2/dof)


def fit_satlog():
    """chi = chi_inf + c / ln(L)  (bounded saturating log; check only)."""
    XL = np.log(L)
    X = np.column_stack([np.ones(n), 1.0/XL])
    WX = w[:, None] * X
    A = X.T @ WX
    b = X.T @ (w * chi)
    theta = np.linalg.solve(A, b)
    resid = chi - X @ theta
    chi2, dof, cov = chi2_and_cov(resid, X, w)
    perr = np.sqrt(np.diag(cov))
    return dict(chi_inf=theta[0], chi_inf_err=perr[0], c=theta[1], c_err=perr[1],
                chi2=chi2, dof=dof, chi2_dof=chi2/dof)


def fit_pure_power():
    """chi = A L^q  (chi_inf = 0 exactly; pure power, no finite limit).
    Fit in log space weighting by relative error (se/chi)."""
    y  = np.log(chi)
    sy = se / chi
    wl = 1.0 / sy**2
    X  = np.column_stack([np.ones(n), np.log(L)])
    A  = X.T @ (wl[:, None] * X)
    b  = X.T @ (wl * y)
    theta = np.linalg.solve(A, b)
    resid = y - X @ theta
    chi2 = float(np.sum(wl * resid**2))
    dof = n - 2
    cov = np.linalg.inv(A)
    perr = np.sqrt(np.diag(cov))
    return dict(A=np.exp(theta[0]), q=theta[1], q_err=perr[1],
                chi2=chi2, dof=dof, chi2_dof=chi2/dof)


def plateau():
    """chi = const (pure plateau = weighted mean)."""
    W = w.sum()
    mu = np.sum(w*chi)/W
    mu_err = 1.0/np.sqrt(W)
    r = chi - mu
    chi2 = float(np.sum(w*r**2))
    dof = n - 1
    return dict(chi_inf=mu, chi_inf_err=mu_err, chi2=chi2, dof=dof,
                chi2_dof=chi2/dof)


def fmt(x, e, nd=3):
    return f"{x:.{nd}f} ± {e:.{nd}f}"


if __name__ == "__main__":
    pdg = np.sqrt(0.0)
    print("="*78)
    print("WG finite-size fit  --  chi_O(L) = chi_inf + c/L^p   (weighted LS, w=1/se^2)")
    print("="*78)
    print("data:")
    for li, ci, si in zip(L, chi, se):
        print(f"  L={int(li):2d}   chi={ci:.4f} ± {si:.4f}   (pull from flat {((ci-np.mean(chi))/si):+.2f})")
    print()

    print("-"*78)
    print("FIXED-EXPONENT fits")
    print("-"*78)
    res = {}
    for p in (2.0, 1.0):
        r = fit_fixed_p(p)
        res[f"p={p:.0f}"] = r
        print(f"  p={p:.0f} fixed:  chi_inf = {fmt(r['chi_inf'], r['chi_inf_err'])}"
              f"   c = {r['c']:+.3f} ± {r['c_err']:.3f}"
              f"   chi2/dof = {r['chi2']:.2f}/{r['dof']} = {r['chi2_dof']:.2f}")
    # literal free-massless growth chi = chi_inf + c L^2  (p = -2)
    rfree = fit_fixed_p(-2.0)
    res["p=-2"] = rfree
    print(f"  p=-2 (free L^2 growth): chi_inf = {fmt(rfree['chi_inf'], rfree['chi_inf_err'])}"
          f"   c = {rfree['c']:+.4f} ± {rfree['c_err']:.4f}"
          f"   chi2/dof = {rfree['chi2']:.2f}/{rfree['dof']} = {rfree['chi2_dof']:.2f}")
    rp = plateau()
    res["p=0"] = rp
    print(f"  p=0 (plateau): chi_inf = {fmt(rp['chi_inf'], rp['chi_inf_err'])}"
          f"                               chi2/dof = {rp['chi2']:.2f}/{rp['dof']} = {rp['chi2_dof']:.2f}")
    print()

    print("-"*78)
    print("FREE-EXPONENT fit")
    print("-"*78)
    rf = fit_free_p(0.05)
    res["free"] = rf
    print(f"  p free (grid >=0.05): p = {fmt(rf['p'], rf['p_err'])}   "
          f"chi_inf = {fmt(rf['chi_inf'], rf['chi_inf_err'])}"
          f"   c = {rf['c']:+.3f} ± {rf['c_err']:.3f}"
          f"   chi2/dof = {rf['chi2']:.2f}/{rf['dof']} = {rf['chi2_dof']:.2f}")
    print("        [DEGENERATE: minimum sits at the p->0 boundary; as p->0 the")
    print("         shift c/L^p -> c is ~constant and chi_inf is unidentifiable")
    print("         (chi_inf -> infinity). Do NOT read this chi_inf as physical.]")
    rb = fit_free_p(0.40)
    res["free_bounded"] = rb
    print(f"  p free (grid >=0.40): p = {fmt(rb['p'], rb['p_err'])}   "
          f"chi_inf = {fmt(rb['chi_inf'], rb['chi_inf_err'])}"
          f"   c = {rb['c']:+.3f} ± {rb['c_err']:.3f}"
          f"   chi2/dof = {rb['chi2']:.2f}/{rb['dof']} = {rb['chi2_dof']:.2f}")
    print()

    print("-"*78)
    print("PURE-POWER fit  chi = A L^q  (chi_inf = 0; no finite limit)")
    print("-"*78)
    rpp = fit_pure_power()
    res["pure_power"] = rpp
    print(f"  A = {rpp['A']:.3f}   q = {fmt(rpp['q'], rpp['q_err'])}"
          f"   chi2/dof = {rpp['chi2']:.2f}/{rpp['dof']} = {rpp['chi2_dof']:.2f}")
    print("        (matches the log-law chi2/dof: both say 'no finite plateau needed')")
    print()

    print("-"*78)
    print("LOG LAWS")
    print("-"*78)
    rl = fit_log()
    res["log"] = rl
    print(f"  chi = A + B ln L :  A = {fmt(rl['A'], rl['A_err'])}   "
          f"B = {rl['B']:+.3f} ± {rl['B_err']:.3f}   chi2/dof = {rl['chi2']:.2f}/{rl['dof']} = {rl['chi2_dof']:.2f}")
    print(f"        -> no finite limit (B>0 diverges); at L=100 chi = {rl['A']+rl['B']*np.log(100.0):.2f}")
    rs = fit_satlog()
    res["satlog"] = rs
    print(f"  chi = chi_inf + c/ln L : chi_inf = {fmt(rs['chi_inf'], rs['chi_inf_err'])}"
          f"   c = {rs['c']:+.2f} ± {rs['c_err']:.2f}"
          f"   chi2/dof = {rs['chi2']:.2f}/{rs['dof']} = {rs['chi2_dof']:.2f}")
    print()

    print("-"*78)
    print("SUMMARY  (chi2/dof ; chi_inf)")
    print("-"*78)
    summary = [
        ("free p",      rf['chi2_dof'],  fmt(rf['chi_inf'], rf['chi_inf_err'])+"*", f"p={fmt(rf['p'], rf['p_err'],2)}"),
        ("free p>=0.4", rb['chi2_dof'],  fmt(rb['chi_inf'], rb['chi_inf_err']), f"p={fmt(rb['p'], rb['p_err'],2)}"),
        ("p=2 fixed",   res['p=2']['chi2_dof'], fmt(res['p=2']['chi_inf'], res['p=2']['chi_inf_err']), ""),
        ("p=-2 freeL2", res['p=-2']['chi2_dof'], fmt(res['p=-2']['chi_inf'], res['p=-2']['chi_inf_err']), ""),
        ("p=1 fixed",   res['p=1']['chi2_dof'], fmt(res['p=1']['chi_inf'], res['p=1']['chi_inf_err']), ""),
        ("p=0 plateau", rp['chi2_dof'], fmt(rp['chi_inf'], rp['chi_inf_err']), ""),
        ("log A+BlnL",  rl['chi2_dof'],  "undefined", f"B={rl['B']:+.3f}±{rl['B_err']:.3f}"),
        ("pure power",  rpp['chi2_dof'], "0 (exact)", f"q={fmt(rpp['q'],rpp['q_err'],2)}"),
        ("sat-log",     rs['chi2_dof'],  fmt(rs['chi_inf'], rs['chi_inf_err']), ""),
    ]
    for name, cd, cinf, extra in summary:
        print(f"  {name:12s}  chi2/dof={cd:7.2f}   chi_inf={cinf:22s} {extra}")
    print("  (* degenerate: at the p->0 boundary, chi_inf is not identifiable)")
    print()
    print("Significance  chi_inf > 0  (fixed-exponent and bounded-free fits):")
    for key, nm in (("p=2","p=2"), ("p=1","p=1"), ("free_bounded","free p>=0.4")):
        r = res[key]
        sig = r['chi_inf']/r['chi_inf_err']
        print(f"  {nm:12s}:  chi_inf/se = {sig:5.2f}  ({'SUPPORTED' if abs(sig)>2 else 'NOT significant'})")
    print()
    print("MODEL COMPARISON (smaller chi2/dof = better):")
    best = min(summary, key=lambda t: t[1])
    print(f"  best chi2/dof: {best[0]} ({best[1]:.2f})")
    print("  plateau (chi_inf const) is EXCLUDED at chi2/dof=7.28 (vs dof=4).")
    print("  log-law and pure-power (both chi_inf=0 / no finite limit) fit as well")
    print("  as any saturating ansatz: the 5 points CANNOT yet distinguish a finite")
    print("  thermodynamic chi_inf from unbounded slow growth. chi_inf>0 is model-")
    print("  dependent (true for p=1,2; absent by construction for log/power).")
    print()
    print("NOTE: no saturating fit has chi2/dof ~ 1; inverse-Hessian errors are formal.")
    print("      With 5 points and 2-3 params, p and chi_inf are poorly constrained.")
    print("      This is the honest reading; L=24/32 will tighten p and chi_inf.")
    print()
    print("-"*78)
    print("FALSIFIABLE PREDICTIONS at L=24,32 (update when data land)")
    print("-"*78)
    for Lp in (24.0, 32.0):
        v_p1 = res['p=1']['chi_inf'] + res['p=1']['c']/Lp
        v_p2 = res['p=2']['chi_inf'] + res['p=2']['c']/Lp**2
        v_pp = rpp['A'] * Lp**rpp['q']
        v_lg = rl['A'] + rl['B']*np.log(Lp)
        v_sl = rs['chi_inf'] + rs['c']/np.log(Lp)
        print(f"  L={int(Lp)}: p=1 sat {v_p1:.3f} | p=2 sat {v_p2:.3f} | sat-log {v_sl:.3f}"
              f" | pure power {v_pp:.3f} | log {v_lg:.3f}")
    print("  saturating family ~0.48-0.53 ; chi_inf=0 family ~0.59-0.70 -> separable at SE~0.05")
    print("="*78)
