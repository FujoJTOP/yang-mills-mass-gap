#!/usr/bin/env python3
"""
sw_marginal.py -- set-wise (ball/doubling/small-ball) geometry of the EXACT
one-plaquette SU(2) marginal, and four refutations of the chain the set-wise
species would need:   (D, alpha, rho''>=kappa)  ==>  Var <= C/kappa .

MODEL.  One plaquette, B=1, Wilson:  S(U) = beta(1 - (1/2)tr U),  dmu = e^{-S} dU.
        Class angle phi in [0,pi]: tr U = 2cos(phi), Haar class measure
        (2/pi) sin^2(phi) dphi.   x := 1 - (1/2)tr U = 1 - cos(phi) in [0,2]:
            rho_beta(x) = Z^{-1} e^{-beta x} sqrt(x(2-x)),   Z = e^{-beta} pi I_1(beta)/beta.
        Substitution x = 1 + sin(u) makes sqrt(x(2-x))dx = cos^2(u)du (smooth), so
        all quadratures below are spectrally accurate.

COMPUTED (per beta):  chi = Var(x) (the one-plaquette G'', NOT the 4D KL moment);
  doubling constant D(r) = sup_x mu(B(x,2r))/mu(B(x,r)); small-ball exponents at
  the mode and at the edge; rho''(mode).  Then the refutations
  (A) SCALE INVARIANCE: mu -> lambda_*mu keeps (D,alpha) fixed, multiplies Var by lambda^2
  (B) TAIL PLACEMENT:   mu -> (1-eps)mu + eps delta_R keeps all local data, Var ~ eps R^2
  (C) CURVATURE SIGN:   rho''(mode) < 0 for EVERY beta  => the hypothesis rho''>=kappa>0 is empty
  (D) FREE GAUSSIAN:    N(0,1/m): Var = 1/m -> inf while D and alpha are m-free.

Run:  py -3.12 sw_marginal.py      (writes sw_marginal.json next to it)
"""
import json
import numpy as np
from scipy.special import iv
from scipy.stats import norm

NU = 40001
U = np.linspace(-np.pi / 2, np.pi / 2, NU)
X = 1.0 + np.sin(U)
COS2 = np.cos(U) ** 2                       # sqrt(x(2-x)) dx = cos^2 u du
rs_sb = np.array([2e-4, 5e-4, 1e-3, 2e-3, 5e-3])


def stats(beta):
    w = np.exp(-beta * X) * COS2
    Z = np.trapezoid(w, U)                  # = int e^{-beta x} sqrt(x(2-x)) dx
    p = w / Z
    E = lambda f: np.trapezoid(p * f, U)
    mean = E(X)
    var = E(X ** 2) - mean ** 2
    cdf = np.concatenate(([0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(U))))
    cdf /= cdf[-1]
    return dict(beta=beta, Z=Z, mean=mean, var=var, x=X, cdf=cdf)


def cdf_at(st, xs):
    return np.interp(xs, st["x"], st["cdf"])


def mass(st, a, b):
    return cdf_at(st, np.clip(b, 0.0, 2.0)) - cdf_at(st, np.clip(a, 0.0, 2.0))


def doubling(st, rs=(0.002, 0.005, 0.01, 0.02, 0.05)):
    xs = np.linspace(0.0, 2.0, 4001)
    out = {}
    for r in rs:
        num = mass(st, xs - 2 * r, xs + 2 * r)
        den = mass(st, xs - r, xs + r)
        ok = den > 1e-14
        out[f"{r:g}"] = float(np.max(num[ok] / den[ok]))
    return out


def sb_exponent(st, x0, rs=rs_sb):
    vals = np.array([mass(st, x0 - r, x0 + r) for r in rs])
    A = np.vstack([np.log(rs), np.ones_like(rs)]).T
    return float(np.linalg.lstsq(A, np.log(vals), rcond=None)[0][0])


def mode_of(beta):
    if beta == 0.0:
        return 1.0
    disc = (2 * beta + 1) ** 2 - 4 * beta
    return ((2 * beta + 1) - np.sqrt(disc)) / (2 * beta)


def rho_norm(st, x):
    b = st["beta"]
    return np.exp(-b * x) * np.sqrt(max(x * (2 - x), 0.0)) / st["Z"]


def main():
    out = {}
    for beta in (0.0, 0.5, 1.0, 2.5, 4.0):
        st = stats(beta)
        x0 = mode_of(beta)
        logrho2 = -0.5 / x0 ** 2 - 0.5 / (2 - x0) ** 2      # (log rho)''(x0)
        # Bessel cross-check of Var = (log Z_bessel)'' , Z_bessel = e^{-b} pi I1(b)/b
        b = max(beta, 1e-5)
        f = lambda t: np.log(iv(1, t) / t) - t
        h = 1e-4
        d2 = (f(b + h) - 2 * f(b) + f(b - h)) / h ** 2
        out[f"beta={beta}"] = dict(
            var=float(st["var"]), var_bessel_xcheck=float(d2), mean=float(st["mean"]),
            mode=float(x0), doubling=doubling(st),
            alpha_centre=sb_exponent(st, x0), alpha_edge=sb_exponent(st, 1e-9),
            rho_mode=float(rho_norm(st, x0)),
            rho2_over_rho=float(logrho2),
            rho2_mode=float(rho_norm(st, x0) * logrho2))

    # --------------------------------------------------- (A) scale invariance
    st = stats(2.5)
    out["A_scale_invariance"] = dict(
        note="mu_lambda := pushforward of mu by x->lambda x.  D(mu_lambda)=D(mu) and "
             "alpha(mu_lambda)=alpha(mu) identically (both dimensionless); Var(mu_lambda)=lambda^2 Var(mu).",
        lam=3.0, D_2p5=float(np.mean(list(doubling(st).values()))),
        Var_mu=float(st["var"]), Var_scaled=float(9 * st["var"]),
        alpha_mu=sb_exponent(st, mode_of(2.5)), alpha_scaled=sb_exponent(st, 3 * mode_of(2.5)))
    # concrete uniform family [0,2] vs [0,6]: same D at every r, Var 0.333 vs 3.0
    def Dr_uniform(lam, r):
        mu = lambda a, b: max(0.0, min(b, 2 * lam) - max(a, 0.0)) / (2 * lam)
        xs = np.linspace(0, 2 * lam, 4001)
        return float(np.max([mu(x - 2 * r, x + 2 * r) / mu(x - r, x + r) for x in xs]))
    out["A_scale_invariance"]["uniform_check"] = {
        f"r={r:g}": [Dr_uniform(1.0, r), Dr_uniform(3.0, r)] for r in (0.005, 0.02, 0.1)}
    out["A_scale_invariance"]["uniform_var"] = [1 / 3, 3.0]

    # ------------------------------------------------- (B) tail placement
    st = stats(2.5)
    var0, m = st["var"], st["mean"]
    rows = []
    for eps in (1e-3, 1e-2, 1e-1):
        for R in (5.0, 20.0, 100.0):
            rows.append(dict(eps=eps, R=R,
                             var=float((1 - eps) * var0 + eps * R ** 2
                                       + eps * (1 - eps) * (R - m) ** 2)))
    out["B_tail_placement"] = dict(
        note="for r<=r0 with R-r0>2 the local data (doubling at scale r<=r0, alpha at the mode, "
             "shape of rho near the mode) are unchanged up to the explicit factor (1-eps); "
             "Var grows like eps*R^2.",
        base_var=float(var0), base_mean=float(m), rows=rows)

    # ---------------------------------------------- (D) free Gaussian family
    xs = np.linspace(-8, 8, 3201)
    gD = {}
    for r in (0.01, 0.05, 0.2):
        num = norm.cdf(xs + 2 * r) - norm.cdf(xs - 2 * r)
        den = norm.cdf(xs + r) - norm.cdf(xs - r)
        gD[f"{r:g}"] = float(np.max(num[den > 1e-15] / den[den > 1e-15]))
    out["D_free_gaussian"] = dict(
        doubling=gD, alpha_centre=1.0,
        note="N(0,sigma^2) is self-similar: D and alpha are sigma-free, while Var=sigma^2=1/m. "
             "A free MASSIVE field has m>0, Var=1/m finite, f1^2=0; m->0 sends Var->inf with "
             "D,alpha fixed.",
        family=[dict(m=mm, sigma2=1 / mm, rho2_at_mode=float(-norm.pdf(0.0) / (1 / mm) ** 1.5))
                for mm in (4.0, 1.0, 0.25, 0.0625, 0.01)])

    # ---------------------------------------------- (C) curvature sign summary
    out["C_curvature_sign"] = dict(
        note="(log rho)''(x) = -1/(2x^2) - 1/(2(2-x)^2) < 0 on all of (0,2); equal to -1 at x=1. "
             "Hence rho''(mode) < 0 for every beta, and for every Gaussian: the hypothesis "
             "'rho'' >= kappa > 0 at the centre' is empty in every exactly solvable marginal. "
             "Where rho''>0 (near the edges) the value is fixed by the Haar Jacobian, not by beta.",
        logrho2_at_1=-1.0,
        logrho2_by_beta={str(b): float(-0.5 / mode_of(b) ** 2 - 0.5 / (2 - mode_of(b)) ** 2)
                         for b in (0.0, 0.5, 1.0, 2.5, 4.0)})

    print(json.dumps(out, indent=1))
    with open("sw_marginal.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
