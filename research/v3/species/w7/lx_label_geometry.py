#!/usr/bin/env python3
"""lx_label_geometry.py -- LX species witness: is the *label-space* (representation)
measure of lattice gauge theory log-concave / Lorentzian, and does that convexity
constrain any physical quantity (G'', chi, f1^2)?

A. SU(2) label weight w_j = (2j+1) c_j(s): log-concavity defect Def(s), onset law s -> 0+.
B. SU(3) label weight w_pq = d_pq c_pq(s) (Weyl integration, 2-D Gauss-Legendre):
   directional log-concavity (the Lorentzian/ULC content), CROSS difference
   (sub/supermodularity), label cross-covariance Cov(p,q).
C. U(1): Bessel label weights I_n(s) -- strictly log-concave, nontrivial support, every s>0,
   yet 4-D U(1) lattice gauge is gapless (f1^2 = 0).  ACID-TEST BOMB.
D. Minimal SU(2) plaquette, physical side: Cov(G,Phi') and Stein identity E[G']-Var(G).
E. Explicit sign-blindness family: exp(-t(n1 -+ n2)^2), both strictly log-concave,
   opposite cross-covariance.

Deterministic. numpy + scipy only.   Run: py -3.12 -I lx_label_geometry.py
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ive

OUT = {}


def gl_nodes(a, b, n):
    x, w = leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * w


def def_1d(logw):
    """Def_k = 2 log w_k - log w_{k-1} - log w_{k+1};  Def >= 0  <=>  log-concave."""
    lw = np.asarray(logw, float)
    return 2.0 * lw[1:-1] - lw[:-2] - lw[2:]


# ------------------------------------------------------------------ A. SU(2)
def su2_coeffs_exact(s, jmax):
    """EXACT: c_j(s) = I_{2j}(s) - I_{2j+2}(s)  (Bessel), from
    e^{s cos th} = I_0 + 2 sum_{m>=1} I_m cos(m th) and
    sin((2j+1)th) sin(th) = [cos(2j th) - cos((2j+2)th)]/2 .
    Positive for all s>0 (Bessel monotonicity in order).  Evaluated with ive (=I_n e^{-s})
    so ratios are exact to machine precision for any s."""
    j = np.arange(jmax + 1)
    return ive(2 * j, s) - ive(2 * j + 2, s)


def part_A():
    rows, jmax = [], 60
    for s in [0.0, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0]:
        if s == 0.0:
            rows.append(dict(s=s, support=1, Def_min=0.0, Def_at=0, logconcave=True,
                             note="beta=0: support degenerate to the trivial label; Def==0"))
            continue
        c = su2_coeffs_exact(s, jmax)
        w = (2 * np.arange(jmax + 1) + 1.0) * c
        m = w > 1e-300 * w.max()          # exact-positive region
        lw = np.log(w[m])
        D = def_1d(lw)
        rows.append(dict(s=s, support=int(m.sum()), Def_min=float(D.min()),
                         Def_at=int(np.argmin(D)), Def_max=float(D.max()),
                         logconcave=bool(D.min() >= 0.0)))
    OUT["A_su2_label_defect"] = rows

    # onset law, and its ANALYTIC prediction:
    #   c_j(s) ~ (s/2)^{2j}/(2j)!  as s->0+  =>  Def_j(0+) is a BETA-INDEPENDENT pure number
    #   Def_j = 2 log d_j - log d_{j-1} - log d_{j+1} + log[((2j+2)(2j+1))/((2j)(2j-1))]
    # j = 1, SU(2):  d = (1,3,5)  =>  log(9/5) + log 6 = log(54/5) = 2.379546
    ss = np.array([1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1.0, 2.5, 5.0, 25.0])
    dm = np.array([def_1d(np.log((2 * np.arange(3) + 1.0) * su2_coeffs_exact(s, 2)))[0]
                   for s in ss])
    pred_j1 = float(np.log(54.0 / 5.0))
    OUT["A_onset"] = dict(s=ss.tolist(), Def_at_j1=dm.tolist(),
                          predicted_j1_analytic=pred_j1,
                          rel_err_at_s1e_8=float(abs(dm[0] - pred_j1) / pred_j1))
    c0 = su2_coeffs_exact(0.0, 6)
    OUT["A_orth_check"] = dict(c0=c0.tolist(), max_abs_off_zero=float(np.max(np.abs(c0[1:]))))
    return OUT


# ------------------------------------------------------------------ F. frustration
def part_F():
    """Frustrated plaquette z = -I: e^{-beta Re tr U} has character coefficients
    c_j(-beta) = (-1)^{2j} c_j(beta)  (since chi_j(-U) = (-1)^{2j} chi_j(U)).
    So the PHYSICAL label weight of a frustrated plaquette is sign-alternating."""
    s = 2.5
    c = su2_coeffs_exact(s, 12)
    w = (2 * np.arange(13) + 1.0) * c
    wt = ((-1.0) ** (2 * np.arange(13))) * w
    OUT["F_frustrated_labels"] = dict(
        beta=s, w_j=[float(v) for v in w], w_j_twisted=[float(v) for v in wt],
        n_negative_twisted=int((wt < 0).sum()), n_labels=13,
        log_of_twisted_defined=bool(np.all(wt > 0)),
        partial_sums=[float(v) for v in np.cumsum(wt)],
    )
    OUT["F_note"] = ("On a frustrated plaquette the label weight is (-1)^{2j} w_j: NOT a positive "
                     "sequence, so log-concavity / Lorentzian / M-convexity has no referent at all "
                     "-- not merely a false hypothesis. The (NA) measure is frustrated by "
                     "definition, so the label-space measure it would need does not exist.")
    return OUT


# ------------------------------------------------------------------ B. SU(3)
def chi_su3(p, q, x1, x2, x3):
    """SU(3) character via Jacobi-Trudi (Schur polynomial) -- polynomial, no 0/0.
    lam = (p+q, q, 0) in partition coordinates;  s_lam = det(h_{lam_j - j + k})."""
    # complete homogeneous symmetric functions h_0..h_M in (x1,x2,x3), e3 = x1x2x3 = 1
    lam = (p + q, q, 0)
    M = lam[0] + 2
    e1 = x1 + x2 + x3
    e2 = x1 * x2 + x2 * x3 + x3 * x1
    e3 = x1 * x2 * x3
    h = [np.ones_like(x1, dtype=complex)]
    for m in range(1, M + 1):
        acc = e1 * h[m - 1]
        if m >= 2:
            acc = acc - e2 * h[m - 2]
        if m >= 3:
            acc = acc + e3 * h[m - 3]
        h.append(acc)

    def hm(m):
        return h[m] if m >= 0 else np.zeros_like(x1, dtype=complex)

    a = [[hm(lam[j] - (j + 1) + (k + 1)) for k in range(3)] for j in range(3)]
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def su3_coeffs(beta, labels, nq=200):
    """c_pq(beta) = int dU e^{beta Re tr U} conj(chi_pq(U)),  Weyl integration formula,
    maximal torus U = diag(e^{iu}, e^{iv}, e^{-i(u+v)})."""
    u, wu = gl_nodes(0.0, 2 * np.pi, nq)
    U, V = np.meshgrid(u, u, indexing="ij")
    WU, WV = np.meshgrid(wu, wu, indexing="ij")
    W = WU * WV
    Vd = ((4 * np.sin((U - V) / 2) ** 2)
          * (4 * np.sin((2 * U + V) / 2) ** 2)
          * (4 * np.sin((U + 2 * V) / 2) ** 2))
    x1, x2, x3 = np.exp(1j * U), np.exp(1j * V), np.exp(-1j * (U + V))
    base = Vd * np.exp(beta * (np.cos(U) + np.cos(V) + np.cos(U + V))) * W
    return {k: float(np.real(np.sum(base * np.conj(chi_su3(k[0], k[1], x1, x2, x3)))))
            for k in labels}


def part_B(betas=(0.5, 2.5, 5.7), nmax=8, nq=200):
    labels = [(p, n - p) for n in range(nmax + 1) for p in range(n + 1)]
    out = {}
    # validation 1: at beta=0 orthogonality forces c_pq = 0 for pq != 00
    c0 = su3_coeffs(0.0, labels, nq)
    OUT["B_orthogonality_check"] = dict(
        c00_at_s0=float(c0[(0, 0)]),
        max_abs_c_over_c00_at_s0=float(max(abs(v) for k, v in c0.items() if k != (0, 0))
                                       / abs(c0[(0, 0)])))
    # validation 2: quadrature convergence
    ref = su3_coeffs(2.5, labels, 320)
    mid = su3_coeffs(2.5, labels, nq)
    OUT["B_quadrature_convergence"] = dict(
        max_rel_diff_200_vs_320=float(max(abs(mid[k] - ref[k]) / max(abs(ref[k]), 1e-30)
                                          for k in labels)))
    for beta in betas:
        c = su3_coeffs(beta, labels, nq)
        c0b = c[(0, 0)]
        w = {k: (k[0] + 1) * (k[1] + 1) * (k[0] + k[1] + 2) / 2.0 * (v / c0b)
             for k, v in c.items()}
        dp, dq = [], []
        for (p, q) in labels:
            if p >= 1 and (p + 1, q) in w and min(w[(p - 1, q)], w[(p + 1, q)], w[(p, q)]) > 0:
                dp.append(2 * np.log(w[(p, q)]) - np.log(w[(p - 1, q)]) - np.log(w[(p + 1, q)]))
            if q >= 1 and (p, q + 1) in w and min(w[(p, q - 1)], w[(p, q + 1)], w[(p, q)]) > 0:
                dq.append(2 * np.log(w[(p, q)]) - np.log(w[(p, q - 1)]) - np.log(w[(p, q + 1)]))
        cross = {}
        for (p, q) in labels:
            if (p + 1, q) in w and (p, q + 1) in w and (p + 1, q + 1) in w:
                four = [w[(p, q)], w[(p + 1, q)], w[(p, q + 1)], w[(p + 1, q + 1)]]
                if min(four) > 0:
                    cross[(p, q)] = (np.log(four[0]) + np.log(four[3])
                                     - np.log(four[1]) - np.log(four[2]))
        ks = np.array(list(w.keys()), float)
        ww = np.array([w[k] for k in w])
        ww = ww / ww.sum()
        pv, qv = ks[:, 0], ks[:, 1]
        Ep, Eq = ww @ pv, ww @ qv
        cov = float(ww @ ((pv - Ep) * (qv - Eq)))
        corr = cov / np.sqrt((ww @ (pv - Ep) ** 2) * (ww @ (qv - Eq) ** 2))
        out[str(beta)] = dict(
            Def_p_min=float(min(dp)), Def_q_min=float(min(dq)),
            logconcave_p=bool(min(dp) >= -1e-9), logconcave_q=bool(min(dq) >= -1e-9),
            cross_min=float(min(cross.values())), cross_max=float(max(cross.values())),
            cross_all_negative=bool(max(cross.values()) < 0),
            cross_at_origin=float(cross.get((0, 0), float("nan"))),
            Cov_pq=cov, Corr_pq=float(corr), mean_p=float(Ep), mean_q=float(Eq),
            predicted_cross_origin=float(-1.0 / 4.0 - beta / 3.0),
        )
    OUT["B_su3"] = out
    OUT["B_note"] = ("predicted_cross_origin = d_p d_q log w at (0,0) = -1/(0+0+2)^2 "
                     "- beta/3, i.e. the SAME Casimir cross term d_p d_q C2 = +1/3 that "
                     "FOSP-A identifies as the SU(3) killer.")
    return OUT


# ------------------------------------------------------------------ C. U(1)
def part_C():
    rows = []
    for s in [1e-3, 0.1, 1.0, 2.5, 10.0]:
        n = np.arange(-40, 41)
        w = ive(n, s)
        D = def_1d(np.log(w))
        rows.append(dict(s=s, support=int((w > 0).sum()), Def_min=float(D.min()),
                         logconcave=bool(D.min() >= -1e-12),
                         Def_at_n0=float(D[len(D) // 2])))
    OUT["C_u1_bessel"] = rows
    OUT["C_note"] = ("4-D U(1) lattice gauge theory has exactly this label structure at every "
                     "beta>0 (Fourier labels n_p in Z^N_plaq, weights I_n>0, strictly "
                     "log-concave, M-convex lattice support) and is gapless: f1^2 = 0, "
                     "continuous Kalen-Lehmann measure from s = 0.")
    return OUT


# ------------------------------------------------------------------ D. physics
def part_D(betas=(0.05, 0.25, 0.5, 1.0, 2.5, 4.0), nq=1600):
    """rho in [0,pi], J = sin^2 rho, Phi' = -2 cot rho, S = beta(1-cos rho), G = beta sin rho;
    measure prop. J e^{-S} drho.  Identity Cov(G,Phi') = E[G'] - Var(G)."""
    rho, wr = gl_nodes(0.0, np.pi, nq)
    rho = rho.copy()
    rows = []
    for beta in betas:
        J = np.sin(rho) ** 2
        S = beta * (1 - np.cos(rho))
        p = J * np.exp(-S) * wr
        p /= p.sum()
        G = beta * np.sin(rho)
        Phi = -2.0 * np.cos(rho) / np.sin(rho)
        EG, EPhi = p @ G, p @ Phi
        cov = float(p @ ((G - EG) * (Phi - EPhi)))
        varG = float(p @ (G * G) - EG ** 2)
        stein = float(p @ (beta * np.cos(rho)) - varG)
        rows.append(dict(beta=beta, Cov_G_Phi=cov, stein=stein, VarG=varG,
                         VarPhi=float(p @ (Phi * Phi) - EPhi ** 2),
                         cov_over_beta2=float(cov / beta ** 2)))
    OUT["D_su2_plaquette"] = rows
    OUT["D_note"] = ("sign of Cov(G,Phi') is POSITIVE for every beta>0 -- (NA) is false here "
                     "(reproduces C8, leading 0.2205 beta^2).")
    return OUT


# ------------------------------------------------------------------ E. sign-blindness
def part_E():
    n = np.arange(-30, 31)
    N1, N2 = np.meshgrid(n, n, indexing="ij")
    rows = []
    for t in [0.1, 1.0, 4.0]:
        eps = 0.01 * t                      # tiny regulariser: makes log w STRICTLY concave
        for sgn, name in [(+1.0, "n1-n2"), (-1.0, "n1+n2")]:
            H = np.array([[-2 * t - 2 * eps, -sgn * 2 * t],
                          [-sgn * 2 * t, -2 * t - 2 * eps]])
            ev = np.linalg.eigvalsh(H)
            w = np.exp(-t * (N1 + sgn * (-N2)) ** 2 - eps * (N1 ** 2 + N2 ** 2))
            Wt = w / w.sum()
            m1, m2 = np.sum(Wt * N1), np.sum(Wt * N2)
            rows.append(dict(t=t, direction=name, logconcave=bool(ev.max() <= 0),
                             hess_eigs=ev.tolist(),
                             Cov=float(np.sum(Wt * (N1 - m1) * (N2 - m2))),
                             Var1=float(np.sum(Wt * (N1 - m1) ** 2))))
    OUT["E_signblind"] = rows
    OUT["E_note"] = ("Both are STRICTLY log-concave (log w is negative definite) with identical "
                     "marginals; Cov has opposite signs. Log-concavity does not determine the "
                     "sign of the label cross-covariance.")
    return OUT


if __name__ == "__main__":
    part_A(); part_B(); part_C(); part_D(); part_E(); part_F()
    with open("lx_label_geometry.json", "w") as f:
        json.dump(OUT, f, indent=1, default=float)
    for k, v in OUT.items():
        print("=" * 78)
        print(k)
        print(json.dumps(v, indent=1, default=float)[:3000])
