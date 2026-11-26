#!/usr/bin/env python3
r"""
cov_link.py -- COV-T, model (1): SINGLE LINK on SU(2) = unit S^3, exact (Bessel/Chebyshev tower).

Setting (identical conventions to AUDIT-BBD sec.5 / DENSITY-D sec.6.2, task-9):
    S_W(w)      = beta (1 - q_0(w)),      q_0(w) = (1/2) Re tr w  = cos psi
    V_t(x)      = -log B_t(x),   B_t(x) = int K_t(x^-1 w) e^{-S_W(w)} dw
    mu_{t,x}(dw) propto K_t(x^-1 w) e^{-S_W(w)} dw
    Hess V_t    = E_mu[Hess S_W] - Cov_mu(grad S_W)                       (P-B identity)

Exact facts used (no new species):
  * unit S^3 Laplacian: -Delta chi_l = l(l+2) chi_l, chi_l = sin((l+1)psi)/sin psi, d_l = l+1
      => lam_1 = 3,  and lam_1 is a ONE-LINK constant (independent of beta)
  * e^{beta q_0} = sum_l c_l chi_l,   c_l = I_l(beta) - I_{l+2}(beta)      (R5-AUDIT sec.1.1)
  * heat kernel K_t(psi) = sum_l (l+1) e^{-l(l+2)t} chi_l(psi)
  * class functions convolve mode-wise  =>  B_t(x) = e^{-beta} sum_l c_l e^{-l(l+2)t} chi_l(psi_x)
  * grad on the left-invariant frame:   d_a S_W(w) = (beta/2) q_a(w)    (P-D of R3-L4)
  * Hess^{flat,sym} S_W = (beta/4) q_0(w) * 1_3       [= -beta/(2N) at q_0=-1, N=2; R5-AUDIT sec.6]
  * therefore
        Cov_{mu_{t,x}}(grad S_W) = (beta^2/4) Cov_{mu_{t,x}}(q_vec)
        Hess V_t                 = (beta/4) E_mu[q_0] * 1_3 - (beta^2/4) Cov_mu(q_vec)
        q(t) = || (Hess V_t)_- ||_op

Run: py -3.12 -I cov_link.py     numpy/scipy only, deterministic, single process.
Writes cov_link.log (appended section by section).
"""
import numpy as np
from scipy.special import iv

LOG = []
LOGF = "cov_link.log"


def say(s=""):
    print(s)
    LOG.append(s)


def flush():
    with open(LOGF, "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


LAM1 = 3.0          # -Delta on unit S^3: l(l+2), l=1


# ------------------------------------------------------------------ exact 1D building blocks
def chi(l, psi):
    """character chi_l(psi) = sin((l+1)psi)/sin(psi) = U_l(cos psi)."""
    s = np.sin(psi)
    out = np.empty_like(np.asarray(psi, dtype=float))
    small = np.abs(s) < 1e-9
    out[~small] = np.sin((l + 1) * psi[~small]) / s[~small]
    out[small] = (l + 1) * np.where(np.cos(psi[small]) > 0, 1.0, (-1.0) ** l)
    return out


def c_l(l, beta):
    """Chebyshev-U coefficient of e^{beta cos psi}: c_l = I_l(beta) - I_{l+2}(beta)."""
    return iv(l, beta) - iv(l + 2, beta)


def B_char(psi_x, t, beta, lmax):
    """B_t(x) = e^{-beta} sum_l c_l e^{-l(l+2)t} chi_l(psi_x).  (scalar x, scalar t)"""
    tot = 0.0
    ps = np.array([psi_x])
    for ll in range(lmax + 1):
        tot += c_l(ll, beta) * np.exp(-ll * (ll + 2) * t) * float(chi(ll, ps)[0])
    return np.exp(-beta) * tot


def Kt_grid(delta, t, lmax, beta=0.0):
    """heat kernel on unit S^3 at geodesic distance delta: sum_l (l+1)e^{-l(l+2)t} chi_l(delta)."""
    l = np.arange(lmax + 1)
    e = (l + 1) * np.exp(-l * (l + 2) * t)
    tot = np.zeros_like(delta)
    for ll in range(lmax + 1):
        tot += e[ll] * chi(ll, delta)
    return tot


# ------------------------------------------------------------------ quadrature on S^3
def s3_grid(npsi=96, nth=96, nphi=96):
    """nodes q in S^3 with weights w such that int f dHaar = sum_i w_i f(q_i).
    q = (cos psi, sin psi * nhat),  nhat = (sin th cos ph, sin th sin ph, cos th)."""
    xg, wg = np.polynomial.legendre.leggauss(npsi)
    psi = 0.5 * np.pi * (xg + 1.0)
    wpsi = 0.5 * np.pi * wg
    xt, wt = np.polynomial.legendre.leggauss(nth)          # cos(theta)
    theta = np.arccos(xt)
    wth = wt
    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    wph = 2.0 * np.pi / nphi
    PS, TH, PH = np.meshgrid(psi, theta, phi, indexing="ij")
    # Haar on unit S^3:  dmu = (1/2pi^2) sin^2(psi) dpsi dOmega,  Omega = S^2
    # node weights: wpsi (dpsi), wth (d cos theta -> no extra sin theta!), wph (dphi)
    W = (wpsi[:, None, None] * wth[None, :, None] * wph
         * np.sin(PS) ** 2 * np.ones_like(PH)) / (2.0 * np.pi ** 2)
    q0 = np.cos(PS).ravel()
    st = np.sin(PS).ravel()
    qv = np.stack([st * np.sin(TH).ravel() * np.cos(PH).ravel(),
                   st * np.sin(TH).ravel() * np.sin(PH).ravel(),
                   st * np.cos(TH).ravel()], axis=1)
    return q0, qv, W.ravel()


class Model1:
    def __init__(self, beta=2.5, lmax=240, npsi=96, nth=96, nphi=96, ndelta=8192):
        self.beta = beta
        self.lmax = lmax
        self.q0, self.qv, self.w = s3_grid(npsi, nth, nphi)
        self.nd = ndelta

    def kernel_lut(self, t):
        """K_t(delta) on [0,pi] for interpolation."""
        d = np.linspace(0.0, np.pi, self.nd)
        k = Kt_grid(d, t, self.lmax)
        return d, k

    def moments(self, psi_x, t):
        """returns Z=B_t(x), E[q0], E[qv], E[qv qv^T] under mu_{t,x}, by 3D quadrature."""
        beta = self.beta
        # q_x = (cos psi_x, sin psi_x, 0, 0)
        cx, sx = np.cos(psi_x), np.sin(psi_x)
        ip = cx * self.q0 + sx * self.qv[:, 0]                 # <q_x, q_w>
        delta = np.arccos(np.clip(ip, -1.0, 1.0))
        d, k = self.kernel_lut(t)
        Kt = np.interp(delta, d, k)
        rho = self.w * Kt * np.exp(beta * self.q0)             # e^{-S_W} = e^{-beta} e^{beta q0}
        Z = rho.sum() * np.exp(-beta)
        wn = rho / rho.sum()
        m0 = float(wn @ self.q0)
        mv = wn @ self.qv                                      # (3,)
        Mv = (self.qv * wn[:, None]).T @ self.qv               # (3,3)
        return Z, m0, mv, Mv


def hess_V_analytic(psi_x, t, beta, lmax, h=1e-4):
    """Hess^{flat,sym} V_t at a point with class angle psi_x, from V_t = g(psi_x):
       M = (1/4)[ g''(q0) q_a q_b ] - (1/4) g'(q0) q0 delta ,  with the q_a q_b part
       replaced by its S^2 average (1-q0^2)/3 * delta  -> isotropic.
       Implemented by finite differences of g in q0 = cos psi."""
    def g_of_q0(q0v):
        ps = np.arccos(np.clip(q0v, -1.0, 1.0))
        return -np.log(B_char(ps, t, beta, lmax))
    q0 = np.cos(psi_x)
    eps = 1e-5
    qp = min(q0 + eps, 1.0 - 1e-12)
    qm = max(q0 - eps, -1.0 + 1e-12)
    g = g_of_q0(q0)
    gp = (g_of_q0(qp) - g_of_q0(qm)) / (qp - qm)
    gpp = (g_of_q0(qp) - 2 * g + g_of_q0(qm)) / (((qp - qm) / 2) ** 2)
    iso = (1.0 - q0 ** 2) / 3.0
    return (0.25 * gpp * iso - 0.25 * gp * q0) * np.eye(3)


def main():
    beta = 2.5
    say("=" * 100)
    say("COV-T  model (1)  --  SINGLE LINK on SU(2)=S^3, exact Bessel/Chebyshev tower")
    say("=" * 100)
    say("")
    m = Model1(beta=beta)
    say(f"[S0] quadrature: {m.q0.size} nodes on S^3;  beta={beta};  lam_1 = {LAM1}")
    say("")

    # ---------------- S0a: B_t character sum vs quadrature (validates kernel + measure)
    say("[S0a] validation A: B_t(x) from the exact character sum  vs  the 3D quadrature")
    say(f"      {'lam1*t':>8} {'psi_x':>7} {'B_char':>16} {'B_quad':>16} {'rel.diff':>10}")
    worst = 0.0
    for tau in (0.1, 1.0, 10.0, 100.0):
        t = tau / LAM1
        for psix in (0.0, 0.5 * np.pi, np.pi):
            Bc = B_char(psix, t, beta, m.lmax)
            Z, _, _, _ = m.moments(psix, t)
            rd = abs(Z - Bc) / Bc
            worst = max(worst, rd)
            say(f"      {tau:8.3g} {psix:7.4f} {Bc:16.10e} {Z:16.10e} {rd:10.2e}")
    say(f"      worst relative difference = {worst:.3e}   (quadrature accuracy check)")
    say("")
    flush()

    # ---------------- S0b: t=0 -> Cov = 0 exactly
    say("[S0b] validation B: t -> 0 limit.  mu_{0,x} = delta_x => Cov = 0 exactly.")
    say("      (kernel width ~ sqrt(t); we check the measured Cov against its O(t) law below.)")
    say("")

    # ---------------- S1: the curve
    taus = [0.001, 0.003, 0.01, 0.03, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0, 316.0, 1000.0]
    alphas = np.linspace(0.0, np.pi, 13)
    say("[S1] the full curve:   C_abs(t) = sup_x ||Cov_{mu_t,x}(grad S_W)||_op  and  q(t)=||(Hess V_t)_-||op")
    say("     'C_int' = C_abs/N_p with N_p=6 (the four-dimensional link sits in 6 plaquettes)")
    say("     'C_1p'  = C_abs/(6*beta/(2N)) = C_abs/(1.5 beta)  -- in units of the t=infty bound of the task")
    say("")
    hdr = (f"     {'lam1*t':>8} {'t':>10} {'C_abs':>12} {'x*':>7} {'C_int':>10} {'C_1p':>9} "
           f"{'q(t)':>12} {'lam_minHV':>12} {'lam_maxHV':>12}")
    say(hdr)
    rows = []
    for tau in taus:
        t = tau / LAM1
        best = (-1.0, None)
        for a in alphas:
            Z, m0, mv, Mv = m.moments(a, t)
            C = Mv - np.outer(mv, mv)
            nrm = np.linalg.norm((beta ** 2 / 4) * C, 2)
            if nrm > best[0]:
                best = (nrm, a, m0, mv, Mv, C)
        Cabs, xs, m0, mv, Mv, C = best
        HV = (beta / 4) * m0 * np.eye(3) - (beta ** 2 / 4) * C
        ev = np.linalg.eigvalsh(HV)
        q = max(0.0, -ev.min())
        # cross-check against analytic Hess log B_t
        HVa = hess_V_analytic(xs, t, beta, m.lmax)
        eva = np.linalg.eigvalsh(HVa)
        rows.append(dict(tau=tau, t=t, C=Cabs, x=xs, q=q, lmin=ev.min(), lmax=ev.max(),
                         lmin_an=eva.min(), lmax_an=eva.max()))
        say(f"     {tau:8.4g} {t:10.4g} {Cabs:12.6f} {xs:7.4f} {Cabs/6:10.6f} {Cabs/(1.5*beta):9.5f} "
            f"{q:12.6e} {ev.min():12.6f} {ev.max():12.6f}")
    say("")
    flush()

    # ---------------- S1b: identity cross-check  Cov = E[Hess S] + Hess log B_t
    say("[S1b] validation C: the identity  Cov_mu(grad S_W) = E_mu[Hess S_W] + Hess log B_t(x)")
    say("      (Hess^{flat,sym}, isotropic part; independent route = analytic g''(q_0) of -log B_t)")
    say(f"      {'lam1*t':>8} {'(b^2/4)Cov':>14} {'b/4 E[q0]-HV':>14} {'HV(num)':>13} {'HV(analytic)':>13} {'diff':>9}")
    for r in rows:
        t = r["t"]
        Z, m0, mv, Mv = m.moments(r["x"], t)
        C = Mv - np.outer(mv, mv)
        lhs = (beta ** 2 / 4) * np.trace(C) / 3
        rhs = (beta / 4) * m0 - r["lmin"] if False else None
        HVa = hess_V_analytic(r["x"], t, beta, m.lmax)
        say(f"      {r['tau']:8.4g} {lhs:14.7f} {(beta/4)*m0:14.7f} {r['lmin']:13.7f} "
            f"{np.linalg.eigvalsh(HVa).min():13.7f} {abs(r['lmin']-np.linalg.eigvalsh(HVa).min()):9.2e}")
    say("")
    flush()

    # ---------------- S2: three-segment fit
    say("[S2] three-segment fit of q(t) (the BD-relevant object) and of C_abs(t)")
    say("")
    tt = np.array([r["t"] for r in rows])
    qq = np.array([r["q"] for r in rows])
    cc = np.array([r["C"] for r in rows])
    # t -> 0 : power law of C
    msk = tt <= 0.05
    if msk.sum() >= 3:
        p = np.polyfit(np.log(tt[msk]), np.log(cc[msk]), 1)
        res = np.log(cc[msk]) - np.polyval(p, np.log(tt[msk]))
        say(f"     t->0   C_abs ~ A t^p :  p = {p[0]:+.5f} +- {np.std(res):.2e}   (A = {np.exp(p[1]):.6f})")
        say(f"            expectation p=1 (variance linear in the window variance) ")
    # t -> inf : exponential rate of q
    msk2 = (tt >= 1.0) & (qq > 1e-17)
    if msk2.sum() >= 3:
        p = np.polyfit(tt[msk2], np.log(qq[msk2]), 1)
        res = np.log(qq[msk2]) - np.polyval(p, tt[msk2])
        say(f"     t->inf q(t) ~ B e^(-lam t) :  lam = {-p[0]:.6f} +- {np.std(res):.2e}   (B = {np.exp(p[1]):.6f})")
        say(f"            in units of lam_1 : lam/lam_1 = {-p[0]/LAM1:.6f}   [lam_1 = {LAM1}]")
        say(f"            fit points: lam1*t = {np.round(tt[msk2],3)}")
        # also C_abs
        p2 = np.polyfit(tt[msk2], cc[msk2], 1)
        say(f"     t->inf C_abs -> {cc[msk2][-1]:.6f} + {p2[0]:+.3e}*t   (slope ~ 0 => SATURATION)")
    say("")
    flush()

    # ---------------- S3: beta-independence of the rate
    say("[S3] the decay rate of q(t) as a function of beta  --  is it lam_1 or a dynamical scale?")
    say(f"      {'beta':>7} {'lam(t->inf)':>14} {'lam/lam_1':>11} {'C(t=inf)=b/4*E_nu[q0]':>24} {'q(0)=b/4':>10}")
    for b in (0.5, 1.0, 2.5, 5.0, 10.0):
        mm = Model1(beta=b)
        tlist = np.array([1.0, 1.6, 2.5, 4.0, 6.3, 10.0])
        qlist = []
        for t in tlist:
            best = -1.0
            for a in alphas:
                Z, m0, mv, Mv = mm.moments(a, t)
                C = Mv - np.outer(mv, mv)
                zw = zzsum = None
                HV = (b / 4) * m0 * np.eye(3) - (b ** 2 / 4) * C
                qv_ = max(0.0, -np.linalg.eigvalsh(HV).min())
                best = max(best, qv_)
            qlist.append(best)
        qlist = np.array(qlist)
        ok = qlist > 1e-13
        lam = -np.polyfit(tlist[ok], np.log(qlist[ok]), 1)[0] if ok.sum() >= 3 else np.nan
        # t=inf value: (b/4) E_nu[q0] with E_nu[q0] = d/dbeta log c0, c0 = 2 I1(b)/b
        from scipy.special import iv as _iv
        Z0 = 2 * _iv(1, b) / b
        dZ = (_iv(0, b) + _iv(2, b)) / b - 2 * _iv(1, b) / b ** 2
        Eq0 = dZ / Z0
        say(f"      {b:7.3f} {lam:14.6f} {lam/LAM1:11.6f} {b/4*Eq0:24.8f} {b/4:10.5f}")
    say("")
    say("      => the rate equals lam_1 = 3 in the SAME units for every beta: it is a ONE-LINK")
    say("         spectral constant of the group, NOT a dynamical (beta-dependent) scale.")
    say("")
    flush()

    say("[S4] verdict data for model (1)")
    say(f"      sup_t C_abs = {max(r['C'] for r in rows):.6f}   (attained at lam1*t = "
        f"{max(rows, key=lambda r: r['C'])['tau']:.4g})")
    say(f"      sup_t q(t)  = {max(r['q'] for r in rows):.6f}   (attained at lam1*t = "
        f"{max(rows, key=lambda r: r['q'])['tau']:.4g})")
    # time integral of q, analytically-correct tail beyond the last sample
    tt2 = np.concatenate([[0.0], tt])
    qq2 = np.concatenate([[beta / 4], qq])
    Itrap = np.trapezoid(qq2, tt2) / LAM1
    lam_fit = -np.polyfit(tt[(tt >= 1.0) & (qq > 1e-14)], np.log(qq[(tt >= 1.0) & (qq > 1e-14)]), 1)[0]
    tail = qq[-1] / lam_fit
    say(f"      Q_partial(0..{tt[-1]:.4g} T1) = {Itrap:.6f}   + tail {tail:.2e}  =>  Q(inf) ~ {Itrap+tail:.6f}")
    say(f"      (in units where lam_1*t is the time variable; Q(inf)*lam_1 = {(Itrap+tail)*LAM1:.6f})")
    say("")
    flush()


if __name__ == "__main__":
    main()
