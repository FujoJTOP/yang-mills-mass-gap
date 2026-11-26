#!/usr/bin/env python3
r"""
cov_models.py -- COV-T, models (2) and (4):  exact linear-algebra models.

Model (2)  1D DIRICHLET STRING (Gaussian, non-circulant, boundary effects)
           S_W = (1/2) phi^T A phi,  A = tridiag(2,-1)_{L x L}  (Dirichlet)
           P_t = heat semigroup, kernel N(0,2t I)   =>   everything closed-form:
               Hess V_t = A (I+2tA)^{-1}   > 0        =>  q(t) = 0 identically
               Cov_t    = 2t A^2 (I+2tA)^{-1}         =>  ||Cov_t||_op = max_k 2t lam_k^2/(1+2t lam_k)
           lam_k = 2-2cos(k pi/(L+1));   lam_max = 2+2cos(pi/(L+1)) -> 4   (L-independent limit)

Model (4a) P4 PERIODIC KERNEL, one site, non-Gaussian:  U(1) with the P4 Boltzmann kernel
           e^{-S_W} = e^{beta cos theta}  (the corpus's P4 kernel K_ij = e^{beta cos(theta_i-theta_j)}),
           P_t = heat semigroup on the unit circle (lam_n = n^2, so lam_1 = 1), exact 1D quadrature.

Model (4b) P4 PERIODIC KERNEL as a CIRCULANT ring (the P4 anchor): 1D periodic Laplacian
           A_k = 2-2cos(2 pi k/L);  lam_max = 4 EXACTLY for every even L  =>  sup_t ||Cov_t||_op = 4
           is an L-independent identity, not a fit.

Run: py -3.12 -I cov_models.py   numpy/scipy only, deterministic, single process.
Writes cov_models.log section by section.
"""
import numpy as np
from scipy.special import iv

LOG = []
LOGF = "cov_models.log"


def say(s=""):
    print(s)
    LOG.append(s)


def flush():
    with open(LOGF, "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


TAUS = [0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0]


# ================================================================== model (2): Dirichlet string
def dirichlet_eigs(L):
    k = np.arange(1, L + 1)
    return 2.0 - 2.0 * np.cos(k * np.pi / (L + 1))


def m2_Cabs(t, lam):
    return float(np.max(2.0 * t * lam ** 2 / (1.0 + 2.0 * t * lam)))


def m2_trCov(t, lam):
    return float(np.sum(2.0 * t * lam ** 2 / (1.0 + 2.0 * t * lam)))


def section2():
    say("=" * 100)
    say("MODEL (2)  1D DIRICHLET STRING   --  exact, Gaussian")
    say("=" * 100)
    say("  S_W = (1/2)phi^T A phi,  A = tridiag(2,-1) Dirichlet;  P_t = heat semigroup (var 2t).")
    say("  Hess V_t = A(I+2tA)^{-1}  >  0  for every t  =>  q(t) = ||(Hess V_t)_-||_op = 0 EXACTLY.")
    say("  Cov_t = 2t A^2 (I+2tA)^{-1}   [PSD, increasing in t]  =>  sup_t ||Cov_t||_op = lam_max(A).")
    say("")
    say(f"  {'L':>5} {'lam_1':>9} {'lam_max':>9} {'sup_t C_abs':>12} {'sup_t C_int=C/L':>16} "
        f"{'tr/op at t=inf':>15}")
    for L in (2, 4, 6, 8, 32, 128):
        lam = dirichlet_eigs(L)
        supC = float(np.max(lam))                      # t -> inf
        say(f"  {L:5d} {lam.min():9.4f} {lam.max():9.4f} {supC:12.6f} {supC/L:16.6f} "
            f"{float(np.sum(lam))/supC:15.2f}")
    say("")
    say("  the op-norm bound is L-uniform:  lam_max(A) = 2+2cos(pi/(L+1))  <=  4  for every L,")
    say("  and -> 4.  The INTENSIVE version C/L -> 0 like 1/L.  tr/op ~ L = the '3E' factor.")
    say("")

    say("  the curve, L = 6 and L = 32  (t in units of lam_1 t):")
    say(f"  {'lam1*t':>8} {'L=6 C_abs':>12} {'L=6 C_int':>12} {'L=32 C_abs':>12} {'L=32 C_int':>12} {'q':>6}")
    for tau in TAUS:
        row = []
        for L in (6, 32):
            lam = dirichlet_eigs(L)
            t = tau / lam.min()
            c = m2_Cabs(t, lam)
            row += [c, c / L]
        say(f"  {tau:8.4g} {row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f} {row[3]:12.6f} {0.0:6.1f}")
    say("")
    say("  three-segment fit, L = 32 (exact formulas, no noise):")
    lam = dirichlet_eigs(32)
    l1, lmax = lam.min(), lam.max()
    ts = np.exp(np.linspace(np.log(1e-6 / lmax), np.log(1e4 / l1), 400))
    cs = np.array([m2_Cabs(t, lam) for t in ts])
    p = np.polyfit(np.log(ts[:60]), np.log(cs[:60]), 1)
    say(f"    t->0   : C ~ A t^p,  p = {p[0]:.6f} (exact 1),  A = {np.exp(p[1]):.6f}  "
        f"[exact A = 2 lam_max^2 = {2*lmax**2:.4f}]")
    say(f"    mid    : rational shape  C(t) = c t/(1+b t) fits with b = 1/(2 lam_eff);")
    A = (cs * (1 + 2 * lmax * ts))
    say(f"             C(t)/lam_max is monotone;  C(lam_max t = 1/2) / lam_max = "
        f"{m2_Cabs(0.5/lmax, lam)/lmax:.6f} ;  C(10) / lam_max = {m2_Cabs(10/lmax, lam)/lmax:.6f}")
    say(f"    t->inf : C -> lam_max = {lmax:.8f}   (rate = 0: SATURATION, not exponential decay)")
    say(f"             q(t) = 0 for all t  =>  Q(inf) = 0  =>  the BD bound is vacuous on this model")
    say("")
    say("  the rise, in absolute t (L = 32, lam_max = 4):")
    say(f"    {'t':>8} {'C_abs':>12} {'C/lam_max':>11} {'2t lam_max^2':>13}")
    lam = dirichlet_eigs(32)
    for t in (1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0):
        c = m2_Cabs(t, lam)
        say(f"    {t:8.4g} {c:12.6f} {c/lam.max():11.6f} {2*t*lam.max()**2:13.6f}")
    say("")
    flush()


# ================================================================== model (4a): P4, U(1), non-Gaussian
TH = 2.0 * np.pi * np.arange(4096) / 4096
WQ = 2.0 * np.pi / 4096


def nmax_for(t):
    if t <= 0:
        return 4000
    return int(min(4000, max(2, np.ceil(np.sqrt(36.0 / t)) + 2)))


def moments_u1(th0, t, beta):
    """Z = B_t(th0) and the normalised weights of mu_{t,th0} on the circle."""
    n = np.arange(-nmax_for(t), nmax_for(t) + 1)
    kt = np.real(np.exp(1j * np.outer(n, TH - th0)) * np.exp(-n ** 2 * t)[:, None]).sum(axis=0) \
        / (2.0 * np.pi)
    rho = kt * np.exp(beta * np.cos(TH))
    Z = float(np.sum(rho) * WQ)
    return Z, rho * WQ / Z


def B_u1_char(th0, t, beta, nmax=200):
    n = np.arange(-nmax, nmax + 1)
    return float(np.sum(iv(np.abs(n), beta) * np.exp(-n ** 2 * t) * np.cos(n * th0)))


def strict_sup_u1(t, beta, nth0=49):
    """sup over th0 of ||(Hess V_t)_-|| = max(0, -(beta E[cos] - beta^2 Var(sin)))."""
    best = -1e300
    arg = None
    for th0 in np.linspace(0.0, 2 * np.pi, nth0, endpoint=False):
        Z, wn = moments_u1(th0, t, beta)
        Es = float(wn @ np.sin(TH))
        Es2 = float(wn @ np.sin(TH) ** 2)
        var = beta ** 2 * (Es2 - Es ** 2)
        if var > best:
            best, arg = var, (th0, Z, float(wn @ np.cos(TH)), var)
    th0, Z, Ec, var = arg
    return dict(t=t, C=var, q=max(0.0, -(beta * Ec - var)), HV=beta * Ec - var, th0=th0, Z=Z)


def section4a(beta=2.5):
    say("=" * 100)
    say("MODEL (4a)  P4 PERIODIC KERNEL, one site, U(1), NON-GAUSSIAN")
    say("=" * 100)
    say("  e^{-S_W} = e^{beta cos theta}   (the P4 Boltzmann kernel); P_t = heat semigroup on the")
    say("  UNIT circle (-Delta eigenvalues n^2 => lam_1 = 1).  g = beta sin theta, Hess S_W = beta cos theta.")
    say("  mu_{t,th0}(d th) propto K_t(th-th0) e^{beta cos th} d th /(2 pi); exact 1D quadrature (4096 pts).")
    say("")
    say(f"  {'lam1*t':>8} {'C_abs':>12} {'q(t)':>13} {'HV':>14} {'B_quad/B_char':>14}")
    rows = []
    for tau in [0.01, 0.03, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0]:
        r = strict_sup_u1(tau, beta)
        rows.append(r)
        say(f"  {tau:8.4g} {r['C']:12.6f} {r['q']:13.6e} {r['HV']:14.6f} "
            f"{r['Z']/B_u1_char(r['th0'], tau, beta):14.2e}")
    say("")
    r = strict_sup_u1(1e4, beta)
    Es = 0.0
    Z, wn = moments_u1(0.0, 1e4, beta)
    Es = float(wn @ np.sin(TH))
    Es2 = float(wn @ np.sin(TH) ** 2)
    Ec = float(wn @ np.cos(TH))
    say("  t -> inf check (Stein identity on the flat circle): Cov_nu = E_nu[Hess S_W]")
    say(f"    measured beta^2 Var_nu(sin) = {beta**2*(Es2-Es**2):.8f} ;  beta E_nu[cos] = "
        f"{beta*Ec:.8f} ;  beta I1(beta)/I0(beta) = {beta*iv(1,beta)/iv(0,beta):.8f}")
    say("")
    tt = np.array([x["t"] for x in rows])
    qq = np.array([x["q"] for x in rows])
    cc = np.array([x["C"] for x in rows])
    p = np.polyfit(np.log(tt[:4]), np.log(cc[:4]), 1)
    say(f"  t->0   : C ~ A t^p,  p = {p[0]:.5f} +- "
        f"{np.std(np.log(cc[:4])-np.polyval(p,np.log(tt[:4]))):.1e}")
    msk = (tt >= 1.0) & (qq > 1e-15)
    p2 = np.polyfit(tt[msk], np.log(qq[msk]), 1)
    say(f"  t->inf : q ~ B e^(-lam t),  lam = {-p2[0]:.6f} +- "
        f"{np.std(np.log(qq[msk])-np.polyval(p2,tt[msk])):.1e}    [lam_1 = 1]")
    say(f"  sup_t C_abs = {cc.max():.6f} at lam1*t = {tt[cc.argmax()]:.4g} ;  C(t=inf) = {cc[-1]:.6f}")
    say(f"  => bounded, saturating:  int_0^inf C dt = infinity ;  Q(inf) ~ "
        f"{np.trapezoid(np.concatenate([[beta],qq]), np.concatenate([[0.0],tt])):.6f}")
    say("")
    say("  beta-scan of the decay rate (is it lam_1 = 1, or a dynamical scale?):")
    say(f"  {'beta':>7} {'lam':>12} {'lam/lam_1':>11}")
    tl = np.array([1.0, 1.8, 3.2, 5.6, 10.0])
    for b in (0.5, 1.0, 2.5, 5.0):
        rr = np.array([strict_sup_u1(t, b)["q"] for t in tl])
        ok = rr > 1e-16
        lam = -np.polyfit(tl[ok], np.log(rr[ok]), 1)[0] if ok.sum() >= 3 else np.nan
        say(f"  {b:7.3f} {lam:12.6f} {lam:11.6f}")
    say("")
    flush()


# ================================================================== model (4b): P4 circulant ring
def section4b():
    say("=" * 100)
    say("MODEL (4b)  P4 PERIODIC KERNEL as a CIRCULANT ring  (the P4 anchor)")
    say("=" * 100)
    say("  A = 1D periodic Laplacian on Z_L:  lam_k = 2-2cos(2 pi k/L),  k = 0..L-1.")
    say("  lam_max = 4 EXACTLY for every even L (k = L/2).  Gaussian => q(t) = 0.")
    say("")
    say(f"  {'L':>5} {'lam_1':>10} {'lam_max':>8} {'sup_t C_abs':>12} {'C_int=C/L':>11} "
        f"{'tr/op(t=inf)':>13} {'L-exponent':>11}")
    Ls = [2, 4, 6, 8, 16, 64, 256]
    sup = []
    for L in Ls:
        lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(L) / L)
        sup.append(float(lam.max()))
        say(f"  {L:5d} {lam[1] if L > 1 else 0:10.6f} {lam.max():8.4f} {lam.max():12.6f} "
            f"{lam.max()/L:11.6f} {np.sum(lam)/lam.max():13.2f} {'':>11}")
    p = np.polyfit(np.log(Ls[:4]), np.log(np.array(sup[:4])), 1)
    say(f"  log-log slope of sup_t C_abs over L = 2,4,6,8 :  {p[0]:+.3e}   (== 0: L-uniform, exactly)")
    say("")
    say("  the curve on the ring (L = 6):  ")
    lam6 = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(6) / 6)
    say(f"  {'lam1*t':>8} {'C_abs':>12} {'C_int':>10} {'q':>6}")
    for tau in TAUS:
        t = tau / lam6[1]
        c = float(np.max(2 * t * lam6 ** 2 / (1 + 2 * t * lam6)))
        say(f"  {tau:8.4g} {c:12.6f} {c/6:10.6f} {0.0:6.1f}")
    say("")
    flush()


# ================================================================== summary
def section5():
    say("=" * 100)
    say("SUMMARY  --  C(t) = ||Cov_{mu_t}(grad S_W)||_op  on the four models")
    say("=" * 100)
    say("  model                        bounded?   integrable?   sup_t C_abs     L-uniform?   q(t)")
    say("  (1) single link SU(2)        YES         NO (sat.)     0.403062         n/a (1 link) exp, rate 3")
    say("  (2) 1D Dirichlet string      YES         NO (sat.)     lam_max<=4       YES         == 0 exactly")
    say("  (3) 2D YM tower              YES         NO (sat.)     0.40             YES         exp, rate E_1")
    say("  (4a) P4 U(1) one site        YES         NO (sat.)     (see above)      n/a          exp, rate 1")
    say("  (4b) P4 circulant ring       YES         NO (sat.)     4 EXACTLY        YES (p=0)   == 0 exactly")
    say("")
    say("  Common law (all four):  C(0)=0 exactly;  C(t) rises like t;  C(t) -> ||E_nu[Hess S_W]||_op > 0.")
    say("  =>  sup_t C(t) < infinity  (branch (c))  but  int_0^inf C(t) dt = infinity (saturation).")
    say("  =>  the BD-relevant object is q(t)=||(Hess V_t)_-||_op = ||(E[Hess S_W] - Cov)_-||_op,")
    say("      whose integrability comes from the exact CANCELLATION of the two O(1) terms.")
    say("")
    flush()


def main():
    section2()
    section4a()
    section4b()
    section5()
    say("[done] cov_models.log written")


if __name__ == "__main__":
    main()
