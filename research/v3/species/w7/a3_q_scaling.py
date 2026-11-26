#!/usr/bin/env python3
r"""
a3_q_scaling.py -- A3-Q-SCALING: what happens to Q(inf) = int_0^inf q(t) dt
                    under a -> 0 (fixed physical quantities)?

Question (task):  q(t) = sup_phi || (Hess V_t)_- ||,  Q(inf) = int_0^inf q dt.
                  Is Q constant / -> 0 / divergent as the lattice spacing a -> 0 ?

Three jobs, in this order:
  (J1) UNIT CONFIRMATION.  In which unit is the measured t?  The model-1 tower has
       lam_l = l(l+2) on the UNIT round S^3 <=> lam_1 = 3.  Recompute q(tau) for
       beta = 2.5 on the COV-T table-T2 grid and settle, by arithmetic, whether
       Q = 0.109856 is int q dt (raw group-metric time) or int q d(tau) (tau = lam_1 t).
  (J2) SCALING.  t -> a^2 t (heat-kernel convention) and t -> a t, with q FROZEN
       (naive) and with q TRANSPORTED with the metric (matched pair).  Fit Q(a) = Q(1) a^p.
  (J3) DIMENSION BOOKKEEPING + BETA SCAN.  S_W is dimensionless, Hess S_W (group
       coordinates) is dimensionless, q is dimensionless, so Q carries one power of
       [t] => e^{-2Q} needs a fixed reference time t_*.  Identify t_* (a-independent?).
       Then measure Q(beta): if Q ~ c beta, the continuum trajectory (beta -> inf as
       a -> 0 at fixed physics) is where Q actually moves.

Model (identical to cov_link.py = COV-T sec.3, exact, no new species):
    S_W(w) = beta (1 - q_0),   q_0 = cos psi,   w in SU(2) = unit S^3
    mu_{t}(dw) propto K_t(w) e^{-S_W(w)} dw   (x = identity; everything is a class function)
    Hess V_t = (beta/4) E_mu[q_0] 1_3 - (beta^2/4) Cov_mu(q_vec)
             = [ (beta/4) E_mu[q_0] - (beta^2/12) E_mu[sin^2 psi] ] 1_3      (isotropy)
    q(t) = max(0, -that scalar)          <- exact 1-D reduction, no 96^3 grid needed

Class-function route (exact, cheap).  With chi_l orthonormal in L^2(Haar),
    e^{-S_W} = sum_l d_l chi_l,   d_l = e^{-beta} (I_l(beta) - I_{l+2}(beta))
    cos psi e^{-S_W} = sum_l p_l chi_l,   sin^2 psi e^{-S_W} = sum_l s_l chi_l
    d_l = e^{-beta} a_l,  p_l = e^{-beta} a_l'(beta),  s_l = d_l - e^{-beta} a_l''(beta)
    a_l(beta) = int e^{beta cos psi} chi_l dU = 2(l+1) I_{l+1}(beta)/beta
    B_t = sum_l d_l (l+1) e^{-lam_l t},   E_mu[f] = sum_l f_l (l+1) e^{-lam_l t} / B_t   (x = id)
Verified below against (i) the closed form B_t = e^{-beta} sum c_l e^{-lam t} chi_l(psi_x),
(ii) direct psi-quadrature, (iii) COV-T table T2.

Run: py -3.12 -I a3_q_scaling.py     (numpy/scipy, deterministic, single process, < 30 s)
Writes a3_q_scaling.log.  No git.
"""
import numpy as np
from scipy.special import iv

LOG = []
LOGF = "a3_q_scaling.log"


def say(s=""):
    print(s)
    LOG.append(str(s))


def flush():
    with open(LOGF, "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


LAM1 = 3.0                     # -Delta on the UNIT round S^3, l = 1
LMAX = 1400                    # enough for t_min = 1e-4  (l ~ sqrt(30/t) ~ 950)


# ---------------------------------------------------------------- exact coefficients
def chi(l, psi):
    """chi_l(psi) = sin((l+1)psi)/sin(psi), vectorised in psi."""
    psi = np.asarray(psi, dtype=float)
    s = np.sin(psi)
    out = np.empty_like(psi)
    small = np.abs(s) < 1e-12
    out[~small] = np.sin((l + 1) * psi[~small]) / s[~small]
    ll = np.where(np.cos(psi[small]) > 0, 1.0, (-1.0) ** l)
    out[small] = (l + 1) * ll
    return out


def a_of(beta, ls):
    """a_l(beta) = 2(l+1) I_{l+1}/beta : character coefficients of e^{beta cos psi}."""
    ls = np.asarray(ls)
    return 2.0 * (ls + 1.0) * iv(ls + 1.0, beta) / beta


def coeffs(beta, L=LMAX, x_is_antipode=True):
    """chi-coefficients of  w,  cos psi w,  sin^2 psi w   with  w = e^{-S_W(x z)}.

    x = identity  :  w = e^{-beta(1-cos psi)}
    x = antipode  :  w = e^{-beta(1+cos psi)}   (COV-T table T2 is taken here: q(0)=beta/4)
    Closed form via a_l(beta) = 2(l+1) I_{l+1}(beta)/beta  (character coefficient of
    e^{beta cos psi});  for the antipode use beta -> -beta, i.e. a_l(-beta) = (-1)^l a_l.
    """
    ls = np.arange(L + 1, dtype=float)
    h = 1e-4 * max(1.0, beta)
    a0 = a_of(beta, ls)
    ap = a_of(beta + h, ls)
    am = a_of(beta - h, ls)
    a1 = (ap - am) / (2.0 * h)
    a2 = (ap - 2.0 * a0 + am) / (h * h)
    eb = np.exp(-beta)
    sgn = (-1.0) ** ls if x_is_antipode else np.ones_like(ls)
    d = eb * sgn * a0                      # coefficient of w
    e = -eb * sgn * a1                     # coefficient of cos psi * w
    f = eb * sgn * (a0 - a2)               # coefficient of sin^2 psi * w
    return d, e, f


def lam(ls):
    return ls * (ls + 2.0)


def lcut(t, L=LMAX):
    """modes needed for exp(-lam_l t) < 1e-30."""
    if t <= 0:
        return L
    return int(min(L, np.ceil(np.sqrt(35.0 / t)) + 20))


def moments(t, beta, coef):
    """(B_t, E_mu[q_0], E_mu[sin^2 psi]) at the chosen x, exact class-function sums."""
    d, e, f = coef
    L = lcut(t)
    ls = np.arange(L + 1, dtype=float)
    w = (ls + 1.0) * np.exp(-lam(ls) * t)
    B = float(np.dot(d[:L + 1], w))
    if not np.isfinite(B) or B <= 0.0:
        return np.nan, np.nan, np.nan          # alternating Bessel sum lost precision (large beta)
    eq0 = float(np.dot(e[:L + 1], w)) / B
    es2 = float(np.dot(f[:L + 1], w)) / B
    return B, eq0, es2


def q_of_t(t, beta, coef):
    """q(t) = max(0, -lambda_min Hess V_t), t in RAW group-metric units (lam_1 = 3).

    x = antipode, so q_0(xz) = -q_0(z) and
      Hess V_t = (beta/4) E_mu[q_0(xz)] 1_3 - (beta^2/12) E_mu[sin^2 psi] 1_3
               = -[ (beta/4) E_mu[q_0(z)] + (beta^2/12) E_mu[sin^2 psi] ] 1_3
    => q(t) = max(0, (beta/4) E_mu[q_0(z)] + (beta^2/12) E_mu[sin^2 psi]).
    Sign pinned by  q(0+) = beta/4  (E_mu = delta_{z=1}: E[q_0(z)] = 1, E[sin^2] = 0).
    """
    _, eq0, es2 = moments(t, beta, coef)
    neg = 0.25 * beta * eq0 + (beta * beta / 12.0) * es2      # = -lambda_min Hess V_t
    return max(0.0, neg), neg


def q_of_tau(tau, beta, coef):
    """same, parametrised by tau = lam_1 t (dimensionless group-spectral time)."""
    return q_of_t(tau / LAM1, beta, coef)


# ---------------------------------------------------------------- checks
def sect(title):
    say()
    say("=" * 78)
    say(title)
    say("=" * 78)


def check_closed_form(beta, coef, taus=(0.001, 0.1, 1.0, 10.0)):
    """B_t from the l-sum vs COV-T's closed form (x = antipode, chi_l(pi) = (-1)^l (l+1))
    vs a direct 1-D psi-quadrature of  int K_t(psi) e^{-beta(1+cos psi)} dHaar."""
    say("  [check A] B_t:  l-sum  vs  COV-T closed form  vs  1-D psi-quadrature")
    xg, wg = np.polynomial.legendre.leggauss(2400)
    psi = 0.5 * np.pi * (xg + 1.0)
    wpsi = 0.5 * np.pi * wg * (2.0 / np.pi) * np.sin(psi) ** 2   # Haar: (2/pi) sin^2 psi dpsi
    ls = np.arange(1201, dtype=float)
    c = iv(ls, beta) - iv(ls + 2.0, beta)                        # = 2(l+1) I_{l+1}/beta
    say(f"    recurrence check  max|c_l - 2(l+1)I_(l+1)/beta| = "
        f"{np.max(np.abs(c - 2*(ls+1)*iv(ls+1,beta)/beta)):.2e}")
    for tau in taus:
        t = tau / LAM1
        B_sum, _, _ = moments(t, beta, coef)
        B_cl = np.exp(-beta) * np.sum(c * (-1.0) ** ls * np.exp(-lam(ls) * t) * (ls + 1.0))
        K = np.zeros_like(psi)
        Lq = lcut(t)
        for l in range(Lq + 1):
            K += (l + 1.0) * np.exp(-l * (l + 2.0) * t) * chi(l, psi)
        B_q = float(np.dot(wpsi, K * np.exp(-beta * (1.0 + np.cos(psi)))))
        say(f"    tau={tau:g}:  l-sum={B_sum:.15e}  closed={B_cl:.15e}  quad={B_q:.15e}"
            f"   relerr(sum,quad)={abs(B_sum-B_q)/B_q:.2e}")


def check_table_t2(beta=2.5):
    """COV-T table T2, beta = 2.5:  q(tau) at the tabulated taus."""
    coef = coeffs(beta)
    tab = [(0.001, 0.62542), (0.01, 0.62913), (0.1, 0.66105), (0.316, 0.66686),
           (1.0, 0.31493), (3.16, 0.02326), (10.0, 2.30e-05)]
    say(f"  [check B] q(tau) vs COV-T table T2 (beta={beta}), and q(0) = beta/4 = {beta/4:g}")
    say("    tau        q(this)        q(T2)      ratio")
    for tau, qtab in tab:
        q, _ = q_of_tau(tau, beta, coef)
        say(f"    {tau:<9g}  {q:.8e}   {qtab:.8e}   {q/qtab:.4f}")
    say(f"    q(tau->0) = {q_of_tau(1e-9, beta, coef)[0]:.12f}   (beta/4 = {beta/4:.12f})")
    return coef


def q_curve(beta, n=900, tau_min=1e-4, tau_max=40.0):
    """q on a log grid in tau (dimensionless), tau = lam_1 t."""
    coef = coeffs(beta)
    taus = np.exp(np.linspace(np.log(tau_min), np.log(tau_max), n))
    q = np.array([q_of_tau(x, beta, coef)[0] for x in taus])
    return taus, q, coef


def Q_of_beta(beta, **kw):
    """Q = int_0^inf q dt  (RAW t units) and  A = int q dtau = lam_1 Q  (dimensionless)."""
    taus, q, coef = q_curve(beta, **kw)
    if not np.all(np.isfinite(q)):
        return np.nan, np.nan, q, taus
    A = float(np.trapezoid(q, taus))          # tau units
    q0, _ = q_of_tau(1e-12, beta, coef)
    A += q0 * taus[0]                         # boundary layer tau < tau_min:  q(0) * tau_min
    Qraw = A / LAM1
    return Qraw, A, q, taus


def reproduce_cov_t_Q(beta=2.5):
    """Re-run cov_link.py S4 verbatim (13 log-spaced tau nodes, trapezoid over raw t, then
    DIVIDED by lam_1) to explain where its 0.109856 comes from."""
    taus13 = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6,
                       100.0, 316.0, 1000.0])
    coef = coeffs(beta)
    tt = taus13 / LAM1
    qq = np.array([q_of_tau(x, beta, coef)[0] for x in taus13])
    tt2 = np.concatenate([[0.0], tt])
    qq2 = np.concatenate([[beta / 4], qq])
    Itrap = float(np.trapezoid(qq2, tt2))
    say("  [reproduction of cov_link.py S4, beta = 2.5]")
    say(f"    trapz(q, t) on the 13 log nodes        = {Itrap:.6f}   (raw t units)")
    say(f"    cov_link.py then divides by lam_1      = {Itrap/LAM1:.6f}   <-- its 'Q(inf)'")
    say(f"    correct: divide not needed            = {Itrap:.6f}")
    say(f"    correct dimensionless  lam_1 * intqdt = {Itrap*LAM1:.6f}")
    say(f"    dense grid (900 pts): int q dt = {Q_of_beta(beta)[0]:.6f}"
        f"   ({Q_of_beta(beta)[0]/Itrap:.3f}x the 13-node trapezoid)")
    return Itrap / LAM1


# ---------------------------------------------------------------- J1: units
def job1():
    sect("J1  UNIT CONFIRMATION (which unit is the measured t?)")
    beta = 2.5
    coef = check_table_t2(beta)
    check_closed_form(beta, coef)

    Qraw, A, q, taus = Q_of_beta(beta)
    say()
    say("  Integrals of the SAME curve, two unit conventions (beta = 2.5):")
    say(f"    q(0)                                  = {q[0]:.8f}")
    say(f"    int q d(tau),  tau = lam_1 t          = {A:.8f}        <- dimensionless")
    say(f"    int q dt,      t raw (lam_1 = 3)      = {Qraw:.8f}        <- carries [t]")
    say(f"    lam_1 * (int q dt)                    = {LAM1*Qraw:.8f}")
    say()
    say("    COV-T sec.3 reports Q(inf) = 0.109856 '(lam_1 t time unit)' and lam_1 Q = 0.32957.")
    say("    Recompute:  int q dtau = %.6f ;  int q dt = %.6f ;  lam_1*int q dtau = %.6f"
        % (A, Qraw, LAM1 * A))
    say("    => neither 0.109856 nor 0.32957 is reproduced by this curve. See E1 below.")
    reproduce_cov_t_Q(beta)
    # local rate in raw t units:  -dlog q/dt  at the largest usable tau
    for tau in (8.0, 10.0, 12.0):
        eps = 1e-3
        q1, _ = q_of_tau(tau - eps, beta, coef)
        q2, _ = q_of_tau(tau + eps, beta, coef)
        rate_t = -(np.log(q2) - np.log(q1)) / (2 * eps / LAM1)
        say(f"    local rate -dlog q/dt at tau={tau:g}: {rate_t:.6f}   (lam_1 = 3, lam/lam_1 ="
            f" {rate_t/LAM1:.6f})")
    return Qraw, A


# ---------------------------------------------------------------- J2: scaling
def job2(Q1_raw, A1):
    sect("J2  SCALING  t -> a^2 t  and  t -> a t ;  fit Q(a) = Q(1) a^p")
    say("  Notation: g = bi-invariant metric on the link group, unit round S^3 (lam_1 = 3).")
    say("            q_g = ||(Hess_g V)_-||  (dimensionless in group coordinates)")
    say("            t_g = heat-kernel time of g     => Q_g = int q_g dt_g  carries [t_g]")
    say()
    say("  (i) NAIVE (q frozen): the ONLY thing that changes is the time unit.")
    a = np.array([2.0 ** -k for k in range(1, 13)])          # a = 1/2 ... 1/4096
    say("      a            Q_a2 = a^2 Q_g        Q_a1 = a Q_g          (Q_g = %.8f)"
        % Q1_raw)
    for ai in a:
        say(f"      {ai:.8e}   {ai**2*Q1_raw:.8e}      {ai*Q1_raw:.8e}")
    p2 = np.polyfit(np.log(a), np.log(a ** 2 * Q1_raw), 1)[0]
    p1 = np.polyfit(np.log(a), np.log(a * Q1_raw), 1)[0]
    say(f"      log-log fit:  p(t ~ a^2 t) = {p2:.12f}        p(t ~ a t) = {p1:.12f}")
    say()
    say("  (ii) MATCHED (q transported with the metric):  g -> g' = a^2 g  gives")
    say("       lam_l -> lam_l/a^2,  Hess_{g'} = a^{-2} Hess_g  =>  q_{g'} = a^{-2} q_g,")
    say("       and the time unit must move the other way: t_{g'} = a^2 t_g  (so that t Delta")
    say("       is invariant).  Product:")
    say("         int q_{g'} dt_{g'} = int (a^{-2} q_g)(a^2 dt_g) = int q_g dt_g   =>  p = 0.")
    # numerical no-op verification of the matched pair on the a^2 grid
    say("       numerical check (substitute tau -> a^2 tau inside the integral):")
    cf25 = coeffs(2.5)
    for ai in (0.5, 0.25, 1.0 / 64):
        # rescaled model:  lam_l -> lam_l/a^2,  q_g'(t) = a^{-2} q_g(t/a^2)
        # (t = the rescaled model's own raw time; its support sits at t ~ a^2/lam_1)
        tt = np.exp(np.linspace(np.log(1e-9), np.log(16.0), 900))
        qv = np.array([(1.0 / ai ** 2) * q_of_t(x / ai ** 2, 2.5, cf25)[0] for x in tt])
        Qm = float(np.trapezoid(qv, tt))
        say(f"         a={ai:<10g}  int q_g' dt_g' = {Qm:.8f}   (Q_g = {Q1_raw:.8f},"
            f" rel diff {abs(Qm-Q1_raw)/Q1_raw:.1e})")
    say()
    say("  => p = +2 (t ~ a^2, frozen q) ; p = +1 (t ~ a, frozen q) ; p = 0 in BOTH")
    say("     conventions once the metric that defines Hess is transported with t.")
    return p2, p1


# ---------------------------------------------------------------- J3: dimension + beta
def job3():
    sect("J3  DIMENSION BOOKKEEPING AND THE ONE PARAMETER THAT ACTUALLY MOVES Q")
    say("  * S_W = beta sum_p (1 - Re tr U_p/2) is a PURE NUMBER: no a, no length.")
    say("  * group coordinates (psi, nhat) are dimensionless => Hess S_W is dimensionless")
    say("    => q is dimensionless => Q = int q dt carries exactly ONE power of [t].")
    say("  * therefore e^{-2Q} needs Q/t_* with t_* FIXED.  In the model:")
    say(f"      lam_1 = 3 exactly (unit round S^3, beta-independent)  =>  t_* = 1/lam_1 = {1/LAM1:.6f}")
    say("      (BD's own normalisation Ric = (N/4) g rescales g by 8/3: lam_1' = 9/8,")
    say("       t_*' = 8/9.  Either way t_* is a PURE NUMBER, a-independent.)")
    say("  * a enters NOWHERE in q(t): q = q(beta, t) only.  Checked structurally: the")
    say("    entire computation of Hess V_t uses (beta, lam_l, chi_l) and nothing else.")
    say()
    say("  beta scan: does Q move with beta?  (a -> 0 at fixed physics <=> beta -> inf)")
    say("    beta      q(0)=beta/4    A=int q dtau      Q=int q dt      A/beta     lam_loc(t)")
    betas = [0.5, 1.0, 2.5, 5.0, 10.0, 20.0, 50.0]
    rows = []
    for b in betas:
        Qraw, A, q, taus = Q_of_beta(b, n=700)
        coef = coeffs(b)
        e = 1e-2
        q1, _ = q_of_tau(9.0 - e, b, coef)
        q2, _ = q_of_tau(9.0 + e, b, coef)
        lamloc = -(np.log(q2) - np.log(q1)) / (2 * e) if (q1 > 0 and q2 > 0) else np.nan
        q_check, _ = q_of_tau(1e-4, b, coef)              # must equal beta/4 at tau -> 0
        bad = ((not np.isfinite(A)) or abs(A) > 20.0
               or abs(q_check - b / 4.0) > 0.05 * b / 4.0)
        if bad:                                           # float64 break-down of the l-sum
            A = Qraw = np.nan
        rows.append((b, b / 4.0, A, Qraw, A / b if np.isfinite(A) else np.nan, lamloc))
        say(f"    {b:<9g} {b/4:.6f}        {A:.8f}      {Qraw:.8f}    {A/b:.6f}   {lamloc:.6f}")
    ok = [r for r in rows if np.isfinite(r[2])]
    A_over_beta = np.array([r[4] for r in ok])
    bb = np.array([r[0] for r in ok])
    AA = np.array([r[2] for r in ok])
    pbeta = np.polyfit(np.log(bb), np.log(AA), 1)[0] if ok.__len__() >= 3 else np.nan
    say(f"    A/beta over the reliable window (beta <= {ok[-1][0]:g}): "
        f"min={A_over_beta.min():.6f} max={A_over_beta.max():.6f}  (RISING with beta)")
    say(f"    log-log fit over that window:  A = int q dtau ~ beta^{pbeta:.3f}"
        f"   (linear would be 1.000; >1 = mildly superlinear)")
    say("    (beta >= 20: the alternating Bessel sum for B_t cancels below float64 -> nan)")
    say("    => Q is ~linear in beta with a slowly RISING coefficient: Q = O(beta), not O(1).")
    say()
    say("  CONSEQUENCE for the continuum trajectory (fixed physical quantities):")
    say("    a -> 0 forces beta = 2N/g0^2 -> inf  =>  Q(beta) ~ c beta -> inf  =>")
    say("    e^{-2Q} -> 0  and  gamma >= (N/4) e^{-2Q} -> 0 : the criterion empties out.")
    say("    But the cause is BETA, not a.  At fixed beta, Q is exactly a-independent.")
    return rows


# ---------------------------------------------------------------- main
def main():
    say("A3-Q-SCALING  --  Q(inf) = int_0^inf q(t) dt under a -> 0")
    say("model 1 (single link SU(2), exact): S_W = beta(1-cos psi), q = ||(Hess V_t)_-||")
    say("python: numpy %s / scipy -- single process, deterministic" % np.__version__)
    say()
    say("Reading of the object: q(t) is DIMENSIONLESS (pure-number action, dimensionless")
    say("group coordinates).  t is the group heat-kernel time in the unit-round-S^3 metric.")
    say("Hence Q carries one power of [t] and e^{-2Q} presupposes a normalisation t_* = 1/lam_1.")
    Qraw, A = job1()
    p2, p1 = job2(Qraw, A)
    rows = job3()
    sect("VERDICT (one line)")
    say("  Q is a-INDEPENDENT (p = 0 in every metric-consistent convention; the frozen-q")
    say("  exponents +2 / +1 are the artefacts of rescaling t without the metric that")
    say("  defines Hess).  The normalisation scale is t_* = 1/lam_1 = 1/3 (pure number,")
    say("  fixed by the link-group metric, a-independent).  So the answer to (a)/(b)/(c) is")
    say("  (b): Q constant in a.  The residual beta-dependence Q ~ c beta is the real")
    say("  killer along the fixed-physics trajectory, and it is a DIFFERENT failure mode")
    say("  from (a).")
    flush()


if __name__ == "__main__":
    main()
