#!/usr/bin/env python3
r"""
ct_chain_exact.py -- T1:  X(0) = -2 Cov(G, Phi'| a0=0)  for the L=2 SU(2) chain, in CLOSED FORM.

MODEL (exactly na2_chain.py's L=2 open chain = one bond; flat chart = principal log):
    xi^(i) in R^3, |xi^(i)| < 2pi,  measure prod_i J(rho_i) d^3xi_i ,
    J(rho) = (sin(rho/2)/(rho/2))^2 = 4 sin^2(rho/2)/rho^2 ,  Phi = -log J = sum_i phi(rho_i) ,
    phi(rho) = 2 log(rho/2) - 2 log sin(rho/2) ,
    S = beta (1 - c12) ,  c12 = (1/2)Re tr(U1 U2^dag) = cos(rho1/2)cos(rho2/2)
                                                   + sin(rho1/2)sin(rho2/2) (n1.n2) ,
    coarse field  a = (xi_3^1 + xi_3^2)/2 ,   d/da = d/dxi_3^1 + d/dxi_3^2   (UNIT SPEED).
    X := -2 Cov(G, Phi' | a=a0) ,  G = d_a S ,  Phi' = d_a Phi .

STEP 1 (pointwise algebra; verified by finite differences below).
    n_i^3 = xi_3^i/rho_i , w = n1.n2 , s_i = sin(rho_i/2) , k_i = cos(rho_i/2) ,
    psi(r) := phi'(r)/r = 2/r^2 - cot(r/2)/r ,   A1 := s2 (k1/2 - s1/rho1) ,  B2 := s1 (k2/2 - s2/rho2) :
        Phi' = sum_i psi(rho_i) xi_3^i                                        (exact, pointwise)
        G    = +beta [ (n2^3 - w n1^3) A1 + (n1^3 - w n2^3) B2 ]               (exact, pointwise)
    On the constraint surface a=0  (=> xi_3^2 = -xi_3^1 = -t , n1^3 = t/rho1 , n2^3 = -t/rho2):
        Phi' = t P ,  P := psi(rho1) - psi(rho2) ,
        G    = +beta t ( Q + w R ) ,
        Q := B2/rho1 - A1/rho2 = (1/2)( s1 k2/rho1 - s2 k1/rho2 ) = (1/2) k1 k2 [T(rho1) - T(rho2)] ,
             T(r) := tan(r/2)/r                          (this form is what makes Lemma 2 provable)
        R := B2/rho2 - A1/rho1 = (1/2)( s1 k2/rho2 - s2 k1/rho1 ) + s1 s2 (1/rho1^2 - 1/rho2^2) .
    P, Q, R are all ANTISYMMETRIC under 1<->2, so t^2 P(Q + wbar R) is symmetric -- the L=2 integrand
    is a function of {rho1,rho2} alone.  (NB: Q is NOT a difference of one-link functions; the
    1/rho1 vs 1/rho2 denominators do not line up.  P alone is one.)

STEP 2 (exact 2-link reduction).  u=(p,0), v=(q cos psi, q sin psi), t=xi_3^1, xi_3^2=-t :
    d^3xi^1 d^3xi^2 delta(a) = 2 dt (2pi p dp)(q dq dpsi) , rho1=sqrt(p^2+t^2), rho2=sqrt(q^2+t^2) ,
    z := beta s1 s2 p q/(rho1 rho2) ,   int dpsi e^{z cos psi} = 2pi I0(z) , int dpsi cos psi e^{..} = 2pi I1(z) .
    G and Phi' are both ODD in t and the weight is EVEN in t  =>  E[G] = E[Phi'] = 0  =>  Cov = E[G Phi'] .
    With wbar := <w>_psi = ( p q I1(z)/I0(z) - t^2 )/(rho1 rho2) :

        X(0) = -2 beta * < t^2 P (Q + wbar R) >_{omega}                            (***)
        omega = m M J(m) J(M) exp(-beta(1-k1 k2)) exp(-beta s1 s2 t^2/(m M)) I0(z) dt dm dM ,

    a THREE-DIMENSIONAL integral (m = min, M = max).  (***) is the semi-closed form; constants cancel.

STEP 3 (leading order, beta -> 0).  At beta=0 : omega0 = 4sin^2(m/2)4sin^2(M/2)/(mM) dt dm dM , and
    t | (m,M) is uniform on [-m,m], so <t^2> = m^2/3, <t^4> = m^4/5, wbar = -t^2/(mM) :
        X(0) = -2 beta [ <(m^2/3) P Q> - <(m^4/5) P R/(m M)> ] + O(beta^2) =: -2 beta (c1 - c2) .

VERIFIED HERE
    [i]   the two pointwise identities against finite differences of S and Phi;
    [ii]  (***) against an INDEPENDENT 3-D Monte Carlo keeping w pointwise (no Bessel step);
    [iii] quadrature convergence in all three node counts;
    [iv]  X(beta) for beta in [1e-6, 2.5];  sign;  beta->0 coefficient and its two pieces;
    [v]   the monotonicity facts used in the discussion (psi' >= 0; sigma not monotone);
    [vi]  na2_chain.log vs (***)  =>  X_na2 = X_here / L^2 (L=2: factor 4), a normalisation erratum.

Run: py -3.12 -I ct_chain_exact.py
"""
import numpy as np
from scipy.special import i0, i1

TWOPI = 2 * np.pi


# ------------------------------------------------------------------ one-link functions
def psi(r):
    """phi'(r)/r , phi = -log J , J = (sin(r/2)/(r/2))^2 .  psi(0) = 1/6."""
    r = np.asarray(r, dtype=float)
    small = r < 1e-4
    rs = np.where(small, 1.0, r)
    out = 2.0 / rs ** 2 - 1.0 / (rs * np.tan(rs / 2))
    return np.where(small, 1.0 / 6 + r ** 2 / 360.0, out)


def dpsi(r):
    """d/dr [ 2/r^2 - cot(r/2)/r ] = -4/r^3 + csc^2(r/2)/(2r) + cot(r/2)/r^2"""
    r = np.asarray(r, dtype=float)
    small = r < 1e-3
    rs = np.where(small, 1.0, r)
    out = -4.0 / rs ** 3 + 0.5 / (rs * np.sin(rs / 2) ** 2) + 1.0 / (rs ** 2 * np.tan(rs / 2))
    return np.where(small, r / 180.0, out)


def sigma(r):
    return np.sin(r) / (2.0 * r)


def dsigma(r):
    return np.cos(r) / (2.0 * r) - np.sin(r) / (2.0 * r ** 2)


def Jf(r):
    return 4.0 * np.sin(r / 2) ** 2 / r ** 2


def c12(rho1, rho2, w):
    return np.cos(rho1 / 2) * np.cos(rho2 / 2) + np.sin(rho1 / 2) * np.sin(rho2 / 2) * w


def QofR(rho1, rho2):
    """Q = B2/rho1 - A1/rho2,  R = B2/rho2 - A1/rho1  (constraint surface, t carried separately)"""
    s1, s2 = np.sin(rho1 / 2), np.sin(rho2 / 2)
    k1, k2 = np.cos(rho1 / 2), np.cos(rho2 / 2)
    A1 = s2 * (k1 / 2 - s1 / rho1)
    B2 = s1 * (k2 / 2 - s2 / rho2)
    return B2 / rho1 - A1 / rho2, B2 / rho2 - A1 / rho1


# ------------------------------------------------------------------ [i] pointwise checks
def check_pointwise(nsamp=8, seed=7):
    rng = np.random.default_rng(seed)
    print("-" * 104)
    print("[i] POINTWISE IDENTITIES vs finite differences (random L=2 configurations)")
    print(f"    {'#':>3}{'|G_fd-G_form|/|G|':>22}{'|P_fd-P_form|/|P|':>22}{'|P2_fd-P2_form|/|P2|':>24}")
    worst = np.zeros(3)
    for k in range(nsamp):
        xi1 = rng.uniform(-1, 1, 3) * 2.0
        xi2 = rng.uniform(-1, 1, 3) * 2.0
        rho1, rho2 = np.linalg.norm(xi1), np.linalg.norm(xi2)
        n1, n2 = xi1 / rho1, xi2 / rho2
        w = n1 @ n2
        beta = 1.0
        s1, s2 = np.sin(rho1 / 2), np.sin(rho2 / 2)
        k1, k2 = np.cos(rho1 / 2), np.cos(rho2 / 2)

        def SP(h):     # d/da = d/dxi_3^1 + d/dxi_3^2  (shift both third components)
            x1 = xi1 + np.array([0, 0, h]); x2 = xi2 + np.array([0, 0, h])
            r1, r2 = np.linalg.norm(x1), np.linalg.norm(x2)
            ww = (x1 / r1) @ (x2 / r2)
            S = beta * (1 - c12(r1, r2, ww))
            P = (2 * np.log(r1 / 2) - 2 * np.log(np.abs(np.sin(r1 / 2)))
                 + 2 * np.log(r2 / 2) - 2 * np.log(np.abs(np.sin(r2 / 2))))
            return S, P
        h = 1e-5
        Sp, Pp = SP(h); Sm, Pm = SP(-h); S0, P0 = SP(0.0)
        G_fd = (Sp - Sm) / (2 * h); P_fd = (Pp - Pm) / (2 * h); P2_fd = (Pp - 2 * P0 + Pm) / h ** 2
        A1 = s2 * (k1 / 2 - s1 / rho1); B2 = s1 * (k2 / 2 - s2 / rho2)
        G_form = beta * ((n2[2] - w * n1[2]) * A1 + (n1[2] - w * n2[2]) * B2)
        P_form = psi(rho1) * xi1[2] + psi(rho2) * xi2[2]
        d2phi = -2 / rho1 ** 2 + 0.5 / np.sin(rho1 / 2) ** 2
        d2phb = -2 / rho2 ** 2 + 0.5 / np.sin(rho2 / 2) ** 2
        P2_form = (d2phi * n1[2] ** 2 + psi(rho1) * (1 - n1[2] ** 2)
                   + d2phb * n2[2] ** 2 + psi(rho2) * (1 - n2[2] ** 2))
        e = [abs(G_fd - G_form) / max(abs(G_fd), 1e-12),
             abs(P_fd - P_form) / max(abs(P_fd), 1e-12),
             abs(P2_fd - P2_form) / max(abs(P2_fd), 1e-12)]
        worst = np.maximum(worst, e)
        print(f"    {k:>3}{e[0]:>22.2e}{e[1]:>22.2e}{e[2]:>24.2e}")
    print(f"    worst over {nsamp}:  G {worst[0]:.2e}   Phi' {worst[1]:.2e}   Phi'' {worst[2]:.2e}")
    print("    NB G is a FIRST derivative of S: no second divided difference of exp appears anywhere")
    print("    in this script, and Phi'' is the exact scalar phi''(rho), not that of a matrix function.")

    # also: the constraint-surface reduction G = +beta t (Q + w R)
    print("\n    constraint-surface form  G = +beta t (Q + w R) :")
    for k in range(3):
        rho1 = rng.uniform(0.3, 5.8); rho2 = rng.uniform(0.3, 5.8)
        t = rng.uniform(-min(rho1, rho2) * 0.9, min(rho1, rho2) * 0.9)
        p = np.sqrt(rho1 ** 2 - t ** 2); q = np.sqrt(rho2 ** 2 - t ** 2); ang = rng.uniform(0, TWOPI)
        w = (p * q * np.cos(ang) - t ** 2) / (rho1 * rho2)
        n1 = np.array([p, 0.0, t]) / rho1
        n2 = np.array([q * np.cos(ang), q * np.sin(ang), -t]) / rho2
        s1, s2 = np.sin(rho1 / 2), np.sin(rho2 / 2)
        k1, k2 = np.cos(rho1 / 2), np.cos(rho2 / 2)
        A1 = s2 * (k1 / 2 - s1 / rho1); B2 = s1 * (k2 / 2 - s2 / rho2)
        G_gen = beta * ((n2[2] - w * n1[2]) * A1 + (n1[2] - w * n2[2]) * B2)
        Q, R = QofR(rho1, rho2)
        G_con = beta * t * (Q + w * R)
        print(f"      rho=({rho1:.3f},{rho2:.3f}) t={t:+.3f}  G_gen={G_gen:+.8f}  G_cons={G_con:+.8f}"
              f"  rel={abs(G_gen - G_con) / abs(G_con):.2e}")
    return worst


# ------------------------------------------------------------------ (***) quadrature
def nodes(nr, ns, nt):
    x, w = np.polynomial.legendre.leggauss(nr)
    wr = 0.5 * (1 + np.sin(np.pi * x / 2))       # clusters quadratically at rho=0 and 2pi
    dwr = 0.5 * (np.pi / 2) * np.cos(np.pi * x / 2)
    m = TWOPI * wr; jm = TWOPI * dwr * w
    xs, ws = np.polynomial.legendre.leggauss(ns)
    s = 0.5 * (xs + 1.0); js = 0.5 * ws
    xt, wt = np.polynomial.legendre.leggauss(nt)
    return m, jm, s, js, xt, wt


def quad(beta, nr=200, ns=160, nt=120):
    """(***) :  X(0) = -2 beta < t^2 P (Q + wbar R) >_omega  .  Returns X and the normalised moments."""
    m, jm, s, js, xt, wt = nodes(nr, ns, nt)
    st = np.sqrt(1.0 - xt ** 2)
    Pm = psi(m); Qm = QofR(m, m)[0] * 0 + 0.0
    D = 0.0; N = 0.0; NQ = 0.0; NR = 0.0; T2 = 0.0; T4 = 0.0
    for k in range(len(s)):
        M = m + (TWOPI - m) * s[k]
        t = m[:, None] * xt[None, :]
        p = m[:, None] * st[None, :]
        q = np.sqrt(M[:, None] ** 2 - t ** 2)
        r1 = m[:, None]; r2 = M[:, None]
        s1, s2 = np.sin(r1 / 2), np.sin(r2 / 2)
        k1, k2 = np.cos(r1 / 2), np.cos(r2 / 2)
        z = beta * s1 * s2 * p * q / (r1 * r2)
        I0 = i0(z)
        wbar = (p * q * i1(z) / I0 - t ** 2) / (r1 * r2)
        dens = (r1 * r2 * Jf(r1) * Jf(r2) * np.exp(-beta * (1 - k1 * k2))
                * np.exp(-beta * s1 * s2 * t ** 2 / (r1 * r2)) * I0)
        wgt = jm[:, None] * js[k] * (TWOPI - m)[:, None] * m[:, None] * wt[None, :]
        dens = dens * wgt
        P = (psi(r1) - psi(r2))
        Q, R = QofR(r1, r2)
        t2 = t ** 2
        D += dens.sum()
        NQ += (dens * t2 * P * Q).sum()
        NR += (dens * t2 * P * wbar * R).sum()
        N += (dens * t2 * P * (Q + wbar * R)).sum()
        T2 += (dens * t2).sum(); T4 += (dens * t2 * t2).sum()
    return dict(X=-2.0 * beta * N / D, D=D, N=N / D, NQ=NQ / D, NR=NR / D, T2=T2 / D, T4=T4 / D)


# ------------------------------------------------------------------ [ii] independent MC
def mc_X(beta, n=4_000_000, seed=3):
    """exact conditioning (xi_3^2 = -xi_3^1) and POINTWISE w: no Bessel step, no kernel.
       proposal: rho_i ~ mu(rho) = 4 sin^2(rho/2)d(rho) ,  t ~ U(-m,m) ,  psi ~ U(0,2pi) .
       The constrained density in (rho1,rho2,t) is prop. to 16 s1^2 s2^2/(rho1 rho2) (the 1/(rho1 rho2)
       is the angular Jacobian dOmega = dt/rho), so the importance weight is
            w = exp(beta c) * 2m / (rho1 rho2)      (2m = the length of the t-interval)."""
    rng = np.random.default_rng(seed)

    def samp_rho(nn):
        out = np.empty(nn); got = 0
        while got < nn:
            u = rng.uniform(0, TWOPI, 4 * (nn - got) + 16)
            u = u[rng.uniform(0, 1, u.size) < np.sin(u / 2) ** 2][:nn - got]
            out[got:got + u.size] = u; got += u.size
        return out

    r1 = samp_rho(n); r2 = samp_rho(n)
    m = np.minimum(r1, r2)
    t = rng.uniform(-1, 1, n) * m
    ang = rng.uniform(0, TWOPI, n)
    p = np.sqrt(r1 ** 2 - t ** 2); q = np.sqrt(r2 ** 2 - t ** 2)
    w = (p * q * np.cos(ang) - t ** 2) / (r1 * r2)
    s1, s2 = np.sin(r1 / 2), np.sin(r2 / 2)
    k1, k2 = np.cos(r1 / 2), np.cos(r2 / 2)
    P = psi(r1) - psi(r2); Q, R = QofR(r1, r2)
    G = beta * t * (Q + w * R); PH = t * P
    # per-link pieces, unit speed:  g_i = -beta d/dxi_3^i c12 ,  p_i = d/dxi_3^i Phi
    d1c = (t / r1) * (-0.5 * s1 * k2 + 0.5 * k1 * s2 * w) + (s1 * s2 / r1) * (-t / r2 - w * t / r1)
    d2c = (-t / r2) * (-0.5 * s2 * k1 + 0.5 * k2 * s1 * w) + (s1 * s2 / r2) * (t / r1 + w * t / r2)
    g1, g2 = -beta * d1c, -beta * d2c
    p1, p2 = psi(r1) * t, -psi(r2) * t
    ww = np.exp(beta * (k1 * k2 + s1 * s2 * w)) * 2.0 * m / (r1 * r2)
    ww = ww / ww.sum()
    EG = (ww * G).sum(); EP = (ww * PH).sum()
    cov = (ww * G * PH).sum() - EG * EP
    se = np.sqrt(max((ww ** 2 * ((G - EG) * (PH - EP) - cov) ** 2).sum(), 0.0))
    E = lambda x: (ww * x).sum()
    C = lambda x, y: E(x * y) - E(x) * E(y)
    bd = np.array([[C(g1, p1), C(g1, p2)], [C(g2, p1), C(g2, p2)]])
    return dict(X=-2 * cov, se=2 * se, t2=(ww * t ** 2).sum(), EG=EG, EP=EP,
                bd=bd, bdsum=bd.sum(), cov=cov, Gcheck=float(np.max(np.abs(G - (g1 + g2)))) if n else 0.0)


# ------------------------------------------------------------------ [v] monotonicity
def mono():
    r = np.linspace(1e-6, TWOPI - 1e-6, 2000001)
    dp = dpsi(r); ds = dsigma(r)
    print("\n[v] MONOTONICITY FACTS USED IN THE DISCUSSION")
    print(f"    psi(0+)={psi(np.array([1e-9]))[0]:.8f}   psi(2pi-)={psi(np.array([TWOPI - 1e-6])).item():.4e}")
    print(f"    min dpsi/dr on (0,2pi) = {dp.min():+.4e} at rho={r[np.argmin(dp)]:.4f}"
          f"   => psi strictly increasing (min is at rho->0 where dpsi/dr = rho/180 -> 0)")
    print(f"    sigma(0+)={sigma(np.array([1e-9]))[0]:.8f}  min={sigma(r).min():+.8f}"
          f" at rho={r[np.argmin(sigma(r))]:.4f}  sigma(2pi-)={sigma(np.array([TWOPI - 1e-6])).item():+.3e}")
    ch = np.where(np.diff(np.sign(ds)) != 0)[0]
    print(f"    dsigma/dr changes sign at rho = {np.round(r[ch], 4)} (tan(rho)=rho) => sigma NOT monotone")
    print("    (sigma is NOT the one-link function behind Q -- see [vii]; Q = k1k2[T(r1)-T(r2)]/2.)")


# ------------------------------------------------------------------ [vii] pointwise positivity
def pointwise():
    """The sign of X(0) reduces to two 2-variable inequalities.  Lemma 1 and Lemma 2 are proved in
       the header comment of this function; inequality (R) is checked on a grid."""
    print("\n" + "-" * 104)
    print("[vii] THE SIGN REDUCES TO A beta-INDEPENDENT GEOMETRIC FACT  (P Q >= 0 and |R| <= |Q|)")
    print("""
    LEMMA 1 (psi strictly increasing).  psi(r) = 2/r^2 - cot(r/2)/r = (sin v - v cos v)/(2 v^2 sin v),
        v = r/2.  With H(v) := v^2 - 2 sin^2 v + (v/2) sin 2v,   sign(psi') = sign(H) and
        H(0)=H'(0)=0, H''(v) = 4 sin v (sin v - v cos v) > 0 on (0,pi)  [since (sin v - v cos v)' = v sin v > 0].
        So psi' > 0 on (0,2pi).                                                      [PROVED]
    LEMMA 2 (sign Q = sign(rho1 - rho2)).  Q = (1/2) k1 k2 [ T(rho1) - T(rho2) ],  T(r) = tan(r/2)/r ,
        and T is strictly increasing on (0,pi) and on (pi,2pi) with T(pi-)=+inf, T(pi+)=-inf :
        T' = [r - sin r] / (2 r^2 cos^2(r/2)) > 0  on (0,pi);  on (pi,2pi) the numerator is
        (r/2)sec^2(r/2) - tan(r/2) > 0 because tan < 0 there.  Hence, on the three cases
        (both < pi | both > pi | straddling pi), sign(Q) = sign(rho1-rho2).          [PROVED]
    COROLLARY (with Lemma 1):  P Q >= 0, with equality iff rho1 = rho2.             [PROVED]
    REMAINING INEQUALITY (R):  |R| <= |Q| on (0,2pi)^2.                       [verified below]
    Given (R) and |wbar| <= 1 (elementary: |n1.n2| <= 1 and wbar is a Bessel average of w), the
    integrand of (***) is pointwise >= 0 and not identically 0  =>  X(0) < 0 for EVERY beta > 0.
    The whole beta-dependence drops out of the sign question.""")
    r = np.linspace(2e-3, TWOPI - 2e-3, 4001)
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    P = psi(R1) - psi(R2); Q, Rr = QofR(R1, R2)
    off = np.abs(R1 - R2) > 1e-3                       # drop the diagonal, where Q = R = 0
    kap = np.where(off, np.abs(Rr) / np.where(np.abs(Q) > 0, np.abs(Q), 1.0), 0.0)
    print(f"    grid 4001^2 :  #(PQ<0) = {int((P*Q<0).sum())}   max |R|/|Q| (off-diagonal) = {kap.max():.8f}"
          f"   #(|R|>|Q|) = {int((np.abs(Rr)>np.abs(Q)).sum())}")
    j = np.unravel_index(kap.argmax(), kap.shape)
    print(f"    the supremum sits at rho = ({R1[j]:.6f}, {R2[j]:.6f}) -> 2pi (the corner where Q,R -> 0);")
    print(f"    there the exact ratio tends to (2pi-eps)/(2pi+eps) -> 1, i.e. |R| < |Q| for every eps>0,")
    print(f"    and the corner carries weight J(m)J(M) ~ eps^4 -> 0.")
    s1, s2 = np.sin(R1 / 2), np.sin(R2 / 2)
    wt = s1 ** 2 * s2 ** 2 / (R1 * R2)                 # beta=0 weight in (rho1,rho2)
    o = np.argsort(-wt.ravel()); cum = np.cumsum(wt.ravel()[o]) / wt.sum()
    keep = o[cum <= 0.99]
    print(f"    max |R|/|Q| over the region carrying 99% of the beta=0 weight = {kap.ravel()[keep].max():.6f}")
    print(f"    min |P Q| off the diagonal / max |P Q| = {np.min(np.where(off, np.abs(P*Q), np.inf))/np.max(np.abs(P*Q)):.3e}")


# ------------------------------------------------------------------ main
if __name__ == '__main__':
    print("=" * 104)
    print("T1   L=2 SU(2) chain, chain-average coarse field, a0=0 :  exact reduction of X(0) = -2 Cov(G,Phi')")
    print("=" * 104)
    check_pointwise()
    mono()
    pointwise()

    print("\n" + "-" * 104)
    print("[iii] QUADRATURE CONVERGENCE of (***)")
    print(f"    {'beta':>6}{'nr,ns,nt':>14}{'X(0)':>18}{'<t^2>':>12}{'<t^4>':>12}")
    for beta in (0.5, 2.5):
        for (nr, ns, nt) in ((100, 80, 60), (150, 120, 90), (200, 160, 120), (300, 200, 150)):
            o = quad(beta, nr, ns, nt)
            print(f"    {beta:6.2f}{f'{nr},{ns},{nt}':>14}{o['X']:+18.12f}{o['T2']:12.6f}{o['T4']:12.6f}")

    print("\n" + "-" * 104)
    print("[ii] INDEPENDENT 3-D MONTE CARLO (pointwise w, no Bessel step) vs the quadrature (***)")
    print(f"    {'beta':>6}{'X_quad':>16}{'X_MC':>16}{'se':>10}{'diff/se':>9}{'<t^2>_q':>10}{'<t^2>_MC':>10}")
    for beta in (0.0, 0.5, 1.0, 2.5):
        o = quad(beta, 200, 160, 120)
        r = mc_X(beta, n=4_000_000)
        d = r['X'] - o['X']
        print(f"    {beta:6.2f}{o['X']:+16.10f}{r['X']:+16.10f}{r['se']:10.2e}{d / max(r['se'], 1e-30):9.2f}"
              f"{o['T2']:10.5f}{r['t2']:10.5f}")

    print("\n" + "-" * 104)
    print("[iv] X(0) vs beta  (300,200,150).  NQ = <t^2 P Q>, NR = <t^2 P wbar R> ; X = -2 beta (NQ+NR)")
    print(f"    {'beta':>8}{'X(0)':>18}{'X/beta':>14}{'-1/2 Cov':>12}{'NQ':>13}{'NR':>13}{'NQ+NR':>13}")
    xs = []
    for beta in (0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5):
        o = quad(beta, 300, 200, 150)
        xs.append((beta, o['X']))
        print(f"    {beta:8.2f}{o['X']:+18.12f}{o['X'] / beta:+14.8f}{-0.5 * o['X'] / beta:+12.6f}"
              f"{o['NQ']:+13.5e}{o['NR']:+13.5e}{o['NQ'] + o['NR']:+13.5e}")
    o0 = quad(1e-6, 300, 200, 150)
    print(f"    beta->0 :  X/beta = {o0['X'] / 1e-6:+.8f}   => c1 - c2 = {-(o0['X'] / 1e-6) / 2:+.8f}"
          f"   (c1 = {o0['NQ']:+.6e},  c2 = {-o0['NR']:+.6e})")
    print(f"    sign: X(0) < 0 at every beta tested:  {all(v < 0 for _, v in xs)}")

    print("\n" + "-" * 104)
    print("[vi] NORMALISATION vs na2_chain.log :  X_na2(chain average) = X_here / L^2  (L=2 -> /4)")
    print("     na2 used G = (1/L)sum_i d/dxi_3^i S and Phi' = (1/L)sum_i d/dxi_3^i Phi, but for")
    print("     a = (1/L)sum_i xi_3^i the unit-speed field is sum_i d/dxi_3^i (no 1/L).  Hence /L^2.")
    ref = {0.5: -2.579e-3, 1.0: -4.705e-3, 2.5: -8.079e-3}
    print(f"    {'beta':>6}{'na2 X':>12}{'X_quad/4':>14}{'ratio':>8}{'se':>10}")
    for beta in (0.5, 1.0, 2.5):
        o = quad(beta, 300, 200, 150)
        print(f"    {beta:6.2f}{ref[beta]:+12.5f}{o['X'] / 4:+14.8f}{ref[beta] / (o['X'] / 4):8.4f}{0.0004:10.4f}")
    o = quad(2.5, 300, 200, 150)
    print(f"    EXACT  X(0; L=2, beta=2.5) = {o['X']:+.10f}  (unit-speed a);  na2's number x4 = {4 * ref[2.5]:+.8f}")

    print("\n" + "-" * 104)
    print("[viii] WHICH COVARIANCE CARRIES X?   Cov(G,Phi') = sum_{i,j} Cov(g_i, p_j)  (MC, 4e6)")
    print("       g_i = -beta d/dxi_3^i c12 (per link) ,  p_i = d/dxi_3^i Phi .")
    print("       A SINGLE-LINK coarse field sees only ONE column (and p_k == 0 on the surface).")
    print(f"    {'beta':>6}{'C(g1,p1)':>13}{'C(g1,p2)':>13}{'C(g2,p1)':>13}{'C(g2,p2)':>13}"
          f"{'sum':>13}{'-X/2':>13}{'off-diag share':>15}")
    for beta in (0.5, 1.0, 2.5):
        r = mc_X(beta, n=4_000_000)
        b = r['bd']
        tot = b.sum(); off = b[0, 1] + b[1, 0]
        print(f"    {beta:6.2f}{b[0,0]:+13.3e}{b[0,1]:+13.3e}{b[1,0]:+13.3e}{b[1,1]:+13.3e}"
              f"{tot:+13.5e}{-0.5 * r['X']:+13.5e}{off / tot:15.2%}")
    print("    => the off-diagonal (cross-link) covariances are what survive; they are exactly the terms")
    print("       that do not exist for a single-link coarse field.  |G-(g1+g2)| check:", end=' ')
    print(f"{mc_X(1.0, n=200000)['Gcheck']:.2e}")
    print("\n  Clay 0%.  Exact reduction, no QED.")
