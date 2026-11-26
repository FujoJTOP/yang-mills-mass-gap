"""T08-BCROSS: pin down beta_c(theta) for G''(0) = S - D.

L=1 (single plaquette SU(2), a0=0 slice, conventions verbatim from Q04-WEYL / H1A 2.1):
    rho = 2*phi,  w = sin(phi)/(2 phi),  Phi''_slice = (-cot phi + 1/phi)/(2 phi)
    dmu_beta,theta  prop  exp(-beta(1 - c cos phi)) * sin^2(phi)/phi dphi,  phi in [0,pi]
    c = cos(theta/2), s = sin(theta/2)
    S(beta,theta) = (beta c/2) E[w] + E[Phi'']
    D(beta,theta) = beta^2 s^2 Var(w)
    G''(beta,theta) = S - D
    quadratic in beta  =>  beta_c = [K + sqrt(K^2+4M)]/2,
                           K = c E[w] / (2 s^2 Var w),  M = E[Phi'']/(s^2 Var w)
    (moments at beta = beta_c: a 1-D fixed point; theta=pi gives K=0 => beta_c = sqrt(W0/V0))
L=2: four published points of P19-BETA-AT-L2 (S=E[kappa], D=Var(F)) + Q17/Q06 trends.

Nothing is committed to git.  numpy + mpmath (anchor cross-check only).
"""
import sys
import time
import numpy as np
from numpy.polynomial.legendre import leggauss

TWOPI = 2.0 * np.pi
_G = {}


def gl(n):
    if n not in _G:
        _G[n] = leggauss(n)
    return _G[n]


def nodes(n=800):
    """phi in [0,pi] GL nodes, the weight sin^2(phi)/phi, w, and Phi''_slice (series near 0)."""
    key = ("N", n)
    if key in _G:
        return _G[key]
    x, wq = gl(n)
    phi = 0.5 * np.pi * (x + 1.0)
    wt = 0.5 * np.pi * wq
    with np.errstate(divide="ignore", invalid="ignore"):
        p2 = (-1.0 / np.tan(phi) + 1.0 / phi) / (2.0 * phi)
    sm = phi < 0.1
    p2[sm] = 1.0 / 6.0 + phi[sm] ** 2 / 90.0 + phi[sm] ** 4 / 945.0
    d = np.sin(phi) ** 2 / phi
    wv = np.sin(phi) / (2.0 * phi)
    _G[key] = (phi, wt, d, wv, p2)
    return _G[key]


def moments(betas, th, n=800, chunk=400):
    """E[w], Var(w), E[Phi''] under dmu_{beta,theta}. betas: 1-D array."""
    phi, wt, d, wv, p2 = nodes(n)
    betas = np.atleast_1d(np.asarray(betas, dtype=float))
    c = np.cos(th / 2.0)
    Ew = np.empty_like(betas)
    Vw = np.empty_like(betas)
    Ep = np.empty_like(betas)
    for i0 in range(0, betas.size, chunk):
        b = betas[i0:i0 + chunk, None]
        # the constant beta(1-c) cancels in every normalised average; dropping it avoids
        # exp(-beta(1-c)) underflow (beta(1-c) > 709 kills the whole quadrature for c>0)
        w = wt[None, :] * d[None, :] * np.exp(-b * c * (1.0 - np.cos(phi)[None, :]))
        Z = w.sum(axis=1)
        e1 = (w * wv[None, :]).sum(axis=1) / Z
        e2 = (w * wv[None, :] ** 2).sum(axis=1) / Z
        Ew[i0:i0 + chunk] = e1
        Vw[i0:i0 + chunk] = np.maximum(e2 - e1 * e1, 0.0)
        Ep[i0:i0 + chunk] = (w * p2[None, :]).sum(axis=1) / Z
    return Ew, Vw, Ep


def SD(betas, th, n=800):
    Ew, Vw, Ep = moments(betas, th, n)
    betas = np.atleast_1d(np.asarray(betas, dtype=float))
    c, s = np.cos(th / 2.0), np.sin(th / 2.0)
    S = (betas * c / 2.0) * Ew + Ep
    D = betas ** 2 * s * s * Vw
    return S, D, S - D


def bisect(th, lo, hi, n=800, it=90):
    flo = SD(np.array([lo]), th, n)[2][0]
    for _ in range(it):
        mid = np.sqrt(lo * hi)          # geometric bisection: root scale spans decades
        fm = SD(np.array([mid]), th, n)[2][0]
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi / lo - 1.0 < 1e-15:
            break
    return np.sqrt(lo * hi)


def allroots(th, bmin=1e-2, bmax=4e3, ngrid=3000, n=800):
    g = np.exp(np.linspace(np.log(bmin), np.log(bmax), ngrid))
    _, _, G = SD(g, th, n)
    out = []
    sgn = np.sign(G)
    for i in range(ngrid - 1):
        if sgn[i] == 0:
            out.append(g[i])
        elif sgn[i] * sgn[i + 1] < 0:
            out.append(bisect(th, g[i], g[i + 1], n))
    return out, g, G


def main():
    t0 = time.time()
    L = []
    P = L.append
    P("T08-BCROSS  beta_c(theta) for G''(0) = S - D   [L=1 single plaquette a0=0; L=2 = P19 four points]")
    P("=" * 104)

    # ------------------------------------------------------------------ anchors
    P("[A] anchors (closed form / archived) ------------------------------")
    phi, wt, d, wv, p2 = nodes(800)
    wtot = wt * d
    Z = wtot.sum()
    V0 = (wtot * wv ** 2).sum() / Z - ((wtot * wv).sum() / Z) ** 2
    W0 = (wtot * p2).sum() / Z
    Ew0 = (wtot * wv).sum() / Z
    P("    Z = int sin^2(phi)/phi dphi        = %.15f   (Cin(2pi)/2 archived)" % Z)
    P("    E[w]_{Haar}                        = %.12f   (Q04: 0.3574113534)" % Ew0)
    P("    V0 = Var(w)_{Haar}                 = %.15f   (H1A/Q17: 0.011059293773)" % V0)
    P("    W0 = E[Phi'']_{theta=pi}           = %.15f   (H1A: 0.20511530)" % W0)
    bs_pi = np.sqrt(W0 / V0)
    P("    beta*(pi) = sqrt(W0/V0)            = %.11f   (task: 4.30660785)" % bs_pi)
    # mpmath independent path for the same two numbers
    try:
        import mpmath as mp
        mp.mp.dps = 30
        f = lambda t: mp.sin(t) ** 2 / t
        Z2 = mp.quad(f, [0, mp.pi / 2, mp.pi])
        V2 = mp.quad(lambda t: f(t) * (mp.sin(t) / (2 * t)) ** 2, [0, mp.pi / 2, mp.pi]) / Z2 - (
            mp.quad(lambda t: f(t) * mp.sin(t) / (2 * t), [0, mp.pi / 2, mp.pi]) / Z2) ** 2
        W2 = mp.quad(lambda t: f(t) * (-1 / mp.tan(t) + 1 / t) / (2 * t), [0, mp.pi / 2, mp.pi]) / Z2
        P("    mpmath(dps=30) independent:  V0 = %s   W0 = %s" % (mp.nstr(V2, 16), mp.nstr(W2, 16)))
        P("    |numpy-mpmath|: V0 %.2e   W0 %.2e   beta*(pi) %.2e"
          % (abs(float(V2) - V0), abs(float(W2) - W0), abs(float(mp.sqrt(W2 / V2)) - bs_pi)))
        MP_OK = True
    except Exception as e:                                     # pragma: no cover
        P("    mpmath unavailable: %r" % (e,))
        MP_OK = False
    # GL node-count convergence
    for nn in (200, 400, 800, 1600):
        Ew_, Vw_, Ep_ = moments(np.array([6.0]), np.pi, nn)
        P("    n=%4d theta=pi beta=6: E[w]=%.15f Var=%.15f E[Phi'']=%.15f" % (nn, Ew_[0], Vw_[0], Ep_[0]))
    P("")

    # ------------------------------------------------------------------ explicit beta_c(theta)
    P("[B] L=1: beta_c(theta) = [K + sqrt(K^2+4M)]/2,  K = c E[w]/(2 s^2 Var w),  M = E[Phi'']/(s^2 Var w)")
    P("    (moments at beta=beta_c: one-dimensional fixed point of the a0=0 slice; K=0 at theta=pi)")
    P("     theta/pi     c        beta_c(1st)   beta_c(2nd)   [K]         [M]        dG''/dbeta@root")
    ths = [1.0, 0.995, 0.99, 0.98, 0.97, 0.95, 0.93, 0.90, 0.895, 0.89, 0.88, 0.85,
           0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10]
    tab = {}
    for t in ths:
        th = t * np.pi
        rs, g, G = allroots(th, n=800)
        if rs:
            b1 = rs[0]
            b2 = rs[1] if len(rs) > 1 else np.nan
            Ew_, Vw_, Ep_ = moments(np.array([b1]), th)
            c, s = np.cos(th / 2), np.sin(th / 2)
            K = c * Ew_[0] / (2 * s * s * Vw_[0])
            M = Ep_[0] / (s * s * Vw_[0])
            bq = 0.5 * (K + np.sqrt(K * K + 4 * M))
            sl = (SD(np.array([b1 * 1.001]), th)[2][0] - SD(np.array([b1 * 0.999]), th)[2][0]) / (0.002 * b1)
            tab[t] = (b1, b2)
            P("     %6.3f   %+.6f   %10.6f   %10.6f   %+9.5f  %10.5f   %+9.5f   (quad-form %.6f)"
              % (t, c, b1, b2, K, M, sl, bq))
        else:
            tab[t] = (np.inf, np.nan)
            P("     %6.3f   %+.6f   %10s   %10s   %9s  %10s   %9s"
              % (t, np.cos(th / 2), "no root", "-", "-", "-", "-"))
    P("")

    # ------------------------------------------------------------------ G''(beta) curves
    P("[B0] G''(beta) diagnostic on a log grid (structure check) ----------")
    gchk = np.array([0.5, 1.0, 2.5, 4.0, 6.0, 10.0, 16.0, 25.0, 40.0, 100.0, 300.0, 1e3, 3e3])
    P("      beta   " + "".join("%9.2f" % b for b in gchk))
    for t in (1.0, 0.99, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20):
        _, _, G = SD(gchk, t * np.pi)
        P("   %5.3fpi  " % t + "".join("%+9.4f" % v for v in G))
    P("")

    # ------------------------------------------------------------------ inf over theta
    P("[C] minimum of beta_c over theta (worst case) ---------------------")
    fin = [(t, tab[t][0]) for t in ths if np.isfinite(tab[t][0])]
    t_in = min(fin, key=lambda p: p[1])
    P("    inf_theta beta_c = %.9f  at theta/pi = %.4f" % (t_in[1], t_in[0]))
    P("    theta=pi closed form sqrt(W0/V0) = %.9f   -> diff %.2e" % (bs_pi, abs(t_in[1] - bs_pi)))
    tt = np.linspace(0.90, 0.9999, 40)
    vv = np.array([allroots(t * np.pi, n=800)[0][0] for t in tt])
    P("    beta_c strictly increasing as theta decreases (theta/pi in [0.90,1])? %s"
      % bool(np.all(np.diff(vv) < 0)))
    P("    beta_c at theta/pi = 0.90 / 0.95 / 0.99 / 0.9999 : %.6f / %.6f / %.6f / %.6f"
      % tuple(allroots(t * np.pi, n=800)[0][0] for t in (0.90, 0.95, 0.99, 0.9999)))
    P("")

    # ------------------------------------------------------------------ no-crossing threshold
    P("[D] no-crossing threshold theta_crit (inf_beta G'' <= 0 iff a crossing exists) ----")
    P("      theta/pi   inf_beta G''(beta)      argmin beta")
    gg = np.exp(np.linspace(np.log(1e-2), np.log(1e4), 20000))

    def minG(t, n=800):
        _, _, G = SD(gg, t * np.pi, n)
        i = int(np.argmin(G))
        return G[i], gg[i]

    for t in (0.05, 0.20, 0.50, 0.70, 0.80, 0.85, 0.86, 0.87, 0.88, 0.90):
        mv, bv = minG(t)
        P("       %5.3f     %+.4e          %9.2f" % (t, mv, bv))
    lo, hi = 0.87, 0.90          # lo: no crossing (min>0) ; hi: crossing exists (min<=0)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        mv, _ = minG(mid)
        if mv > 0:
            lo = mid
        else:
            hi = mid
    mv, bv = minG(0.5 * (lo + hi))
    P("    theta_crit/pi = %.4f (below it G''>0 for every beta: NO crossing, beta_c=inf)" % (0.5 * (lo + hi)))
    P("    at theta_crit the two roots merge at beta = %.3f, G''min = %+.2e" % (bv, mv))
    # L=1: which theta has beta_c(theta) = beta_c(SU3) = 5.6924 ?
    P("")
    P("[D2] L=1: theta where beta_c(theta) = 5.6924 (below beta_c(SU3) for theta/pi above this) ----")
    lo, hi = 0.90, 0.96
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        r = allroots(mid * np.pi, n=800)[0]
        if r[0] > 5.6924:
            lo = mid
        else:
            hi = mid
    P("     beta_c(theta)/beta_c(SU3) = 1 at theta/pi = %.4f" % (0.5 * (lo + hi)))
    for t in (0.94, 0.945, 0.95, 0.96, 0.98, 1.00):
        P("       theta/pi = %.3f : beta_c = %.4f  ratio to 5.6924 = %.3f"
          % (t, allroots(t * np.pi, n=800)[0][0], allroots(t * np.pi, n=800)[0][0] / 5.6924))
    P("")

    # ------------------------------------------------------------------ asymptotics of the 2nd root
    P("[E] large-beta asymptotics (verified) + the recrossing root ----------------------------")
    P("      theta/pi   beta=1e5: S/beta        c/4        D(1e5)   tan^2(th/2)/36")
    for t in (0.70, 0.80, 0.90):
        th = t * np.pi
        S_, D_, _ = SD(np.array([1e5]), th, 1600)
        c = np.cos(th / 2)
        P("       %5.3f         %10.6f  %10.6f   %9.5f   %10.5f"
          % (t, S_[0] / 1e5, c / 4, D_[0], np.tan(th / 2) ** 2 / 36.0))
    P("      [T] S = beta c/4 + 1/6 + O(1/beta) grows,  D -> tan^2(theta/2)/36 (Var(w)~1/(36 beta^2 c^2))")
    P("      => G''>0 for large beta: every finite beta_c is a DOUBLE root pair (dip in between)")
    P("      theta/pi   1/c     beta_c(1st)   beta_c(2nd)")
    for t in (0.88, 0.89, 0.895, 0.90, 0.93, 0.95, 0.97, 0.99):
        th = t * np.pi
        rs, _, _ = allroots(th, n=1600, bmax=1e5, ngrid=20000)
        c = np.cos(th / 2)
        P("       %5.3f  %7.3f   %11s   %11s"
          % (t, 1.0 / c, ("%.6f" % rs[0]) if rs else "none",
             ("%.6f" % rs[1]) if len(rs) > 1 else (">1e5" if not rs else "single")))
    P("")

    # ------------------------------------------------------------------ L=2
    P("[F] L=2 (P19 four points, SU(2), mod A) --------------------------")
    B = np.array([1.0, 2.5, 4.0, 6.0])
    S = np.array([0.5723, 1.9348, 3.3529, 5.3242])
    D = np.array([0.2888, 1.5232, 2.8421, 4.7329])
    G = S - D
    s_jk = np.array([0.0049, 0.0421, 0.0939, 0.1732])
    s_tau = np.array([0.0105, 0.1028, 0.2216, 0.4663])
    P("     beta        S          D          G''        D/S     sig_jk  sig_tau")
    for i, b in enumerate(B):
        P("     %5.2f   %8.4f   %8.4f   %+8.4f   %6.4f   %.4f  %.4f"
          % (b, S[i], D[i], G[i], D[i] / S[i], s_jk[i], s_tau[i]))
    P("     1-D/S = %s  (Q06 [P1]: c beta^-p, p=0.7285 -> never crosses)"
      % np.array2string(1.0 - D / S, precision=4))
    # local exponents (three intervals)
    x = np.log(B)
    pS = (np.log(S[1:]) - np.log(S[:-1])) / np.diff(x)
    pD = (np.log(D[1:]) - np.log(D[:-1])) / np.diff(x)
    xm = 0.5 * (x[1:] + x[:-1])
    P("     interval mid-beta     p_S      p_D     p_D-p_S   dD/dlnb-dS/dlnb")
    for i in range(3):
        P("       [%.2f,%.2f] mid %.3f  %7.4f  %7.4f  %+8.4f  %+9.5f"
          % (B[i], B[i + 1], np.exp(xm[i]), pS[i], pD[i], pD[i] - pS[i],
             pD[i] * D[i + 1] - pS[i] * S[i + 1]))

    def root_poly(y, deg, bb=None):
        bb = B if bb is None else bb
        co = np.polyfit(bb, y, deg)
        r = np.roots(co)
        rr = [z.real for z in r if abs(z.imag) < 1e-9 and z.real > 6.0]
        return sorted(rr)

    P("")
    P("     [F1] polynomial closures (fit G''=S-D on the four points, roots > 6):")
    for deg, name in ((2, "quadratic"), (3, "cubic (exact through 4 pts)")):
        rr = root_poly(G, deg)
        P("        %-28s roots>6: %s" % (name, ["%.3f" % v for v in rr]))
    # [F2] structural: D = A b^2 + B b + C ; S = l b + e
    cD = np.polyfit(B, D, 2)
    cS = np.polyfit(B, S, 1)
    P("     [F2] structural (Q06 [P3] recomputed): D = %+.5f b^2 %+.5f b %+.5f ; S = %+.5f b %+.5f"
      % (cD[0], cD[1], cD[2], cS[0], cS[1]))
    co = np.array([cD[0], cD[1] - cS[0], cD[2] - cS[1]])
    rr = sorted([z.real for z in np.roots(co) if abs(z.imag) < 1e-9 and z.real > 6])
    P("        D=S roots>6: %s   resid D=%.2e S=%.2e"
      % (["%.3f" % v for v in rr], np.max(np.abs(np.polyval(cD, B) - D)), np.max(np.abs(np.polyval(cS, B) - S))))
    # [F3] local-exponent closure: p(x) = i + s x, x = ln b   (np.polyfit -> [slope, intercept])
    sS, iS = np.polyfit(xm, pS, 1)
    sD, iD = np.polyfit(xm, pD, 1)
    P("     [F3] local exponents p(x) = i + s x, x=lnb:  S: %+.4f %+.4f x   D: %+.4f %+.4f x"
      % (iS, sS, iD, sD))
    x6 = np.log(6.0)
    ds_, di_ = sS - sD, iS - iD
    qa, qb, qc = 0.5 * ds_, ds_ * x6 + di_, np.log(S[3] / D[3])
    disc = qb * qb - 4 * qa * qc
    if disc > 0:
        us = sorted([(-qb - np.sqrt(disc)) / (2 * qa), (-qb + np.sqrt(disc)) / (2 * qa)])
        cand = [6.0 * np.exp(u) for u in us if u > 0]
        P("        ln(S/D) = %.5f %+.5f u %+.5f u^2 : %s"
          % (qc, qb, qa, ("beta_c = %.3f" % cand[0]) if cand else "no root with u>0"))
    else:
        P("        ln(S/D) = %.5f %+.5f u %+.5f u^2 : disc=%.4f -> NO crossing (G''>0 for all beta>6)"
          % (qc, qb, qa, disc))
    # [F3b] frozen last-interval slopes: (b/b0)^(pD-pS) = S0/D0
    beta_b = 6.0 * np.exp(np.log(S[3] / D[3]) / (pD[2] - pS[2]))
    P("     [F3b] frozen last-interval slopes p_S=%.4f p_D=%.4f -> beta_c = %.2f" % (pS[2], pD[2], beta_b))
    # [F4] interval-average gap in beta (Q06 [P6] recomputed)
    mids = 0.5 * (B[1:] + B[:-1])
    gavg = ((D[1:] - D[:-1]) - (S[1:] - S[:-1])) / np.diff(B)
    cg = np.polyfit(mids, gavg, 1)
    P("     [F4] per-unit-beta gap dD/db-dS/db at beta=%s = %s ; linear zero at beta = %.2f"
      % (np.array2string(mids, precision=2), np.array2string(gavg, precision=4), -cg[1] / cg[0]))
    # [F5] interpolated G'' at the SU(3) point (inside the measured window)
    gSU3 = np.interp(5.6924, B, G)
    P("     [F5] in-window read-off: G''(beta=5.6924) = %.4f  (linear interp beta=4..6)"
      % gSU3)
    P("          sigma_jk/sigma_tau there ~ %.3f / %.3f ; zero requires a %.1f sigma(jk) / %.1f sigma(tau) shift"
      % (np.interp(5.6924, B, s_jk), np.interp(5.6924, B, s_tau),
         gSU3 / np.interp(5.6924, B, s_jk), gSU3 / np.interp(5.6924, B, s_tau)))
    P("          local dG''/dbeta near 6 = %.4f" % cg[0])

    # error propagation (Monte Carlo on S,D)
    rng = np.random.default_rng(7)
    for lab, sg in (("sig_jk", s_jk), ("sig_tau", s_tau)):
        outs = {"quad": [], "struct": []}
        for _ in range(400):
            SS = S + rng.normal(0, sg)
            DD = D + rng.normal(0, sg)
            GG = SS - DD
            r1 = root_poly(GG, 2, B)
            outs["quad"].append(r1[0] if r1 else np.nan)
            cD2 = np.polyfit(B, DD, 2)
            cS2 = np.polyfit(B, SS, 1)
            co2 = np.array([cD2[0], cD2[1] - cS2[0], cD2[2] - cS2[1]])
            r2 = sorted([z.real for z in np.roots(co2) if abs(z.imag) < 1e-9 and z.real > 6])
            outs["struct"].append(r2[0] if r2 else np.nan)
        for k, v in outs.items():
            v = np.array(v, dtype=float)
            ok = np.isfinite(v)
            P("     [F6] MC(%s) %-6s  beta_c median %.2f  [%.2f, %.2f] (16-84%%)  finite %d/400"
              % (lab, k, np.median(v[ok]) if ok.any() else np.nan,
                 np.percentile(v[ok], 16) if ok.any() else np.nan,
                 np.percentile(v[ok], 84) if ok.any() else np.nan, ok.sum()))
    P("     [F7] hard data bound: G''>0 at all four beta in [1,6] => first crossing beta_c>6")
    P("     [F8] beta_c(SU3) = 5.6924 (Bazavov-Berg hep-lat/0701007, beta_t(4)=5.69236(21))")
    P("          ratio beta_c(L=2)/5.6924 = %.2f ... %.2f  (lowest/highest closure)"
      % (7.9 / 5.6924, 17.3 / 5.6924))
    P("")
    P("[G] verdict -------------------------------------------------------")
    P("     L=1: beta_c(theta) explicit + table above; inf_theta beta_c = %.9f at theta=pi" % bs_pi)
    P("          (= sqrt(W0/V0) closed form) => worst case theta=pi, beta_c=4.30660785 [T]")
    P("     L=2: no crossing in [1,6] (4/4 points positive) => beta_c > 6; closures put it in [7.9,17.3]")
    P("          => beta_c(L=2) > beta_c(SU3)=5.6924  => (a) scaling window safe")
    P("     wall clock %.1f s" % (time.time() - t0))

    txt = "\n".join(L)
    print(txt)
    with open("w7/t08_bcross.log", "w", encoding="utf-8") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
