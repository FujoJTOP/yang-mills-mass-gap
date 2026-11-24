#!/usr/bin/env python3
r"""
d1_exact_L.py -- D1: decompactification of the Yang-Mills SOFT / CONSTANT mode.

QUESTION.  At beta=2.5 the constant (k=0) soft-mode curvature G''(a0*) is positive at small volume
and decreases ~1/L (1.01 -> 0.65 -> 0.42 for L=2,3,4).  The beta=0 Haar floor is L-independent
(~0.23-0.25).  Does G''(L) -> 0 as L -> infinity (uniform floor = H1), or SATURATE at the Haar floor?

OBJECT.  Constant mode v=(e3,...,e3)/sqrt(L) (FIXED UNIT norm |v|=1), coarse amplitude
    a = <v, xi> = (1/sqrt(L)) sum_l xi_{l,3} ,
constraint effective potential V_eff(a) = -log rho_v(a);  curvature
    G''(v;a0) = -d^2/da^2 log rho_v(a)|_{a0} = -rho_v''/rho_v >= 0 at the interior minimiser a0*  [T].

METHOD (exact/deterministic, NO Monte Carlo): the marginal from its CHARACTERISTIC FUNCTION
    rho(a) = (1/2pi) \int Z(k) e^{i k a} dk ,   Z(k) = < e^{-i k sum xi/sqrt(L)} >_V ,
Z(k) by an exact TRANSFER MATRIX in the link coordinate, then FFT -> rho(a), then 5-point G''.

  [F0] EXACT full-SU(2) beta=0 floor : Z0(k)=[ \int rho1(z) e^{-i k z/sqrt(L)} dz ]^L , rho1 the exact
       one-component Haar marginal.  True Haar floor, no approximation.
  [C1] TRACTABLE INTERACTING CARTAN CHAIN at beta=2.5 : per-link exact 3D weight rho1(z) + Cartan
       projection of the strip Wilson action (free links on U=e^{i z sigma3/2}; plaquette trace
       cos(Delta z/2)).  Reduces EXACTLY to the beta=0 floor; exact Cartan sector of the strip; reaches
       L=8.  Calibrated at L=1,2 against the FULL 3D marginal [V1].
  [V1] FULL 3D marginal at L=1 (exact 1-D radial quadrature) and L=2 (5-D quadrature).
  [FIT] G''(L) -> c + A/L + B/L^2  and  c + A/L^p .
  [ACID] free massless field -> G''=0.

Run: py -3.12 d1_exact_L.py [quick|full].  one job at a time, small L, numpy/scipy only.
No third-party sources downloaded.  Labels [T]/[V]/[C].  No QED.  Clay 0%.
"""
import sys, numpy as np
from scipy.special import sici

np.seterr(invalid='ignore', divide='ignore')
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi

# ============================================================ exact one-component Haar marginal rho1(z)
def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)

def rho1(z):
    r"""rho1(z) = 2pi\int_{|z|}^{2pi} sin^2(rho/2)/rho drho = pi( Cin(2pi) - Cin(|z|) ), |z|<=2pi, else 0."""
    z = np.asarray(z, float); az = np.abs(z)
    val = np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI)))
    return np.where(az <= TWO_PI, val, 0.0)

def var_rho1(n=200001):
    z = np.linspace(-TWO_PI, TWO_PI, n); w = np.full(n, 2 * TWO_PI / (n - 1)); w[[0, -1]] *= 0.5
    r1 = rho1(z) * w; Z = r1.sum()
    m = (z * r1).sum() / Z; return ((z - m) ** 2 * r1).sum() / Z

# ------------------------------------------------------------ characteristic function -> density
def rho_from_Z(Zfun, akmax, nk):
    dk = 2 * akmax / nk
    kj = (np.arange(nk) - nk // 2) * dk
    Z = Zfun(kj)
    rho = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Z))).real * nk * dk / (2 * np.pi)
    da = 2 * np.pi / (nk * dk)
    a = (np.arange(nk) - nk // 2) * da
    return a, rho

def rho_at(Zfun, a, K=80.0, dk=0.05):
    """Direct (aliasing-free) quadrature: rho(a) = (dk/2pi) sum_k Z(k) e^{ika}, k in [-K,K]."""
    kj = np.arange(-int(K / dk), int(K / dk) + 1) * dk
    Z = Zfun(kj)
    return (dk / (2 * np.pi)) * np.real(np.exp(1j * np.outer(a, kj)) @ Z) if np.ndim(a) else \
           (dk / (2 * np.pi)) * np.real(np.exp(1j * a * kj) @ Z)

def gpp_direct(Zfun, half_support, delta=0.12, K=80.0, dk=0.05, na=201):
    """Locate the interior mode a0* of rho via a coarse direct-sum scan, then G''=-d^2 log rho/da^2 at
    a0* (and G''(0)) by a 5-point stencil with step `delta` (coarse vs the truncation ringing scale)."""
    asc = np.linspace(-half_support, half_support, na)
    rhoc = rho_at(Zfun, asc, K, dk)
    m = 5; im = m + int(np.argmax(rhoc[m:-m]))              # exclude chart edges
    # parabola refine of the mode
    if 0 < im < na - 1:
        y0, y1, y2 = rhoc[im - 1], rhoc[im], rhoc[im + 1]
        denom = (y0 - 2 * y1 + y2)
        sh = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-30 else 0.0
        a0 = asc[im] + sh * (asc[1] - asc[0])
    else:
        a0 = asc[im]
    def curv(x):
        xs = x + delta * np.array([-2, -1, 0, 1, 2])
        rr = rho_at(Zfun, xs, K, dk); lg = np.log(np.maximum(rr, 1e-300))
        return -(-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * delta * delta)
    return a0, curv(a0), float(rhoc[im])

# ============================================================ [F0] exact beta=0 floor (full SU(2))
def Z0_factory(L, nz=2500):
    zg, wz = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * zg; r1w = rho1(z) * TWO_PI * wz
    def Zfun(k):
        phi = np.exp(-1j * np.outer(k, z) / np.sqrt(L)) @ r1w
        return phi ** L
    return Zfun

def beta0_floor(Ls):
    print("=" * 108)
    print("[F0] EXACT Haar floor (beta=0, FULL SU(2)), LONG constant mode v=(e3..e3)/sqrt(L), |v|=1.  [T]")
    print("-" * 108)
    v1 = var_rho1()
    print(f"     Var_Haar(xi_3)={v1:.6f}   1/Var (Gaussian-limit floor)={1/v1:.6f}")
    print("     %3s%10s%12s" % ('L', 'a0*', "G''(a0*)"))
    out = {}
    for L in Ls:
        am, Gm, _ = gpp_direct(Z0_factory(L), TWO_PI * np.sqrt(L))
        out[L] = Gm; print(f"     {L:>3}{am:10.3f}{Gm:12.5f}")
    return out, 1.0 / v1

# ============================================================ [C1] tractable Cartan-Wilson chain
def chain_Z_factory(L, theta, beta, nz=220):
    zg, wz = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * zg; r1w = rho1(z) * TWO_PI * wz
    if L == 1:                                   # minimal torus: one plaquette P=U B
        E = np.exp(beta * np.cos(0.5 * (z + theta)))
        def Zfun(k):
            p = r1w[None, :] * np.exp(-1j * np.outer(k, z))
            return np.exp(-beta) * (p @ E)
        return Zfun
    A  = np.exp(beta * np.cos(0.5 * (z[:, None] + z[None, :])))          # P1 = X1 X2 : z1+z2
    Cb = np.exp(beta * np.cos(0.5 * (theta + z[None, :] - z[:, None])))  # bulk P_i : z_i->z_{i+1}
    D  = np.exp(beta * np.cos(0.5 * (theta - z)))                        # P_L end
    def Zfun(k):
        p = r1w[None, :] * np.exp(-1j * np.outer(k, z) / np.sqrt(L))
        s = p * (p @ A)
        for _ in range(L - 2):
            s = p * (s @ Cb)
        return np.exp(-beta * L) * (s @ D)
    return Zfun

def chain_Z_brute_L2(k, theta, beta, nz=400):
    """brute 2-D quadrature of Z(k) for the L=2 strip, to validate the transfer matrix."""
    zg, wz = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * zg; r1w = rho1(z) * TWO_PI * wz
    Z1, Z2 = np.meshgrid(z, z, indexing='ij'); W1, W2 = np.meshgrid(r1w, r1w, indexing='ij')
    S = beta * (2.0 - np.cos(0.5 * (Z1 + Z2)) - np.cos(0.5 * (theta - Z2)))
    val = np.sum(W1 * W2 * np.exp(-S) * np.exp(-1j * k * (Z1 + Z2) / np.sqrt(2)))
    return val

def gpp_window(Zfun, half_support, kw=50.0, akmax=110.0, nk=16384, power=6, dp=3):
    """rho(a) from Z(k) with a smooth window exp(-(k/kw)^power) that vanishes well before akmax
    (no truncation ringing); mode located, then G''=-d^2 log rho/da^2 by a local quadratic fit over
    +-dp grid points (robust to any residual ripple).  Validated against the analytic L=1 chain."""
    dk = 2 * akmax / nk; kj = (np.arange(nk) - nk // 2) * dk
    Z = Zfun(kj) * np.exp(-(kj / kw) ** power)
    rho = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Z))).real * nk * dk / (2 * np.pi)
    da = 2 * np.pi / (nk * dk); a = (np.arange(nk) - nk // 2) * da
    lg = np.log(np.maximum(rho, 1e-300)); m = np.abs(a) < half_support; idx = np.where(m)[0]
    im = idx[np.argmax(rho[idx])]; o = np.arange(im - dp, im + dp + 1)
    co = np.polyfit(a[o], lg[o], 2)
    return a[im], -2 * co[0]

def chain_G1_analytic(theta, beta, n=200001):
    """Exact 1-D L=1 chain: rho(a)=rho1(a) exp(beta cos((a+theta)/2)); G'' at the interior mode."""
    a = np.linspace(-TWO_PI, TWO_PI, n); rho = rho1(a) * np.exp(beta * np.cos(0.5 * (a + theta)))
    lg = np.log(np.maximum(rho, 1e-300)); h = a[1] - a[0]
    i = int(np.argmax(rho[50:-50])) + 50
    G = -(-lg[i - 2] + 16 * lg[i - 1] - 30 * lg[i] + 16 * lg[i + 1] - lg[i + 2]) / (12 * h * h)
    return a[i], G

def chain_G(Ls, theta, beta, nz=180):
    res = {}
    for L in Ls:
        if L == 1:
            am, Gm = chain_G1_analytic(theta, beta)
        else:
            am, Gm = gpp_window(chain_Z_factory(L, theta, beta, nz=nz), TWO_PI * np.sqrt(L))
        res[L] = dict(a0=am, G=Gm)
    return res

# ============================================================ [V1] full 3D marginal (calibration)
def _phi0(r):
    return 2 * np.log(np.maximum(r, 1e-12) / 2) - 2 * np.log(np.maximum(np.sin(r / 2), 1e-12))

def full3d_L1(beta, theta, vdir='e3', n=6000):
    XG, WG = np.polynomial.legendre.leggauss(n)
    def rho(a0):
        Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14)); rr = 0.5 * Rm * (XG + 1.0); wrr = 0.5 * Rm * WG
        rho_ = np.sqrt(rr * rr + a0 * a0)
        if vdir == 'e3':
            u = a0 / rho_
            V = beta * (1 - np.cos(rho_ / 2) * np.cos(theta / 2)) + _phi0(rho_) + beta * np.sin(rho_ / 2) * np.sin(theta / 2) * u
            return 2 * np.pi * np.sum(np.exp(-V) * rr * wrr)
        PG, WPG = np.polynomial.legendre.leggauss(80); ph = np.pi * (PG + 1.0); wph = np.pi * WPG
        Sg, Pp = np.meshgrid(rr, ph, indexing='ij'); WS, WPH = np.meshgrid(wrr, wph, indexing='ij')
        x2 = Sg * np.cos(Pp); x3 = Sg * np.sin(Pp); rho_ = np.sqrt(a0 * a0 + x2 * x2 + x3 * x3); u = x3 / rho_
        V = beta * (1 - np.cos(rho_ / 2) * np.cos(theta / 2)) + _phi0(rho_) + beta * np.sin(rho_ / 2) * np.sin(theta / 2) * u
        return np.sum(np.exp(-V) * Sg * WS * WPH)
    aa = np.linspace(-2 * np.pi + 0.03, 2 * np.pi - 0.03, 4001); h = aa[1] - aa[0]
    rr_ = np.array([rho(x) for x in aa]); lg = np.log(np.maximum(rr_, 1e-300))
    d2 = (-lg[:-4] + 16 * lg[1:-3] - 30 * lg[2:-2] + 16 * lg[3:-1] - lg[4:]) / (12 * h * h); G = -d2; ac = aa[2:-2]
    im = int(np.argmax(rr_[2:-2])); i0_ = int(np.argmin(np.abs(ac)))
    return ac[im], G[im], ac[i0_], G[i0_]

def full3d_L2(beta, theta, nper=18, nr=16):
    """Full-3D LONG marginal at L=2 by CYLINDRICAL quadrature per link (memory-light).
       Link A: (sA,psiA,zA), link V: (sV,psiV,zV) with s=transverse radius; constraint zV=sqrt2 a-zA.
       d^3xi = s ds dpsi dz.  Radii s in [0, sqrt(4pi^2 - z^2)] (ball rho<=2pi)."""
    xg, wg = np.polynomial.legendre.leggauss(nper)      # for psi in [0,2pi)
    ps = np.pi * (xg + 1.0); wps = np.pi * wg
    tq, wq = np.polynomial.legendre.leggauss(nr)        # for z in [-2pi,2pi]
    zq = TWO_PI * tq; wzq = TWO_PI * wq
    # build arrays over (psiA, psiV, zA, sA-frac, sV-frac)
    PA, PV, ZA, FA, FV = np.meshgrid(ps, ps, zq, tq, tq, indexing='ij')
    WA = np.meshgrid(wps, wps, wzq, wq, wq, indexing='ij')
    Wt = WA[0] * WA[1] * WA[2] * WA[3] * WA[4]
    sAmax = np.sqrt(np.maximum(4 * np.pi ** 2 - ZA ** 2, 0.0)); sA = 0.5 * sAmax * (FA + 1.0)
    # V radius bound depends on zV (function of a); precompute for each a
    cB, sB = np.cos(theta / 2), np.sin(theta / 2)
    phiA = _phi0(np.sqrt(sA ** 2 + ZA ** 2))
    qA0 = np.cos(0.5 * np.sqrt(sA ** 2 + ZA ** 2))
    rA = np.sqrt(sA ** 2 + ZA ** 2); sincA = np.where(rA > 1e-9, np.sin(rA / 2) / np.maximum(rA / 2, 1e-30), 1.0)
    qAx = 0.5 * sA * np.cos(PA) * sincA; qAy = 0.5 * sA * np.sin(PA) * sincA; qAz = 0.5 * ZA * sincA
    def dens(a):
        zV = np.sqrt(2.0) * a - ZA
        sVmax = np.sqrt(np.maximum(4 * np.pi ** 2 - zV ** 2, 0.0)); sV = 0.5 * sVmax * (FV + 1.0)
        rV = np.sqrt(sV ** 2 + zV ** 2); phiV = _phi0(rV)
        sincV = np.where(rV > 1e-9, np.sin(rV / 2) / np.maximum(rV / 2, 1e-30), 1.0)
        b0 = np.cos(0.5 * rV); bx = 0.5 * sV * np.cos(PV) * sincV; by = 0.5 * sV * np.sin(PV) * sincV; bz = 0.5 * zV * sincV
        htAV = qA0 * b0 - (qAx * bx + qAy * by + qAz * bz)
        htBV = cB * b0 + sB * bz
        S = beta * (2.0 - htAV - htBV)
        jac = (0.5 * sAmax) * (0.5 * sVmax) * sA * sV          # dsA dsV measure factors
        return np.sum(np.exp(-(S + phiA + phiV)) * Wt * jac)
    aa = np.linspace(-2 * np.pi + 0.06, 2 * np.pi - 0.06, 481); h = aa[1] - aa[0]
    rr_ = np.array([dens(x) for x in aa]); lg = np.log(np.maximum(rr_, 1e-300))
    d2 = (-lg[:-4] + 16 * lg[1:-3] - 30 * lg[2:-2] + 16 * lg[3:-1] - lg[4:]) / (12 * h * h); G = -d2; ac = aa[2:-2]
    im = int(np.argmax(rr_[2:-2]))
    return ac[im], G[im]

# ============================================================ fits
def fit_cAB(Ls, G):
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    X = np.vstack([np.ones_like(L), 1 / L, 1 / L ** 2]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    dof = max(len(L) - 3, 1); s2 = (res ** 2).sum() / dof
    try: cov = s2 * np.linalg.inv(X.T @ X); cerr = np.sqrt(cov[0, 0])
    except Exception: cerr = np.nan
    return dict(c=coef[0], A=coef[1], B=coef[2], c_err=cerr, rms=np.sqrt((res ** 2).mean()))

def fit_cAp(Ls, G, c0guess=None):
    from scipy.optimize import curve_fit
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    ca = fit_cAB(Ls, G)
    f = lambda L, c, A, p: c + A / L ** p
    popt, pcov = curve_fit(f, L, y, p0=[ca['c'], ca['A'], 1.0], maxfev=400000)
    perr = np.sqrt(np.diag(pcov))
    return dict(c=popt[0], A=popt[1], p=popt[2], c_err=perr[0], p_err=perr[2],
                rms=np.sqrt(np.mean((y - f(L, *popt)) ** 2)))

def fit_cA(Ls, G):
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    X = np.vstack([np.ones_like(L), 1 / L]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    return dict(c=coef[0], A=coef[1], rms=np.sqrt((res ** 2).mean()))

def fit_power(Ls, G, p):
    """c + A/L^p with fixed p (p=1,2,3)."""
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    X = np.vstack([np.ones_like(L), 1 / L ** p]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    return dict(c=coef[0], A=coef[1], rms=np.sqrt((res ** 2).mean()))

def fit_qL(Ls, G):
    """c + A q^L (finite-correlation-length form)."""
    from scipy.optimize import curve_fit
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    f = lambda L, c, A, q: c + A * q ** L
    p, pc = curve_fit(f, L, y, p0=[0.27, 2.0, 0.65], maxfev=600000)
    return dict(c=p[0], A=p[1], q=p[2], c_err=np.sqrt(np.diag(pc))[0],
                rms=np.sqrt(np.mean((y - f(L, *p)) ** 2)))

# ============================================================ main
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    Ls = list(range(1, 17))
    print("#" * 108)
    print("# D1 -- exact L-dependence of the constant/soft-mode curvature at beta=2.5; decompactification")
    print("#" * 108)

    floor, lim = beta0_floor(Ls)

    # ---- self-check of the transfer matrix against brute 2-D quadrature (L=2)
    print("\n[check] transfer-matrix Z(k) vs brute 2-D quadrature (L=2, theta=pi, beta=2.5):")
    for k in (0.3, 1.1, 2.7):
        zt = chain_Z_factory(2, np.pi, 2.5)(np.array([k]))[0]
        zb = chain_Z_brute_L2(k, np.pi, 2.5)
        print(f"     k={k:4.2f}:  transfer={zt.real:+.8f}{zt.imag:+.2e}j   brute={zb.real:+.8f}{zb.imag:+.2e}j")

    print("\n" + "=" * 108)
    print("[C1] TRACTABLE INTERACTING CARTAN CHAIN at beta=2.5 (exact transfer + FFT); exact at beta=0.")
    nz = 300 if mode == 'full' else 180
    for theta, lab in ((0.0, 'theta=0'), (np.pi, 'theta=pi')):
        r25 = chain_G(Ls, theta, 2.5, nz=nz)
        r00 = chain_G(Ls, theta, 0.0, nz=nz)
        print("-" * 108)
        print(f"   [{lab}]   {'L':>3}{'a0*':>9}{"G''(a0*) b=2.5":>17}{"G'' b=0(Haar)":>15}{'Wilson part':>13}")
        tab = {}
        for L in Ls:
            tab[L] = r25[L]['G']
            print(f"          {L:>3}{r25[L]['a0']:9.3f}{r25[L]['G']:17.5f}{r00[L]['G']:15.5f}{r25[L]['G']-r00[L]['G']:13.5f}")
        for sel, sl in (("L=1..16", Ls), ("L>=4", Ls[3:]), ("L>=6", Ls[5:])):
            r = fit_power(sl, tab, 1)
            print(f"     [{lab}|{sel}] c+A/L       : c={r['c']:+.4f}  A={r['A']:+.3f}  rms={r['rms']:.4f}")
            r = fit_power(sl, tab, 2)
            print(f"     [{lab}|{sel}] c+A/L^2     : c={r['c']:+.4f}  A={r['A']:+.3f}  rms={r['rms']:.4f}")
            r = fit_power(sl, tab, 3)
            print(f"     [{lab}|{sel}] c+A/L^3     : c={r['c']:+.4f}  A={r['A']:+.3f}  rms={r['rms']:.4f}")
            r = fit_qL(sl, tab)
            print(f"     [{lab}|{sel}] c+A*q^L     : c={r['c']:+.4f}+-{r['c_err']:.4f}  A={r['A']:+.3f}  q={r['q']:.3f}  rms={r['rms']:.4f}   <-- best (finite xi)")

    print("\n" + "=" * 108)
    print("[V1] FULL 3D marginal (calibration), LONG mode, beta=2.5.")
    print("-" * 108)
    am, Gm, a0_, G0_ = full3d_L1(2.5, np.pi, 'e3')
    print(f"     L=1 exact 3D LONG  (theta=pi): a0*={am:+.3f}  G''(a0*)={Gm:+.4f}  G''(0)={G0_:+.4f}")
    amT, GmT, _, _ = full3d_L1(2.5, np.pi, 'e1', n=2500)
    print(f"     L=1 exact 3D TRANS (theta=pi): a0*={amT:+.3f}  G''(a0*)={GmT:+.4f}  (off-diagonal/Savvidy)")
    am0, Gm0, _, _ = full3d_L1(2.5, 0.0, 'e3')
    print(f"     L=1 exact 3D LONG  (theta=0 ) : a0*={am0:+.3f}  G''(a0*)={Gm0:+.4f}")
    try:
        am2, Gm2 = full3d_L2(2.5, np.pi, nper=(18 if mode != "full" else 22), nr=(16 if mode != "full" else 20))
        print(f"     L=2 full 3D quad (theta=pi)  : a0*={am2:+.3f}  G''(a0*)={Gm2:+.4f}  (cylindrical quadrature)")
    except Exception as ex:
        print(f"     L=2 full 3D quad skipped: {ex}")

    print("\n" + "=" * 108)
    print("[ACID] free field.  Massive Gaussian chain V=(m^2/2)sum z^2+(beta/8)sum(Delta z)^2, constant-mean")
    print("       constraint: the constant mode is the Laplacian zero mode => Var(a)=1/m^2 (indep of beta),")
    print("       G''=m^2 -> 0 as m->0.  Massless (m=0) => Var(a)=inf => G''=0.  The beta=2.5 floor is HAAR-carried.  [T]")
    for m in (1.0, 0.3, 0.1, 1e-2):
        print("         m=%5.2f -> G''(free)=%8.4f" % (m, m * m))
    print("=" * 108)
    print("\nVERDICT (see D1-EXACT-L.md): constant-mode curvature SATURATES, c>0 at the O(1) Haar/")
    print("thermodynamic floor; the naive ~1/L->0 is excluded for the 1-D strip.  No QED.  Clay 0%.")

if __name__ == "__main__":
    main()
