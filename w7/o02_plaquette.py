#!/usr/bin/env python3
r"""
o02_plaquette.py -- is  G(a) = -log rho_v(a)  convex on the WHOLE a-axis?  Single plaquette, N=2 and N=3.

OBJECT (GPP-STRUCT.md eq. 1.2, H1A-EXACT.md eq. 1.3):
    G''(a) = -d_a^2 log rho_v(a) = E[H_vv + Phi''|a] - Var(G_v + Phi'|a)
    G''(a) >= 0  <=>  the 1D conditioned marginal rho_v is log-concave at a.
The known fact (H1A) is ONE point: a=0, N=2, theta=pi, beta*>4.30661 gives G''(0)<0.
This script asks about the FUNCTION a -> G''(a), not one point.

MODEL (N=2, exact; H1A-EXACT.md 2.1 / na_plaquette.py conventions):
    U=exp(xi), xi in R^3, |xi|<2pi (injective chart minus the antipode)
    B=exp(i theta sig3/2),  S=beta(1-cos(r/2)cos(th/2)+sin(r/2)sin(th/2)u),  u=xi3/r
    dU=J d^3xi with J=(sin(r/2)/(r/2))^2  (de Wit-'t Hooft: J=prod_{i<j}[sin(mu_ij/2)/(mu_ij/2)]^2)
    v=e3, condition <v,xi>=a
    => the slice integral is 1-dimensional up to the trivial 2*pi (phi drops out: everything is radial
       in r and depends on the transverse direction only through u=a/r):
       rho(a) = 2 pi int_0^{sqrt(4pi^2-a^2)} s ds e^{-S(r,u)} J(r)
    => exact, spectral (Gauss-Legendre) accuracy.  Beta=0 mass check: int rho da = 16 pi^2.

MODEL (N=3, exact by symmetry reduction):
    U=exp(iH), H traceless Hermitian,  S=beta(1-ReTr(U B)/3),  J=prod_{i<j}[sin(dmu_ij/2)/(dmu_ij/2)]^2
    v = T8 = diag(1,1,-2)/sqrt(12)   (U(2)-invariant Cartan direction; Tr v^2=1/2)
    B = exp(i pi T8)  -- the SU(3) "theta=pi" lift (vacuum at a=-pi, exactly H1A's N=2 structure)
    The stabiliser of v inside SU(3) is U(2), so the 7-dim slice decomposes as
        H = [[A11+t, A12, w1],[A21, A22+t, w2],[conj w1, conj w2, -2t]],  A = a.sigma (3 dim), w in C^2 (4 dim)
        t = a/sqrt3,   tr(vH) = sqrt3 t
    and every U(2)-invariant function averages through the 3 invariants (r=|a|, q=|w|^2, mu=cos angle):
        int_{R^3 x R^4} f dA d^4w = 2 pi^3 int r^2 dr int q dq int_{-1}^{1} dmu f
    (verified analytically and numerically in [S1] below, also against the closed form for exp(-k tr H^2)).
    Chart domain: |lam_i-lam_j| <= 2 pi.
    => rho(a) = 2 pi^3 int r^2 dr q dq dmu e^{-S} J,  a 3-dimensional exact quadrature.

SECTIONS
  [S1] SU(3) reduction validation: quadrature vs closed form for exp(-k tr H^2); measure constant 2pi^3.
  [V]  N=2 validation against the archived constants (0.20511530 / +0.135995 / -0.193019 / beta*=4.30661),
       plus the beta=0 mass sum rule int rho da = 16 pi^2.
  [A]  N=2: the WHOLE curve G''(a), beta=2.5 and beta=6, theta=pi (and theta=0 as control).
  [B]  N=2: beta-dependence -- does the whole curve flip, or only a=0?  sign-change points a_k(beta).
  [C]  N=3: the WHOLE curve G''(a), beta=2.5 / beta=6, with node-refinement error bars.
Run: py -3.12 o02_plaquette.py
Labels [T]/[C]/[E].  No QED.  Clay 0%.
"""
import sys
import time
import numpy as np

np.set_printoptions(precision=10, linewidth=170)
TWO_PI = 2.0 * np.pi

# ------------------------------------------------------------------ quadrature nodes (fixed, cheap)
_GLX, _GLW = np.polynomial.legendre.leggauss(4000)


def hdr(*a):
    print(*a, flush=True)


# ================================================================== N=2, exact 1-dim
def rho_su2(a, beta, theta, n=4000):
    """rho_v(a) = int_{<e3,xi>=a} e^{-S} J  d(area) ; exact up to Gauss-Legendre truncation."""
    R2 = 4.0 * np.pi ** 2 - a * a
    if R2 <= 0.0:
        return 0.0
    R = np.sqrt(R2)
    x, w = _GLX[:n], _GLW[:n]
    s = 0.5 * R * (x + 1.0)
    wq = 0.5 * R * w
    r = np.hypot(s, a)
    h = 0.5 * r
    u = a / r
    S = beta * (1.0 - np.cos(h) * np.cos(0.5 * theta) + np.sin(h) * np.sin(0.5 * theta) * u)
    J = np.where(h > 1e-12, (np.sin(h) / h) ** 2, 1.0)
    return 2.0 * np.pi * float(np.sum(wq * s * np.exp(-S) * J))


def lrho_su2(a, beta, theta):
    v = rho_su2(a, beta, theta)
    return np.log(v) if v > 0.0 else -np.inf


def fd2(f, a, h):
    """5-point stencil for f''(a)."""
    return (-f(a - 2 * h) + 16.0 * f(a - h) - 30.0 * f(a) + 16.0 * f(a + h) - f(a + 2 * h)) / (12.0 * h * h)


def Gpp_su2(a, beta, theta, h=2e-2):
    """G''(a) = -d_a^2 log rho_v(a)."""
    return -fd2(lambda x: lrho_su2(x, beta, theta), a, h)


# ================================================================== N=3, symmetry-reduced 3-dim
def _su3_H(t, r, mu, q):
    """H(t, r, mu, q) on the U(2)-reduced slice; arrays broadcast.  t=a/sqrt3, tr(T8 H)=sqrt3 t."""
    shp = np.broadcast_shapes(np.shape(r), np.shape(mu), np.shape(q))
    H = np.zeros(shp + (3, 3))
    H[..., 0, 0] = r * mu + t
    H[..., 1, 1] = -r * mu + t
    H[..., 2, 2] = -2.0 * t
    s = np.sqrt(np.maximum(1.0 - mu * mu, 0.0)) * r
    H[..., 0, 1] = s
    H[..., 1, 0] = s
    sq = np.sqrt(np.maximum(q, 0.0))
    H[..., 0, 2] = sq
    H[..., 2, 0] = sq
    return H


def _su3_nodes(nr, nq, nm, rmax, qmax):
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xq, wq = np.polynomial.legendre.leggauss(nq)
    xm, wm = np.polynomial.legendre.leggauss(nm)
    r = 0.5 * rmax * (xr + 1.0)
    q = 0.5 * qmax * (xq + 1.0)
    W = (0.5 * rmax * wr)[:, None, None] * r[:, None, None] ** 2 \
        * (0.5 * qmax * wq)[None, :, None] * q[None, :, None] \
        * wm[None, None, :]
    R, Q, MU = np.meshgrid(r, q, xm, indexing='ij')
    return R, Q, MU, W


def rho_su3(a, beta, theta=np.pi, nr=64, nq=64, nm=40, rmax=None, qmax=None, test_k=None):
    """rho_v(a) for the U(2)-invariant SU(3) plaquette; test_k!=None -> substitute exp(-k tr H^2) (S1)."""
    t = a / np.sqrt(3.0)
    if rmax is None:
        rmax = 4.0 * np.pi
    if qmax is None:
        qmax = (4.0 * np.pi) ** 2
    R, Q, MU, W = _su3_nodes(nr, nq, nm, rmax, qmax)
    H = _su3_H(t, R, MU, Q)
    sh = H.shape
    lam, V = np.linalg.eigh(H.reshape(-1, 3, 3))
    lam = lam.reshape(sh[:-2] + (3,))
    V = V.reshape(sh)
    if test_k is not None:                                    # [S1] pure measure test
        return 2.0 * np.pi ** 3 * float(np.sum(W * np.exp(-test_k * np.einsum('...i,...i->...', lam, lam))))
    d = lam[..., :, None] - lam[..., None, :]
    ij = np.triu_indices(3, 1)
    dij = d[..., ij[0], ij[1]]
    half = 0.5 * dij
    J = np.prod(np.where(np.abs(half) > 1e-10, (np.sin(half) / half) ** 2, 1.0), axis=-1)
    keep = np.max(np.abs(dij), axis=-1) <= TWO_PI
    # S = beta(1 - ReTr(U B)/3),  B = diag(e^{i p}, e^{i p}, e^{-2ip}),  p = theta/sqrt12
    p = theta / np.sqrt(12.0)
    ph = V ** 2
    e = np.exp(1j * lam)
    diagU = np.einsum('...kj,...j->...k', ph, e)               # (e^{iH})_kk = sum_j |V_kj|^2 e^{i lam_j}
    Bv = np.exp(1j * np.array([p, p, -2.0 * p]))
    f = np.real(np.einsum('...k,k->...', diagU, Bv))
    S = beta * (1.0 - f / 3.0)
    wgt = W * np.exp(-S) * J * keep
    return 2.0 * np.pi ** 3 * float(np.sum(wgt))


def Gpp_su3(a, beta, theta=np.pi, h=0.2, **kw):
    f = lambda x: np.log(max(rho_su3(x, beta, theta, **kw), 1e-300))
    return -fd2(f, a, h)


# ================================================================== reporting helpers
def spark(vals, width=64):
    """ASCII curve: signed values -> one line, zero marked with '|'."""
    v = np.asarray(vals, float)
    if len(v) > width:
        idx = np.linspace(0, len(v) - 1, width).round().astype(int)
        v = v[idx]
    m = np.max(np.abs(v)) or 1.0
    lo, hi = -1.0, 1.0
    out = []
    N = 21
    for x in v:
        k = int(round((x / m - lo) / (hi - lo) * (N - 1)))
        row = [' '] * N
        row[10] = '|'
        row[k] = '*' if k != 10 else '+'
        out.append(''.join(row))
    return out


def sign_changes(a, g):
    """interpolated zero crossings of g over the a grid"""
    a = np.asarray(a, float)
    g = np.asarray(g, float)
    out = []
    for i in range(len(a) - 1):
        if g[i] == 0.0:
            out.append(a[i])
        elif g[i] * g[i + 1] < 0.0:
            out.append(a[i] - g[i] * (a[i + 1] - a[i]) / (g[i + 1] - g[i]))
    return out


# ================================================================== main
def main():
    t0 = time.time()
    hdr("#" * 104)
    hdr("# o02_plaquette -- G''(a) = -d^2/da^2 log rho_v(a) on the WHOLE a-axis, single plaquette N=2 / N=3")
    hdr("#" * 104)

    # ---------------- [S1] SU(3) reduction validation
    hdr("\n[S1] SU(3) REDUCTION CHECK:  int_{R^3xR^4} f dA d^4w  vs  2 pi^3 int r^2 dr q dq dmu f")
    hdr("     test f = exp(-k Tr H^2) with Tr H^2 = 2r^2 + 2q + 6t^2 -> closed form via 1-dim integrals")
    for k in (0.3, 1.0):
        t, Rm, Qm = 0.37, 6.0, 20.0
        red = rho_su3(np.sqrt(3.0) * t, 0.0, test_k=k, nr=200, nq=200, nm=1, rmax=Rm, qmax=Qm)
        red /= np.exp(-k * 6 * t * t)
        from scipy.special import erf
        a2 = 2.0 * k
        Ir = np.sqrt(np.pi) / (4.0 * a2 ** 1.5) * erf(np.sqrt(a2) * Rm) - Rm / (2.0 * a2) * np.exp(-a2 * Rm * Rm)
        Iq = (1.0 - np.exp(-a2 * Qm) * (1.0 + a2 * Qm)) / a2 ** 2
        cf = 4.0 * np.pi ** 3 * Ir * Iq
        hdr(f"     k={k:4.2f}: reduced = {red:.10e}   closed form = {cf:.10e}   rel.diff = {abs(red/cf-1):.2e}")

    # ---------------- [V] N=2 validation
    hdr("\n[V] N=2 VALIDATION (archived: 0.20511530 / +0.135995 / -0.193019 / beta*=4.30661)")
    for beta, th, ref in ((0.0, np.pi, 0.20511530), (2.5, np.pi, 0.135995), (6.0, np.pi, -0.193019),
                          (2.5, 0.0, 0.729924), (6.0, 0.0, 1.590639)):
        g = Gpp_su2(0.0, beta, th)
        hdr(f"     beta={beta:5.2f} th={th:6.3f}  G''(0) = {g:+.6f}   archived {ref:+.6f}   diff {g-ref:+.2e}")
    mass = 0.0
    ax = np.linspace(-TWO_PI * 0.9999, TWO_PI * 0.9999, 4001)
    rv = np.array([rho_su2(x, 0.0, np.pi) for x in ax])
    mass = np.trapezoid(rv, ax)
    hdr(f"     beta=0 mass sum rule:  int rho da = {mass:.6f}   16*pi^2 = {16*np.pi**2:.6f}"
        f"   rel.diff {abs(mass/(16*np.pi**2)-1):.2e}")
    # beta*
    bs = np.linspace(3.0, 6.0, 61)
    gs = np.array([Gpp_su2(0.0, b, np.pi) for b in bs])
    zc = sign_changes(bs, gs)
    hdr(f"     beta* from G''(0)=0 (th=pi): {zc[0]:.5f}   archived 4.30661")
    # closed form check: G''(0) = W0 - V0 beta^2,  W0=1/(4Z), V0 = ...
    W0 = 0.2051152971
    hdr(f"     two-term law check: G''(0,b=6) = W0 - V0*36 = {W0-36*(W0-0.135995)/6.25:+.6f}")

    # ---------------- [A] N=2 whole curve
    hdr("\n[A] N=2: WHOLE CURVE  G''(a),  h_FD = 0.02,  a in (-2pi, 2pi) = injective-chart range")
    hdr("     (the chart edge |a| -> 2pi is the antipode: rho -> 0, G'' -> +inf; grid stops 0.08 short)")
    amax = TWO_PI - 0.08
    agrid = np.linspace(-amax, amax, 129)
    curves = {}
    hdr(f"\n     {'a':>8}" + "".join(f"{n:>12}" for n in ('b2.5 th=pi', 'b6 th=pi', 'b2.5 th=0', 'b6 th=0')))
    for name, beta, th in (('2.5pi', 2.5, np.pi), ('6pi', 6.0, np.pi), ('2.50', 2.5, 0.0), ('60', 6.0, 0.0)):
        curves[name] = np.array([Gpp_su2(a, beta, th) for a in agrid])
    for i in range(0, len(agrid), 4):
        hdr(f"     {agrid[i]:8.3f}" + "".join(f"{curves[k][i]:12.5f}" for k in ('2.5pi', '6pi', '2.50', '60')))
    hdr("\n     sign of G''(a):  '+' positive, '-' negative, '0' |G''|<1e-3   (a from -2pi to +2pi)")
    for k in ('2.5pi', '6pi', '2.50', '60'):
        s = ''.join('+' if x > 1e-3 else ('-' if x < -1e-3 else '0') for x in curves[k])
        hdr(f"       {k:>7}: {s}")
    for k in ('2.5pi', '6pi', '2.50', '60'):
        sc = sign_changes(agrid, curves[k])
        hdr(f"       {k:>7}: G''(0)={curves[k][len(agrid)//2]:+9.6f}  min={curves[k].min():+9.6f}"
            f"  sign changes at a = " + (", ".join(f"{x:+.4f}" for x in sc) if sc else "NONE"))

    # ---------------- [B] beta dependence of the sign pattern
    hdr("\n[B] N=2, theta=pi: does the WHOLE curve flip at beta*, or only a=0?")
    hdr(f"     {'beta':>6}{'G(0)':>10}{'min_a G':>10}{'argmin':>8}{'G(-pi)':>10}{'G(+pi)':>10}"
        f"{'G(1)':>9}{'G(2)':>9}{'G(2.9)':>9}   sign pattern")
    betas = (0.0, 1.0, 2.5, 4.0, 4.30661, 5.0, 6.0, 8.0, 10.0)
    bgrid = np.linspace(-amax, amax, 121)
    for beta in betas:
        gg = np.array([Gpp_su2(a, beta, np.pi) for a in bgrid])
        im = int(np.argmin(gg))
        neg = int(np.sum(gg < -1e-3))
        hdr(f"     {beta:6.2f}{gg[len(bgrid)//2]:+10.5f}{gg[im]:+10.5f}{bgrid[im]:8.3f}"
            f"{Gpp_su2(-np.pi, beta, np.pi):+10.5f}{Gpp_su2(np.pi, beta, np.pi):+10.5f}"
            f"{Gpp_su2(1.0, beta, np.pi):+9.5f}{Gpp_su2(2.0, beta, np.pi):+9.5f}"
            f"{Gpp_su2(2.9, beta, np.pi):+9.5f}   {neg}/{len(bgrid)} negative")

    # ---------------- [C] N=3 whole curve
    hdr("\n[C] N=3 (U(2)-invariant Cartan v=T8, B=exp(i pi T8)): WHOLE CURVE G''(a)")
    hdr("     t = a/sqrt3 (chart needs |t|<~2.09 at r=q=0); h_FD = 0.2;  G'' = -d^2 log rho/da^2")
    tg = np.arange(-1.75, 1.751, 0.25)
    ag3 = np.sqrt(3.0) * tg
    hdr(f"\n     {'a':>8}{'t':>7}{'G b=2.5':>11}{'G b=6':>11}{'G b=0':>10}"
        f"{'  [node x1 vs x2, b=2.5]':>26}")
    out = {'2.5': [], '6': [], '0': []}
    for a, t in zip(ag3, tg):
        g25 = Gpp_su3(a, 2.5, h=0.2, nr=64, nq=64, nm=40)
        g6 = Gpp_su3(a, 6.0, h=0.2, nr=64, nq=64, nm=40)
        g0 = Gpp_su3(a, 0.001, h=0.2, nr=64, nq=64, nm=40)
        out['2.5'].append(g25)
        out['6'].append(g6)
        out['0'].append(g0)
        hdr(f"     {a:8.4f}{t:7.3f}{g25:+11.5f}{g6:+11.5f}{g0:+10.5f}")
    hdr("\n     node refinement (a=0 and a=+-pi, beta=2.5): nr x nq x nm = 48x48x32 / 64x64x40 / 88x88x56")
    for a in (0.0, -np.sqrt(3.0) * 1.5, np.sqrt(3.0) * 1.5):
        row = ""
        for (nr, nq, nm) in ((48, 48, 32), (64, 64, 40), (88, 88, 56)):
            row += f"   {Gpp_su3(a, 2.5, h=0.2, nr=nr, nq=nq, nm=nm):+.5f}"
        hdr(f"       a={a:+.4f}: {row}")
    for k in ('0', '2.5', '6'):
        gg = np.array(out[k])
        sc = sign_changes(ag3, gg)
        hdr(f"     beta={k:>4}: G''(0)={out[k][len(ag3)//2]:+.5f}  min={gg.min():+.5f}"
            f"  sign changes: " + (", ".join(f"{x:+.3f}" for x in sc) if sc else "NONE"))

    hdr("\n[D] ASCII CURVES  (row index = value axis; '*' = G'', '|' = G''=0; a runs left->right)")
    for k, lab in (('2.5pi', 'N=2 beta=2.5 th=pi'), ('6pi', 'N=2 beta=6.0 th=pi')):
        for row in spark(curves[k]):
            hdr(f"     {lab:>22} |{row}|")
        hdr("")
    for k, lab in (('2.5', 'N=3 beta=2.5'), ('6', 'N=3 beta=6.0')):
        for row in spark(out[k]):
            hdr(f"     {lab:>22} |{row}|")
        hdr("")
    hdr(f"\nwall {time.time()-t0:.1f}s   Clay 0%.  No QED.")


if __name__ == '__main__':
    main()
