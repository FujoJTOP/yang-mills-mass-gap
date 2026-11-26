#!/usr/bin/env python3
r"""
o02_p4_2dym.py -- G''(a) = -d_a^2 log rho_v(a) on the WHOLE a-axis for the two other exactly solvable
models: (i) the P4 PERIODIC constant-mode transfer kernel, (ii) 2D Yang-Mills (Migdal).

(i) P4 PERIODIC (research/v3/w7/p4_periodic_transfer.py, P4-PERIODIC-TRANSFER.md):
    cyclic chain of L SU(2) links, flat Cartan chart z_l, exact one-link Haar marginal
        rho1(z) = pi (Cin(2pi) - Cin|z|),  |z|<2pi
    S_W = beta sum_i (1-cos((z_{i+1}-z_i)/2)),  constant mode a = (1/sqrt L) sum_l z_l
    Z(k) = <e^{-i k a}> = Tr[M_k^L],  M_k(z,z') = sqrt(rho1 rho1') e^{-i k(z+z')/(2 sqrt L)} e^{-beta(1-cos((z'-z)/2))}
    M_k is similar to the HERMITIAN D^{1/2} M_0 D^{1/2} (D=diag(e^{-ikz/(2 sqrt L)})), so Z(k) is real:
        Z(k) = sum_n lambda_n(k)^L ,  lambda_n real.
    rho_v(a) = (1/2pi) int Z(k) e^{i k a} dk ,   G''(a) = -d_a^2 log rho_v(a).
    (The archived number is G''(0); here we get the whole function.)

(ii) 2D YM (Migdal/heat kernel, exact):  rho_t(phi) = K_t(phi) sin^2 phi,
    K_t = sum_j (2j+1) chi_j(phi) e^{-t j(j+1)},  chi_j = sin((2j+1)phi)/sin phi,   t = g^2 A/2.
    G''(phi) = -(K''/K - (K'/K)^2) + 2/sin^2 phi     (analytic derivatives of the sum, no FD).
    phi=0 is the singular point (J=sin^2 phi -> 0): the conditioning point is NOT phi=0 here.

SECTIONS
  [P0] P4 machinery check: Z(k) vs the beta=0 closed form; nz convergence of G''(0) vs archived numbers.
  [P1] P4 whole curve G''(a), beta=2.5, L=4/6/8 + sign structure.
  [P2] 2D YM whole curve G''(phi), several t; min over phi; t->inf limit.
  [P3] ASCII curves.
Run: py -3.12 o02_p4_2dym.py [quick|full]
Labels [T]/[C].  No QED.  Clay 0%.
"""
import sys
import time
import numpy as np
from scipy.special import sici

np.seterr(invalid='ignore', divide='ignore')
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi


def hdr(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- exact one-link Haar marginal
def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)


def rho1(z):
    z = np.asarray(z, float)
    az = np.abs(z)
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)


def nodes(nz):
    xg, wg = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * xg
    return z, rho1(z) * TWO_PI * wg


def Zk_periodic(L, beta, kgrid, z, r1w):
    """Z(k)=Tr[M_k^L].  M_k = D M_0 D (D=diag(e^{-ikz/(2 sqrtL)})) is complex SYMMETRIC, NOT Hermitian
    (D^{1/2}M_0D^{1/2} is only *similar* to it; it is not Hermitian either) -> general eigvals.
    Z(k) is real by the z -> -z symmetry; we keep .real."""
    sqL = np.sqrt(L)
    s = np.sqrt(np.maximum(r1w, 0.0))
    K0 = np.exp(beta * np.cos(0.5 * (z[None, :] - z[:, None])))
    M0 = (s[:, None] * K0) * s[None, :]
    M0 = 0.5 * (M0 + M0.T)
    pref = np.exp(-beta * L)
    Z = np.empty(len(kgrid))
    for i, k in enumerate(kgrid):
        ph = np.exp(-1j * k * z / (2.0 * sqL))
        Mk = M0 * (ph[:, None] * ph[None, :])
        ev = np.linalg.eigvals(Mk)
        Z[i] = pref * np.sum(ev ** L).real
    return Z


def fd2(f, a, h):
    return (-f(a - 2 * h) + 16.0 * f(a - h) - 30.0 * f(a) + 16.0 * f(a + h) - f(a + 2 * h)) / (12.0 * h * h)


def sign_changes(a, g):
    a = np.asarray(a, float)
    g = np.asarray(g, float)
    out = []
    for i in range(len(a) - 1):
        if g[i] * g[i + 1] < 0.0:
            out.append(a[i] - g[i] * (a[i + 1] - a[i]) / (g[i + 1] - g[i]))
    return out


def spark(vals, width=64):
    v = np.asarray(vals, float)
    if len(v) > width:
        idx = np.linspace(0, len(v) - 1, width).round().astype(int)
        v = v[idx]
    m = np.max(np.abs(v)) or 1.0
    out = []
    N = 21
    for x in v:
        k = int(round(((x / m) + 1.0) / 2.0 * (N - 1)))
        row = [' '] * N
        row[10] = '|'
        row[k] = '+' if k == 10 else '*'
        out.append(''.join(row))
    return out


# ---------------------------------------------------------------- 2D YM exact
def gpp_2dym(phi, t, JMAX=140):
    r"""G''(phi) = -d^2/dphi^2 log[K_t(phi) sin^2 phi], analytic derivatives of the character sum.
    f_j = A sin(A phi)/sin(phi), A=2j+1:
      f_j'' = A[ -A^2 sin(Aphi)/S - 2 A cos(Aphi)cos(phi)/S^2 + sin(Aphi)sin(phi)/S^2
                + 2 sin(Aphi)cos^2(phi)/S^3 ]
    and  -2 d^2/dphi^2 log sin(phi) = +2/sin^2(phi)."""
    phi = np.asarray(phi, float)
    S, C = np.sin(phi), np.cos(phi)
    j = np.arange(JMAX + 1)
    A = 2.0 * j + 1.0
    E = np.exp(-t * j * (j + 1.0))
    s = np.sin(np.outer(phi, A))
    cs = np.cos(np.outer(phi, A))
    Sm, Cm = S[:, None], C[:, None]
    c0 = A * s / Sm * E
    c1 = (-A * A * cs / Sm + A * s * Cm / Sm ** 2) * E
    c2 = (A * (-A ** 2 * s / Sm - 2.0 * A * cs * Cm / Sm ** 2
               + s * Sm / Sm ** 2 + 2.0 * s * Cm ** 2 / Sm ** 3)) * E
    K0 = c0.sum(1)
    K1 = c1.sum(1)
    K2 = c2.sum(1)
    return -(K2 / K0 - (K1 / K0) ** 2) + 2.0 / S ** 2


def main():
    t0 = time.time()
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    hdr("#" * 104)
    hdr("# o02_p4_2dym -- G''(a) = -d^2/da^2 log rho_v(a) on the WHOLE a-axis: P4 periodic kernel, 2D YM")
    hdr("#" * 104)

    # ---------------- [P0]
    hdr("\n[P0] P4 MACHINERY CHECK")
    nz = 160 if mode != 'full' else 260
    z, r1w = nodes(nz)
    hdr(f"     nz = {nz} Gauss nodes on |z|<2pi   (archived runs used 320/460)")
    dk = 0.0125 if mode != 'full' else 0.008
    K = 8.0
    nh = int(round(K / dk))
    kgrid = np.arange(-nh, nh + 1) * dk
    hdr(f"     k-grid: {kgrid[0]:+.2f}..{kgrid[-1]:+.2f}, dk={dk}  ({len(kgrid)} points)")

    def Z_beta0_closed(L, kk):
        sqL = np.sqrt(L)
        return (np.exp(-1j * np.outer(kk, z) / sqL) @ r1w) ** L

    for L in (4,):
        Zc = Zk_periodic(L, 0.0, np.array([0.3, 1.1, 2.5]), z, r1w)
        Zb = Z_beta0_closed(L, np.array([0.3, 1.1, 2.5]))
        hdr(f"     beta=0, L={L}, Z(k) eigvals vs closed form: " +
            "  ".join(f"k={k:3.1f}: {a:+.8f}/{b:+.8f}" for k, a, b in zip((0.3, 1.1, 2.5), Zc, Zb)))

    # convergence of G''(0) (moment ratio) against archived 0.22196 (L->inf) / 0.2165 (L=8)
    hdr("\n     G''(0) via moment ratio int k^2 Z / int k Z  (archived: L=4 -> 0.2115, L=8 -> 0.2165, fit c=0.22196)")
    Zs = {}
    for L in (4, 8):
        Z = Zk_periodic(L, 2.5, kgrid, z, r1w)
        Zs[L] = Z
        g = np.trapezoid(kgrid ** 2 * Z, kgrid) / np.trapezoid(Z, kgrid)
        hdr(f"       L={L}: G''(0) = {g:.5f}")

    # ---------------- [P1] whole curve
    hdr("\n[P1] P4: WHOLE CURVE G''(a), beta=2.5.  rho via inverse FT; G'' by 5-point FD in a (h=0.02)")
    hdr("     (a -> +-2 pi sqrt(L) is the support edge; the curve is even in a by z -> -z symmetry)")
    hstep = 0.02
    agrid = np.arange(-4.0, 4.0001, 0.05)
    curves = {}
    for L in (4, 8):
        Z = Zs[L]

        def lrho(x, Z=Z):
            v = (dk / TWO_PI) * np.real(np.exp(1j * x * kgrid) @ Z)
            return np.log(v) if v > 0 else -np.inf

        curves[L] = np.array([-fd2(lrho, a, hstep) for a in agrid])
        sc = sign_changes(agrid, curves[L])
        hdr(f"     L={L}: G''(0)={curves[L][len(agrid)//2]:+.6f}  min={curves[L].min():+.6f}"
            f"  max={curves[L].max():+.6f}  sign changes: " + (", ".join(f"{x:+.3f}" for x in sc) if sc else "NONE"))
    hdr(f"\n     {'a':>8}{'L=4':>12}{'L=8':>12}")
    for i in range(0, len(agrid), 8):
        hdr(f"     {agrid[i]:8.3f}{curves[4][i]:12.6f}{curves[8][i]:12.6f}")
    for L in (4, 8):
        s = ''.join('+' if x > 1e-3 else ('-' if x < -1e-3 else '0') for x in curves[L])
        hdr(f"     L={L} sign pattern (a=-4..4): {s}")
    # a few off-origin checks with a coarser FD step (smoothing test)
    for L in (8,):
        Z = Zk_periodic(L, 2.5, kgrid, z, r1w)
        for a in (0.5, 1.5, 2.5, 3.5):
            row = ""
            for h in (0.01, 0.02, 0.05, 0.1):
                def lrho(x, Z=Z):
                    v = (dk / TWO_PI) * np.real(np.exp(1j * x * kgrid) @ Z)
                    return np.log(max(v, 1e-300))
                row += f"  h={h:.2f}:{-fd2(lrho, a, h):+.5f}"
            hdr(f"     L=8 a={a:+.2f} FD-step stability:{row}")

    # ---------------- [P2] 2D YM
    hdr("\n[P2] 2D YM: WHOLE CURVE G''(phi), analytic (no FD).  phi=0 is singular: G''-> +inf")
    phig = np.linspace(0.02, np.pi - 0.02, 100)
    hdr(f"\n     {'t':>7}{'G(.3)':>10}{'G(.6)':>10}{'G(1.0)':>10}{'G(1.5)':>10}{'G(2.0)':>10}"
        f"{'G(2.6)':>10}{'min':>9}{'argmin':>8}{'all>0?':>8}")
    for t in (0.5, 1.0, 2.0, 5.0, 40.0, 200.0):
        vals = gpp_2dym(np.array([0.3, 0.6, 1.0, 1.5, 2.0, 2.6]), t)
        full = gpp_2dym(phig, t)
        hdr(f"     {t:7.1f}" + "".join(f"{v:10.4f}" for v in vals) +
            f"{full.min():9.4f}{phig[int(np.argmin(full))]:8.3f}{str(bool(np.all(full > 0))):>8}")
    hdr("\n     (archived H1A table for cross-check: t=40 -> 22.9011 6.2731 2.8246 2.0101 2.4189 7.5261)")
    hdr(f"     minimum over phi in (0,pi), scanned in t:")
    hdr(f"       {'t':>8}{'min_phi G':>12}{'argmin':>9}{'  any negative in phi?':>24}")
    for t in (0.1, 0.5, 1, 2, 5, 10, 20, 40, 100, 400):
        full = gpp_2dym(phig, t)
        hdr(f"       {t:8.1f}{full.min():12.5f}{phig[int(np.argmin(full))]:9.3f}"
            f"{('  NONE (all > 0)' if full.min() > 0 else '  WINDOW with G<0'):>24}")
    # convergence in jmax
    for t in (40.0,):
        for JMAX in (60, 100, 140):
            v = gpp_2dym(np.array([1.5]), t, JMAX=JMAX)[0]
            hdr(f"     jmax check t=40 phi=1.5: jmax={JMAX} -> {v:.6f}")

    # ---------------- [P3]
    hdr("\n[P3] ASCII CURVES ('*'=G'', '|'=zero)")
    for L in (8,):
        for row in spark(curves[L]):
            hdr(f"     {'P4 L=8 b=2.5':>16} |{row}|")
    for t in (2.0, 40.0):
        for row in spark(gpp_2dym(phig, t)):
            hdr(f"     {'2DYM t='+str(t):>16} |{row}|")
    hdr(f"\nwall {time.time()-t0:.1f}s   Clay 0%.  No QED.")


if __name__ == '__main__':
    main()
