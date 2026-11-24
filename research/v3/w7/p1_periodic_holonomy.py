#!/usr/bin/env python3
r"""
p1_periodic_holonomy.py -- P1: the EXACT global constant-holonomy anchor on the PERIODIC torus.

QUESTION (periodic torus, no twist, no fixed boundary): does the constant-mode curvature at
beta=2.5 saturate at a positive L-independent value, or decay?

TASK.  The exact global constant-holonomy anchor.
 [A] Define the global constant link configuration; compute every plaquette; show S_W == 0 exactly.
 [B] The marginal is therefore pure Haar on the single modulus U.  Exact curvature
         G''_anch(0) = -rho_Haar''(0)/rho_Haar(0)
     for SU(2) in CLOSED FORM, and for SU(3) by direct Haar sampling.  Var_Haar for both.
     *** The exact object is -rho''/rho, NOT 1/Var.  They agree only in the GAUSSIAN limit. ***
 [C] Object B (k=0 FOURIER mode of the fluctuations) vs object A (global holonomy): the k=0
     marginal is the L-fold convolution and Gaussianises -> 1/Var_Haar.  A != B.
 [D] ACID TEST (free field): U(1) with flat Haar (J=1) gives G''=0; 1/Var = 3/pi^2 is a SPURIOUS
     floor -> proves the exact object is -rho''/rho, not 1/Var.
 [E] SU(3): Var_Haar(xi3) and G''_anch(0), cross-validated on SU(2) where the exact value is known.

Run: py -3.12 p1_periodic_holonomy.py [quick|full].  numpy only; one job at a time; small L.
No third-party sources.  Labels [T]/[V]/[C].  No QED.  Clay 0%.
"""
import sys
import numpy as np
from scipy.special import sici

np.seterr(all='ignore')
MODE = sys.argv[1] if len(sys.argv) > 1 else 'quick'
NMC = {'quick': 1_500_000, 'full': 8_000_000}[MODE]
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi

# ============================================================ SU(2) Haar marginal (exact)
def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)

def rho1(z):
    r"""exact one-Cartesian-component Haar marginal of U=exp(i xi.sigma/2):
       rho1(z) = 2pi \int_{|z|}^{2pi} sin^2(rho/2)/rho drho = pi ( Cin(2pi) - Cin(|z|) ),  |z|<=2pi."""
    z = np.asarray(z, float); az = np.abs(z)
    val = np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI)))
    return np.where(az <= TWO_PI, val, 0.0)

def var_rho1(n=400001):
    z = np.linspace(-TWO_PI, TWO_PI, n)
    w = np.full(n, 2 * TWO_PI / (n - 1)); w[[0, -1]] *= 0.5
    r1 = rho1(z) * w
    Z = r1.sum(); m = (z * r1).sum() / Z
    return ((z - m) ** 2 * r1).sum() / Z

# ============================================================ [A] plaquette / action check
def su2_from_vec(vec):
    """SU(2) from a 3-vector vec (angle = |vec|), U = exp(i vec.sigma/2)."""
    th = np.linalg.norm(vec)
    if th < 1e-12:
        return np.eye(2, dtype=complex)
    n = vec / th
    return np.cos(th / 2) * np.eye(2) + 1j * np.sin(th / 2) * (n[0] * np.array([[0, 1], [1, 0]]) +
           n[1] * np.array([[0, -1j], [1j, 0]]) + n[2] * np.array([[1, 0], [0, -1]]))

def plaquette_check(L=2, N=2):
    """Periodic L^4 torus, all links in direction 3 equal U, others identity.  All plaquettes = I,
    S_W = 0.  Also the all-links-equal-U config.  Returns max plaquette deviation from I for both."""
    rng = np.random.default_rng(0)
    U = su2_from_vec(rng.uniform(-np.pi, np.pi, 3))  # a generic SU(2) element
    I2 = np.eye(2, dtype=complex)
    def build(all_equal):
        links = {}
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        for mu in range(4):
                            if all_equal:
                                links[(x0, x1, x2, x3, mu)] = U
                            else:
                                links[(x0, x1, x2, x3, mu)] = U if mu == 3 else I2
        return links
    def sh(x, mu):
        y = list(x); y[mu] = (y[mu] + 1) % L; return tuple(y)
    def s_w(links):
        tot = 0.0; maxdev = 0.0
        site = [(a, b, c, d) for a in range(L) for b in range(L) for c in range(L) for d in range(L)]
        for x in site:
            for mu in range(4):
                for nu in range(mu + 1, 4):
                    xmu = sh(x, mu); xnu = sh(x, nu)
                    # P = U_mu(x) U_nu(x+mu) U_mu(x+nu)^-1 U_nu(x)^-1
                    P = (links[x + (mu,)] @ links[xmu + (nu,)]
                         @ links[xnu + (mu,)].conj().T @ links[x + (nu,)].conj().T)
                    maxdev = max(maxdev, np.max(np.abs(P - I2)))
                    tot += 1.0 - 0.5 * np.real(np.trace(P))   # S_W/ beta, per plaquette
        return tot, maxdev
    S_a, dev_a = s_w(build(False))
    S_b, dev_b = s_w(build(True))
    return dict(dir3=(S_a, dev_a), all_equal=(S_b, dev_b), U=U)

# ============================================================ [C] k=0 Fourier mode (object B)
def rho1_hat(q, nz=4000):
    r"""normalized characteristic function of rho1: \hat rho1(q) = < e^{-i q z} >_rho1."""
    zg, wz = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * zg; w = rho1(z) * TWO_PI * wz
    return (np.exp(-1j * np.outer(np.atleast_1d(q), z)) @ w) / w.sum()

def Z0(k, L):
    r"""k=0 mode a=(1/sqrt L)sum xi_l:  Z_a(k) = [ \hat rho1(k/sqrt L) ]^L."""
    return rho1_hat(np.atleast_1d(k) / np.sqrt(L)) ** L

def rho_from_Z(Zfun, half, nk=1 << 15, kw=40.0):
    dk = 2 * kw / nk; kj = (np.arange(nk) - nk // 2) * dk
    Z = Zfun(kj) * np.exp(-(kj / kw) ** 8)          # smooth window kills truncation ringing
    rho = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Z))).real * nk * dk / (2 * np.pi)
    da = 2 * np.pi / (nk * dk); a = (np.arange(nk) - nk // 2) * da
    return a, rho, da

def Gpp0_from_Z(Zfun, half, dp=4):
    """G'' at a=0 (the symmetric mode) of the k=0 marginal, local quadratic fit over +-dp points."""
    a, rho, da = rho_from_Z(Zfun, half)
    i0 = np.argmin(np.abs(a)); o = np.arange(i0 - dp, i0 + dp + 1)
    lg = np.log(np.maximum(rho[o], 1e-300))
    co = np.polyfit(a[o], lg, 2)
    return a[i0], -2 * co[0]

# ============================================================ [E] Haar samplers (flat-chart Cartan)
def haar_su2(n, rng):
    g = rng.standard_normal((n, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
    q0, q1, q2, q3 = g[:, 0], g[:, 1], g[:, 2], g[:, 3]
    U = np.empty((n, 2, 2), complex)
    U[:, 0, 0] = q0 + 1j * q3; U[:, 0, 1] = q2 + 1j * q1
    U[:, 1, 0] = -q2 + 1j * q1; U[:, 1, 1] = q0 - 1j * q3
    return U

def log_su2_xi(U):
    """principal xi in U=exp(i xi.sigma/2), returns (n,3)."""
    tr = np.clip(0.5 * np.real(U[:, 0, 0] + U[:, 1, 1]), -1, 1)
    q1 = np.imag(U[:, 0, 1]); q2 = np.real(U[:, 0, 1]); q3 = np.imag(U[:, 0, 0])
    s = np.sqrt(np.maximum(q1 * q1 + q2 * q2 + q3 * q3, 1e-300)); ang = 2 * np.arccos(tr)
    return ang[:, None] * np.stack([q1, q2, q3], axis=1) / s[:, None]

def haar_su3(n, rng):
    Z = (rng.standard_normal((n, 3, 3)) + 1j * rng.standard_normal((n, 3, 3))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    return Q * (d / np.abs(d))[:, None, :]

_L3 = np.diag([1, -1, 0]).astype(complex)
def log_su3_xi3(U):
    """Hm = -i log U (Hermitian, principal); xi^3 = Re tr(Hm lambda_3)."""
    w, V = np.linalg.eig(U); phi = np.angle(w)
    Hm = V @ (phi[:, :, None] * np.conjugate(np.swapaxes(V, 1, 2)))
    return np.real(np.einsum('nij,ji->n', Hm, _L3))

def sample_cartan(group, Ntot, rng, chunk=200_000):
    """single-link (global holonomy) Cartan amplitude samples over Haar."""
    out = []
    done = 0
    while done < Ntot:
        n = min(chunk, Ntot - done); done += n
        if group == 2:
            out.append(log_su2_xi(haar_su2(n, rng))[:, 2])
        else:
            out.append(log_su3_xi3(haar_su3(n, rng)))
    return np.concatenate(out)

def kde_curv0(a, h):
    """Gaussian-KDE estimate of G''(0) = -rho''(0)/rho(0)  (mode at 0 by symmetry)."""
    x = a / h
    phi = np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
    rho0 = phi.sum() / (len(a) * h)
    rho0pp = ((x * x - 1) * phi).sum() / (len(a) * h ** 3)
    return -rho0pp / rho0

# ============================================================ main
def main():
    print("#" * 100)
    print("# P1 -- exact GLOBAL CONSTANT-HOLONOMY anchor on the PERIODIC torus; object A vs k=0 mode B")
    print("#" * 100)

    print("\n[A] PLAQUETTE / ACTION CHECK (periodic L=2 4D torus, generic SU(2) element U):")
    pc = plaquette_check(L=2)
    print(f"    dir-3 = U, rest = I :  S_W/beta = {pc['dir3'][0]:+.3e}   max|P-I| = {pc['dir3'][1]:.3e}")
    print(f"    all links = U       :  S_W/beta = {pc['all_equal'][0]:+.3e}   max|P-I| = {pc['all_equal'][1]:.3e}")
    print("    => every plaquette is the identity, S_W == 0 EXACTLY, for every L, beta and U.  [T]")

    print("\n" + "=" * 100)
    print("[B] SU(2) EXACT GLOBAL-HOLONOMY CURVATURE (single-link / single-modulus Haar marginal):")
    print("-" * 100)
    # exact closed form G''_anch(0) = 1/(2 Cin(2pi))
    G_exact = 1.0 / (2.0 * Cin(TWO_PI))
    # independent numeric check by fine finite differences on rho1
    z = np.linspace(-0.02, 0.02, 9)
    lg = np.log(rho1(z)); h = z[1] - z[0]
    G_fd = -(-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * h * h)
    v1 = var_rho1()
    print(f"    exact  G''_anch(0) = 1/(2 Cin(2pi))            = {G_exact:.6f}")
    print(f"    numeric 5-pt FD of -rho1''/rho1 at 0           = {G_fd:.6f}")
    print(f"    Var_Haar(xi_3)  (exact (4pi^2-6)/9)            = {v1:.6f}   [(4pi^2-6)/9 = {(4*np.pi**2-6)/9:.6f}]")
    print(f"    1/Var_Haar  (Gaussian limit)                   = {1.0/v1:.6f}   (= 9/(4pi^2-6))")
    print(f"    *** G''_anch = -rho''/rho = {G_exact:.4f}  !=  1/Var = {1/v1:.4f}  (gap {1/v1-G_exact:+.4f}): the marginal is NON-Gaussian ***")

    print("\n" + "=" * 100)
    print("[C] OBJECT B -- the k=0 FOURIER mode a=(1/sqrt L) sum_l xi_{l,3} (L-fold convolution, beta=0):")
    print("-" * 100)
    print("     %3s%12s" % ('L', "G''_L(0)"))
    for L in (1, 2, 3, 4, 6, 8, 12, 16, 32):
        _, G = Gpp0_from_Z(lambda k, L=L: Z0(k, L), TWO_PI * np.sqrt(L))
        print(f"     {L:>3}{G:12.5f}")
    print(f"     L->inf (CLT)  ->  1/Var_Haar = {1.0/v1:.5f}   (object B decompactification floor, beta=0)")

    print("\n" + "=" * 100)
    print("[D] ACID TEST -- free field U(1), J=1 (flat Haar, uniform on [-pi,pi]):")
    print("-" * 100)
    # amplitude = theta uniform => rho const => rho''=0 => G''=0 ; Var=pi^2/3 => 1/Var=3/pi^2
    print(f"    exact G''_anch(0) = -rho''/rho = 0           (uniform density; the field is FREE)")
    print(f"    1/Var_Haar(theta) = 3/pi^2 = {3/np.pi**2:.4f}   <- would be a SPURIOUS floor; fails it")
    print("    => -rho''/rho PASSES the acid test (0); 1/Var FAILS it (0.304).  Exact object = -rho''/rho.  [T]")

    print("\n" + "=" * 100)
    print("[E] SU(3) (and SU(2) cross-check) via direct Haar sampling; flat-chart Cartan amplitude:")
    print("-" * 100)
    rng = np.random.default_rng(12345)
    for Ngrp, name in ((2, 'SU(2)'), (3, 'SU(3)')):
        a = sample_cartan(Ngrp, NMC, rng)
        var = float(np.var(a))
        hs = [0.10, 0.15, 0.20, 0.30]
        Gs = [kde_curv0(a, h) for h in hs]
        print(f"    {name}: Var_Haar(xi_3) = {var:.5f}   1/Var = {1.0/var:.5f}")
        print(f"          G''_anch(0) [KDE, h={hs}] = " + "  ".join(f"{g:.4f}" for g in Gs))
        if Ngrp == 2:
            print(f"          (exact SU(2) reference: Var={v1:.5f}, G''={G_exact:.5f})")
        else:
            print(f"          SU(3): G''_anch(0) ~ {np.median(Gs):.4f} (h-stable), 1/Var={1.0/var:.4f}")

    print("\n" + "=" * 100)
    print("VERDICT: the periodic global constant-holonomy anchor is EXACT and L-INDEPENDENT (S_W==0)")
    print("=> marginal pure Haar => G''_anch(0) = -rho_Haar''/rho_Haar  [T], every L; SU(2)=1/(2Cin2pi).")
    print("It is NOT 1/Var_Haar (that is the Gaussian / k=0 L->inf limit); acid test U(1) proves which")
    print("is exact.  The PHYSICAL decompactification soft mode is object B (the k=0 fluctuation mode),")
    print("whose beta=0 L->inf floor IS 1/Var_Haar; object A (global holonomy) is the flat-connection")
    print("moduli measure, an exact solvable anchor but NOT the soft-mode marginal.  No QED.  Clay 0%.")
    print("=" * 100)

if __name__ == '__main__':
    main()
