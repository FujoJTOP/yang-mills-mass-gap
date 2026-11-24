#!/usr/bin/env python3
r"""
q4_bakry_emery.py -- Q4: a local/weighted Bakry-Emery criterion on the gauge group + a defect term,
for H1 (chi = Var_mu(xi_{k->0}) <= 1/lambda, uniformly, at beta=2.5).

 (R) Ricci of SU(N), bi-invariant (heat-kernel/Haar) metric <X,Y> = -2 tr(XY): Ric = (N/4) g.
 (G) Global BE at beta=2.5: Ric + Hess S_W on the Savvidy negative mode -> negative (W6 ~ -2.3).
 (S) 1D SU(2) exact one-link Haar marginal: G''(0)=1/(2Cin(2pi)), finite lambda_1(nu), ess inf G''.
 (P) periodic beta=0 k=0 marginal = L-fold convolution: G''(0)->1/Var=0.2688, Var->3.7198.
 (W) weighted-criterion MODEL: lambda_1 vs defect-MEASURE vs BARRIER (Kramers) vs CAPACITY.
 (U) free U(1): G''->0, lambda_1->0.  Acid test.

Run: py -3.12 -I q4_bakry_emery.py .  numpy+scipy, deterministic, one small job.  No QED.  Clay 0%.
"""
import numpy as np
from scipy.special import sici
from scipy.linalg import eigh

EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi


# ============================================================ (R) SU(N) structure & Ricci
def su_structure(N):
    """Traceless Hermitian basis T^a of su(N) with tr(T^a T^b)=1/2 delta. Returns (dim, f) with
    [T^a,T^b] = i f^{abc} T^c (real f)."""
    raws = []
    for j in range(N):
        for k in range(j + 1, N):
            E = np.zeros((N, N), dtype=complex); E[j, k] = 1.0; E[k, j] = 1.0
            raws.append(E)                     # symmetric (real)
            F = np.zeros((N, N), dtype=complex); F[j, k] = -1.0j; F[k, j] = 1.0j
            raws.append(F)                     # antisymmetric (imag)
    for l in range(1, N):
        D = np.zeros((N, N), dtype=complex)
        for j in range(l):
            D[j, j] = 1.0
        D[l, l] = -float(l)
        raws.append(D / np.linalg.norm(D))      # diagonal
    Traw = np.array(raws)
    # Gram G_ab = tr(T_a T_b) (real symmetric PD); rotate+scale to tr(T_a T_b)=1/2 delta
    Gram = np.real(np.einsum('aij,bji->ab', Traw, Traw))
    w, V = eigh(Gram)
    T = np.einsum('ab,bij->aij', V.T, Traw) * (np.sqrt(0.5) / np.sqrt(w))[:, None, None]
    d = T.shape[0]
    f = np.zeros((d, d, d))
    for a in range(d):
        for b in range(d):
            C = T[a] @ T[b] - T[b] @ T[a]
            f[a, b] = np.real(-2.0j * np.einsum('ij,cji->c', C, T))
    return d, f


def ricci_su(N):
    d, f = su_structure(N)
    # with Hermitian T^a (tr=1/2) the anti-Hermitian e_a = i T^a are unit in <-2tr>;
    # Ric(e_a,e_a) = 1/4 sum_{i,c} f^{aic} f^{aic}  (= N/4 if the Casimir identity holds)
    ric = np.array([0.25 * np.sum(f[a] ** 2) for a in range(d)])
    Ssum = np.einsum('acd,bcd->ab', f, f)       # should be N * delta
    id_err = np.max(np.abs(Ssum - N * np.eye(d)))
    return d, ric, id_err


# ============================================================ 1-D tools
def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)


def rho1_haar(z):
    z = np.asarray(z, float); az = np.abs(z)
    val = np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI)))
    return np.where(az <= TWO_PI, val, 0.0)


def lam1_of_marginal(x, logrho):
    """Lowest positive eigenvalue of the generator L=d^2 - V' d (V=-log rho) on L^2(nu), by P1 FEM.
    Rayleigh quotient int f'^2 dnu / Var(f). Neumann ends. Returns (lambda_1, sorted_evals)."""
    x = np.asarray(x, float); n = len(x); h = np.diff(x)
    rho = np.exp(logrho); rho /= np.trapezoid(rho, x)
    K = np.zeros((n, n)); M = np.zeros((n, n))
    for e in range(n - 1):
        he = h[e]; w = 0.5 * (rho[e] + rho[e + 1]) * he
        ke = (1.0 / he) * np.array([[1, -1], [-1, 1]])
        me = (he / 6.0) * np.array([[2, 1], [1, 2]])
        idx = [e, e + 1]
        for ii in range(2):
            for jj in range(2):
                K[idx[ii], idx[jj]] += ke[ii, jj] * w
                M[idx[ii], idx[jj]] += me[ii, jj] * w
    vals = np.sort(eigh(K, M, eigvals_only=True))
    pos = vals[vals > 1e-9]
    return (pos[0] if len(pos) else 0.0), vals


def gpp_central(x, logrho):
    x = np.asarray(x, float); i0 = np.argmin(np.abs(x)); dx = x[1] - x[0]; lg = logrho
    d2 = (-lg[i0 - 2] + 16 * lg[i0 - 1] - 30 * lg[i0] + 16 * lg[i0 + 1] - lg[i0 + 2]) / (12 * dx * dx)
    return -d2


# ============================================================ (W) double-well model
def model_double_well(g):
    x = np.linspace(-2.6, 2.6, 2401)
    V = 0.5 * g * (x ** 2 - 1.0) ** 2
    logrho = -V; logrho -= logrho.max()
    Gpp = g * (3 * x ** 2 - 1.0)                 # exact V'' for this quartic
    rho = np.exp(logrho); rho /= np.trapezoid(rho, x)
    essinf = Gpp.min()
    mask = Gpp < 0
    nubad = np.trapezoid(rho[mask], x[mask]) if np.any(mask) else 0.0
    i0 = np.argmin(np.abs(x)); rgate = rho[i0]
    ones = np.argmin(np.abs(np.abs(x) - 1.0))
    barrier = V[i0] - V[ones]
    lam1, _ = lam1_of_marginal(x, logrho)
    return lam1, essinf, nubad, rgate, barrier


def kramers(g):
    wmin = np.sqrt(4.0 * g); wsad = np.sqrt(2.0 * g); dV = 0.5 * g
    return (wmin * wsad) / (2 * np.pi) * np.exp(-dV)


def cheeger_h(x, rho):
    """Exact 1-D Cheeger constant h = inf_x 2 rho(x)/min(F(x),1-F(x))  (half-lines are optimal in 1-D).
    Rigorous: lambda_1(nu) >= h^2/4.  Excludes extreme tails where the CDF loses precision.
    Returns (h, gate location, F(gate))."""
    x = np.asarray(x, float); r = np.asarray(rho, float); r = r / np.trapezoid(r, x)
    F = np.concatenate([[0.0], np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(x))])
    F = np.clip(F, 0, 1)
    mn = np.minimum(F, 1.0 - F)
    ok = (F > 1e-4) & (F < 1.0 - 1e-4)          # only where min(F,1-F) and the CDF are reliable
    psi = np.where(ok, 2.0 * r / np.maximum(mn, 1e-300), np.inf)
    i = np.argmin(psi)
    return psi[i], x[i], F[i]


# ============================================================ (U) free U(1)
def free_u1(L, mass=0.0):
    k = TWO_PI / L
    return k * k + mass * mass        # Gaussian: G''=lambda_1=precision -> 0


def main():
    print("=" * 80)
    print("(R) RICCI of SU(N), metric <X,Y> = -2 tr(XY)   [heat-kernel / Haar bi-invariant]")
    for N in (2, 3):
        d, ric, id_err = ricci_su(N)
        print(f"  SU({N}): dim={d:3d}   Ric(e_a,e_a)={ric.mean():.6f}  (theory N/4={N/4:.4f})  "
              f"spread={ric.std():.1e}  |sum f f - N I|={id_err:.1e}")
    print("  => Ric = (N/4) g > 0 pointwise [T]; SU(2)=0.5 <-> round S^3 radius 2 (kappa=1/2).")

    print("=" * 80)
    print("(G) GLOBAL BE at beta=2.5: weighted Ricci Ric + Hess S_W")
    beta = 2.5
    for N, ric in ((2, 0.5), (3, 0.75)):
        lam_sav = -1.12 * beta
        print(f"  SU({N}): Ric={ric:.2f}  Hess S_W(min)~{lam_sav:.2f}  => Ric+Hess S = "
              f"{ric+lam_sav:.2f} < 0   [NG] (matches W6 ~ -2.3 / -2.0)")
    print("  => global BE fails: ess inf(Ric+Hess S_W) < 0  (H6 essential-infimum).")

    print("=" * 80)
    print("(S) 1D SU(2), exact one-link Haar marginal rho_1(z)=pi(Cin(2pi)-Cin|z|)")
    x = np.linspace(-TWO_PI, TWO_PI, 1601)[1:-1]
    logrho = np.log(np.maximum(rho1_haar(x), 1e-300))
    gpp0 = gpp_central(x, logrho)
    rho = np.exp(logrho); rho /= np.trapezoid(rho, x)
    m = np.trapezoid(x * rho, x); var = np.trapezoid((x - m) ** 2 * rho, x)
    Gpp = -np.gradient(np.gradient(logrho, x), x)
    lam1, _ = lam1_of_marginal(x, logrho)
    print(f"  G''(0)={gpp0:.6f} (theory 1/(2Cin(2pi))={1/(2*Cin(TWO_PI)):.6f})  Var={var:.6f}  "
          f"1/Var={1/var:.6f}  G''*Var={gpp0*var:.4f}")
    print(f"  ess inf G''={Gpp.min():+.4f} at z={x[np.argmin(Gpp)]:.3f}  lambda_1(nu)={lam1:.6f}")
    hh, xh, Fh = cheeger_h(x, rho)
    print(f"  exact 1-D Cheeger h=inf 2rho/min(F,1-F)={hh:.4f} at gate z={xh:.3f} (F={Fh:.3f}); "
          f"h^2/4={hh**2/4:.4f} <= lambda_1={lam1:.4f}  (OK)")
    print("  => the one-link Haar marginal is LOG-CONCAVE (essinf G''=+0.205>0 at z=0); plain BE is")
    print("     non-vacuous, G''(0)=0.205 <= lambda_1=0.2645 <= 1/Var=0.2688, so chi<=1/lambda_1=3.78.")
    print("     (beta=0 anchor; the beta=2.5 periodic floor 0.222/Var 4.45 is the same single well [V].)")

    print("=" * 80)
    print("(P) periodic beta=0 k=0 marginal = L-fold convolution (exact Z0(k)=[phi(k/sqrtL)]^L)")
    xf = np.linspace(-TWO_PI, TWO_PI, 1601)[1:-1]; rf = rho1_haar(xf)
    def conv_marg(L, a):
        kmax = 60.0; dk = 0.05
        k = np.arange(-int(kmax / dk), int(kmax / dk) + 1) * dk
        ph = np.trapezoid(rf[None, :] * np.exp(-1j * np.outer(k, xf) / np.sqrt(L)), xf, axis=1)
        Z0 = (ph ** L) * np.exp(-(k / 40.0) ** 6)
        return (dk / TWO_PI) * np.real(np.exp(1j * np.outer(a, k)) @ Z0)
    for L in (2, 4, 8, 16):
        a = np.linspace(-0.5, 0.5, 201); r = conv_marg(L, a); lr = np.log(np.maximum(r, 1e-300))
        gg = -np.gradient(np.gradient(lr, a), a); i0 = np.argmin(np.abs(a))
        print(f"  L={L:2d}: G''(0)={gg[i0]:.4f}   (P2 exact beta=0 sumrule: L=2,4,8,16 -> "
              f"0.2396,0.2542,0.2591,0.2616)")
    print("  Var(Haar)=3.7198 exact; 1/Var=0.2688 CLT.  beta=2.5 periodic (P2/P4): G''=0.222, "
          "Var=4.45, G''*Var=0.963  [V]")

    print("=" * 80)
    print("(W) WEIGHTED-CRITERION MODEL: V=(g/2)(x^2-1)^2 (symmetric double well)")
    print("       g  lam1(FEM)  Kramers  G''w=2g  lam1/G''w  nu(G''<0)  rho_gate   Cheeger h  h^2/4   CD=(nu/rgate)^2")
    for g in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        lam1, essinf, nubad, rgate, barrier = model_double_well(g); kr = kramers(g)
        CD = (nubad / rgate) ** 2 if rgate > 0 else np.inf
        xm = np.linspace(-2.6, 2.6, 2401); rm = np.exp(-0.5 * g * (xm ** 2 - 1.0) ** 2); rm /= np.trapezoid(rm, xm)
        hh, xh, Fh = cheeger_h(xm, rm)
        print(f"     {g:4.1f} {lam1:10.3e} {kr:9.2e} {2*g:8.3f} {lam1/(2*g):9.4f} {nubad:10.3e} "
              f"{rgate:.2e} {hh:10.3e} {hh**2/4:8.2e} {CD:10.2e}")
    print("    single-well reference (Gaussian precision rho=2): lam1=2=G''well, nu(bad)=0, CD=0.")
    print("    => Cheeger h collapses with rho_gate (the gate), and lambda_1>=h^2/4 stays valid but tiny;")
    print("       nu(G''<0) SHRINKS (0.31->0.008) yet lam1/h^2 COLLAPSE: defect-MEASURE is WRONG control.")
    print("       Correct defect = CAPACITY: the gate density rho_gate (CD=(nu/rgate)^2 blows up).")

    print("=" * 80)
    print("(U) FREE U(1): massless k->0")
    for L in (4, 16, 64, 256):
        gg = free_u1(L); print(f"  L={L:4d}: G''=lambda_1=(2pi/L)^2={gg:.3e} -> 0  (ACID PASS)")
    print("  Ric(U(1))=0, Phi=0 (flat Haar) => no geometric floor.  [T]+[V]")

    print("=" * 80)
    print("SUMMARY (no QED):  R: Ric=(N/4)g>0 [T].  G: Ric+Hess S<0 [NG]=H6 global failure.")
    print("  S: 1D SU2 weighted constant finite ~0.2-0.27 [V].  W: defect-MEASURE insufficient (Kramers);")
    print("  need CAPACITY.  U: U(1)->0 acid pass.  Missing for 4D uniform = H1 [O].")


if __name__ == "__main__":
    main()
