#!/usr/bin/env python3
r"""
p5_periodic.py -- P5: ADVERSARIAL universality test for the PERIODIC constant mode.

QUESTION.  For the 1xL chain with fixed boundary holonomy, D5 found (i) the Z2 twist theta=pi
SATURATES at a positive L-independent floor and (ii) the trivial fixed boundary b=0 makes Var(a)
GROW (G'' decays).  The PERIODIC chain (no boundary, links wrap around) is neither: the constant
configuration has S_W = 0 EXACTLY for every L (no end plaquettes to pin the uniform mode).  Is the
periodic constant-mode curvature POSITIVE and L-INDEPENDENT (Route B, saturation at 1/chi) or does
Var(a) GROW / G'' DECAY (Route A, equipartition / free zero mode)?

ROUTES.
  (A) measure concentration / equipartition: the uniform link mode is an EXACT zero mode of the
      periodic Wilson difference operator; a = (1/sqrt L) sum_l z_l = sqrt(L)*c for a uniform value c,
      so Var(a) = L*Var(c) -> grows -> G'' ~ 1/L -> 0.         [the free-field / zero-mode route]
  (B) Haar / Prekopa saturation: the uniform mode is pinned by the local Haar measure (compactness /
      Jacobian Phi=-log J), chi = sum_d C(d) < infty => Var(a) -> chi finite => G'' -> 1/chi > 0.
      At beta=0 this is EXACT: Var(a)=Var_Haar(xi3)=3.71981, G''->0.26883, L-independent.

MODEL.  Pure-Cartan-sector chain (exact diagonal sector at every L; the beta=0 measure part is EXACT
full-3D).  Per-link EXACT one-component Haar marginal rho1(z) (transverse components integrated out).
  open strip (D1/D5):  P1=X1X2, P_i=B X_i^-1 X_{i+1}, P_L=B X_L^-1, z_l in [-2pi,2pi], Haar weight rho1.
  periodic:            P_i = X_i^-1 X_{i+1},  z_{L+1}=z_1  (trace of the transfer matrix).
U(1) (flat chart, no Haar Jacobian) is the abelian/free DISCRIMINATOR (acid test).

METHOD (deterministic, NO Monte Carlo).
  Var(a) EXACT for the periodic chain: Z(k)=tr(M_k^L), M_k = D K_k D, D=diag(sqrt(rho1(z_a) w_a)),
  K_k[a,b]=exp(-beta(1-cos((z_b-z_a)/2))) exp(-i k (z_a+z_b)/(2 sqrt L));  then Var(a)=-Z''(0)/Z(0),
      Z''(0) = L tr(M^{L-1} M'') + L(L-1) tr(M^{L-2} M'^2),
      M' = (-i/(2 sqrt L)) D H D,  H[a,b]=(g_a+g_b) B[a,b],  g_a=z_a,
      M'' = (-1/(4 L)) D [(g_a+g_b)^2 B] D .   (Z'(0)=0 by z->-z symmetry.)
  G''(a0*) from rho(a)=(1/2pi) int Z(k) e^{ika} dk by a modest direct quadrature (L<=5).
  G''~1/Var once the marginal Gaussianises (D5 established R=G''*Var->1); 1/Var is reported for every L.

Run: py -3.12 -I p5_periodic.py [quick|full].  numpy/scipy only, one job at a time.  No downloads.
Labels [T]/[V]/[C]/[NG]/[O].  No QED.  Clay 0%.
"""
import sys, time
import numpy as np
from scipy.special import sici

np.seterr(invalid='ignore', divide='ignore')
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi
MODE = sys.argv[1] if len(sys.argv) > 1 else 'quick'
NZ = {'quick': 80, 'full': 160}[MODE]

def Cin(x):
    x = np.asarray(x, float); _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)

def rho1(z):
    z = np.asarray(z, float); az = np.abs(z)
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)

def var_rho1(n=400001):
    z = np.linspace(-TWO_PI, TWO_PI, n); w = np.full(n, 2 * TWO_PI / (n - 1)); w[[0, -1]] *= 0.5
    r1 = rho1(z) * w; Z = r1.sum(); m = (z * r1).sum() / Z
    return ((z - m) ** 2 * r1).sum() / Z

# ---------------- periodic trace: nodes, kernel, exact Var(a) -----------------------------------
def periodic_nodes(kind):
    xg, wg = np.polynomial.legendre.leggauss(NZ)
    if kind == 'su2':
        z = TWO_PI * xg; w = TWO_PI * wg; rho = rho1(z)
    else:
        z = np.pi * xg; w = np.pi * wg; rho = np.ones_like(z)
    return z, w, rho

def periodic_DB(kind, z, w, rho, beta):
    zp, zm = z[:, None], z[None, :]
    d = 0.5 * (zm - zp) if kind == 'su2' else (zm - zp)
    B = np.exp(-beta * (1.0 - np.cos(d)))
    D = np.diag(np.sqrt(rho * w))
    return B, D, z

def periodic_var_a(kind, L, beta):
    z, w, rho = periodic_nodes(kind)
    if L == 1:                                   # one plaquette X^-1 X = I -> pure Haar
        p = rho * w; Z0 = p.sum(); m = (z * p).sum() / Z0
        return float(((z - m) ** 2 * p).sum() / Z0)
    B, D, g = periodic_DB(kind, z, w, rho, beta)
    M = D @ B @ D
    Mp = (-1j / (2 * np.sqrt(L))) * (D @ ((g[:, None] + g[None, :]) * B) @ D)
    Mpp = (-1.0 / (4 * L)) * (D @ (((g[:, None] + g[None, :]) ** 2) * B) @ D)
    ML1 = np.linalg.matrix_power(M, L - 1)
    ML2 = np.linalg.matrix_power(M, L - 2)
    Z0 = np.trace(np.linalg.matrix_power(M, L))
    Zpp = L * np.trace(ML1 @ Mpp) + L * (L - 1) * np.trace(ML2 @ (Mp @ Mp))
    return float(np.real(-(Zpp / Z0)))

def periodic_Zfun(kind, L, beta):
    z, w, rho = periodic_nodes(kind)
    B, D, g = periodic_DB(kind, z, w, rho, beta)
    M0 = D @ B @ D
    # eigen-decompose only if real-symmetric path; else matrix power per k
    def Zfun(k):
        k = np.atleast_1d(np.asarray(k, float)); out = np.empty(k.shape[0], complex)
        for i, kk in enumerate(k):
            K = B * np.exp(-1j * kk * (g[:, None] + g[None, :]) / (2 * np.sqrt(L)))
            out[i] = np.trace(np.linalg.matrix_power(D @ K @ D, L))
        return out
    return Zfun

def gpp_periodic(kind, L, beta, K=40.0, dk=0.05, na=1601):
    """rho(a) by direct quad; G'' at the mode by 5-pt stencil."""
    z, w, rho = periodic_nodes(kind)
    B, D, g = periodic_DB(kind, z, w, rho, beta)
    kj = np.arange(-int(K / dk), int(K / dk) + 1) * dk
    Zv = np.empty(kj.shape[0], complex)
    for i, kk in enumerate(kj):
        Kk = B * np.exp(-1j * kk * (g[:, None] + g[None, :]) / (2 * np.sqrt(L)))
        Zv[i] = np.trace(np.linalg.matrix_power(D @ Kk @ D, L))
    a = np.linspace(-TWO_PI * np.sqrt(L), TWO_PI * np.sqrt(L), na)
    rho_a = (dk / TWO_PI) * np.real(np.exp(1j * np.outer(a, kj)) @ Zv)
    rho_a = np.maximum(rho_a, 1e-300)
    lg = np.log(rho_a); h = a[1] - a[0]; m = 4
    im = m + int(np.argmax(rho_a[m:-m]))
    G = -(-lg[im-2] + 16*lg[im-1] - 30*lg[im] + 16*lg[im+1] - lg[im+2]) / (12*h*h)
    return a[im], G

# ---------------- exact U(1) open Var (D5 product quadrature, small L) --------------------------
def u1_open_var_exact(L, beta, th, nx):
    XG, WG = np.polynomial.legendre.leggauss(nx)
    pts = np.pi * (XG + 1.0); wts = np.pi * WG
    grids = [g.ravel() for g in np.meshgrid(*[pts] * L, indexing='ij')]
    Wg = np.meshgrid(*[wts] * L, indexing='ij'); W = np.ones_like(Wg[0])
    for k in range(L): W = W * Wg[k]
    W = W.ravel()
    if L == 1:
        pl = [grids[0] + th]
    else:
        pl = [grids[0] + grids[1]]
        for j in range(2, L): pl.append(th - grids[j - 1] + grids[j])
        pl.append(th - grids[L - 1])
    S = np.zeros(grids[0].size)
    for p in pl: S = S + beta * (1 - np.cos(p))
    a = np.sum(grids, axis=0) / np.sqrt(L); w = W * np.exp(-S)
    m = np.sum(w * a) / np.sum(w)
    return np.sum(w * a * a) / np.sum(w) - m * m

def main():
    t0 = time.time()
    print("#" * 100)
    print("# P5 -- periodic constant mode: Route A (decay) vs Route B (saturation); adversarial test")
    print("#" * 100)
    vh = var_rho1()
    print(f"[ref] Var_Haar(xi3)={vh:.5f}   1/Var_Haar={1/vh:.5f}  (the beta=0 periodic floor)  [T]")
    Ls = [1, 2, 3, 4, 5, 6, 8, 10]

    print("\n" + "=" * 100)
    print("[U1] U(1) (flat chart, NO Haar Jacobian) -- the abelian/free DISCRIMINATOR ; Var(a)")
    for beta in (0.0, 1.0, 2.5):
        nx = 60 if True else 40
        row_o0 = [u1_open_var_exact(L, beta, 0.0, nx if L <= 3 else 32) if L <= 4 else np.nan for L in Ls]
        row_op = [u1_open_var_exact(L, beta, np.pi, nx if L <= 3 else 32) if L <= 4 else np.nan for L in Ls]
        row_pe = [periodic_var_a('u1', L, beta) for L in Ls]
        fmt = lambda r: "".join(("%9.3f" % x) if np.isfinite(x) else "      ---" for x in r)
        print(f"   beta={beta}:   L      ={''.join('%9d'%L for L in Ls)}")
        print(f"               open b=0  ={fmt(row_o0)}")
        print(f"               open b=pi ={fmt(row_op)}")
        print(f"               periodic  ={fmt(row_pe)}   <-- P5 TARGET")

    print("\n" + "=" * 100)
    print("[SU2] SU(2) pure-Cartan chain -- periodic Var(a) EXACT (analytic -Z''/Z)")
    for beta in (0.0, 1.0, 2.5):
        row_pe = [periodic_var_a('su2', L, beta) for L in Ls]
        fmt = lambda r: "".join("%9.4f" % x for x in r)
        print(f"   beta={beta}:   L      ={''.join('%9d'%L for L in Ls)}")
        print(f"               periodic  ={fmt(row_pe)}   <-- P5 TARGET  (plateau => Route B)")

    print("\n" + "=" * 100)
    print("[GPP] periodic SU(2) G''(a0*) by direct-quad rho(a), L=2..5, beta=0 and 2.5  [V]")
    for beta in (0.0, 2.5):
        print(f"   beta={beta}:")
        for L in (2, 3, 4, 5):
            try:
                a0, G = gpp_periodic('su2', L, beta)
                va = periodic_var_a('su2', L, beta)
                print(f"      L={L}: a0*={a0:+.3f}  G''(a0*)={G:+.4f}   Var(a)={va:.4f}   "
                      f"1/Var={1/va:.4f}   R=G''*Var={G*va:+.3f}")
            except Exception as ex:
                print(f"      L={L}: skipped ({ex})")

    print("\n" + "=" * 100)
    print("[ACID] discriminator:")
    print("   free/abelian: uniform mode = zero mode, Var(a)~L^2/grows -> G''->0.  U(1) periodic shows it.")
    print("   SU(2): Haar Jacobian Phi=-log J pins the uniform mode.  Periodic Var(a) plateau => Route B;")
    print("          growth => Route A.  The SU(2) table above decides.")
    print(f"\nP5 periodic done in {time.time()-t0:.0f}s.  No QED.  Clay 0%.")

if __name__ == '__main__':
    main()
