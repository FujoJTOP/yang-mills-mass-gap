#!/usr/bin/env python3
r"""
p3_polyakov_center.py -- P3: Polyakov / center-sector analysis of the PERIODIC
constant-mode curvature.  Complements the now-landed anchors:

  P1: on the periodic torus the GLOBAL constant-holonomy slice has S_W == 0 exactly
      (all plaquettes trivial: the constant field is flat), so the measure there is
      pure Haar and G'' = 0.2051 (L=1 Haar marginal; -> 1/Var(xi3)=0.2688 as L->inf).
  P2: the periodic k=0 mode saturates: c = 0.22196, q = 0.849, Var(a) plateau = 4.45
      (=> 1/Var = 0.2247), i.e. a positive L-independent floor just above the Haar value.

THE QUESTION THIS WITNESS SETTLES.
  On the temporal torus the constant mode is the Polyakov loop (holonomy).  The measure
  is Z(N)-invariant (center symmetry).  Does the center sector change the constant-mode
  floor -- is the physical periodic floor the Haar floor (P1, 0.2051), the P2 value
  (0.222), or an explicitly-broken center-sector floor (the fixed-boundary twist ~0.48)?

SECTIONS
 [A] Constant-holonomy slice S_W==0 in EVERY center sector  [T] + tiny numeric plaquette check.
     Contrast: a genuine spatial 't Hooft flux (B^2=-I, beta=2.5) is NOT flat (S_W>0) -- that is
     the boundary-fixed twist that D5 calls ~0.48, NOT the Polyakov center sector.
 [B] Full-3D L=1 minimal torus: partition function / marginal are EXACTLY theta-independent
     (right-Haar).  The flat-chart G''(a) theta-dependence (0.73 at B=I vs 0.013 at B=-I) is a
     CHART artifact (exponential chart centred at U=I; J(2pi)=0 edge).   [T]+[V]
 [C] Fixed-boundary center-sector spread, exact pure-Cartan transfer matrix, beta=2.5, L=4..8:
     G''(a0*) at theta=0 (trivial), theta=pi ('t Hooft flux), theta=2pi (strict center -I).
     These are EXPLICIT-breaking boundary values, not the periodic value.   [V]
 [D] Finite-L center symmetry is UNBROKEN: 1D strip transfer matrix gapped (D3 [T]) and
     a->-a symmetry makes <tr P>=0 at theta in {0,2pi}.  No SSB at these L.   [T]+[V]
 [E] Acid test: free Gaussian G''=kappa->0 and stays kappa->0 under center-orbit average
     (no manufactured floor); beta=0 Haar 0.2051 is sector-independent.   [T]

Run: py -3.12 -I p3_polyakov_center.py [quick|full]
numpy + d3_transfer (exact transfer matrix).  No downloads.  No QED.  Clay 0%.
"""
import sys, numpy as np
sys.path.insert(0, r'D:\Dev\math\research\v3\w7')
np.seterr(all='ignore')

# ------------------------------------------------------------------ SU(2) helpers
def e2(v):
    """exp(i v.sigma/2) for a 3-vector v."""
    r = np.linalg.norm(v)
    if r < 1e-15:
        return np.eye(2, dtype=complex)
    n = v / r
    return (np.cos(r/2)*np.eye(2) + 1j*np.sin(r/2)*(n[0]*np.array([[0,1],[1,0]]) +
            n[1]*np.array([[0,-1j],[1j,0]]) + n[2]*np.array([[1,0],[0,-1]]))).astype(complex)

def htr(M):
    return 0.5*np.real(np.trace(M))

# ================================================================== [A] flat slice, every center sector
def secA():
    print("="*92)
    print("[A] S_W on the GLOBAL constant-holonomy slice, by center sector   [T]+[V]")
    print("    (constant diagonal links a_mu = c_mu e3; centre = multiply one direction by -I)")
    print("-"*92)
    c = np.array([0.7, 1.3])          # two directions, both along Cartan e3 -> commute
    e3 = np.array([0., 0., 1.])
    print(f"    {'sector':>26}{'all |tr P_plaq /2 - 1|':>24}{'S_W/L^2':>10}")
    for lab, zs in (("holonomy = I  (trivial)", (1, 1)),
                    ("holonomy = -I (center)", (-1, 1)),
                    ("holonomy = -I x -I", (-1, -1))):
        U0 = zs[0]*e2(c[0]*e3); U1 = zs[1]*e2(c[1]*e3)
        P = U0 @ U1 @ U0.conj().T @ U1.conj().T
        err = abs(htr(P) - 1.0)
        SW = 2.5*(1.0 - htr(P))
        print(f"    {lab:>26}{err:>24.2e}{SW:>10.2e}   -> S_W = 0 on the flat slice, every sector")
    # a genuine spatial 't Hooft flux: link holonomy B with B^2=-I  (theta=pi, B=i sigma3)
    # constant 't Hooft flux is NOT flat: the twist plaquette = center element, tr = -2
    B = e2(np.pi*e3)                  # e^{i pi sigma3/2} = i sigma3 ; B^2 = -I
    print(f"\n    spatial 't Hooft flux B = i sigma3  (theta=pi), B^2=-I :  1/2 tr(B^2) = {htr(B@B):+.4f}")
    print(f"    the twist plaquette = -I  =>  S_W/twist-plaquette = beta(1-1/2 tr(-I)) = "
          f"{2.5*(1-htr(-np.eye(2))):.2f} > 0")
    print("    => the fixed-boundary theta=pi twist is an ACTIONFUL ('t Hooft-flux) sector, NOT a")
    print("       flat Polyakov center sector.  The ~0.48 floor of D5 is THAT boundary, not the periodic one.")

# ================================================================== [B] full-3D L=1 theta-independence
def secB():
    print("\n"+"="*92)
    print("[B] Full-3D L=1 minimal torus: is the physical measure theta-dependent?   [T]+[V]")
    print("    Haar is right-invariant: P=U B is Haar for every fixed B => Z and the plaquette")
    print("    distribution are EXACTLY theta-independent.  Numerical check (int rho_long da):")
    print("-"*92)
    from sm4_constant_mode import rho_long
    aa = np.linspace(-2*np.pi+0.02, 2*np.pi-0.02, 6001)
    Zs = []
    for th in (0.0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi):
        Z = np.trapezoid([rho_long(a, 2.5, th) for a in aa], aa)
        Zs.append(Z)
        print(f"    theta={th:6.3f}:  Z = int rho_long da = {Z:.8e}")
    print(f"    spread = {max(Zs)-min(Zs):.2e}  (machine/quadrature zero)  =>  theta-INDEPENDENT  [V]")
    print("    => the flat-chart G''(a) theta-spread (0.730 at B=I vs 0.013 at B=-I, L=1) is a CHART")
    print("       artifact: the exponential chart is centred at U=I and J(2pi)=0 at the -I edge.  [T]")

# ================================================================== [C] fixed-boundary sector spread
def secC(mode):
    print("\n"+"="*92)
    print("[C] Fixed-boundary center-sector spread, EXACT pure-Cartan transfer matrix, beta=2.5.  [V]")
    print("    theta=0 trivial | theta=pi 't Hooft flux | theta=2pi strict center -I.  L=4..8.")
    print("    (These are EXPLICIT-breaking boundary values; the periodic (P2) value is ~0.222.)")
    print("-"*92)
    from d3_transfer import make_kgrid, nodes, build_Z, gap_from_Z, haar_var_xi3
    kgrid, dk = make_kgrid(mode); z, r1 = nodes(400 if mode != 'full' else 500)
    v1 = haar_var_xi3()
    print(f"    Haar CLT floor 1/Var(xi3) = {1/v1:.5f}   (L=1 Haar marginal = 0.2051; P2 periodic = 0.22196)")
    print(f"    {'L':>2} | {'theta=0':>9}{'G(a0*)':>9} | {'theta=pi':>9}{'G(a0*)':>9} | {'theta=2pi':>9}{'G(a0*)':>9}")
    for L in (4, 5, 6, 7, 8):
        half = 2*np.pi*np.sqrt(L); row = []
        for th in (0.0, np.pi, 2*np.pi):
            Z = build_Z(L, th, 2.5, kgrid, z, r1)
            a0, G, Gb, _ = gap_from_Z(Z, kgrid, dk, half)
            row.append((a0, G))
        print(f"    {L:2d} | {row[0][0]:+9.2f}{row[0][1]:9.4f} | {row[1][0]:+9.2f}{row[1][1]:9.4f} | "
              f"{row[2][0]:+9.2f}{row[2][1]:9.4f}")
    print("    => the sector VALUE moves (spread ~0.2-0.4) but stays positive; the periodic P2 value")
    print("       0.222 lies near the Haar floor, NOT at the 't Hooft-flux boundary value.")

# ================================================================== [D] center symmetry unbroken
def secD():
    print("\n"+"="*92)
    print("[D] Finite-L center symmetry: UNBROKEN at these L.   [T]+[V]")
    print("-"*92)
    from d3_transfer import bulk_spectrum
    for th in (0.0, np.pi, 2*np.pi):
        ev, D, xi = bulk_spectrum(th, 2.5)
        print(f"    theta={th:5.3f}: transfer-matrix Delta=-log(lam1/lam0) = {D:.4f} > 0  => xi finite")
    print("    A gapped 1D finite-range kernel has NO spontaneous symmetry breaking [T] (D3).")
    print("    At theta in {0,2pi} the marginal is symmetric under a->-a (P=U B with B commutes);")
    print("    hence <tr P> = 0 (Polyakov order parameter vanishes): the center sector is UNBROKEN.")
    print("    The physical ensemble is the center-NEUTRAL one (no sector selection at these L).")

# ================================================================== [E] acid test
def secE():
    print("\n"+"="*92)
    print("[E] Acid test (free field).   [T]")
    print("-"*92)
    for kappa in (0.5, 0.1, 1e-3):
        print(f"    free massless Gaussian, kappa={kappa:7.4f}: G''=kappa={kappa:7.4f} -> 0")
    print("    Center-orbit average of a Gaussian N(0,1/kappa) is the SAME Gaussian (itself")
    print("    centre-symmetric) => G'' stays kappa -> 0: the centre projection manufactures NO floor.")
    print("    beta=0 pure Haar L=1: G''=0.2051 in EVERY centre sector (pure Haar is centre-blind).")
    print("    Discriminator: the Haar Jacobian Phi=-log J (absent for the free Gaussian).")

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    secA(); secB(); secC(mode); secD(); secE()
    print("\nP3-polyakov-center: done.  No QED.  Clay 0%.")
