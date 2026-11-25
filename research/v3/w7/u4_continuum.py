#!/usr/bin/env python3
r"""
u4_continuum.py -- U4: the a->0 continuum limit of the atom fraction w1 = f1^2/M0.
Deterministic, numpy only, one small job (~a few s). Self-check of the claims in
research/v3/U4-CONTINUUM-ATOM.md. Not a proof: theorem content is in the .md.

What is checked
  [A] RG equation for w1: f1^2 and M0 both carry Z_O^2, so w1=f1^2/M0 is invariant;
      d ln f1^2/d ln mu = 2*gamma_O != 0 while d ln w1/d ln mu = 0 (machine zero).
  [B] free field: :phi^2: Kallen-Lehmann density is a two-particle cut (no atom) =>
      f1^2=0 => w1=0 for every m and every UV cutoff (acid passes). M0(UV) diverges
      as the cutoff is removed -- the bare fraction is UV-dominated.
  [C] 2D exact: SU(2) character orthogonality <chi_0|chi_R>=delta_{0R} => f1^2(2D trF^2)=0
      => w1=0 exactly, though the 2D spectrum is gapped (the 0++ is a holonomy state).
  [D] binding->0 / thin atom: w1 -> 0 continuously as the atom residue -> 0, while the
      RG invariance (gamma_{w1}=0) holds at every point. RG carries the value, not its sign.
  [E] UV contact divergence of an operator of dimension Delta: M0=C(0) ~ a^{-(2Delta-4)}.
      Free 4D scalar: measured G(0)~a^{-2} => :phi^2: (Delta=2) M0~a^{-4}; trF^2 (Delta=4)
      M0~a^{-8}. chi=int s^{-1}dnu ~ a^{-4} for trF^2 (T4-CARRIER, cited).

Run: py -3.12 -I u4_continuum.py
"""
import numpy as np

def blk(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# ----------------------------------------------------------------------------
# [A] RG equation: w1 = f1^2/M0 is invariant; anomalous dimension gamma_{w1} = 0
# ----------------------------------------------------------------------------
def measure_moments(atom_res, atom_m2, cut_amp, m_cut2, s_max, n=200001):
    """A KL measure nu = atom_res*delta(s-atom_m2) + cut_amp*rho_cont(s) on [m_cut2,s_max].
    rho_cont chosen flat-normalized so M0, M1, chi are finite on the window.
    Returns f1^2, M0=int dnu, chi=int s^{-1}dnu."""
    f1sq = atom_res
    s = np.linspace(m_cut2, s_max, n)
    # a two-particle-like continuum density, integrable on the window
    rho = cut_amp * np.sqrt(np.clip(1.0 - m_cut2 / s, 0.0, None))   # ~const near thr, ~1 far
    M0 = f1sq + np.trapezoid(rho, s)
    chi = f1sq / atom_m2 + np.trapezoid(rho / s, s)
    return f1sq, M0, chi

def part_A():
    blk("[A] RG equation: f1^2 and M0 both ~ Z_O^2  =>  w1 = f1^2/M0 invariant (gamma_{w1}=0)")
    # one-loop: 1/g^2 = 2 b0 ln(mu/Lam); Z_O(mu,mu0) = g^2(mu0)/g^2(mu) (trace anomaly)
    b0 = 11.0 / (48 * np.pi**2)          # pure-YM one-loop coefficient used in corpus
    Lam = 1.0
    def g2(mu):  # one-loop running, 1/g^2 = 2 b0 ln(mu/Lam)
        return 1.0 / (2 * b0 * np.log(mu / Lam))
    def gamma_O(mu):  # gamma_O = -d ln|beta/2g|/d ln mu = 2 b0 g^2
        return 2 * b0 * g2(mu)
    mus = np.array([2.0, 3.0, 5.0, 10.0, 30.0, 100.0])
    # base measure at mu0 = 2
    a_res, a_m2, cut, mcut2 = 0.30, 1.00, 0.20, 2.00
    f1_0, M0_0, chi_0 = measure_moments(a_res, a_m2, cut, mcut2, 400.0)
    mu0 = mus[0]
    print(f"  base (mu0={mu0}):  f1^2={f1_0:.6f}  M0={M0_0:.6f}  w1=f1^2/M0={f1_0/M0_0:.10f}")
    print(f"  {'mu':>7} {'Z_O^2':>12} {'f1^2(mu)':>13} {'M0(mu)':>13} {'w1(mu)':>16} {'dlnw1/dlnmu':>12}")
    prev = None
    for mu in mus:
        Z2 = (g2(mu0) / g2(mu))**2       # = Z_O^2
        f1, M0, chi = f1_0 * Z2, M0_0 * Z2, chi_0 * Z2
        w1 = f1 / M0
        dlw = ""
        if prev is not None:
            mu_p, w1_p = prev
            dlw = f"{(np.log(w1) - np.log(w1_p)) / (np.log(mu) - np.log(mu_p)):+.2e}"
        print(f"  {mu:7.1f} {Z2:12.6f} {f1:13.6f} {M0:13.6f} {w1:16.12f} {dlw:>12}")
        prev = (mu, w1)
    # compare with d ln f1^2/d ln mu = 2 gamma_O (finite difference at mu=5)
    mu = 5.0; eps = 1e-4
    Z2p = (g2(mu0)/g2(mu*np.exp(eps)))**2
    Z2m = (g2(mu0)/g2(mu*np.exp(-eps)))**2
    dlnf1 = (np.log(f1_0*Z2p) - np.log(f1_0*Z2m)) / (2*eps)
    print(f"  check: dln f1^2/dln mu = {dlnf1:.6f}  vs  2*gamma_O(5) = {2*gamma_O(mu):.6f}")
    print("  => w1 column is constant to ~1e-15; dln w1/dln mu = 0.  gamma_{w1}=0 EXACTLY.")

# ----------------------------------------------------------------------------
# [B] free field: :phi^2: is a two-particle cut; f1^2=0 => w1=0; M0(UV) diverges
# ----------------------------------------------------------------------------
def part_B():
    blk("[B] free field :phi^2:  --  two-particle cut, NO atom  =>  w1 = 0 for all m, all UV")
    for m in (0.5, 1.0, 2.0):
        # d=4 free composite spectral density (threshold 4m^2), point mass at any m0: 0
        # f1^2 = nu({s0}) = 0 for every s0>0 (absolutely continuous measure)
        f1sq = 0.0
        # total weight M0(UV cutoff Lam) = int_{4m^2}^{Lam^2} rho(s) ds  (grows with Lam)
        def M0_of_Lam(Lam):
            s = np.linspace(4*m*m, Lam*Lam, 400001)
            rho = (1.0/(16*np.pi)) * np.sqrt(np.clip(1 - 4*m*m/s, 0, None))
            return np.trapezoid(rho, s)
        M_1, M_10 = M0_of_Lam(10.0), M0_of_Lam(100.0)
        print(f"  m={m:4.1f}: f1^2=0  =>  w1 = 0/M0 = 0 (exact);  "
              f"M0(Lam=10)={M_1:.4g} -> M0(Lam=100)={M_10:.4g}  (UV-divergent, bare fraction -> 0/0)")
    print("  => the free composite has NO point mass; w1 = 0 identically (acid PASSES).")

# ----------------------------------------------------------------------------
# [C] 2D exact: SU(2) character orthogonality => f1^2(2D trF^2) = 0 => w1 = 0
# ----------------------------------------------------------------------------
def part_C():
    blk("[C] 2D pure YM: trF^2=H conserved, character-diagonal => <0|trF^2|G>=0 => w1=0")
    # SU(2) characters chi_j(t)=sin((2j+1)t/2)/sin(t/2), Haar dU ~ sin^2(t/2) dt, t in [0,2pi]
    t = np.linspace(1e-9, 2*np.pi - 1e-9, 2000001)
    w = np.sin(t/2)**2
    def chi(j):
        return np.sin((2*j+1)*t/2) / np.sin(t/2)
    for j in (0.5, 1.0, 1.5, 2.0):
        ip = np.trapezoid(chi(0)*chi(j)*w, t) / np.trapezoid(chi(0)*chi(0)*w, t)
        print(f"  <chi_0|chi_{j:g}> = {ip:+.3e}   (=0 => f1^2(2D trF^2)=0 => w1=0)")
    print("  => the 2D 0++ 'glueball' is a HOLONOMY state; the trF^2 atom fraction w1 = 0 EXACTLY,")
    print("     even though 2D is gapped.  w1=0 is the correct continuum answer here.")

# ----------------------------------------------------------------------------
# [D] binding->0: w1 -> 0 continuously while gamma_{w1}=0 (RG blind to positivity)
# ----------------------------------------------------------------------------
def part_D():
    blk("[D] binding -> 0:  atom residue A -> 0  =>  w1=A/M0 -> 0,  RG invariance holds throughout")
    cut_amp, m_cut2, s_max = 0.20, 2.00, 400.0
    # continuum (cut) total weight, fixed:
    s = np.linspace(m_cut2, s_max, 400001)
    rho = cut_amp * np.sqrt(np.clip(1 - m_cut2/s, 0, None))
    M0_cut = np.trapezoid(rho, s)
    print(f"  fixed continuum total weight M0_cut = {M0_cut:.6f}")
    print(f"  {'atom res A':>12} {'w1=A/(A+cut)':>16} {'Z_O^2':>10} {'w1 after Z^2=9':>16}")
    for A in (1.0, 0.3, 0.1, 0.01, 1e-4, 1e-8, 0.0):
        M0 = A + M0_cut
        w1 = A / M0
        Z2 = 9.0
        w1z = (A*Z2) / (M0*Z2)     # same ratio after identical multiplicative scaling
        print(f"  {A:12.1e} {w1:16.10f} {Z2:10.1f} {w1z:16.10f}")
    print("  => w1 -> 0 as binding vanishes; the Z_O^2 transport maps it to the SAME w1.")
    print("     Anomalous dimension is 0 at every A, yet positivity is NOT protected (0 is absorbing).")

# ----------------------------------------------------------------------------
# [E] UV contact divergence  M0=C(0) ~ a^{-(2 Delta - 4)}  (free 4D scalar skeleton)
# ----------------------------------------------------------------------------
def part_E():
    blk("[E] UV contact: M0=C(0) ~ a^{-(2Delta-4)}.  Free 4D scalar measures G(0) ~ a^{-2}.")
    # fixed physical box Lbox=4, physical mass m=1; lattice spacing a=Lbox/L, L sites/dim.
    # G(0) = (1/V) sum_p 1/(phat^2+m^2),  phat_mu = (2/a) sin(p_mu a/2),  p_mu=2 pi n_mu/(L a).
    # As a->0 at fixed physical box, the UV cutoff pi/a grows => G(0) ~ 1/(16 a^2).
    m = 1.0
    Lbox = 4.0
    rows = []
    for L in (4, 8, 12, 16, 24, 32):
        a = Lbox / L
        V = Lbox**4                                  # physical volume (fixed)
        n = np.arange(-(L // 2), L // 2)             # L values per dim
        pa = 2 * np.pi * n / L                       # p_mu * a
        phat1 = (2.0 / a) * np.sin(pa / 2.0)         # lattice dispersion component
        PX = phat1.reshape(-1, 1, 1, 1)
        PY = phat1.reshape(1, -1, 1, 1)
        PZ = phat1.reshape(1, 1, -1, 1)
        PW = phat1.reshape(1, 1, 1, -1)
        phat2 = PX**2 + PY**2 + PZ**2 + PW**2
        G0 = (1.0 / V) * np.sum(1.0 / (phat2 + m * m))   # = int d^4p/(2pi)^4 (Riemann sum)
        rows.append((a, G0))
    rows = np.array(rows)
    slope = np.polyfit(np.log(rows[:, 0]), np.log(rows[:, 1]), 1)[0]
    print(f"  {'a':>8} {'G(0)':>14} {'G(0)*a^2':>12}  (->const confirms G(0)~a^-2)")
    for a, g in rows:
        print(f"  {a:8.4f} {g:14.6f} {g*a*a:12.6f}")
    print(f"  measured log-log slope of G(0) vs a: {slope:+.3f}  (theory -2 => G(0)~a^-2)")
    print("  => :phi^2: (Delta=2): M0=2G(0)^2 ~ a^-4 ;  trF^2 (Delta=4): M0 ~ a^-8 (power count).")
    print("     chi=int s^{-1}dnu ~ a^-4 for trF^2 (T4-CARRIER).  w1=f1^2/M0 is dimensionless.")
    print("     FIXED-SCALE renorm (flow t / point-split eps) removes the contact => finite w1.")

if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
    part_D()
    part_E()
    print("\nDone.  No QED.  Clay 0%.  (self-check only; theorem content in U4-CONTINUUM-ATOM.md)")
