"""s05_acid.py — acid test of (star):  G''_j = E[H_j] - Var_j + W_j >= delta_W > 0
Four models: free field / 2D YM / M6 / single plaquette.
Deterministic, no sampling.  py -3.12, numpy only.
Definitions (corpus, verbatim):
  G'' = -d^2/da^2 log rho_v(a)|_0 = Delta_dyn + W_full ,  W_full = W_geo + X
  W_geo = E[Phi''] - Var(Phi'),  Phi = -log J            (T03-TOY:24 ; N04-2D:26-30)
"""
import numpy as np

def chi_pp(phi, w):
    """d^2/dphi^2 [ sin(w phi)/sin(phi) ]  (w = 2j+1 odd integer)."""
    s, c = np.sin(phi), np.cos(phi)
    sw, cw = np.sin(w * phi), np.cos(w * phi)
    return (s**3 * sw * (1.0 - w * w) - 2.0 * s * c * (w * cw * s - sw * c)) / s**4

def Kt(phi, t, jmax=60):
    """Migdal kernel K_t(phi) with analytic 1st/2nd phi-derivatives (2D SU(2) YM)."""
    K = Kp = Kpp = 0.0
    s, c = np.sin(phi), np.cos(phi)
    for j in range(jmax + 1):
        w = 2 * j + 1
        wt = (2 * j + 1) * np.exp(-t * j * (j + 1))
        sw, cw = np.sin(w * phi), np.cos(w * phi)
        K += wt * sw / s
        # d/dphi [ sin(w phi)/sin phi ] = (w cos(w phi) sin phi - sin(w phi) cos phi)/sin^2 phi
        Kp += wt * (w * cw * s - sw * c) / s**2
        Kpp += (2 * j + 1) * np.exp(-t * j * (j + 1)) * chi_pp(phi, w)
    return K, Kp, Kpp

def two_d(phi, t):
    K, Kp, Kpp = Kt(phi, t)
    W_geo = 2.0 / np.sin(phi)**2          # -d^2 log sin^2(phi)
    D_dyn = -(Kpp / K) + (Kp / K)**2      # -d^2 log K_t
    return W_geo + D_dyn, W_geo, D_dyn

print("== [1] FREE FIELD : G'' = Delta_dyn + W_full ==")
print("  periodic continuum / flat chart: Phi = 0, S quadratic, Stein  E[H] = Var(G)")
print("     => Delta_dyn = 0 ,  W_full = 0  =>  G'' = 0   (exact)", )
print("  finite Dirichlet box (the '0.6169' column): G'' = pi^2/L^2 -> 0")
for L in (2, 4, 8, 16, 64):
    print(f"     L={L:3d}  pi^2/L^2 = {np.pi**2/L**2:.6f}   4pi^2/L^2 = {4*np.pi**2/L**2:.6f}")
print("  => free field is FALSE for any fixed delta_W>0 (delta_W may -> 0 with L, or G''=0)")

print("\n== [2] 2D YM (Migdal exact, space circle) : G'' = W_geo + Delta_dyn ==")
print("     phi     W_geo    D_dyn(t=.5)   G''(t=.5)   D_dyn(t=5)   G''(t=5)   D_dyn(t=40)  G''(t=40)")
for phi in (0.3, 0.6, 1.0, 1.5, 2.0, 2.6):
    g05, W, d05 = two_d(phi, 0.5)
    g5, _, d5 = two_d(phi, 5.0)
    g40, _, d40 = two_d(phi, 40.0)
    print(f"   {phi:5.2f} {W:9.5f} {d05:12.6f} {g05:11.6f} {d5:12.3e} {g5:10.6f} {d40:11.2e} {g40:10.6f}")
phi_grid = np.linspace(0.05, np.pi - 0.05, 401)
for t in (0.5, 1.0, 2.0, 5.0, 10.0, 40.0):
    vals = np.array([two_d(p, t)[0] for p in phi_grid])
    Ws = np.array([2.0 / np.sin(p)**2 for p in phi_grid])
    print(f"   t={t:5.1f}  min_phi G'' = {vals.min():+9.5f}   min_phi W_geo = {Ws.min():6.4f}"
          f"   #(G''<2) = {int((vals < 2).sum())}/401")
print("  identity residual  G'' - (W_geo + Delta_dyn) at phi=1.0,t=5 :",
      abs(two_d(1.0, 5.0)[0] - (two_d(1.0, 5.0)[1] + two_d(1.0, 5.0)[2])))
print("  => 2D: G'' = 2/sin^2(phi) + O(e^{-2t}) >= 2  (kinematic, t-independent floor) TRUE")

print("\n== [3] M6  (oscillator (x) free flat direction) : G'' = 2*omega on a-sector ==")
for om in (0.5, 1.0, 2.0):
    print(f"   omega={om:4.2f}   G''_a = 2*omega = {2*om:6.3f}   ;  m_L = min(omega, 2pi^2/L^2) -> 0")
print("   b-sector: Lebesgue (S=0, Phi=0, J=1) => Delta_dyn=0, W=0 => G''_b = 0  EXACT")
print("   => pointwise-in-a-mode reading: (star) TRUE on M6 ; inf-over-sectors reading: 0 -> FALSE")

print("\n== [4] SINGLE PLAQUETTE (0-dimensional, no transfer direction) ==")
# exact anchor carried from T03-TOY:50 / q11 (not recomputed here):  Z = Cin(2pi)/2
print("   G'' = S - D = +0.135994711041988  [E] (T03-TOY:50, q11 anchor, |diff| 1.1e-16)")
print("   m : UNDEFINED (0 dims, no transfer matrix, no E1).  no scale ladder j => 'uniform in j' empty")

print("\n== [5] SIGNATURE of (star) vs the four-model touchstone (False/True/False/undef) ==")
print("   reading (P) pointwise in the coarse amplitude, a-sector only:")
print("       free 0 (F) | 2D >=2 (T) | M6 2omega>0 (T) | plaquette 0.136 but m undef -> (F,T,T,undef)  => FAIL")
print("   reading (Q) inf over all coarse sectors:")
print("       free 0 (F) | 2D >=2 (T) | M6 min(2omega,0)=0 (F) | plaquette undef -> (F,T,F,undef)  => PASS")
print("   BUT (Q) is carried on 2D by W_geo = 2/sin^2(phi): pure Haar/Jacobian, t-independent, kinematics.")
print("   => (Q) is a 'is the coarse coordinate non-flat?' detector, not a gap detector.")
print("   witness: at t=40 D_dyn ~ 1e-34 (dynamics switched off), G'' >= 2 unchanged.")
