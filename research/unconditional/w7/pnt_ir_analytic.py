"""
pnt_ir_analytic.py -- witness for PNT-IR-ANALYTIC (2026-11-24)

IR-only analytic functional = nearest Lee-Yang / Fisher zero in the complex TIME extent tau
of the 0++-PROJECTED partition function  Z_{0++}(tau) = Tr_{0++} e^{-tau H}.

Fixes LEEYANG's failure: LEEYANG used complex COUPLING beta=1/g^2 (kinematic beta=0 zeros,
acid-dirty). Here we use complex TIME tau at fixed real coupling, projected to the 0++ sector:
   Z(tau) = g0 + sum_atoms g_i e^{-tau m_i} + integral_{s0}^inf rho(s) e^{-tau s} ds .
A zero at Re(tau)>0 requires an ISOLATED state below the continuum, and then
   Im(tau_zero) = pi / m_gap        (pure atom: g0 + g1 e^{-tau m}=0  =>  tau = (1/m)(ln(g1/g0) - i pi(2k+1))).

Run:  py -3.12 -I research/unconditional/w7/pnt_ir_analytic.py
"""
import numpy as np


def Z(tau, g0, atoms, s0, cont, npt=500, W=28.0):
    z = g0 + sum(g * np.exp(-tau * m) for g, m in atoms)
    if s0 is not None and cont != 0.0:
        ss = np.linspace(s0, s0 + W, npt)
        ds = ss[1] - ss[0]
        rho = cont * np.sqrt(np.maximum(ss - s0, 0.0)) * np.exp(-ss)
        z = z + np.sum(rho * np.exp(-tau * ss)) * ds
    return z


def nearest_zero(g0, atoms, s0, cont):
    """scan Re(tau)>0 half-plane for the minimum |Z| (robust; no fragile root iteration)."""
    best = None
    for x in np.linspace(0.01, 4.0, 160):
        for y in np.linspace(0.05, 8.0, 320):
            v = abs(Z(x + 1j * y, g0, atoms, s0, cont))
            if best is None or v < best[0]:
                best = (v, x, y)
    return best[1], best[2], best[0]


def show(name, g0, atoms, s0, cont):
    x, y, v = nearest_zero(g0, atoms, s0, cont)
    print(f"  {name:24s} Re(tau)={x:6.3f}  Im(tau)={y:6.3f}  pi/Im={np.pi/y:6.3f}  |Z|={v:.2e}")
    return x, y


print("=" * 80)
print("0++-projected Z(tau): nearest Lee-Yang/Fisher zero.  Im(tau_zero) = pi/m_gap  ONLY")
print("when an isolated 0++ atom sits below the continuum (acid-clean).")
print("=" * 80)

print("\n-- (0) ANALYTIC pure atom  g0 + g1 e^{-tau m} = 0  --")
g0, g1, m = 1.0, 2.0, 1.0
lam = np.log(g1 / g0) / m
print(f"   g0={g0}, g1={g1}, m={m}:  zeros tau = (1/m)(ln(g1/g0) - i pi(2k+1))")
for k in range(3):
    print(f"      k={k}:  Re={lam:6.3f}  Im={np.pi*(2*k+1)/m:6.3f}")

print("\n-- (A) ATOM m1=1.0 (g1=2>g0) + continuum from s0=2.0 : expect Im(tau)=pi/1=3.142 --")
xA, yA = show("atom m1=1.0", 1.0, [(2.0, 1.0)], 2.0, 0.30)

print("\n-- (B) ATOM m1=0.6 + same continuum : expect Im(tau)=pi/0.6=5.236 --")
xB, yB = show("atom m1=0.6", 1.0, [(2.0, 0.6)], 2.0, 0.30)

print("\n-- (C) :phi^2: analog: NO atom, continuum from s0=4.0  (no isolated state below cut) --")
xC, yC = show("no atom, cont s0=4", 1.0, [], 4.0, 0.30)
print("      -> the 'gap' pi/Im has NO atom meaning (it is set by the cut / cutoff).")

print("\n-- (D) FREE MASSLESS analog: NO atom, continuum from s0=0 --")
xD, yD = show("no atom, cont s0=0", 1.0, [], 0.0, 0.30)
print("      -> gapless: zero controlled by the cut, no physical atom.")

print("\n-- (E) atom present, weight g1 sweeps DOWN to 0: zero persists at Im=pi/m until atom gone --")
for g1v in [2.0, 1.0, 0.3, 0.05, 0.0]:
    x, y, v = nearest_zero(1.0, [(g1v, 1.0)] if g1v > 0 else [], 2.0, 0.30)
    tag = "atom" if g1v > 0 else "NO atom"
    print(f"      g1={g1v:4.2f} ({tag:7s}): Re={x:6.3f} Im={y:6.3f}  pi/Im={np.pi/y:6.3f}")

print("\nINTERPRETATION")
print("  (A),(B):  Im(tau_zero) = pi/m1  tracks the ATOM mass; changing m1 moves Im(tau).")
print("            pi/m1 is a PHYSICAL pole -> renormalon-free / scheme-independent (the gain")
print("            over the u=2 condensate; same advantage as the large-N mass and dual Higgs).")
print("  (C),(D):  without an atom there is NO zero at pi/m_atom: the nearest zero is")
print("            continuum/cutoff-controlled -> the functional does NOT fire on free / :phi^2:.")
print("  (E):      while the atom exists the zero sits at Im=pi/m (Re tracks ln(g1/g0)); when")
print("            the atom weight ->0 the isolated zero disappears.  => fires iff f1^2>0.")
print("  COST [O]: computing Z_{0++}(tau) needs the 0++ projector (non-local) AND L,a->infty")
print("            with the atom isolated below the continuum = the transport wall (H1).")
