"""
ROOF-RESURGENCE witness (cheap, one shot, numpy/scipy only).
Toy renormalon = Euler-type series with Borel transform B(t)=1/(1-t), pole at t=1.
With x := 1/(2s) the nonperturbative scale is e^{-1/x} = e^{-2s}  ->  the "u=2" IR renormalon
of the gluon condensate (the pole strength/power is what matters, not the coupling bookkeeping).

What this demonstrates (ALL exact by the residue theorem / resurgent structure):
  (A) the lateral-Borel ambiguity equals 2*pi*i * Res = 2*pi*i * (-e^{-1/x})
      -> a RESURGENCE INVARIANT, computable from large-order coefficients a_n ~ n! ;
  (B) the MEDIAN sum (= PV) is real and canonical, BUT the transseries parameter
      sigma (the condensate) is FREE: F = PV + sigma*e^{-1/x} is real for any real sigma,
      and the SAME perturbative series is consistent with any sigma ;
  (C) reality / perturbation theory do NOT pin sigma. Resurgence pins the OBSTRUCTION,
      not the condensate magnitude.
No claim about Yang-Mills. Illustrative/toy. Clay 0%; no QED.
"""
import numpy as np, math
from scipy.integrate import quad

def pv_borel(x, delta=1e-3):
    """Principal-value Borel sum of the Euler series: PV int_0^inf e^{-t/x}/(1-t) dt."""
    f = lambda t: np.exp(-t/x)/(1.0 - t)
    a = quad(f, 0.0, 1.0-delta, limit=500, epsabs=1e-13, epsrel=1e-11)[0]
    b = quad(f, 1.0+delta, np.inf, limit=500, epsabs=1e-13, epsrel=1e-11)[0]
    return a + b

print("="*74)
print("(A) Stokes discontinuity  S+ - S-  =  2*pi*i*Res   (the resurgence invariant)")
print("="*74)
print(f"{'x':>7} {'s=1/(2x)':>9} {'PV/median':>15} {'Res=-e^-1/x':>15} {'2*pi*|Res|':>14}")
for x in [0.5, 0.2, 0.1, 0.05, 0.02]:
    pv = pv_borel(x)
    res = -np.exp(-1/x)            # Res_{t=1} e^{-t/x}/(1-t)
    disc = 2*np.pi*abs(res)        # |S+ - S-|
    print(f"{x:7.3f} {1/(2*x):9.3f} {pv:15.6e} {res:15.6e} {disc:14.6e}")
print("=> S_+-S_- = 2*pi*i*Res EXACTLY (residue theorem). |disc| = 2*pi*e^{-1/x}.")
print("   The residue is the Stokes constant = a resurgence invariant; it is")
print("   readable from the large-order asymptotics a_n ~ n! (see C). It is the AMBIGUITY.")

print()
print("="*74)
print("(B) the condensate sigma is FREE: reality does not pin it")
print("="*74)
x = 0.1
pv = pv_borel(x); e = np.exp(-1/x)
print("   F(x) = PV(x) + sigma*e^{-1/x};  PV,e real  =>  F real for ANY real sigma.")
for sigma in [0.0, -0.5, 1.7, 10.0]:
    print(f"     sigma={sigma:6.2f}  ->  F = {pv + sigma*e: .6e}   (real)")
print("=> SAME perturbative series + SAME PV, arbitrary condensate sigma: sigma is a")
print("   free transseries parameter (integration constant). The median sum is canonical")
print("   but is NOT the physical value unless sigma is supplied NONPERTURBATIVELY.")

print()
print("="*74)
print("(C) large-order coefficients encode the residue (Stokes), NOT sigma")
print("="*74)
a = np.array([math.factorial(n) for n in range(1, 12)], float)
print("   a_n/n!  =", np.round(a/np.array([math.factorial(n) for n in range(1,12)],float), 8))
print("=> a_n ~ n!  (Borel pole at t=1, residue -1): the invariant is READ off PT.")
print("   sigma does NOT appear in {a_n}. This is the exact, sharpest statement of the wall:")
print("   perturbation theory (with resurgence) knows the ambiguity, never the condensate.")
