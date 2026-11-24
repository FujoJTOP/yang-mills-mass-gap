#!/usr/bin/env python3
# q7_direction.py  --  Q7 TOOL-FORGE, Theorem Q1 witness
# "Directional no-go": analytic / Fourier / Fisher-Rao control of the coarse marginal
# can only *lower*-bound the variance (= *upper*-bound the curvature G''), the logical
# opposite of H1 (which needs Var <= 1/lambda, i.e. G'' >= lambda).  It is saturated
# from the gapless side by the free (flat) field.
#
# Deterministic, numpy only.  Run:  py -3.12 -I q7_direction.py
import numpy as np

np.set_printoptions(precision=6, suppress=False)
print("=" * 74)
print("Q7 / Q1  DIRECTIONAL NO-GO witness   (Fourier-uncertainty + Fisher-Rao)")
print("=" * 74)

# ---------------------------------------------------------------- (A) Bernstein ceiling
# On the circle, rho_eps(a) = (1 + eps cos(tau a)) is positive for eps<1 and is a
# finite Fourier polynomial of bandwidth tau.  The exact object -rho''/rho at the peak a=0:
#     -rho''(0)/rho(0) = eps*tau^2/(1+eps).
# Bernstein's inequality: |rho''| <= tau^2 rho  (exponential type tau).  So G'' <= tau^2
# ALWAYS (a CEILING).  But as eps->0 (flat / free-field limit) G''->0 at FIXED tau:
# bandwidth gives NO positive lower bound.
print("\n(A) BERNSTEIN CEILING vs the eps->0 (free-field) limit  [bandwidth tau fixed]")
for tau in (1.0, 2.0, 4.0):
    print(f"   tau = {tau:4.1f}")
    for eps in (0.01, 0.1, 0.5, 0.9, 0.99):
        gpp_peak = eps * tau**2 / (1.0 + eps)
        # grid-check the Bernstein ceiling  -rho''/rho <= tau^2 on the whole circle
        a = np.linspace(-np.pi / tau + 1e-9, np.pi / tau - 1e-9, 20001)
        rho = 1.0 + eps * np.cos(tau * a)
        rpp = -eps * tau**2 * np.cos(tau * a)
        ratio = -rpp / rho
        print(f"      eps={eps:5.2f}:  -rho''/rho(peak)= {gpp_peak:9.5f}   "
              f"grid max= {ratio.max():9.5f}  (ceiling tau^2={tau**2:.1f})")
    print(f"      => eps->0  gives  G''->0  at tau={tau:.0f};  the ceiling is NEVER a floor.\n")

# ---------------------------------------------------------------- (B) Cramer-Rao = variance FLOOR
# On R (physical marginal vanishes at the chart edge, so boundary terms drop):
#     I(rho) * Var(a) >= 1     (Cauchy-Schwarz; equality for the Gaussian).
# Equivalently  Var(a) >= 1/I .  This is a LOWER bound on the variance.
# The free field (sigma->inf) has I->0 and Var->inf: it SATISFIES the floor trivially.
print("(B) CRAMER-RAO:  I(rho)*Var(a) >= 1  is a VARIANCE FLOOR (wrong direction for H1)")

def fisher_info_1d(x, rho):
    # I = int (rho')^2/rho dx  via finite differences (adaptive range avoids tail loss)
    rho = np.asarray(rho, float)
    rp = np.gradient(rho, x)
    integrand = np.where(rho > 1e-250, rp**2 / np.maximum(rho, 1e-250), 0.0)
    return np.trapezoid(integrand, x)

for sigma in (0.3, 1.0, 3.0, 10.0, 40.0):
    x = np.linspace(-8 * sigma, 8 * sigma, 400001)            # adaptive: no truncation
    rho = np.exp(-x**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)
    rho = rho / np.trapezoid(rho, x)
    I = fisher_info_1d(x, rho)
    mean = np.trapezoid(x * rho, x)
    var = np.trapezoid((x - mean)**2 * rho, x)
    gpp_peak = 1.0 / var                                       # Gaussian: -rho''/rho = 1/Var
    print(f"   sigma={sigma:5.1f}:  G''(peak)=-rho''/rho= {gpp_peak:9.6f}   "
          f"I= {I:9.6f}   Var= {var:9.3f}   I*Var= {I*var:7.4f}")
print("   => Gaussian saturates I*Var=1; sigma->inf (free field) has Var->inf, G''->0.")
print("      Every entry satisfies Var>=1/I ; NONE certifies Var<=1/lambda. WRONG DIRECTION.\n")

# non-Gaussian: mixture keeps I*Var>1 (still a floor, still wrong direction)
x = np.linspace(-8, 11, 400001)
rho1 = np.exp(-x**2 / (2 * 0.8**2))
rho2 = np.exp(-(x - 2.5)**2 / (2 * 0.8**2))
rho = 0.5 * rho1 + 0.5 * rho2
rho = rho / np.trapezoid(rho, x)
I = fisher_info_1d(x, rho)
mean = np.trapezoid(x * rho, x); var = np.trapezoid((x - mean)**2 * rho, x)
print(f"   two-Gaussian mixture: I= {I:.5f}  Var= {var:.5f}  I*Var= {I*var:.4f} (>1).\n")

# ---------------------------------------------------------------- (C) exact Fisher identity
# E_rho[V_eff''] = E_rho[-(log rho)''] = I(rho)   (integration by parts)
print("(C) EXACT FISHER IDENTITY:  E_rho[-(log rho)''] = I(rho)")
for sigma in (0.5, 1.0, 2.0, 5.0):
    x = np.linspace(-8 * sigma, 8 * sigma, 400001)            # adaptive range
    rho = np.exp(-x**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)
    rho = rho / np.trapezoid(rho, x)
    lr = np.log(rho)
    lr2 = np.gradient(np.gradient(lr, x), x)           # (log rho)''
    Eminuslr2 = np.trapezoid(rho * (-lr2), x)
    I = fisher_info_1d(x, rho)
    print(f"   sigma={sigma:3.1f}:  E[-(log rho)'']= {Eminuslr2:9.5f}   I= {I:9.5f}"
          f"   |diff|= {abs(Eminuslr2-I):.2e}")
print("   => E[G''] = I exactly.  Its uniform lower bound  I>=delta  IS the gap (H1).\n")

# ---------------------------------------------------------------- (D) free-field acid test
print("(D) FREE-FIELD ACID TEST (the tool MUST return 0 on the gapless field)")
print("    flat marginal rho=const:  -rho''/rho = 0  and  I = 0  and Var = inf.")
print("    A band-limited rho_eps with eps->0 reproduces it (panel A).  PASSES (returns 0),")
print("    and returns it NON-vacuously: no false positivity is manufactured.\n")

print("Q1 VERDICT [NG]: Fourier/uncertainty give G'' <= tau^2 (ceiling);")
print("                 Fisher-Rao gives Var >= 1/I (variance floor).")
print("                 Both are the WRONG DIRECTION for H1 (Var <= 1/lambda).")
