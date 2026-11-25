# trg_ekr.py -- witness for research/v3/species/TRG-EKR.md
# EKR-type strict tensor RG (Ebel-Kennedy-Rychkov) vs 4D pure YM at beta=2.5.
# Light, single process, deterministic. py -3.12, numpy 2.4.6 + scipy 1.18.1.
#
# Sections
# [A] exact SU(2) character coefficients c_n = I_n - I_{n+2}; completeness sum = e^beta
# [B] EKR "deviation from A_*" proxy: r(beta) = ||nontrivial||/||trivial|| = Phi(0)
#     and the beta at which r equals EKR's proven basin (0.02 per sector, 0.02*sqrt(63) total)
# [C] high-T calibration: EKR proven Ising bound beta<=0.12 vs tanh(0.12);  XY 0.18995 vs beta_c
# [D] irrep truncation j<=J at beta=2.5: single-plaquette weight tail, J for 1% / 0.1%
#     and the *volume-consistent* J for a 4D L^4 torus
# [E] tail asymptotics: factorial (superexponential) vs e^{-cJ} vs polynomial

import numpy as np
from scipy.special import iv

np.set_printoptions(precision=6, suppress=False)
out = []


def log(s):
    out.append(s)
    print(s)


# ----------------------------------------------------------------------------
# [A] exact SU(2) character coefficients
# Wilson action S = beta * sum_p (1 - (1/2) Re tr U_p);  e^{beta cos th} = sum_n c_n chi_n
# c_n(beta) = I_n(beta) - I_{n+2}(beta),  n = 2j,  dim = n+1
# ----------------------------------------------------------------------------
NMAX = 60


def coeffs(beta, nmax=NMAX):
    n = np.arange(nmax + 1)
    return iv(n, beta) - iv(n + 2, beta)


log("[A] exact SU(2) character coefficients  c_n = I_n(beta) - I_{n+2}(beta)")
log("    beta        sum(n+1)c_n        e^beta          rel.err     c_0/e^beta   u=c_1/c_0")
for beta in [1 / 48, 0.12, 0.30, 0.5, 1.0, 2.0, 2.30, 2.5, 4.0, 6.0]:
    c = coeffs(beta)
    n = np.arange(len(c))
    tot = float(np.sum((n + 1) * c))
    log(f"    {beta:<10.5g}  {tot:16.10f}  {np.exp(beta):16.10f}  {abs(tot/np.exp(beta)-1):9.2e}"
        f"  {c[0]/np.exp(beta):11.6f}  {c[1]/c[0]:10.6f}")

# ----------------------------------------------------------------------------
# [B] EKR basin vs the beta=2.5 deviation.
# A = A_* + b with A_* the trivial (all-index-0) tensor, normalized A_0000 = 1.
# In the character expansion the natural relative size of the non-trivial part is
#     r(beta) := Phi(0) = sum_{n>=1} (n+1) c_n / c_0 = (e^beta - c_0)/c_0 ,
# i.e. ||b||/||A_*|| measured on one plaquette (the local weight).
# EKR Theorem 3.1: basin = 0.02 in HS norm per sector, 63 sectors => 0.02*sqrt(63) total.
# ----------------------------------------------------------------------------
BASIN_SECTOR = 0.02
BASIN_TOTAL = 0.02 * np.sqrt(63)
log("")
log(f"[B] EKR proven basin: {BASIN_SECTOR} per sector, {BASIN_TOTAL:.5f} total (HS norm, 63 sectors)")


def r_of_beta(beta):
    c = coeffs(beta)
    return (np.exp(beta) - c[0]) / c[0]


from scipy.optimize import brentq

b_tot = brentq(lambda x: r_of_beta(x) - BASIN_TOTAL, 1e-6, 10.0)
b_sec = brentq(lambda x: r_of_beta(x) - BASIN_SECTOR, 1e-6, 10.0)
log(f"    r(beta) = ||b||/||A_*|| (one plaquette):  beta=1/48 -> {r_of_beta(1/48):.5f}"
    f" ; beta=2.5 -> {r_of_beta(2.5):.5f}")
log(f"    beta where r = basin_total  ({BASIN_TOTAL:.4f}) : beta* = {b_tot:.5f}")
log(f"    beta where r = basin_sector ({BASIN_SECTOR:.4f}) : beta* = {b_sec:.5f}")
log(f"    overshoot at the target  beta=2.5 :  r/basin_total = {r_of_beta(2.5)/BASIN_TOTAL:.2f}x"
    f" , r/basin_sector = {r_of_beta(2.5)/BASIN_SECTOR:.2f}x")
log(f"    SFZ/OS78 strong-coupling radius 1/48 = {1/48:.5f}  =>  beta*(EKR basin)/ (1/48) = {b_tot/(1/48):.1f}x")
log(f"    SU(2) bulk crossover beta_c ~ 2.2986  =>  beta*(EKR basin)/beta_c = {b_tot/2.2986:.3f}")

# ----------------------------------------------------------------------------
# [C] calibration: how far does an EKR-type theorem actually reach?
# proven: 2D Ising beta <= 0.12 (beta_c = 0.4407), 2D XY beta <= 0.18995 (beta_c ~ 1.12)
# ----------------------------------------------------------------------------
bc_ising = 0.5 * np.log(1 + np.sqrt(2))
log("")
log("[C] calibration of the proven reach (fraction of beta_c)")
log(f"    2D Ising : beta_max = 0.12 , beta_c = {bc_ising:.5f} -> {0.12/bc_ising:.3f} beta_c"
    f" ;  tanh(0.12) = {np.tanh(0.12):.4f}  (compare basin_total {BASIN_TOTAL:.4f})")
log(f"    2D XY    : beta_max = 0.18995 , beta_c ~ 1.12 -> {0.18995/1.12:.3f} beta_c")
log(f"    SU(2) YM4: target beta = 2.5 , beta_c ~ 2.2986 -> {2.5/2.2986:.3f} beta_c"
    f"   (i.e. ABOVE the crossover)")
log(f"    same-fraction extrapolation to SU(2) YM4: 0.27 * beta_c = {0.12/bc_ising*2.2986:.3f}"
    f"  < 2.5  (short by {2.5/(0.12/bc_ising*2.2986):.2f}x)")

# ----------------------------------------------------------------------------
# [D] irrep truncation j <= J (n = 2j) at beta = 2.5
# w_n = (n+1) c_n / e^beta is the exact per-plaquette weight of irrep n.
# ----------------------------------------------------------------------------
beta = 2.5
c = coeffs(beta, 200)
n = np.arange(len(c))
w = (n + 1) * c / np.exp(beta)
w = w / w.sum()
log("")
log(f"[D] irrep truncation at beta = {beta}:  w_n = (n+1)c_n/e^beta   (sum = {w.sum():.12f})")
log("    n :  " + "  ".join(f"{k:>9d}" for k in range(0, 13)))
log("    w :  " + "  ".join(f"{w[k]:9.2e}" for k in range(0, 13)))
log(f"    argmax_n w_n = {int(np.argmax(w))}  (j = {int(np.argmax(w))/2:.1f}) , w_max = {w.max():.4f}")
log("    J (j<=J)   tail T_J = sum_{n>2J} w_n     1/T_J")
for J in range(0, 13):
    tail = float(np.sum(w[2 * J + 1:]))
    log(f"    J={J:<2d}       {tail:12.4e}            {1/tail if tail>0 else np.inf:12.4e}")
J1 = next(J for J in range(0, 60) if np.sum(w[2 * J + 1:]) < 1e-2)
J2 = next(J for J in range(0, 60) if np.sum(w[2 * J + 1:]) < 1e-3)
log(f"    single plaquette: J(1%) = {J1} , J(0.1%) = {J2}")

# volume-consistent: 4D L^4 torus, #plaquettes = 6 L^4
log("")
log("    volume-consistent truncation (4D L^4 torus, #plaq = 6L^4, need 6L^4 * T_J < eps):")
for L in [4, 6, 8, 12, 16]:
    nplaq = 6 * L ** 4
    for eps in [1e-2, 1e-3]:
        need = eps / nplaq
        Jn = next((J for J in range(0, 200) if np.sum(w[2 * J + 1:]) < need), None)
        log(f"      L={L:<3d} nplaq={nplaq:8d}  eps={eps:.0e}  need T_J<{need:.2e}  =>  J>={Jn}")

# ----------------------------------------------------------------------------
# [E] tail asymptotics: factorial or exponential?
# c_n ~ I_n(beta) ~ (beta/2)^n / n!  =>  w_n ~ (n+1)(beta/2)^n / n!   (superexponential)
# effective per-unit-J rate  c_eff(J) = -ln(T_J)/J ; if exponential it would be constant.
# ----------------------------------------------------------------------------
log("")
log("[E] tail law at beta=2.5:  c_eff(J) = -ln(T_J)/J   (constant => e^{-cJ}; growing => factorial)")
for J in [1, 2, 3, 4, 5, 6, 8, 10, 12]:
    tail = float(np.sum(w[2 * J + 1:]))
    asym = ((beta / 2) ** (2 * J + 2)) / np.math.factorial(2 * J + 2) if hasattr(np, "math") else None
    log(f"    J={J:<3d} T_J={tail:12.4e}  -ln T_J/J = {-np.log(tail)/J:7.3f}")
# closed-form check of the leading factorial term
import math
for J in [2, 4, 6, 8]:
    nlead = 2 * J + 2
    lead = (nlead + 1) * (beta / 2) ** nlead / math.factorial(nlead) / math.exp(beta)
    tail = float(np.sum(w[2 * J + 1:]))
    log(f"    check J={J}: leading factorial term = {lead:.3e} , exact T_J = {tail:.3e} , ratio = {lead/tail:.2f}")

open("trg_ekr.out", "w").write("\n".join(out) + "\n")
print("\n[written to trg_ekr.out]")
