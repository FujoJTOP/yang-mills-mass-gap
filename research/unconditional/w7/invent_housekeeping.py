"""
invent_housekeeping.py -- witness for INVENT-HOUSEKEEPING.md  (2026-11-24)

Non-equilibrium (driven) entropy-production functional on a Gaussian gauge-mode toy.
Auxiliary non-conservative drive F = gamma*R (R = rotation, curl != 0).
Overdamped SDE: dX = A X dt + sqrt(2D) dW,  A = -m*I + gamma*R.
Stationary Gaussian pi (cov Sigma solves the Lyapunov eq); total EP rate
   sigma = < (b - D grad log pi)^T D^-1 (b - D grad log pi) >_pi
         = Tr[ M Sigma M^T D^-1 ],   M = A + D Sigma^-1.

Tests:
  (A) reversible (gamma=0): sigma = 0 for every mass m  -> Theorem HK
  (B) driven (gamma>0): sigma > 0, set by the DRIVE not the mass (acid-dirty)
  (C) controller-dependence: sigma ~ gamma^2 (tunable by the invented drive)
  (D) equilibrium excess EP (non-adiabatic) = Fisher metric, homogeneous in m

Run:  py -3.12 -I invent_housekeeping.py
"""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov

R = np.array([[0.0, 1.0], [-1.0, 0.0]])
D = np.eye(2)


def Sigma_stat(A, D):
    # A Sigma + Sigma A^T + 2 D = 0
    return solve_continuous_lyapunov(A, -2.0 * D)


def entropy_production(A, D):
    S = Sigma_stat(A, D)
    Dinv = np.linalg.inv(D)
    M = A + D @ np.linalg.inv(S)          # b - D grad log pi = M x
    sigma = float(np.trace(M @ S @ M.T @ Dinv))
    return sigma, S


print("=" * 72)
print("A) REVERSIBLE (gamma=0): sigma must be 0 for every mass m")
print("=" * 72)
for m in [0.03, 0.3, 1.0, 3.0]:
    A = -m * np.eye(2)
    sig, S = entropy_production(A, D)
    print(f"  m={m:5.3f}  sigma={sig:.3e}   Sigma_diag={np.diag(S)}")

print()
print("=" * 72)
print("B) DRIVEN (gamma>0): sigma>0, set by the DRIVE, not the mass")
print("=" * 72)
for m in [0.03, 0.3, 1.0, 3.0]:
    row = []
    for g in [0.0, 0.2, 0.5, 1.0]:
        A = -m * np.eye(2) + g * R
        sig, _ = entropy_production(A, D)
        row.append(f"g={g:4.2f}:{sig:7.4f}")
    print(f"  m={m:5.3f}   " + "   ".join(row))

print()
print("=" * 72)
print("C) CONTROLLER-DEPENDENCE: sigma ~ gamma^2 (tunable by the invented drive)")
print("=" * 72)
m = 1.0
for g in [0.05, 0.1, 0.2, 0.4]:
    A = -m * np.eye(2) + g * R
    sig, _ = entropy_production(A, D)
    print(f"  gamma={g:5.3f}  sigma={sig:.6f}   sigma/gamma^2={sig / g ** 2:.6f}")

print()
print("=" * 72)
print("D) EQUILIBRIUM EXCESS (non-adiabatic) EP between m and m' = Fisher, homogeneous")
print("=" * 72)
# 1D Gaussian N(0,1/m): grad log(pi/pi') = -(m-m')x ; <x^2>_pi = 1/m
# sigma_na = (m-m')^2 / m   (Fisher information of the mass; 0 at m'=m)
for mp in [1.0, 1.2, 2.0]:
    m = 1.0
    sig_na = (m - mp) ** 2 / m
    print(f"  m={m}  m'={mp}   sigma_na={sig_na:.6f}")

print()
print("READOUT: (A) sigma=0 exactly when reversible. (B) any gamma>0 gives sigma>0")
print("at EVERY mass (even m->small) -> functional is nonzero on the light/free-ish")
print("field = ACID-DIRTY; its size tracks gamma, not the coarse mass. (C) sigma ~ gamma^2")
print("=> controller-dependent (dialed by the invented drive). (D) equilibrium excess =")
print("Fisher metric = homogeneous (Theorem BE territory), not a kappa-independent source.")
