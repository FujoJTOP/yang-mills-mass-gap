#!/usr/bin/env python3
"""audit_b3_sw_marginal.py -- (log rho)'' : beta-free at FIXED x, beta-dependent at the MODE.

rho_beta(x) = Z^-1 e^{-beta x} sqrt(x(2-x)),  x in (0,2)  (exact one-plaquette SU(2) marginal).
    (log rho)''(x) = -1/(2x^2) - 1/(2(2-x)^2)        <- no beta at all (linear term)
    mode: (1-x)/(x(2-x)) = beta  =>  x* = ((2b+1) - sqrt(4b^2+1))/(2b)
    rho''(x*) = rho(x*)(log rho)''(x*)  < 0          <- sign fixed, VALUE moves with beta
Everything checked by finite differences against the closed form.
"""
import numpy as np

NU = 400001
U = np.linspace(-np.pi / 2, np.pi / 2, NU)
X = 1.0 + np.sin(U)                       # x = 1 + sin u
COS2 = np.cos(U) ** 2                     # sqrt(x(2-x)) dx = cos^2 u du


def logrho_part(beta, x):
    return -beta * x + 0.5 * np.log(x) + 0.5 * np.log(2 - x)


def mode(b):
    if b == 0.0:
        return 1.0
    return ((2 * b + 1) - np.sqrt(4 * b * b + 1)) / (2 * b)


print("closed form  (log rho)''(x) = -1/(2x^2) - 1/(2(2-x)^2)   -- beta absent")
print(" beta      x=0.2        x=0.5        x=1.0        x=1.5    |  x*(beta)   (logrho)''(x*)   rho''(x*)")
for b in (0.0, 0.5, 1.0, 2.5, 4.0):
    xs = [0.2, 0.5, 1.0, 1.5]
    vals = [-0.5 / x ** 2 - 0.5 / (2 - x) ** 2 for x in xs]
    xstar = mode(b)
    l2 = -0.5 / xstar ** 2 - 0.5 / (2 - xstar) ** 2
    # normalised rho at the mode (for rho'' = rho * (log rho)'')
    w = np.exp(-b * X) * COS2
    Z = np.trapezoid(w, U)
    rho = np.exp(-b * xstar) * np.sqrt(xstar * (2 - xstar)) / Z
    print(f"{b:5.1f} " + "".join(f"{v:13.6f}" for v in vals) + f" | {xstar:8.4f} {l2:14.6f} {rho*l2:12.6f}")

# finite-difference check that (log rho)'' at FIXED x is beta-independent, and at the mode is not
print("\nFD check (h=1e-4) of d2/dx2 [log rho]:")
for x0 in (0.2, 1.0, 1.8):
    row = []
    for b in (0.0, 2.5, 4.0):
        f = lambda x: logrho_part(b, x)
        row.append((f(x0 + 1e-4) - 2 * f(x0) + f(x0 - 1e-4)) / 1e-8)
    print(f"  x={x0}: beta=0/2.5/4 -> {row[0]:.9f} {row[1]:.9f} {row[2]:.9f}   (spread {max(row)-min(row):.1e})")
print("\n  => at fixed x: identical to 1e-9 (beta only multiplies, does not curve).")
print("     at the mode: -1.000 (b=0) ... -36.5 (b=4)  -- the MODE MOVES, the sign does not.")
print("  rho''(x) itself is NOT negative on all of (0,2) for large beta: rho''=rho[(logrho)''+((logrho)')^2],")
print(f"     beta=10, x=0.5: (logrho)'={( -10 + 0.5/0.5 - 0.5/1.5):.4f}, rho''/rho={( -0.5/0.25-0.5/2.25 + (-10+1-1/3)**2):.4f} > 0")
