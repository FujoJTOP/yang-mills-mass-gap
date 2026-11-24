"""route_new_rg.py -- ROUTE-NEW-RG witness.
Shows there is NO small parameter at the crossover: the strong-coupling weight
u = I2/I1 (SU(2)) and the weak-coupling coupling g^2 = 4/beta are BOTH O(1) at
beta = 2.5, whereas each is small in the two rigorous regimes (u small at strong
coupling; g^2 small at weak coupling). This is why no constructive-RG small
parameter exists at the crossover -- the wall is the indefinite (Savvidy) action.
Run: py -3.12 -I route_new_rg.py
"""
import numpy as np
from scipy.special import iv


def u_weight(beta):
    # SU(2): leading strong-coupling character weight u = I2(beta)/I1(beta)
    return iv(2, beta) / iv(1, beta)


def g_sq(beta):
    # SU(2) Wilson: beta = 4/g^2
    return 4.0 / beta


print("beta        u= I2/I1     g^2=4/beta   regime")
for b in (1 / 48, 1.0, 1.6, 2.0, 2.5, 7.0):
    tag = "STRONG (u small)" if u_weight(b) < 0.1 else (
        "WEAK (g^2 small)" if g_sq(b) < 0.6 else "CROSSOVER (no small param)")
    print(f"{b:8.4f}   {u_weight(b):9.4f}   {g_sq(b):9.3f}   {tag}")

u25, g25 = u_weight(2.5), g_sq(2.5)
print("\nAt beta=2.5:  u=%.4f, g^2=%.3f  -> both O(1): NO small parameter."
      % (u25, g25))
print("Strong coupling (beta=1/48): u=%.2e -> small parameter EXISTS."
      % u_weight(1 / 48))
