"""
invent_taildep_witness.py -- INVENT-TAILDEP witness (2026-11-24).

Extreme-value / tail-dependence functional of the 0++ (tr F^2) measure:
  tail index xi (GPD shape), extremal index theta_EV (clustering),
  tail-dependence / extremal coefficient.
Tests the acid test and kappa-homogeneity on Gaussian adversaries:
  free-massless  ~ AR(1) with rho -> 1  (long-range),
  gapped YM-proxy ~ AR(1) small rho    (short-range),
  :phi^2: adversary ~ square of the series (Wick / chi^2 composite).

Run:  py -3.12 -I research/unconditional/w7/invent_taildep_witness.py
"""

import numpy as np
from scipy import stats

rng = np.random.default_rng(20261124)


def ar1(n, rho):
    x = np.empty(n)
    x[0] = rng.standard_normal()
    e = rng.standard_normal(n)
    s = np.sqrt(max(1.0 - rho * rho, 0.0))
    for i in range(1, n):
        x[i] = rho * x[i - 1] + s * e[i]
    return x


def runs_theta(x, q=0.99, run=3):
    """runs-declustering estimate of the extremal index theta_EV."""
    t = np.quantile(x, q)
    ex = np.where(x > t)[0]
    if len(ex) == 0:
        return 1.0, 0
    clusters = 1
    for i in range(1, len(ex)):
        if ex[i] - ex[i - 1] > run:
            clusters += 1
    return clusters / len(ex), len(ex)


def gpd_xi(x, q=0.95):
    """GPD shape (= tail index) from peaks-over-threshold."""
    t = np.quantile(x, q)
    ex = x[x > t] - t
    if len(ex) < 20:
        return np.nan
    xi = stats.genpareto.fit(ex, floc=0)[0]
    return xi


print("=" * 72)
print("W1 | tail index xi and extremal index theta_EV vs correlation strength")
print("     AR(1) Gaussian: rho=0 short-range/gapped; rho->1 long-range/massless")
print("=" * 72)
print(f"{'rho':>6} {'xi(GPD)':>9} {'theta_EV':>9} {'n_exc':>8}")
N = 200000
for rho in [0.0, 0.3, 0.6, 0.9, 0.99]:
    x = ar1(N, rho)
    xi = gpd_xi(x)
    th, nex = runs_theta(x)
    print(f"{rho:>6.2f} {xi:>9.3f} {th:>9.3f} {nex:>8d}")
print()
print("  -> xi ~ 0 (Gumbel, light tail) EVERYWHERE: tail index does NOT discriminate.")
print("  -> theta_EV decreases monotonically as correlation grows (rho->1 = massless):")
print("     theta_EV is a MONOTONE FUNCTION OF CORRELATION STRENGTH = kappa (mass).")

print()
print("=" * 72)
print("W2 | composite adversary :phi^2:  (square the series; Wick / chi^2)")
print("=" * 72)
print(f"{'rho':>6} {'xi':>9} {'theta_EV':>9}")
for rho in [0.0, 0.3, 0.6, 0.9]:
    x = ar1(N, rho)
    y = x * x
    xi = gpd_xi(y)
    th, _ = runs_theta(y)
    print(f"{rho:>6.2f} {xi:>9.3f} {th:>9.3f}")
print()
print("  -> gapped-YM-proxy (small rho) and :phi^2: (gapped composite) give the SAME")
print("     theta_EV (~1): EVT cannot distinguish gapped-YM from gapped-atomless :phi^2:")
print("     => ATOM-BLIND.  Massless (rho->1) fires theta<1 on the WRONG side.")

print()
print("=" * 72)
print("W3 | acid-test sign impossibility")
print("=" * 72)
print("  want F>0 on YM(gapped,atom);  F<=0 on free-massless;  F<=0 on :phi^2:(gapped,atomless)")
print("  EVT only gives:  theta(gapped YM) ~ theta(:phi^2:) ~ 1  >  theta(massless)")
print("  * F = theta-1 : massless<0 OK, :phi^2:=0 OK, but YM=0 (NOT >0)  -> FAIL")
print("  * F = 1-theta : YM=0 (NOT >0), :phi^2:=0, but massless>0       -> FAIL")
print("  => no sign convention passes: the functional cannot be positive ONLY on the atom.")
print("  (B2: zero-T 4D YM is single-phase; EVT/LDP is a transition detector -> nothing to detect.)")
