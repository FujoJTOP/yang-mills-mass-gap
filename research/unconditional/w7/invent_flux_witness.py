"""
invent_flux_witness.py -- INVENT-FLUX acid test.

Species: a center-flux twisted-gap functional Gamma(z): the 0++ spectral gap as a
function of a discrete center flux (twist) threaded through a spatial cycle.

Claim tested: a threaded twist alpha on a spatial circle shifts the momentum
lattice p_n = (2 pi n + alpha)/L. The LOCAL mass gap of a free field is
twist-INDEPENDENT in the thermodynamic limit (the alpha-dependence is O(1/L), a
finite-volume artifact); the twist-sensitive global object is the twisted
partition function (Casimir energy), NOT the local gap.

This is the free-field shadow of the [T] fact that the 0++ operator (tr F^2) is
center-neutral (zero N-ality) and therefore blind to the 1-form center-flux
superselection sector: a neutral operator cannot change the flux sector, so its
correlator -- hence the 0++ mass -- is sector-independent.
"""
import numpy as np


def gap(m, alpha, L):
    """single-particle gap: min_n sqrt(p_n^2 + m^2), p_n = (2 pi n + alpha)/L."""
    n = np.arange(-4000, 4001)
    p = (2.0 * np.pi * n + alpha) / L
    return float(np.min(np.sqrt(p**2 + m**2)))


def logZ(m, alpha, L, beta=1.0):
    """1-loop log of the twisted partition function of a free complex scalar on S^1."""
    n = np.arange(-4000, 4001)
    p = (2.0 * np.pi * n + alpha) / L
    w = np.sqrt(p**2 + m**2)
    x = np.exp(-beta * w)
    return float(-np.sum(np.log1p(-x)))


print("=" * 72)
print("(1) LOCAL GAP vs TWIST alpha, MASSLESS m=0 (acid: must -> 0 as L->inf):")
print("    the twist lifts a zero mode only at finite L; gap = pi*alphabar/L -> 0.")
for L in [4, 8, 16, 32, 64]:
    g0, gq = gap(0.0, 0.0, L), gap(0.0, np.pi / 2, L)
    print(f"  L={L:3d}: gap(alpha=0)={g0:.5f}   gap(alpha=pi/2)={gq:.5f}")

print()
print("(2) LOCAL GAP vs TWIST alpha, MASSIVE m=1 (gapped): must -> m=1 for every alpha.")
for L in [4, 8, 16, 32, 64]:
    g0, gq = gap(1.0, 0.0, L), gap(1.0, np.pi / 2, L)
    print(f"  L={L:3d}: gap(alpha=0)={g0:.5f}   gap(alpha=pi/2)={gq:.5f}")

print()
print("(3) TWISTED PARTITION FUNCTION (Casimir) IS alpha-dependent:")
print("    the flux response lives in logZ, a GLOBAL functional, NOT in the local gap.")
for L in [8, 16]:
    vals = [logZ(1.0, a, L) for a in [0.0, np.pi / 4, np.pi / 2, np.pi]]
    print(f"  L={L}: logZ(alpha=0,pi/4,pi/2,pi) = " + ", ".join(f"{v:.6f}" for v in vals))

print()
print("(4) SECTOR-INDEPENDENCE of a NEUTRAL operator (1-form toy, exact):")
print("    neutral op (0 one-form charge) -> same in every flux sector m;")
print("    charged op (charge q) -> phase e^{2 pi i m q/N}.")
N = 5
for m in range(N):
    phase_neutral = 1.0 + 0.0j
    phase_charged = np.exp(2j * np.pi * m * 1 / N)
    print(f"  flux sector m={m}: neutral-op phase={phase_neutral.real:.3f}   "
          f"charged-op phase={phase_charged.real:+.3f}{phase_charged.imag:+.3f}j")
print("  => the neutral 0++ operator cannot distinguish flux sectors; only 1-form ops can.")
