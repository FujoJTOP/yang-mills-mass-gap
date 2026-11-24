"""
construct_luscher.py -- witness for CONSTRUCT-LUSCHER
Effective-string (Nambu-Goto + closed-string Casimir) prediction for the 0++ glueball
mass and finite-volume scaling. Lightweight numpy only. All masses in sqrt(sigma).

Closed-string spectrum (D=4, D-2=2 transverse bosons), worldsheet spatial circumference l:
    E(N,Nb;l) = sigma*l + (2 pi / l)(N + Nb - (D-2)/12),   level matching N = Nb
    E0(l) = sigma*l - pi(D-2)/(6 l) = sigma*l - pi/(3 l)   [ground state -> 0++ scalar]
    E1(l) = sigma*l + 11 pi/(3 l)                          [first excited: 0++ trace + 2++]
Pure NG: E0'(l) = sigma + pi/(3 l^2) > 0  => NO interior minimum (closed string shrinks).
The stable glueball needs the rigidity kappa:  E(l) = 2 pi sigma l + pi kappa / l,
    l* = sqrt(kappa/(2 sigma)),  m = 2 pi sqrt(2 sigma kappa),  m/sqrt(sigma) = 2 pi sqrt(2 kappa).
Torelon (flux tube wrapping an L-torus): M(L) = sigma L - pi(D-2)/(6 L) = sigma L - pi/(3 L) + O(1/L^3).
Open-string Luscher term:  V(R) = sigma R - pi(D-2)/(24 R) = sigma R - pi/(12 R)  [D=4, famous -pi/12].
"""
import numpy as np

sig = 1.0
D = 4
ntr = D - 2  # transverse bosons = 2

print("=" * 70)
print("(1) Closed-string effective-string spectrum (sigma=1, D=4)")
print("    E0(l) = l - pi/(3l) ;  E1(l) = l + 11pi/(3l)")
print("=" * 70)
for l in [0.3, 0.5, 1.0, 2.0, 3.0]:
    E0 = sig * l - np.pi * ntr / (6 * l)
    E1 = sig * l + 11 * np.pi / (3 * l)
    print(f"  l={l:4.2f}   E0={E0:9.3f}   E1={E1:9.3f}   gap={E1 - E0:8.3f}")

print("\n  dE0/dl = sigma + pi/(3l^2) > 0 for all l>0  => NO interior minimum.")
print("  => pure NG closed string SHRINKS; no stable 0++ glueball without rigidity.")

print("\n" + "=" * 70)
print("(2) Rigid-string completion: E(l) = 2 pi sigma l + pi kappa / l")
print("    l* = sqrt(kappa/(2 sigma));  m = 2 pi sqrt(2 sigma kappa)")
print("    m/sqrt(sigma) = 2 pi sqrt(2 kappa)")
print("=" * 70)


def m_over_sqrtsig(kappa):
    return 2 * np.pi * np.sqrt(2 * kappa)


for target in [3.3, 3.7, 3.9]:
    kappa = (target / (2 * np.pi)) ** 2 / 2
    lstar = np.sqrt(kappa / (2 * sig))
    m = 2 * np.pi * np.sqrt(2 * sig * kappa)
    print(f"  target m/sqrt(sig)={target:.2f} -> kappa={kappa:.4f}, l*={lstar:.3f}, m={m:.3f}")

kmin = (3.3 / (2 * np.pi)) ** 2 / 2
kmax = (3.9 / (2 * np.pi)) ** 2 / 2
print(f"\n  => lattice band 3.3-3.9 gives kappa in [{kmin:.3f}, {kmax:.3f}]; kappa = H1 (rigidity).")

print("\n" + "=" * 70)
print("(3) Torelon finite-volume scaling: M(L)=sigma L - pi(D-2)/(6L)=sigma L - pi/(3L)")
print("    relative correction = -pi/(3 sigma L^2)")
print("=" * 70)
for L in [2, 3, 4, 6, 10]:
    M = sig * L - np.pi * ntr / (6 * L)
    rel = (M - sig * L) / (sig * L)
    print(f"  L={L:3d}   M(L)={M:9.4f}   M/(sigmaL)={M / (sig * L):.5f}   rel.corr={rel:+.5f}")

print("\n" + "=" * 70)
print("(4) Luscher coefficient: open vs closed (D=4)")
print("    open  : V(R) = sigma R - pi(D-2)/(24 R) = sigma R - pi/(12 R)  [famous -pi/12]")
print("    closed: E0(l) = sigma l - pi(D-2)/(6 l)  = sigma l - pi/(3 l)   [4 x open]")
print("=" * 70)
print(f"  open coefficient   -pi(D-2)/24 = {-np.pi * ntr / 24:.6f}")
print(f"  closed coefficient -pi(D-2)/6  = {-np.pi * ntr / 6:.6f}")
print(f"  ratio closed/open  = {(-np.pi * ntr / 6) / (-np.pi * ntr / 24):.1f}")

print("\n" + "=" * 70)
print("VERDICT (no QED): absolute 0++ mass needs kappa (=H1);")
print("effective string fixes the finite-V scaling only (-pi/3 torelon; -pi/12 open).")
print("=" * 70)
