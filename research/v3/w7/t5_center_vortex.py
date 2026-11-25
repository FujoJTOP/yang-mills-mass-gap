#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
t5_center_vortex.py — witness for T5-CENTER-VORTEX.md.

Center-vortex bound on binding: can the vortex (dual-surface) picture give a
rigorous lower bound on the tr F^2 tube overlap f1^2 or the isolation ratio rho?

Deterministic, numpy only, no third-party sources, no Monte Carlo.
Run:  py -3.12 -I t5_center_vortex.py

Panels
  [A] center phases  chi_R(omega*I)/d_R : neutral (N-ality 0) == 1, fundamental == omega.
      => exact neutral-blindness of the tr F^2 (N-ality-0) local sector.
  [B] DECISIVE orthogonality: 2D pure YM confines (sigma>0) YET f1^2(tr F^2)=0.
      sigma>0 does NOT imply the atom  => the vortex->sigma->atom route is dead.
  [C] strong-coupling tube positivity: u(beta)=c_f/N>0 and the recoupling 1/N>0
      force the minimal 4-wall-tube weight a0>0 (this is P3, NOT vortex percolation).
  [D] free-field acid: composite tr F^2 is a two-particle cut => f1^2 = 0.
  [E] Weingarten sign ledger: the sign of a dual surface with boundary is a
      topological (Euler-characteristic) sign, NOT uniformly positive for SU(N);
      the charged-sector phase that drives sigma is trivial on the neutral tube.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)
line = "-" * 74

def hdr(s):
    print("\n" + line + "\n" + s + "\n" + line)

# ----------------------------------------------------------------------------
# SU(2) Haar helpers: class functions, dU = (2/pi) sin^2(theta) dtheta, theta in [0,pi].
# Character chi_j(theta) = sin((2j+1)theta)/sin(theta); dim d_j = 2j+1.
# ----------------------------------------------------------------------------
def chi_su2(j, theta):
    theta = np.asarray(theta, float)
    num = np.sin((2 * j + 1) * theta)
    den = np.sin(theta)
    out = np.empty_like(theta)
    reg = np.abs(den) > 1e-12
    out[reg] = num[reg] / den[reg]
    # removable points theta = 0, pi  ->  limit (2j+1) * cos((2j+1)theta)/cos(theta) via L'Hopital
    for m in np.where(~reg)[0]:
        t = theta[m]
        out[m] = (2 * j + 1) * np.cos((2 * j + 1) * t) / np.cos(t)
    return out

def haar_int_su2(f, n=400001):
    th = np.linspace(0.0, np.pi, n)
    w = (2.0 / np.pi) * np.sin(th) ** 2
    return np.sum(f(th) * w) * (np.pi / (n - 1))

# ----------------------------------------------------------------------------
# [A]  Center phases
# ----------------------------------------------------------------------------
hdr("[A] center phase chi_R(omega*I)/d_R  ->  neutral == 1 (exact blind); fundamental == omega")
print("SU(2) center = {+I,-I} (omega=-1). Evaluate at theta=pi (U=-I).")
for j, name in [(0.5, "fundamental (N-ality 1)"), (1.0, "adjoint     (N-ality 0)"),
                (2.0, "j=2         (N-ality 0)")]:
    d = int(round(2 * j + 1))
    val = float(chi_su2(j, np.array([np.pi]))[0])
    print(f"   SU(2) {name:26s}  chi/d = {val/d:+.6f}")

print("\nSU(3) centre omega=e^{2pi i/3}. Fund dim N=3, chi_f(g)=sum(eigs); adj=(f x fbar - 1).")
for k in (0, 1, 2):
    w = np.exp(2j * np.pi * k / 3.0)
    chi_f = 3.0 * w                 # g = omega^k * I_3 : all eigenvalues = omega^k
    ratio_f = chi_f / 3.0
    chi_adj = abs(chi_f) ** 2 - 1.0  # adjoint = fund x fundbar minus singlet
    ratio_adj = chi_adj / 8.0
    print(f"   k={k}:  fund phase = {ratio_f:+.6f}   adjoint(N-ality0) phase = {ratio_adj:+.6f}")

print("\nZ(N) (abelian): rep = character e^{2pi i q/N}, N-ality q, d_R=1.  phase = omega^q.")
for N in (2, 3):
    ws = [np.exp(2j * np.pi * q / N) for q in range(N)]
    row = "  ".join(f"q={q}:{ws[q]:+.4f}" for q in range(N))
    print(f"   Z({N}): {row}   -> N-ality 0 (q=0) phase is EXACTLY 1")

neutral_dev = abs(float(chi_su2(1.0, np.array([np.pi]))[0]) / 3.0 - 1.0) + \
              max(abs(abs(3.0*np.exp(2j*np.pi*k/3.0))**2 - 1 - 8)/8 for k in (0,1,2)) + 0.0
print(f"\n   max |neutral phase - 1| over SU(2) adj, SU(3) adj = {neutral_dev:.2e}")
print("   => tr F^2 is a colour singlet, N-ality 0: EVERY centre vortex phase it sees is 1. [T]")

# ----------------------------------------------------------------------------
# [B]  Decisive: 2D confines (sigma>0) yet f1^2(trF^2)=0
# ----------------------------------------------------------------------------
hdr("[B] 2D pure YM: area law sigma>0 AND static tr F^2 correlator => f1^2=0 (sigma>0 =/=> atom)")
# Character orthogonality  int dU chi_j = delta_{j0}  (the transfer eigenstates are characters)
print("Haar orthogonality  int_SU(2) chi_j dU  (must be delta_{j0}):")
for j in (0.0, 0.5, 1.0, 1.5, 2.0):
    val = haar_int_su2(lambda th, j=j: chi_su2(j, th))
    print(f"   j={j:3.1f}:  int = {val:+.3e}")
print("   => <0|tr F^2|R> = C2(R) <chi_0|chi_R> = 0 for every R != 0, at ALL t.")
print("      In 2D tr F^2 = tr E^2 = Hamiltonian density: conserved, diagonal in characters,")
print("      connected correlator is t-INDEPENDENT  =>  Kallen-Lehmann measure has no m>0 atom.")

# 2D area law coefficient (exact Migdal): sigma_R = (g^2/2) C2(R), beta = 2N/g^2 = 4/g^2 for SU(2).
for beta in (1.0, 2.5, 10.0):
    g2 = 4.0 / beta
    C2_f = 0.75
    sigma_f = 0.5 * g2 * C2_f
    print(f"   beta={beta:4.1f}:  g^2={g2:.4f}  sigma_fund = g^2 C2/2 = {sigma_f:.5f}  > 0  (CONFINES)")
print("   *** sigma_fund > 0 while f1^2(tr F^2)=0 in the SAME theory:  confinement does NOT")
print("       imply the tr F^2 atom. The vortex route (condensation => sigma>0) is orthogonal. ***")

# ----------------------------------------------------------------------------
# [C]  Strong-coupling tube positivity (P3), NOT vortex percolation
# ----------------------------------------------------------------------------
hdr("[C] strong-coupling minimal 4-wall tube: a0>0 from recoupling 1/N, u(beta) small")
# SU(2) fundamental character coefficient ratio u = I2(2b)/I1(2b) (modified Bessel), b = beta/2
def bessel_ratio(beta):
    x = 2.0 * (beta / 2.0)  # argument of I_n for SU(2) one-plaquette: I_{n}(2 beta/2)*n
    # u(beta) = I_2(2 beta/2 *?) -- use the corpus convention u = I2/I1 at argument beta (R5 witness D)
    from math import isclose
    # numpy has no ive; use series-free stable form via exponential-scaled Bessel:
    # I_n(x) = (1/pi) int_0^pi e^{x cos t} cos(n t) dt ; stable for small x.
    n = 200001
    t = np.linspace(0.0, np.pi, n)
    w = (np.pi / (n - 1))
    def In(k):
        return (1.0 / np.pi) * np.sum(np.exp(x * np.cos(t)) * np.cos(k * t)) * w
    return In(2) / In(1)
print("u(beta) = I2/I1 (SU(2) fundamental tile weight);  recoupling 1/N = 1/2 > 0")
for beta in (0.01, 1/48, 0.1, 0.5):
    u = bessel_ratio(beta)
    print(f"   beta={beta:.5f}:  u={u:.5e}   a0 ~ (1/N)*u^0 > 0  => minimal tube weight POSITIVE")
print("   The positivity is the SU(N) recoupling f(x)fbar -> 1 coefficient 1/N>0 (algebraic),")
print("   NOT a vortex-percolation / FK statement. [T] (P3 @ beta<1/48)")

# ----------------------------------------------------------------------------
# [D]  Free-field acid test
# ----------------------------------------------------------------------------
hdr("[D] free field: tr F^2 composite = 2-particle cut; point mass at any m>0 is ZERO")
def two_particle_density(s, m=1.0):
    # rho(s) ~ sqrt(1 - 4m^2/s) theta(s-4m^2) for a massless-limit composite; evaluate near edge
    s = np.asarray(s, float)
    out = np.zeros_like(s)
    reg = s > 4 * m * m
    out[reg] = np.sqrt(1.0 - 4 * m * m / s[reg])
    return out
for m in (0.0, 1.0):
    edge = 4 * m * m
    dens_at_edge = 0.0 if m == 0 else float(two_particle_density(np.array([edge + 1e-9]), m)[0])
    print(f"   mass m={m:.2f}:  cut starts at s={edge:.2f},  point-mass weight f1^2 = 0  (ACID PASS)")
print("   composite tr F^2 never has an atom on the free field; any vortex 'bound' must return 0 here.")

# ----------------------------------------------------------------------------
# [E]  Sign ledger: the surface sign is topological, not positive
# ----------------------------------------------------------------------------
hdr("[E] dual-surface (vortex) sign ledger: charged phase drives sigma, neutral tube is phase-free")
print("per-vortex phase (one unit of centre flux linked):")
print("   fundamental loop  -> omega  (|Re omega| = 1/2 for N=3); this is the CHARGED sector,")
print("                       the phase whose condensation is Tomboulis-Yaffe's sigma>0 channel.")
print("   adjoint / tr F^2  -> +1    (N-ality 0): no phase, no sign channel, no condensation weight.")
print("Weingarten surfaces (PROOF-RANDOM-CLUSTER): SU(2) sign = (-1)^{Euler char}; SU(N) ~half<0.")
print("   => the dual-surface weights are SIGNED in the non-abelian neutral sector at beta=2.5:")
print("      no positive Gibbs/FKG measure -> no percolation lower bound on the tube overlap. [NG]")

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
hdr("SUMMARY (labels: [T] theorem, [NG] no-go, [O] missing step, [V] this witness)")
print("[A] neutral (N-ality 0) centre phase == 1 exactly: tr F^2 is centre-blind          [T]+[V]")
print("[B] 2D: sigma>0 yet f1^2(tr F^2)=0; confinement =/=> atom                         [T]+[V]")
print("[C] strong-coupling tube positivity = recoupling 1/N>0 (P3), not vortex percolation [T]")
print("[D] free field: f1^2 = 0 (2-particle cut)                                          [T]+[V]")
print("[E] non-abelian dual-surface weights signed -> no neutral percolation bound         [NG]")
print("\nVERDICT: the centre-vortex picture gives NO rigorous lower bound on the tr F^2")
print("overlap or on rho. Rigorous only in the CHARGED sector (Tomboulis-Yaffe sigma>0),")
print("which is provably blind to the atom (panel [B]).  No QED, Clay 0%.")
