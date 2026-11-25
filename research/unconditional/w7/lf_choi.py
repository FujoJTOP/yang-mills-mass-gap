# -*- coding: utf-8 -*-
"""
lf_choi.py -- witness for candidate (2) of the death-cause axis: operator-valued K-theory /
complete positivity / Choi rank.

Question put to the numbers:
  (Q1) What is the Choi rank of the ONLY canonical CP map the theory supplies
       (conjugation by the transfer operator, X -> T X T*) ?
  (Q2) Is any real-valued functional of the Choi operator an IR datum?  Which end of the
       spectrum does it read?
  (Q3) Can a *non-trivial* Choi rank be produced canonically (i.e. as a function of (beta,L,a))?

Conventions: Choi matrix C_Phi = sum_{ij} E_ij (x) Phi(E_ij), block (i,j) = Phi(E_ij).
Choi's theorem: Phi CP <=> C_Phi >= 0.  Choi rank = rank C_Phi = minimal Kraus number.
"""
import numpy as np
from scipy.linalg import expm

np.random.seed(20261125)
OUT = []


def say(s=""):
    print(s)
    OUT.append(str(s))


def choi(Phi, d):
    C = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            E = np.zeros((d, d), dtype=complex)
            E[i, j] = 1.0
            C[i * d:(i + 1) * d, j * d:(j + 1) * d] = Phi(E)
    return C


def crank(C, tol=1e-9):
    return int(np.linalg.matrix_rank(C, tol=tol))


def ad(V):
    return lambda X: V @ X @ V.conj().T


def dephase(lam):
    return lambda X: (1.0 - lam) * X + lam * np.diag(np.diag(X))


def rand_unitary(d):
    A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


# ================================================================ W1  Ad(V): rank 1
say("=" * 78)
say("W1  the canonical CP map:  Phi(X) = V X V*   (V = transfer operator e^{-aH})")
say("    claim: C_Phi = vec(V) vec(V)^*  =>  Choi rank = 1, for every V, every dim")
say("=" * 78)
say(f"{'d':>4} {'rank C_Phi':>11} {'nonzero eig(C)':>18} {'||V||_F^2':>14} {'eig count>tol':>14}")
for d in (2, 3, 4, 8, 16):
    V = rand_unitary(d) @ np.diag(np.exp(-np.linspace(0, 3, d)))
    C = choi(ad(V), d)
    ev = np.sort(np.linalg.eigvalsh(C))[::-1]
    nz = ev[np.abs(ev) > 1e-9]
    say(f"{d:4d} {crank(C):11d} {len(nz):18d} {np.linalg.norm(V,'fro')**2:14.8f} "
        f"{int(np.sum(np.abs(ev)>1e-9)):14d}")
say(f"  -> rank = 1 in every dimension; the Choi operator is a rank-one projector.")
say(f"  -> its ENTIRE spectrum is {{||V||_F^2, 0, ..., 0}}: there is no second eigenvalue,")
say(f"     hence no place for a 'gap' to live.  (Choi's theorem is satisfied trivially.)")

# ================================================================ W2  dephasing
say()
say("=" * 78)
say("W2  the only way to get rank > 1: add a non-unitary ingredient")
say("    Phi_lam(X) = (1-lam) X + lam diag(X)     (Kraus {sqrt(1-lam) I, sqrt(lam)|i><i|})")
say("=" * 78)
d = 4
say(f"  d={d}")
say(f"{'lam':>10} {'rank C':>8}   note")
for lam in (0.0, 1e-12, 1e-6, 0.01, 0.5, 1.0):
    r = crank(choi(dephase(lam), d))
    say(f"{lam:10.0e} {r:8d}   {'rank-1 (conjugation)' if r==1 else 'full rank -> the step'}")
say("  -> a STEP function of lam at lam=0: no threshold-crossing real value, nothing to read")
say("     as '>= c > 0'.  And lam is a MODELLING choice, not a function of (beta,L,a):")
H = np.diag([0.0, 0.7, 1.9, 3.1])
V = expm(-H)
say(f"     same H, same transfer, dephasing 0 vs 0.5:  rank {crank(choi(ad(V),4))} vs "
    f"{crank(choi(dephase(0.5),4))}  -- the dynamics did not change, the rank did.")

# ================================================================ W3  gap closing
say()
say("=" * 78)
say("W3  ACID TEST in numbers: H_beta = 1.5 I + [[0,beta],[beta,0]],  spec = 1.5 +- beta")
say("    the gap 2|1.5-beta| closes at beta=1.5.  Does any Choi datum move at the closure?")
say("=" * 78)
say(f"{'beta':>7} {'gap':>9} {'rank C':>8} {'||V||_F^2 = Z(2a)':>20} {'eig(C) 2nd':>12}")
for beta in (0.0, 0.5, 1.0, 1.4, 1.5, 1.6):
    H = 1.5 * np.eye(2) + np.array([[0.0, beta], [beta, 0.0]])
    w = np.linalg.eigvalsh(H)
    gap = 2.0 * abs(1.5 - beta)
    V = expm(-1.0 * H)
    C = choi(ad(V), 2)
    ev = np.sort(np.linalg.eigvalsh(C).real)[::-1]
    say(f"{beta:7.2f} {gap:9.3f} {crank(C):8d} {np.linalg.norm(V,'fro')**2:20.8f} {ev[1]:12.3e}")
say("  -> the rank is identically 1 ACROSS the closure (identical to the 'gapless' case).")
say("     The one nonzero Choi eigenvalue is ||e^{-aH}||_F^2 = sum_n e^{-2aE_n} = Z(2a):")
say("     the partition function at lattice spacing 2a.  It is extensive/UV-dominated and")
say("     it is SMOOTH through the gap closure  =>  the tool cannot see the IR edge.")
say("     (corpus: PROOF-THERMODYNAMIC-FORMALISM -- pressure = dominant eigenvalue; the gap")
say("      is the SECOND eigenvalue, and no bulk thermodynamic quantity carries it.)")

# ================================================================ W4  what the scalar is
say()
say("=" * 78)
say("W4  is the surviving scalar IR or UV?   decomposing Z(2a)")
say("=" * 78)
spec = np.array([0.0, 0.6444, 1.2, 1.9, 2.7, 3.6, 4.6, 5.7])
for a in (0.25, 0.5, 1.0):
    Z = np.sum(np.exp(-2.0 * a * spec))
    vac_share = np.exp(-2.0 * a * spec[0]) / Z
    uv_share = np.exp(-2.0 * a * spec[-1]) / Z
    say(f"  a={a:4.2f}: Z(2a)={Z:12.8f}   vacuum share={vac_share:.10f}   "
        f"top-level share={uv_share:.3e}")
say("  -> Z(2a) is dominated by the vacuum; its VARIATION with the spectrum is")
say("     sum_n e^{-2aE_n} -- a sum over the WHOLE spectrum, i.e. an extensive quantity")
say("     that changes smoothly with every E_n.  A change of the gap alone moves it by")
say(f"     the amount of the level itself (e^{{-2a*1}} = {np.exp(-2*0.6444):.6f} at a=1),")
say("     indistinguishable from a change of any other level: no threshold, no direction.")

# ================================================================ W5  audit
say()
say("=" * 78)
say("W5  audit of the class prediction made in MG-META (D4 = channel/spin mismatch)")
say("=" * 78)
say("  D4 demands: the tool reads the wrong channel (line/flux/spin-1/N-ality-blind).")
say("  But CP/Choi structure is attached to the ALGEBRA, not to any channel: it is the same")
say("  for every operator in the 0^{++} sector and for every other sector.  A statistic that")
say("  is constant across channels cannot be 'in the wrong channel'.")
say("  => D4 does NOT fire.  The class prediction is wrong.")

say()
say("=" * 78)
say("SUMMARY (candidate 2 / Choi rank)")
say("=" * 78)
say(" - Choi rank of the canonical map (conjugation by T) is 1 identically (W1): a class")
say("   constant, same on the free massless field -- zero discriminating power.")
say(" - the Choi operator of Ad(V) is a RANK-ONE projector: spectrum {||V||_F^2,0,...,0}")
say("   => there is no second eigenvalue, so no coordinate in which a gap could appear.")
say(" - rank > 1 needs a non-unitary ingredient whose strength lam is a modelling choice,")
say("   and the rank is then a STEP function of lam (W2) -- an integer with no threshold value.")
say(" - the only surviving real scalar is ||e^{-aH}||_F^2 = Z(2a), extensive & UV-dominated,")
say("   smooth through the gap closure (W3,W4) -- the corpus's dead thermodynamic reading.")
say(" - verdict:  D16 (integer/rigid invariant vs real infimum) with D10 + D15 + D2 as")
say("   companions.  NO non-circular Choi-rank-type criterion can be constructed.")
