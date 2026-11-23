"""
sdpb_toy.py -- de-risk the conformal-bootstrap SDP machinery on a KNOWN 1D toy.

VALIDATES PLUMBING ONLY: block vectors, crossing, the positivity (nonneg)
cone, the per-operator Z_N phase rotation, and the recovery of a known gap
bound.  It does NOT validate any 4D Yang-Mills block or the defect setup.

Toy: 1D CFT of one real scalar phi (dimension DPHI), Z_2 symmetric.
Known answer: the Generalized Free Field (GFF / mean-field theory),
   lightest exchanged operator at Delta_min = 2*DPHI, OPE data a_n^2 > 0.

Conventions (each cross-checked numerically against the GFF correlator):
  bare 1D block    g_d(u)  = u^d * 2F1(d, d; 2d; u)             [= NUMERICS.md]
  full block       G_d(u)  = u^(d-2 DPHI) * 2F1(d, d; 2d; u)    [external factor * bare]
  identity         G_0(u)  = u^(-2 DPHI)                        [d=0 limit]
  GFF correlator   G(u)    = 1 + u^(-2 DPHI) + (1-u)^(-2 DPHI)  [crossing symmetric]
  OPE / crossing:  G(u) = sum_d a_d^2 G_d(u),  a_d^2 >= 0,  and G(u)=G(1-u).
  crossing vector  R_d(u)  = G_d(u) - G_d(1-u)   (odd about u=1/2; used for phases)

Pipeline / checks:
 [1] blocks + crossing : max|G(u)-G(1-u)| ~ 0
 [2] positivity cone   : fit GFF with a_n^2 >= 0 (NNLS) -> a_n^2 > 0, resid ~ 0
 [3] known gap bound   : scan Delta_gap in the ansatz {Delta_gap + 2n};
                         residual has a sharp zero at Delta_gap = 2*DPHI.
                         Dropping the lightest operator -> residual -> O(1).
 [4] Z_N phase rotation: rotate block columns by omega^(q_n q_test), |omega|=1
                         (ray-positivity cone unchanged); trivial sector real,
                         "untwisted one-form" sector gains Im(cross)=0 equality.

Honest boundary: a naive DUAL-functional LP (Lambda(R_0)=1, Lambda(R_Delta)>=0)
was also tried; its threshold is grid/dof dependent (2.27 with one grid, none
with a finer one) because a point-evaluation LP is NOT a faithful SOS/PSD
positivity test.  The faithful PSD/SOS layer of real SDPB is absent here, so the
quantity we report as the "bound" is the PRIMAL known-answer recovery [3], not a
trustworthy dual bound.  This is the correct scope of what scipy-only de-risks.
"""

import numpy as np
from scipy.special import hyp2f1
from scipy.optimize import nnls

DPHI = 0.5                                   # external scalar dimension


# ----------------------------------------------------------------- blocks
def full_block(delta, u):
    """Full 1D conformal block G_delta(u) (external factor absorbed)."""
    u = np.asarray(u, dtype=float)
    if abs(delta) < 1e-12:                   # identity
        return u ** (-2.0 * DPHI)
    return u ** (delta - 2.0 * DPHI) * hyp2f1(delta, delta, 2.0 * delta, u)


def gff(u):
    """GFF reduced correlator, crossing symmetric by construction."""
    u = np.asarray(u, dtype=float)
    return 1.0 + u ** (-2.0 * DPHI) + (1.0 - u) ** (-2.0 * DPHI)


def cross_vec(delta, u):
    """Crossing vector R_delta(u) = G_delta(u) - G_delta(1-u)."""
    return full_block(delta, u) - full_block(delta, 1.0 - u)


# ------------------------------------------- [1] blocks + crossing check
def check_crossing(nu=200):
    u = np.linspace(0.02, 0.98, nu)
    return float(np.max(np.abs(gff(u) - gff(1.0 - u))))


# ---------------------- [2] positivity cone: recover the known GFF data
def fit_nonneg(deltas, u1=0.8, nu=60):
    """NNLS fit G(u) = 1*G_0(u) + sum_n a_n^2 G_{Delta_n}(u), enforcing a_n^2>=0."""
    ug = np.linspace(0.02, u1, nu)
    A = np.array([full_block(d, ug) for d in deltas]).T          # (nu, N)
    b = gff(ug) - full_block(0.0, ug)                           # identity coeff = 1
    scale = np.linalg.norm(A, axis=0)                           # condition columns
    x, _ = nnls(A / scale, b)
    a2 = x / scale
    resid = float(np.linalg.norm(A @ a2 - b) / np.linalg.norm(b))
    return resid, a2


# --------------------------- [3] known exclusion bound (Delta_min = 2*DPHI)
def gap_scan(gaps, u1=0.8, nu=60, N=25):
    """Residual of the nonneg fit over ansatz {gap + 2n}; min marks the gap."""
    out = []
    for g in gaps:
        r, _ = fit_nonneg([g + 2 * n for n in range(N)], u1=u1, nu=nu)
        out.append(r)
    return np.array(out)


# ---------------------------- [4] Z_N per-operator phase rotation
def phase_rotation_test(a2, N=3, k=1, qtest=1, nu=60):
    """Rotate block column n by omega^(q_n qtest), omega = exp(2 pi i k / N).
         trivial sector  q_n = 0  -> real crossing (standard SDPB)
         one-form sector q_n = n  -> complex crossing; Im(cross)=0 is ONE extra
                                     linear equality (Im term is generally !=0)."""
    ug = np.linspace(0.02, 0.98, nu)
    deltas = np.array([2.0 * DPHI + 2 * n for n in range(len(a2))])
    R = np.array([cross_vec(d, ug) for d in deltas])            # (N, nu)
    omega = np.exp(2j * np.pi * k / N)

    def max_im(q):
        phases = omega ** (q * qtest)                           # per-operator phase
        S = np.sum(a2[:, None] * phases[:, None] * R, axis=0)   # complex crossing
        return float(np.max(np.abs(S.imag))), float(np.max(np.abs(S.real)))

    q_triv = np.zeros(len(a2), dtype=int)
    q_one = np.arange(len(a2)) % N
    im0, re0 = max_im(q_triv)
    im1, re1 = max_im(q_one)
    return abs(omega), im0, re0, im1, re1


# ------------------------------------------------------------------- main
if __name__ == "__main__":
    print("=" * 66)
    print("TOY 1D conformal bootstrap  DPHI=%.2f  -- SDP-plumbing de-risk" % DPHI)
    print("=" * 66)

    # [1]
    print("\n[1] blocks + crossing : max|GFF(u)-GFF(1-u)| = %.2e" % check_crossing())

    # [2]
    resid, a2 = fit_nonneg([2 * DPHI + 2 * n for n in range(20)])
    print("[2] positivity cone   : GFF nonneg fit resid = %.2e" % resid)
    print("    a_n^2[0:4] = %s   (all >= 0 -> GFF is in the cone)"
          % np.array2string(a2[:4], precision=4))
    print("    min a_n^2  = %.2e" % float(a2.min()))

    # [3] known bound: scan Delta_gap
    gaps = np.round(np.arange(0.6, 1.45, 0.05), 3)
    rs = gap_scan(gaps)
    gbest = float(gaps[int(np.argmin(rs))])
    print("[3] known gap bound   : best-fit Delta_gap = %.3f  (resid %.1e);  "
          "known 2*DPHI = %.1f" % (gbest, rs.min(), 2 * DPHI))
    r_drop, _ = fit_nonneg([2 * DPHI + 2 + 2 * n for n in range(20)])  # drop lightest
    print("    drop lightest op (Delta_min=%.1f): resid = %.2e  -> excluded"
          % (2 * DPHI + 2, r_drop))

    # [4] phase rotation
    mod, im0, re0, im1, re1 = phase_rotation_test(a2)
    print("[4] Z_3 phase rotation: |omega| = %.15f  (unit modulus => cone unchanged)" % mod)
    print("    trivial  sector q=0 : max|Im cross| = %.1e  (real; standard SDPB)" % im0)
    print("    one-form sector q=n : max|Im cross| = %.2e  (nonzero => +1 lin. eq.)"
          % im1)
    print("    (Re parts: trivial %.2e, one-form %.2e)" % (re0, re1))
    print("=" * 66)
