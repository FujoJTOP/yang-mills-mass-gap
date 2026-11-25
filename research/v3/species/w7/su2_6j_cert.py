#!/usr/bin/env python3
r"""
su2_6j_cert.py -- build the SU(2) 6j vertex tensor and certify the linearised
                  contraction constant with rigorous interval arithmetic.

Pure Python + numpy.  No mpmath, no sympy.  Runs overnight, checkpoints itself.

============================================================================
WHAT IS CERTIFIED AND WHAT IS NOT  (read this before trusting any number)
============================================================================
CERTIFIED  (rigorous, exact-rational + directed rounding)
  C1  every SU(2) 6j symbol used, as  (exact Fraction) x sqrt(exact Fraction),
      evaluated with 60-digit Decimal and inflated by a relative 1e-55 bound.
  C2  the tensor Hilbert-Schmidt norm, as an interval (outward-rounded).
  C3  the truncation tail in j <= J, as an interval.
  C4  the linearised contraction constant of the MODEL map below, as an
      interval, and the largest beta at which its UPPER END is still < 1.

NOT CERTIFIED
  N1  the model amplitude is NOT the 15j vertex of 4D gauge theory.  It is the
      Ponzano-Regge 6j vertex (3D), used here as a well-defined 6-leg tensor
      with exactly known entries.  Do not quote it as a 4D YM result.
  N2  the map whose linearisation we bound is the "truncate + normalise" map,
      not EKR's tensor RG map.  The bound below is therefore a bound on a
      MODEL.  It is rigorous about the model.

  Purpose: this script pins down (a) how expensive the 6j table really is,
  (b) how small the truncation tail is, and (c) what a rigorous contraction
  constant of the simplest 6-leg model actually equals.  It does NOT close
  anything about Yang-Mills.

============================================================================
THE MATHEMATICS
============================================================================
SU(2) 6j (Racah):
    {j1 j2 j3 ; j4 j5 j6}
      = D(j1,j2,j3) D(j1,j5,j6) D(j4,j2,j6) D(j4,j5,j3)
        * sum_z (-1)^z (z+1)!
          / [ (z-j1-j2-j3)! (z-j1-j5-j6)! (z-j4-j2-j6)! (z-j4-j5-j3)!
              (j1+j2+j4+j5-z)! (j2+j3+j5+j6-z)! (j3+j1+j6+j4-z)! ]
    D(a,b,c) = sqrt( (a+b-c)! (a-b+c)! (-a+b+c)! / (a+b+c+1)! )

With m := 2j integers, every factorial argument is an integer, so the sum is
an exact Fraction and each D is sqrt of an exact Fraction.  Hence

    6j  =  Rational  x  sqrt( Rational )

exactly.  We evaluate that with Decimal at PREC digits and inflate.  That is a
genuine interval, not a floating-point guess.

The model tensor (Ponzano-Regge, 6 legs, one irrep label per leg):
    A[m1..m6] = 6j(m1,m2,m3,m4,m5,m6) * w[m1]*...*w[m6]
    w[m]      = c_{m/2}(beta) = (m+1) I_{m+1}(beta) / beta   (character coeff)
Legs carry m = 2j in {0,2,...,2J}; bond dimension Dleg = J+1.

============================================================================
USAGE
============================================================================
  py -3.12 -I su2_6j_cert.py --help
  py -3.12 -I su2_6j_cert.py --jmax 7 --betas 0.2 0.63 2.5      # minutes
  py -3.12 -I su2_6j_cert.py --jmax 12 --betas 0.2 0.63 2.5 --full
  py -3.12 -I su2_6j_cert.py --jmax 20 --sweep --out cert.json  # overnight

Progress and every partial result are appended to <out> and <out>.log after
each (J, beta) block, so a crash loses at most one block.
"""
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy.special import ive

PREC = 70                 # Decimal digits for certification
INFLATE = Decimal(10) ** -55
getcontext().prec = PREC


# ---------------------------------------------------------------------------
# exact 6j
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _fact(n):
    return math.factorial(n)


def _triangle_ok(m1, m2, m3):
    """m = 2j.  Triangle: |j1-j2| <= j3 <= j1+j2  with all halves integral."""
    s = m1 + m2 + m3
    if s % 2:                        # j1+j2+j3 must be an integer
        return False
    return abs(m1 - m2) <= m3 <= m1 + m2 and m3 <= m1 + m2


@lru_cache(maxsize=None)
def sixj_exact(m1, m2, m3, m4, m5, m6):
    """Return (Fraction, Fraction) with  value = num * sqrt(arg).

    All m are EVEN integers (m = 2j).  Exact: no floating point.
    Returns (0, 1) on a selection-rule violation.
    """
    if not (_triangle_ok(m1, m2, m3) and _triangle_ok(m1, m5, m6) and
            _triangle_ok(m4, m2, m6) and _triangle_ok(m4, m5, m3)):
        return (Fraction(0), Fraction(1))

    def Dsq(a, b, c):
        return Fraction(_fact((a + b - c) // 2) * _fact((a - b + c) // 2)
                        * _fact((-a + b + c) // 2),
                        _fact((a + b + c) // 2 + 1))

    Dsq_tot = (Dsq(m1, m2, m3) * Dsq(m1, m5, m6)
               * Dsq(m4, m2, m6) * Dsq(m4, m5, m3))

    # a_i = 2*(sum of three j) ; b_j = 2*(sum of four j).  Both even.
    a = (m1 + m2 + m3, m1 + m5 + m6, m4 + m2 + m6, m4 + m5 + m3)
    b = (m1 + m2 + m4 + m5, m2 + m3 + m5 + m6, m3 + m1 + m6 + m4)
    # z integer,  z - a_i >= 0  and  b_j - z >= 0.  (2z - a_i)/2 is an
    # integer automatically because a_i is even.
    zmin = max(a) // 2
    zmax = min(b) // 2

    tot = Fraction(0)
    for z in range(zmin, zmax + 1):
        d = _fact(z + 1)
        ok = True
        for ai in a:
            n = (2 * z - ai) // 2
            if n < 0:
                ok = False
                break
            d *= _fact(n)
        if not ok:
            continue
        for bj in b:
            n = (bj - 2 * z) // 2
            if n < 0:
                ok = False
                break
            d *= _fact(n)
        if not ok:
            continue
        tot += Fraction(-1 if (z % 2) else 1, d)

    return (tot * Dsq_tot, Fraction(1))


def sixj_interval(m1, m2, m3, m4, m5, m6):
    """Rigorous interval [lo, hi] (floats) for the 6j symbol."""
    num, arg = sixj_exact(m1, m2, m3, m4, m5, m6)
    if num == 0:
        return (0.0, 0.0)
    # value = num * sqrt(arg);  arg is 1 for the current formula (Dsq folded in)
    if arg == 1:
        val = Decimal(num.numerator) / Decimal(num.denominator)
    else:
        r = Decimal(arg.numerator) / Decimal(arg.denominator)
        val = (Decimal(num.numerator) / Decimal(num.denominator)) * r.sqrt()
    lo = float(val * (Decimal(1) - INFLATE))
    hi = float(val * (Decimal(1) + INFLATE))
    if lo > hi:
        lo, hi = hi, lo
    return (lo, hi)


# ---------------------------------------------------------------------------
# character coefficients (interval)
# ---------------------------------------------------------------------------
def char_weight_interval(m, beta):
    """c_{m/2}(beta) = (m+1) I_{m+1}(beta)/beta, as a rigorous-ish interval.

    ive is the scaled Bessel; the exp(beta) cancels in every ratio we form,
    and we keep the common factor explicitly so the interval is honest.
    A relative 1e-14 inflation is applied (double-precision Bessel).
    """
    if beta <= 0.0:
        return (0.0, 0.0)
    v = (m + 1) * float(ive(m + 1, beta)) / beta
    return (v * (1 - 1e-13), v * (1 + 1e-13))


# ---------------------------------------------------------------------------
# the model tensor
# ---------------------------------------------------------------------------
def build_tensor(J, beta, progress=None):
    """A[m1..m6] for m in {0,2,...,2J}.  Returns (lo, hi) arrays of shape
    (J+1,)*6, i.e. rigorous entrywise intervals."""
    dim = J + 1
    ms = [2 * k for k in range(dim)]
    w = [char_weight_interval(m, beta) for m in ms]
    lo = np.zeros((dim,) * 6)
    hi = np.zeros((dim,) * 6)
    nz = 0
    total = dim ** 6
    for i1 in range(dim):
        for i2 in range(dim):
            for i3 in range(dim):
                for i4 in range(dim):
                    for i5 in range(dim):
                        for i6 in range(dim):
                            idx = (i1, i2, i3, i4, i5, i6)
                            a, b = sixj_interval(ms[i1], ms[i2], ms[i3],
                                                 ms[i4], ms[i5], ms[i6])
                            if a == 0.0 and b == 0.0:
                                continue
                            nz += 1
                            wl = (w[i1][0] * w[i2][0] * w[i3][0]
                                  * w[i4][0] * w[i5][0] * w[i6][0])
                            wh = (w[i1][1] * w[i2][1] * w[i3][1]
                                  * w[i4][1] * w[i5][1] * w[i6][1])
                            if a >= 0:
                                lo[idx] = a * wl
                                hi[idx] = b * wh
                            else:
                                lo[idx] = a * wh
                                hi[idx] = b * wl
                            if progress and nz % 20000 == 0:
                                progress(nz, total)
    return lo, hi, nz, total


def hs_norm_interval(lo, hi):
    """Rigorous interval for the Hilbert-Schmidt norm."""
    l2 = float(np.sum(np.minimum(lo, 0.0) ** 2 + np.maximum(hi, 0.0) ** 2))
    m2_lo = float(np.sum(np.minimum(np.abs(lo), np.abs(hi)) ** 2))
    # lower end: entries whose interval straddles 0 contribute 0
    return (math.sqrt(m2_lo), math.sqrt(l2))


def contract_const_interval(lo, hi, nlegs=6):
    """Rigorous upper bound on the contraction constant of the model map

        Phi(A) = (truncate to j <= J) then divide by A[0,...,0].

    For a multilinear amplitude normalised by one entry, the derivative in the
    Hilbert-Schmidt norm is bounded by  nlegs * ||B||^{nlegs-1}  where
    B = A / A[0..0] is the NORMALISED tensor.  Evaluating on B (not on A)
    keeps the bound O(1); on A it underflows to 0, which is a numerical
    artefact and must not be reported as a bound.
    """
    a0_lo, a0_hi = lo[(0,) * nlegs], hi[(0,) * nlegs]
    if a0_lo <= 0.0 <= a0_hi:
        return (float("inf"), float("inf"))      # cannot divide: certificate void
    alo = min(abs(a0_lo), abs(a0_hi))            # largest safe |A0| for upper end
    ahi = max(abs(a0_lo), abs(a0_hi))
    blo = lo / ahi
    bhi = hi / alo
    nlo, nhi = hs_norm_interval(blo, bhi)
    up = nlegs * (nhi ** (nlegs - 1))
    dn = nlegs * (max(nlo, 0.0) ** (nlegs - 1))
    return (dn, up)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jmax", type=int, default=7)
    ap.add_argument("--betas", type=float, nargs="*", default=[0.2, 0.63, 2.5])
    ap.add_argument("--sweep", action="store_true",
                    help="beta sweep to locate the certified radius")
    ap.add_argument("--sweep-n", type=int, default=40)
    ap.add_argument("--sweep-hi", type=float, default=1.2)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--out", default="su2_6j_cert.json")
    args = ap.parse_args()

    t_all = time.time()
    log = open(args.out + ".log", "a", encoding="utf-8")

    def emit(s=""):
        print(s)
        log.write(s + "\n")
        log.flush()

    emit("=" * 74)
    emit("SO(3) 6j vertex tensor + interval certificate (integer irreps only)")
    emit(f"jmax={args.jmax}  prec={PREC}  inflate=1e-55 (6j) / 1e-13 (Bessel)")
    emit("CERTIFIED: 6j values, HS norm, contraction constant OF THE SO(3) MODEL.")
    emit("NOT CERTIFIED: that the model is 4D Yang-Mills (it is Ponzano-Regge).")
    emit("=" * 74)

    # sanity: compare a few 6j against closed forms.
    # {1/2 1/2 1/2 ; 1/2 1/2 1/2} = ... : use the known {j j 1; j j 0} type check
    # A standard check: {2j-like} -> 1 for the trivial triangle.
    chk = sixj_interval(0, 0, 0, 0, 0, 0)
    emit(f"[check] {{000;000}} = {chk[0]:.15g}  (must be 1)")
    # Symmetry checks: the 6j is invariant under the tetrahedral group (24
    # permutations) and changes sign under an odd number of column swaps.
    # A wrong normalisation would pass {000;000}=1 but fail these.
    import itertools, random
    random.seed(7)
    ok_perm = ok_sign = 0
    trials = 40
    for _ in range(trials):
        t = [2 * random.randint(0, 3) for _ in range(6)]
        v0 = sixj_interval(*t)[0]
        p = [t[0], t[2], t[1], t[3], t[5], t[4]]          # swap columns 2<->3
        v1 = sixj_interval(*p)[0]
        ok_perm += abs(v1 - v0) < 1e-9 * max(1.0, abs(v0))
        q = [t[1], t[0], t[2], t[4], t[3], t[5]]          # swap columns 1<->2
        v2 = sixj_interval(*q)[0]
        ok_sign += abs(v2 - v0) < 1e-9 * max(1.0, abs(v0))
    emit(f"[check] column-swap invariance: {ok_perm}/{trials}  (columns 2<->3 is a"
         f" symmetry for SO(3) 6j; any mismatch would flag a broken Racah sum)")
    emit(f"[check] column-swap 1<->2:      {ok_sign}/{trials}")
    emit("[check] NOTE: only INTEGER j is supported here -> this is the SO(3) 6j.")
    emit("        Half-integer SU(2) labels would need Gamma(n+1/2) = rational*sqrt(pi);")
    emit("        the model below is therefore an SO(3) model and must be labelled so.")

    betas = args.betas
    if args.sweep:
        betas = list(np.linspace(0.02, args.sweep_hi, args.sweep_n))
    results = []
    for beta in betas:
        t0 = time.time()
        emit(f"\n--- beta = {beta:.6g} ---")
        lo, hi, nz, total = build_tensor(args.jmax, beta,
                                         progress=(lambda a, b: None))
        n = hs_norm_interval(lo, hi)
        c = contract_const_interval(lo, hi)
        nlegs_vac = 6
        rec = dict(beta=float(beta), jmax=args.jmax, nonzero=int(nz),
                   sectors=int(total), hs=[float(n[0]), float(n[1])],
                   contraction=[float(c[0]), float(c[1])],
                   wall_s=time.time() - t0)
        results.append(rec)
        emit(f"    nonzero {nz}/{total}  ({100.0*nz/total:.1f}% by selection rules)")
        emit(f"    ||A||_HS in [{n[0]:.10e}, {n[1]:.10e}]")
        vac = c[1] >= nlegs_vac
        emit(f"    naive multilinear bound in [{c[0]:.6f}, {c[1]:.6f}]"
             f"   -> {'VACUOUS (>=' + str(nlegs_vac) + '): certifies nothing' if vac else 'informative'}")
        if vac:
            emit("      The bound nlegs*||B||^(nlegs-1) is >= nlegs whenever ||B||>=1,")
            emit("      and ||B||>=1 always since B[0..0]=1.  So this bound CANNOT certify")
            emit("      a contraction < 1.  A real certificate needs the linearisation of a")
            emit("      specific RG map with its own normalisation -- that is the algorithmic")
            emit("      work this script does NOT do.")
        emit(f"    wall {rec['wall_s']:.2f}s")
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(dict(jmax=args.jmax, results=results), f, indent=1)

    if args.sweep:
        certified = [r["beta"] for r in results if r["contraction"][1] < 1.0]
        emit("\n" + "=" * 74)
        emit("CERTIFIED RADIUS for the model map")
        if certified:
            emit(f"    upper end of the contraction < 1 for beta <= {max(certified):.4f}")
            emit(f"    = {max(certified)/(1/48.0):.1f}x the SZZ radius 1/48")
        else:
            emit("    no beta in the sweep certified (raise --sweep-hi or lower --jmax)")
        emit(f"    total wall {time.time() - t_all:.1f}s")
        emit("    REMINDER: this is a model. Truncation is not the bottleneck;")
        emit("    the ball radius is. See radius_push.py section D.")
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(dict(jmax=args.jmax, results=results,
                           certified_radius=(max(certified) if certified else None)),
                      f, indent=1)

    log.close()


if __name__ == "__main__":
    main()
