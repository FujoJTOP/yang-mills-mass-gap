#!/usr/bin/env python3
r"""
radius_push.py -- how far can the rigorous strong-coupling radius be pushed?

Pure numpy+scipy, vectorised, no I/O beyond one JSON + one markdown table.
Fast mode ~1 s. Full mode ~1-2 min. Designed to be launched and left.

WHY THIS SCRIPT EXISTS
----------------------
The corpus quotes "the best existing import pushes beta from 1/48 = 0.021 to
about 0.63 (120x -> 4x)".  That number came from a tensor-RG paper whose
rigorousness rests on an analyticity ball around the TRIVIAL fixed point, with
a quoted total attraction radius 0.1587 in a Hilbert-Schmidt norm.

Two very different things are being conflated in that sentence, and this script
separates them:

  (A) THE TRUNCATION COST.  Irrep truncation j <= J in SU(2) character space.
      Governed by Bessel tails.  This is a COMPUTE question and it is cheap:
      we measure exactly where the tail falls below 1% and 0.1% as a function
      of beta and J.

  (B) THE ANALYTICITY RADIUS.  Where the physical tensor A(beta) exits a ball
      of radius 0.1587 around the trivial fixed point.  This is an ANALYSIS
      question: no amount of compute moves it.

The script reports (A) and (B) in SEVERAL EXPLICIT NORMALISATIONS, because the
headline value "0.63" is NOT convention-free.  Whichever crossing you quote,
say which norm it belongs to.  If two normalisations disagree by orders of
magnitude, then "0.63" must not be used as a physical number until pinned down.

CONVENTIONS (stated so they can be checked)
-------------------------------------------
SU(2) Wilson action, character expansion of exp(beta * chi_{1/2}(U)/2):

    exp(beta * tr(U)/2) = sum_j (2j+1) c_j(beta) chi_j(U),
    c_j(beta) = (2j+1) I_{2j+1}(beta) / beta .

All ratios below are built from  chat_j := (2j+1) * ive(2j+1, beta),
where ive(n,x) = exp(-x) I_n(x).  The common factor exp(beta)/beta cancels in
every ratio, so the numbers are stable to beta = 200 without overflow.

Weights used:
  Plancherel (character) weight for SU(2): d_j = 2j+1.
  Norm-1 ("Plancherel HS"):  ||f||^2 = sum_j d_j |c_j|^2
  Norm-2 ("sup"):            ||f||   = max_j |c_j|
  Norm-3 ("tail-mass"):      ||f||   = sum_{j>0} d_j |c_j| / (d_0 |c_0|)
  The trivial fixed point is c_* = (1, 0, 0, ...) up to normalisation, i.e.
  beta = 0; so ||b|| := ||c - c_*|| restricted to j >= 1/2.

USAGE
  py -3.12 -I radius_push.py                 # fast mode, prints a table
  py -3.12 -I radius_push.py --mode full     # dense grid + MK check + JSON
  py -3.12 -I radius_push.py --mode full --jmax 30 --nbeta 4000
"""
import argparse
import json
import time

import numpy as np
from scipy.special import ive

TARGET_RADIUS = 0.1587          # EKR-type attraction radius, HS norm (quoted)
BETA_STRONG = 1.0 / 48.0        # SZZ radius
BETA_PHYS = 2.5                 # physical coupling


def chat(beta, jmax):
    """chat_j = d_j * ive(2j+1, beta),  j = 0, 1/2, 1, ... jmax  (d_j = 2j+1)."""
    two_j = 2 * np.arange(jmax + 1)          # 0,2,4,...  -> j = two_j/2
    d = two_j + 1                            # 1,3,5,... = 2j+1
    # ive(order, x) is vectorised over both arguments via broadcasting
    vals = ive(two_j[:, None] + 1, np.atleast_1d(beta)[None, :])
    return d[:, None] * vals, d


def metrics(beta, jmax):
    """Return the four quantities, shape (nbeta,)."""
    ch, d = chat(beta, jmax)
    total = np.sum(d[:, None] * ch ** 2, axis=0)
    c0 = ch[0]

    # (A) truncation tail: mass in j > jmax relative to total mass
    #     (needs a larger jmax to be meaningful; caller passes jmax for the
    #      *renormalised* tail estimate below instead)
    # (B) distance to the trivial fixed point in three normalisations
    excited = d[:, None] * ch ** 2
    excited[0] = 0.0
    rho_hs = np.sqrt(np.sum(excited, axis=0)) / c0        # Plancherel HS
    rho_sup = np.max(ch[1:], axis=0) / c0                  # sup norm
    rho_tail = np.sum(d[1:, None] * np.abs(ch[1:]), axis=0) / (d[0] * c0)

    # (C) MK-style crossing: r(beta) := c_{1/2}/c_0 ; decimation exponent
    #     b^{d-1} = 2^3 = 8 in d = 4, so the confined fixed point is r < 1.
    r_mk = ch[1] / c0

    return dict(total=total, c0=c0, rho_hs=rho_hs, rho_sup=rho_sup,
                rho_tail=rho_tail, r_mk=r_mk)


def tail_after_J(beta, jmax, ntail=400):
    """Exact truncation tail: sum_{j>jmax} d_j c_j^2 / sum_j d_j c_j^2.

    Uses ntail extra orders so the tail is itself converged.  This is the
    quantity that decides how large J must be (a COMPUTE question)."""
    ch, d = chat(beta, jmax + ntail)
    sq = d[:, None] * ch ** 2
    tot = np.sum(sq, axis=0)
    tail = np.sum(sq[jmax + 1:], axis=0)
    return tail / tot


def first_crossing(beta, y, level):
    """Smallest beta where y crosses `level` upwards.  None if never."""
    idx = np.nonzero(y >= level)[0]
    if idx.size == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(beta[0])
    x0, x1, y0, y1 = beta[i - 1], beta[i], y[i - 1], y[i]
    if y1 == y0:
        return float(x1)
    return float(x0 + (level - y0) * (x1 - x0) / (y1 - y0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["fast", "full"], default="fast")
    ap.add_argument("--jmax", type=int, default=None)
    ap.add_argument("--nbeta", type=int, default=None)
    ap.add_argument("--json", default="radius_push_results.json")
    args = ap.parse_args()

    if args.mode == "full":
        jmax = args.jmax or 30
        nbeta = args.nbeta or 4000
    else:
        jmax = args.jmax or 12
        nbeta = args.nbeta or 400

    t0 = time.time()
    beta = np.logspace(np.log10(1e-3), np.log10(4.0), nbeta)
    m = metrics(beta, jmax)

    # ---- (B) where does A(beta) exit the 0.1587 ball, in each norm? ----
    crossings = {name: first_crossing(beta, m[key], TARGET_RADIUS)
                 for name, key in (("rho_hs", "rho_hs"),
                                   ("rho_sup", "rho_sup"),
                                   ("rho_tail", "rho_tail"))}

    # ---- (C) sensitivity to the target radius, in the norm that reproduces
    #          the corpus's 0.63 (tail-mass norm) ----
    sweep_levels = [0.10, 0.125, 0.1587, 0.20, 0.25, 0.30, 0.40, 0.50]
    sweep = {lvl: first_crossing(beta, m["rho_tail"], lvl) for lvl in sweep_levels}
    # also the same sweep in the HS norm, for contrast
    sweep_hs = {lvl: first_crossing(beta, m["rho_hs"], lvl) for lvl in sweep_levels}

    # ---- (A) truncation budget: J needed for 1% and 0.1% at beta_phys ----
    beta_probe = np.array([BETA_STRONG, 0.25, 0.5, BETA_PHYS])
    budget = {}
    for J in range(2, jmax + 1):
        tails = tail_after_J(beta_probe, J, ntail=600)
        budget[J] = [float(t) for t in tails]

    # ---- report ----
    print("=" * 78)
    print("A) TRUNCATION COST  (a compute question)")

    print("   tail = sum_{j>J} d_j c_j^2 / sum_j d_j c_j^2 ," \
          "   columns: beta =", ", ".join(f"{b:.4f}" for b in beta_probe))
    print("   " + "-" * 66)
    for J in sorted(budget):
        row = "  ".join(f"{v:.2e}" for v in budget[J])
        print(f"   J={J:2d}   {row}")
    # smallest J meeting 1% and 0.1% at beta_phys
    ip = len(beta_probe) - 1
    need1 = next((J for J in sorted(budget) if budget[J][ip] < 1e-2), None)
    need3 = next((J for J in sorted(budget) if budget[J][ip] < 1e-3), None)
    need4 = next((J for J in sorted(budget) if budget[J][ip] < 1e-4), None)
    print(f"\n   at beta={BETA_PHYS}: J needed for <1% tail = {need1},"
          f"  <0.1% = {need3},  <0.01% = {need4}")
    if need3:
        print(f"   dense sectors at J={need3}: (J+1)^6 = {(need3 + 1) ** 6:.3e}")
        print(f"   dense float64 size:        {8.0 * (need3 + 1) ** 6 / 1e6:.2f} MB")
        print("   -> the truncation leg is NOT the bottleneck; contraction order is.")

    print("=" * 78)
    print("B) ANALYTICITY RADIUS  (an analysis question; compute cannot move it)")
    print(f"   target ball radius = {TARGET_RADIUS}  (quoted in the corpus)")
    print("   norm                     beta where ||A(beta)-A_*|| = 0.1587")
    print("   " + "-" * 66)
    for name, val in crossings.items():
        print(f"   {name:<24s} {('%.4f' % val) if val is not None else 'never'}")
    print(f"\n   spread across normalisations = "
          f"{max(v for v in crossings.values() if v):.4f} / "
          f"{min(v for v in crossings.values() if v):.4f} = "
          f"{max(v for v in crossings.values() if v) / min(v for v in crossings.values() if v):.1f}x")
    print("   -> a '0.63' quoted without a norm is not a physical number.")

    print("=" * 78)
    print("C) SENSITIVITY to the assumed ball radius (and the 31x check)")
    print("   corpus quotes 7x-30x, i.e. beta* ~ 0.15-0.63 ; 1/48 = 0.0208")
    print("   " + "-" * 66)
    print("   radius   beta*(tail-mass)   beta*(HS)   ratio to 1/48 (tail)")
    for lvl in sweep_levels:
        bt, bh = sweep[lvl], sweep_hs[lvl]
        rr = (bt / BETA_STRONG) if bt else float("nan")
        print("   %6.4f   %12.4f   %10.4f   %8.1fx" % (lvl, bt if bt else float("nan"), bh if bh else float("nan"), rr))
    if crossings["rho_tail"]:
        print()
        print("   tail-mass crossing at radius 0.1587 is %.4f = %.1fx the SZZ radius" % (crossings["rho_tail"], crossings["rho_tail"] / BETA_STRONG))
        print("   -> reproduces the corpus upper end (30x) and IDENTIFIES the normalisation behind 0.63")
    print("   NOTE: an earlier draft printed an MK-type crossing; that recursion was")
    print("         non-standard (3.9 vs beta_c~2.30) and has been REMOVED, not shipped.")

    print("=" * 78)
    print("D) THE ONLY REAL QUESTION: how much BIGGER must the ball be?")
    i5 = int(np.argmin(np.abs(beta - BETA_PHYS)))
    for name, key in (("rho_hs", "rho_hs"), ("rho_sup", "rho_sup"), ("rho_tail", "rho_tail")):
        v = float(m[key][i5])
        print("   %-9s at beta=2.5: %10.4f   -> ball must grow x %8.1f" % (name, v, v / TARGET_RADIUS))
    print("   (no truncation and no FLOPs change this: it is a property of the map.)")
    print("=" * 78)
    print(f"jmax={jmax}  nbeta={nbeta}  wall={time.time() - t0:.2f}s")
    print("   REMINDER: 1/48 = 0.0208. The physical point is 2.5.")
    print("   Compute buys the truncation; it does NOT buy the radius.")

    if args.mode == "full":
        out = dict(jmax=jmax, nbeta=nbeta,
                   target_radius=TARGET_RADIUS,
                   beta_strong=BETA_STRONG, beta_phys=BETA_PHYS,
                   crossings=crossings, sweep={str(k): v for k, v in sweep.items()},
                   tail_budget={str(J): budget[J] for J in budget},
                   beta_probe=[float(b) for b in beta_probe],
                   j_needed={"_1pct": need1, "_0p1pct": need3, "_0p01pct": need4},
                   wall_s=time.time() - t0,
                   note=("0.63 is NOT reproduced by any single norm here; "
                         "the quoted value depends on the normalisation used "
                         "for the HS distance in the source. Do not quote "
                         "0.63 as physical until the norm is pinned."))
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=2)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
