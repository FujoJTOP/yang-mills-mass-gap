#!/usr/bin/env python3
r"""
o02_bscan.py -- the beta-dependence of WHOLE-CURVE convexity (not just G''(0)).

H1A-EXACT.md knows ONE number: G''(0, theta=pi) = W0 - V0 beta^2, zero at beta*=4.30661.
Question here: as beta grows, does the ENTIRE curve G''(a) change sign together with G''(0),
or does a negative window appear somewhere else first?  Output:
  beta_cv  := the beta at which min_a G''(a) first touches 0 (loss of whole-curve convexity)
  a_neg    := the boundaries of {a : G''(a) < 0} as a function of beta
  a_vac    := the vacuum point of the model (H1A: minimum of the action on the slice set),
              for theta=pi this is a=-pi; G''(a_vac) is the curvature AT the physical point.

Same exact 1-dim quadrature as o02_plaquette.py (N=2).  N=3 via the reduced 3-dim one (imported).
Run: py -3.12 o02_bscan.py
"""
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o02_plaquette import (rho_su2, lrho_su2, fd2, Gpp_su2, Gpp_su3, sign_changes,
                           TWO_PI)

AG = np.linspace(-(TWO_PI - 0.08), TWO_PI - 0.08, 97)


def hdr(*a):
    print(*a, flush=True)


def curve(beta, theta, h=2e-2):
    return np.array([Gpp_su2(a, beta, theta, h) for a in AG])


def main():
    t0 = time.time()
    hdr("#" * 100)
    hdr("# o02_bscan -- where does whole-curve convexity die, vs where does G''(0) die?")
    hdr("#" * 100)
    hdr("\n[A] N=2 theta=pi:  beta scan.  min_a G''(a) and G''(0) side by side")
    hdr(f"     {'beta':>7}{'G(0)':>11}{'min_a G':>11}{'argmin':>8}{'G(-pi)=G(a_vac)':>17}"
        f"{'n(a:G<0)':>10}   negative window(s) on a in (-2pi,2pi)")
    prev = None
    bcv = None
    for beta in np.arange(0.0, 8.01, 0.25):
        g = curve(beta, np.pi)
        mn = g.min()
        # negative windows as intervals
        neg = g < 0.0
        wins = []
        i = 0
        while i < len(AG):
            if neg[i]:
                j = i
                while j + 1 < len(AG) and neg[j + 1]:
                    j += 1
                wins.append((AG[i], AG[j]))
                i = j + 1
            else:
                i += 1
        if bcv is None and mn < 0.0:
            bcv = beta
        hdr(f"     {beta:7.2f}{g[len(AG)//2]:+11.5f}{mn:+11.5f}{AG[int(np.argmin(g))]:8.3f}"
            f"{Gpp_su2(-np.pi, beta, np.pi):+17.5f}{int(neg.sum()):10d}   "
            + ("; ".join(f"[{x:+.2f},{y:+.2f}]" for x, y in wins) if wins else "none"))
    hdr(f"\n     beta_cv (whole-curve convexity lost, theta=pi) <= {bcv:.2f}"
        f"   vs   beta* = 4.30661 (G''(0) alone)")

    # bisection on beta_cv
    lo, hi = bcv - 0.25, bcv
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        if curve(mid, np.pi).min() < 0.0:
            hi = mid
        else:
            lo = mid
    hdr(f"     bisection: beta_cv = {0.5*(lo+hi):.6f}   (min_a G'' there = {curve(0.5*(lo+hi), np.pi).min():+.2e})")

    # derivative of min_a G'' at 0 in beta (small-beta behaviour)
    hdr("\n[B] small-beta behaviour of the curve: fitted min_a G''(a) for beta <= 2")
    bb, mm = [], []
    for beta in np.arange(0.0, 2.01, 0.1):
        bb.append(beta)
        mm.append(curve(beta, np.pi).min())
    bb = np.array(bb)
    mm = np.array(mm)
    A = np.vstack([np.ones_like(bb), bb ** 2, bb ** 4]).T
    coef, *_ = np.linalg.lstsq(A, mm, rcond=None)
    hdr(f"     min_a G'' ~ c0 + c2 beta^2 + c4 beta^4 :  c0={coef[0]:+.6f} c2={coef[1]:+.6f} c4={coef[2]:+.6f}")
    hdr("     (c0 = the beta=0 Haar floor of the curve, c2 < 0 = the Savvidy variance term)")

    hdr("\n[C] N=2 theta=0 (vacuum AT the conditioning point): does it stay convex everywhere?")
    hdr(f"     {'beta':>7}{'G(0)':>11}{'min_a G':>11}{'argmin':>8}{'n(a:G<0)':>10}   windows")
    for beta in (0.0, 2.5, 4.0, 6.0, 8.0):
        g = curve(beta, 0.0)
        neg = g < 0.0
        wins = []
        i = 0
        while i < len(AG):
            if neg[i]:
                j = i
                while j + 1 < len(AG) and neg[j + 1]:
                    j += 1
                wins.append((AG[i], AG[j]))
                i = j + 1
            else:
                i += 1
        hdr(f"     {beta:7.2f}{g[len(AG)//2]:+11.5f}{g.min():+11.5f}{AG[int(np.argmin(g))]:8.3f}"
            f"{int(neg.sum()):10d}   " + ("; ".join(f"[{x:+.2f},{y:+.2f}]" for x, y in wins) if wins else "none"))

    hdr("\n[D] N=3 (U(2) Cartan v=T8): beta scan at the conditioning point and at the vacuum a=-pi")
    hdr(f"     {'beta':>7}{'G(0)':>11}{'G(-pi)':>11}{'G(0.43)':>11}{'G(-0.5)':>11}")
    for beta in (0.0, 1.0, 2.5, 4.0, 5.7, 6.0, 8.0):
        b = max(beta, 1e-3)
        hdr(f"     {beta:7.2f}{Gpp_su3(0.0, b, h=0.2):+11.5f}{Gpp_su3(-np.pi, b, h=0.2):+11.5f}"
            f"{Gpp_su3(0.43, b, h=0.2):+11.5f}{Gpp_su3(-0.5, b, h=0.2):+11.5f}")
    hdr(f"\nwall {time.time()-t0:.1f}s   Clay 0%.  No QED.")


if __name__ == '__main__':
    main()
