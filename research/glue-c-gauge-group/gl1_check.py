#!/usr/bin/env python3
"""GL1 executable check: the single-measure dictionary.

The shared construction (gluing §1.2) fixes

    beta = 2 N_c / g^2 ,     S_W = (beta/N_c) sum_p (1 - Re tr U_p),

so as a *probability measure* (the constant term cancels in Z)

    mu_{L,a,beta}  ~  exp( (beta/N_c) sum_p Re tr U_p ) prod dU .

Hence the only parameter that enters the measure is

    gamma := beta / N_c        (coefficient of sum_p Re tr U_p).

Every session whose written weight is exp(-(beta_s/N_s) sum_p (1 - chi))
has gamma_s = beta_s/N_s and is *identical* to the shared measure when
beta_s/N_s = beta_A0/N_c.  A session writing exp(beta_s sum_p chi) has
gamma_s = beta_s, so it must use beta_s = beta_A0/N_c -- a factor N_c
different from the symbol beta_A0.

This script verifies the two claims exactly on a finite group lattice:
  (C1) constant-prefactor lemma: the ratio of the two weights is
       configuration-independent iff the dictionary is applied;
  (C2) normalised expectations agree under the dictionary and disagree
       when the same numeric beta is plugged into both forms.

Run: python gl1_check.py [--group S3] [--lattice 2x1] [--beta 2.0]
"""

from __future__ import annotations
import argparse
import math

from gl6_check import make_z2, make_s3, Lattice, build_tables, decode


def plaquette_chi_sum(g, lat, chi_i, ptab, ds):
    q = g.q
    tot = 0.0
    for p, (a, b, c, d) in enumerate(lat.plaq):
        h = ptab[p][((ds[a] * q + ds[b]) * q + ds[c]) * q + ds[d]]
        tot += chi_i[h]
    return tot


def weights(g, lat, beta, form):
    """Return (unnormalised exponent, list of weights) for a written form."""
    q, L, N = g.q, lat.nlinks, g.q ** lat.nlinks
    chi_i, ptab = build_tables(g, lat)
    ws = []
    for idx in range(N):
        ds = decode(idx, q, L)
        csum = plaquette_chi_sum(g, lat, chi_i, ptab, ds)
        if form == "A0":           # exp(-(beta/N_c) sum_p (1 - chi))
            e = -(beta / g.dim) * (len(lat.plaq) - csum)
        elif form == "os3":        # exp(beta sum_p chi)   (unnormalised trace)
            e = beta * csum
        else:
            raise ValueError(form)
        ws.append(math.exp(e))
    return ws


def expectation(ws, F):
    Z = sum(ws)
    return sum(w * f for w, f in zip(ws, F)) / Z


def run(group_name, Lx, Ly, beta):
    g = make_z2() if group_name == "Z2" else make_s3()
    lat = Lattice(Lx, Ly)
    q, L, N = g.q, lat.nlinks, g.q ** lat.nlinks
    chi_i, ptab = build_tables(g, lat)
    print(f"group={g.name} (dim f={g.dim})  lattice={Lx}x{Ly}  configs={N}  beta_A0={beta}")
    print(f"dictionary: beta_os3 = beta_A0/N_c = {beta/g.dim:.6g}")

    # gauge-invariant observable: chi of the first plaquette (unnormalised)
    F = []
    for idx in range(N):
        ds = decode(idx, q, L)
        a, b, c, d = lat.plaq[0]
        h = ptab[0][((ds[a] * q + ds[b]) * q + ds[c]) * q + ds[d]]
        F.append(chi_i[h])

    wA0 = weights(g, lat, beta, "A0")
    w3_ok = weights(g, lat, beta / g.dim, "os3")
    w3_bad = weights(g, lat, beta, "os3")

    # (C1) ratio configuration-independent?
    ratios = [a / b for a, b in zip(wA0, w3_ok)]
    rmin, rmax = min(ratios), max(ratios)
    dev = (rmax - rmin) / max(abs(rmin), 1e-300)
    print(f"(C1) max relative deviation of w_A0 / w_os3 with beta_os3=beta_A0/N_c: {dev:.3e}")

    EA0 = expectation(wA0, F)
    E3ok = expectation(w3_ok, F)
    E3bad = expectation(w3_bad, F)
    print(f"(C2) E[chi] in A0 form                     : {EA0:.12f}")
    print(f"(C2) E[chi] in os3 form, beta=beta_A0/N_c  : {E3ok:.12f}   diff={abs(E3ok-EA0):.3e}")
    print(f"(C2) E[chi] in os3 form, beta=beta_A0      : {E3bad:.12f}   diff={abs(E3bad-EA0):.3e}")

    ok1 = dev < 1e-12
    ok2 = abs(E3ok - EA0) < 1e-12
    ok3 = abs(E3bad - EA0) > 1e-6
    print(f"(C1) dictionary applied  : {'PASS' if ok1 else 'FAIL'}")
    print(f"(C2) dictionary applied  : {'PASS' if ok2 else 'FAIL'}")
    print(f"(C2) same numeric beta   : {'PASS (differs, as predicted)' if ok3 else 'FAIL'}")

    # SU(3) handoff-point illustration (pure arithmetic, no enumeration)
    print("\nSU(3) illustration: beta_A0 = 2N/g^2 = 6 at g^2 = 1.")
    print("  A0 coefficient gamma = beta_A0/N_c = 2.")
    print("  A session writing exp(beta_their sum Re tr) needs beta_their = 2, not 6.")
    print("  Plugging beta_their = 6 is the g^2 = 1/3 theory: a factor 3 in the exponent.")

    if not (ok1 and ok2 and ok3):
        raise SystemExit("GL1 CHECK FAILED")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="S3", choices=["Z2", "S3"])
    ap.add_argument("--lattice", default="2x1")
    ap.add_argument("--beta", type=float, default=2.0)
    args = ap.parse_args()
    Lx, Ly = (int(s) for s in args.lattice.lower().split("x"))
    run(args.group, Lx, Ly, args.beta)
    print("\nGL1 CHECK: PASS")


if __name__ == "__main__":
    main()
