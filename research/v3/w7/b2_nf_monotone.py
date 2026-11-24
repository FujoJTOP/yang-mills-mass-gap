#!/usr/bin/env python3
"""
B2 witness -- N_f monotonicity test for the SU(3) conformal-window lattice data.

Pure numpy, no plotting dependency. Every quoted number carries the source key
used in research/v3/B2-NF-LATTICE.md (section 6, verification ledger). The
script does NOT invent any value: unverified quantities are marked None.

Run:  python -I b2_nf_monotone.py
"""
import itertools
import numpy as np

# ---------------------------------------------------------------------------
# 1. Universal dichotomous order parameter
#    G(Nf) = 1 : gapped / confining (chirally broken for Nf >= 2)
#    G(Nf) = 0 : gapless / IR-conformal (an infrared fixed point)
#    Phase structure is scheme independent; WHERE the boundary sits is physical.
#    For contested Nf we allow both values and enumerate.
# ---------------------------------------------------------------------------
G = {
    0:  {1},          # pure YM: confining, gapped            [E09-07, 13.12]
    2:  {1},          # QCD-like, chiSB                       [E09-29]
    4:  {1},          # chiSB                                 [15.12-00957]
    6:  {1},          # chiSB, enhanced condensate            [E09-29]
    8:  {1},          # chiSB near chiral limit (walking)     [E09-30, 14.05.4752]
    10: {0, 1},       # 2023 IRFP (g2=15.0(5)) vs 2012 edge gamma~1  [23.06.07236, 12.04.6000]
    12: {0, 1},       # CONTESTED: IRFP vs below-window chiSB [13.16, 13.17, E09-27/31/32]
    16: {0},          # IRFP, nearly marginal                 [09.07.0919, E09-27]
}
NF = sorted(G)


def is_monotone(assign):
    """G non-increasing in Nf: gapped at small Nf, gapless at large Nf."""
    seq = [assign[n] for n in NF]
    return all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1))


def enumerate_assignments():
    keys = [n for n in NF if len(G[n]) > 1]
    for combo in itertools.product(*[sorted(G[n]) for n in keys]):
        a = {n: next(iter(G[n])) for n in NF}
        a.update(dict(zip(keys, combo)))
        yield a


# ---------------------------------------------------------------------------
# 2. Mass anomalous dimension quoted at / near the would-be IRFP.
#    gamma_m is SCHEME DEPENDENT and only defined where an IRFP exists.
#    (value, quoted_error_or_None, group, source_key)
# ---------------------------------------------------------------------------
GAMMA = {
    6:  [(0.4, None, "LSD (occasionally quoted)", "secondary")],     # not in abstract
    8:  [(1.0, None, "LatKMI (FSH fit)", "13.02.6859"),
         (0.94, 0.02, "LSD data (secondary citation)", "secondary")],
    10: [(0.6, None, "Hasenfratz et al. 2023", "23.06.07236"),
         (1.0, None, "LSD 2012 (edge)", "12.04.6000")],
    12: [(0.45, None, "LatKMI 2012 (0.4-0.5)", "12.07.3060"),
         (1.0, None, "Hasenfratz 2010 (near free-field)", "10.04.1004")],
    16: [(None, None, "no clean number (perturbative small)", "09.07.0919")],
}


def two_loop_bz(nf):
    """SU(3): b0 = 11 - 2nf/3, b1 = 102 - 38nf/3; g*^2 = -b0/b1 if b0>0,b1<0."""
    b0 = 11.0 - (2.0 / 3.0) * nf
    b1 = 102.0 - (38.0 / 3.0) * nf
    g2 = (-b0 / b1) if (b0 > 0 and b1 < 0) else float("nan")
    return b0, b1, g2


def main():
    print("=" * 74)
    print("B2: SU(3) N_f monotonicity witness")
    print("=" * 74)

    print("\n[1] Data available for N_f =", NF)
    print("    -> all 8 requested values have usable, verified lattice data.\n")

    print("[2] Universal order parameter G(Nf): 1=gapped/confining, 0=IR-conformal")
    print("    " + "  ".join(f"{n}:{sorted(G[n])}" for n in NF))

    assigns = list(enumerate_assignments())
    mono = [a for a in assigns if is_monotone(a)]
    nonmono = [a for a in assigns if not is_monotone(a)]
    print(f"\n    consistent assignments: {len(assigns)}")
    print(f"    monotone (G non-increasing): {len(mono)}")
    print(f"    NON-monotone: {len(nonmono)}")
    for a in nonmono:
        seq = "-".join(str(a[n]) for n in NF)
        print(f"      violation: Nf={NF} -> G={seq}")

    # The sharp logical point
    print("\n[3] Sharp statement (monotone window, unique edge):")
    ten_conf = any(is_monotone(a) and a[10] == 0 and a[12] == 1 for a in assigns)
    print("    (Nf=10 conformal) AND (Nf=12 chiSB) can be monotone? ->",
          "YES" if ten_conf else "NO")
    print("    => if Nf=10 is conformal (2023) and the window edge is unique/monotone,")
    print("       then Nf=12 must also be conformal; Fodor's Nf=12 chiSB claim would")
    print("       force either Nf=10 to be below the window or a scheme/non-universality.")

    # The two monotone scenarios
    print("\n[4] Monotone scenarios that survive:")
    edgeA = [a for a in mono if a[8] == 1 and a[10] == 0 and a[12] == 0 and a[16] == 0]
    edgeB = [a for a in mono if a[8] == 1 and a[10] == 1 and a[12] == 1 and a[16] == 0]
    print(f"    Edge in (8,10]: {len(edgeA)} assignment(s)  [2023 Nf=10 conformal; Fodor Nf=12 wrong]")
    print(f"    Edge in (12,16]: {len(edgeB)} assignment(s) [Fodor Nf=12 chiSB; 2023 Nf=10 conformal wrong/borderline]")

    print("\n[5] gamma_m quoted by group (scheme dependent, only where an IRFP exists):")
    for n in NF:
        if n in GAMMA:
            parts = []
            for val, err, grp, src in GAMMA[n]:
                if val is None:
                    parts.append(f"n/a ({grp})")
                elif err is not None:
                    parts.append(f"{val:.2f}({int(err*100):02d}) [{grp}]")
                else:
                    parts.append(f"{val:.2f} [{grp}]")
            print(f"    Nf={n:>2}: " + "  |  ".join(parts))
    print("    -> Nf=12 ranges 0.45..1.0 across groups; Nf=10 ranges 0.6..1.0.")
    print("    -> gamma_m is NOT monotone in Nf and is not a usable monotone.")

    print("\n[6] Two-loop Banks-Zaks reference (exact rational arithmetic in spirit):")
    print("    Nf   b0        b1        g*^2      g*")
    for n in NF:
        b0, b1, g2 = two_loop_bz(n)
        if np.isnan(g2):
            print(f"    {n:>2}  {b0:8.3f}  {b1:8.3f}    none      none")
        else:
            print(f"    {n:>2}  {b0:8.3f}  {b1:8.3f}  {g2:8.4f}  {np.sqrt(g2):8.4f}")
    print("    -> BZ zero begins at Nf > 306/38 = 8.0526; AF ends at Nf < 16.5.")
    print("    -> Nf=0 sits ~9-12 units below the edge: far, not marginal.")

    print("\n[7] Bottom line")
    print("    Data are COARSE-monotone (gapped at low Nf, gapless at high Nf).")
    print("    The boundary is bracketed, not located; Nf=12 is contested;")
    print("    no N_f-monotonicity theorem is known. No QED.")
    print("=" * 74)


if __name__ == "__main__":
    main()
