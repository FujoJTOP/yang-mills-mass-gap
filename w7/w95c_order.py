"""W95C-ORDER: numeric backup for research/v3/species/W95C-ORDER.md

Sec 1: u(beta) = I2/I1 (fundamental character coefficient, U01-F1AUDIT.md:111)
       anchor: u(1/48) = 5.208e-3  (T23-TEMPLATE.md:64 [V])
Sec 2: m_c = a_c*ln(1/u) + O(1), a_c(0++)=a_c(2++)=1, others a_c>=2
       -> m2/m0 -> 1 (slowly); (m_c - m0) -> +inf for a_c >= 2
Sec 3: inf spec = min(Delta_0, Delta_2): positive iff BOTH positive.
       Witness: Delta_0 = delta > 0, Delta_2 = 0  =>  inf = 0.
"""
import numpy as np
from scipy.special import iv

def u(beta):
    return iv(2, beta) / iv(1, beta)      # SU(2)-type fundamental character coeff, u ~ beta/4

def hdr(s):
    print("\n" + "=" * 74)
    print(s)
    print("=" * 74)

hdr("W95C-ORDER : mass ordering of 0++, 2++, 0-+, 1--  (numeric backup only)")
print("m_c(beta) = a_c * ln(1/u(beta)) + c_c ,   a_0++ = a_2++ = 1   [PROOF-IN11.md:9,11,23]")
print("c_c = O(1) subleading remainder, unknown sign  [PROOF-IN11.md:23]")
print("MODEL for c_c: symmetric split c_0 = -C, c_2 = +C  (C in {0, 0.5, 1}); illustrative only,")
print("               the remainder is O(1) of UNKNOWN sign [PROOF-IN11.md:23] -> needs m_c>0")

# ---------------------------------------------------------------- Sec 1
hdr("[1] u(beta) and the anchor u(1/48)")
for b in [1/48, 0.05, 0.1, 0.5, 1.0, 2.5]:
    print("    beta = %8.5f   u = %.6e   ln(1/u) = %8.4f" % (b, u(b), np.log(1 / u(b))))
ref = 5.208e-3
print("    anchor: u(1/48) = %.6e  vs corpus 5.208e-3  -> rel.dev = %.2e"
      % (u(1/48), abs(u(1/48) - ref) / ref))

# ---------------------------------------------------------------- Sec 2
hdr("[2] m2/m0 -> 1  and  (m_c - m0) -> +inf")
Cs = (0.0, 0.5, 1.0)
print("     beta        u          ln(1/u)    m2/m0 C=0    C=.5     C=1     (m_c-m0) a_c=2")
for b in [1.0, 0.5, 0.2, 0.1, 0.05, 1/48, 0.01, 0.001]:
    L = np.log(1 / u(b))
    r = [((L + C) / (L - C)) if (L - C) > 0 else np.nan for C in Cs]
    print("   %8.5f  %.4e  %8.4f   %9.6f %8.6f %8.6f   %+10.4f" % (b, u(b), L, r[0], r[1], r[2], L))
print("     (m2/m0 at the corpus window edge beta=1/48: C=1 -> %.2f ; C=.5 -> %.2f)"
      % ((np.log(1 / u(1/48)) + 1) / (np.log(1 / u(1/48)) - 1),
         (np.log(1 / u(1/48)) + 0.5) / (np.log(1 / u(1/48)) - 0.5)))

print("""
    => m2/m0 -> 1 as beta -> 0 (u -> 0), but SLOWLY: ln(1/u) is only O(1)-O(10)
       across the whole strong-coupling window, and the O(1) split can flip the
       sign of (m2 - m0). Degeneracy at leading order != an ordering theorem.
    => a_c >= 2 channels: m_c - m0 >= (a_c-1) ln(1/u) + O(1) -> +inf  [PROOF-IN11.md:13]
       (last column, a_c = 2, C = 0) : strictly growing, so 0-+ never competes.""")

# ---------------------------------------------------------------- Sec 3
hdr("[3] inf spec = min(Delta_0, Delta_2) : m>0 needs BOTH")
print("    certificate accounting: cert = min(cert_0++, cert_2++)  [IN-10-attack.md:7]")
print("    inf_c Delta_c = min(D_0, D_2)                          [PROOF-IN11.md:17]\n")
print("      D_0++    D_2++   min    inf>0 ?   reading")
grid = [(1.0, 1.0), (1.0, 0.5), (0.5, 1.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0), (0.0, 1e-9)]
for d0, d2 in grid:
    m = min(d0, d2)
    ok = m > 0
    why = "gap" if ok else ("0++ certified, gap LOST (2++ floor 0)" if d0 > 0 else
                            ("2++ certified, gap LOST (0++ floor 0)" if d2 > 0 else "no gap"))
    print("   %7.4f  %7.4f  %7.4f   %-6s   %s" % (d0, d2, m, ok, why))

print("""
    min(x,y) > 0  <=>  x > 0  AND  y > 0.
    Witness D_0 = delta > 0, D_2 = 0  =>  inf spec = 0 : covering 0++ alone does
    NOT give m>0. Same with the roles swapped. (Isomorphic to the T6-5(i)
    witness m1 = 0, I = T, G = F  [V03-I2GAP.md:120].)

    Escape hatch 'channel-blind certificate' needs cert_blind <= inf_c Delta_c
    with zero slack = the gap verbatim, not constructible  [C2-bypass-INb.md:13].
    => m>0 requires 0++ AND 2++, plus (b') membership for 0-+/1--/higher spin.""")
hdr("END -- numeric backup only; no theorem is added by this script.")
