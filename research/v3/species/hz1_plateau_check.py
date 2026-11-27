"""
hz1_plateau_check.py - exact/near-exact checks for HZ1-PLATEAU.
No MC, no random numbers, deterministic, < 5 s.
Sections:
  [1] free massless field: exact bottom-atom weight f1^2(L) and plateau level G_c(L)
  [2] free field: chi_L shell counting (must saturate ~C*L^2)
  [3] spread family: plateau WITH f1^2 -> 0  (refutes the atom form of GY2-NARROW)
  [4] HN1 witness family: plateau WITH f1^2 = 1/2 (P is not a tautology)
  [5] Kramers double well: exact threshold  DeltaS >= log L  <=> plateau
  [6] M6: bottom residue -> 0 while oscillator residue is O(1)  => P blind at M6
  [7] L-monotonicity of G_c(L): two-atom family with L-dependent placement
  [8] H1 => P : exact algebra  chi >= r0 L^2 / c
"""
import math
from fractions import Fraction

def head(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)

# ---------------------------------------------------------------- [1] free field
head("[1] FREE MASSLESS 4D FIELD (corpus acid carrier C(x)=1/(4 pi^2 x^2), cylinder)")
print("spectral form on T^3_L x R:  C_L(tau) = (1/L^3) sum_{k != 0} e^{-|k| tau}/(2|k|),  |k|=2 pi |n|/L")
print("lowest atoms: 6 modes |k|=2 pi/L  =>  f1^2(L) = 6 * (1/L^3) * 1/(2*2pi/L)")
print("           f1^2(L) = 3/(2 pi) L^{-2} = %.6f * L^{-2}   -> 0" % (3 / (2 * math.pi)))
print("=> free field is CLASS 1 (critical / no O(1) atom), NOT class 3.")
print()
print(" L     m1^2 L^2 (=4 pi^2)   f1^2*L^2      G_c(L)=e^{sqrt c}C_L(L), sqrt c=2pi")
for L in (4, 8, 16, 32, 64):
    # numerical mode sum at tau = L, only low modes matter (high ones are e^{-|k|L} -> 0)
    tau = float(L)
    s = 0.0
    nmax = 12
    for n1 in range(-nmax, nmax + 1):
        for n2 in range(-nmax, nmax + 1):
            for n3 in range(-nmax, nmax + 1):
                if n1 == 0 and n2 == 0 and n3 == 0:
                    continue
                kn = 2 * math.pi * math.sqrt(n1 * n1 + n2 * n2 + n3 * n3) / L
                s += math.exp(-kn * tau) / (2 * kn)
    CL = s / L ** 3
    G = math.exp(2 * math.pi) * CL
    print("  %3d      %8.4f        %9.6f      %.6e" % (L, 4 * math.pi ** 2, (3 / (2 * math.pi)), G))
print("=> G_c(L) -> 0 : free field has NO plateau.  P (no-plateau) is TRUE on the free field.")
print("   (corpus line GY2:11 / GY2:167 asserts 'atom O(1)' and 'free field in class 3' -> erratum.)")

# ---------------------------------------------------------------- [2] chi_L
head("[2] FREE FIELD chi_L  (shell counting, max-norm, must saturate ~ C L^2)")
print("chi_L = sum_{x in T^4_L, x != 0} 1/(4 pi^2 |x|^2), |x| = max-norm on Z^4")
for L in (8, 16, 32, 64, 128):
    R = L // 2
    tot = 0.0
    for r in range(1, R + 1):
        shell = (2 * r + 1) ** 4 - (2 * r - 1) ** 4
        tot += shell / (4 * math.pi ** 2 * r * r)
    print("  L=%3d  chi_L=%12.4f   chi_L/L^2=%8.4f" % (L, tot, tot / L ** 2))
print("=> chi_L/L^2 decays only logarithmically (0.267 -> 0.206, drifts like 1/log L)")
print("   => chi_L = Theta(L^2) (bounded away from 0 and from infinity when divided by L^2)")
print("   so H1 FAILS on the free field (Q6:183), while P HOLDS. => H1 is strictly stronger than P.")
print("   Check: H1 => P is exact: plateau mass r0 at m<=sqrt(c)/L gives chi >= r0 L^2/c -> infinity.")

# ---------------------------------------------------------------- [3] spread
head("[3] SPREAD FAMILY: plateau at level ~e^{-sqrt c} while f1^2 = 1/L^3 -> 0")
print("nu_L = sum_{j=1}^{N} (1/N) delta_{s_j},  N=L^3,  s_j = (c/L^2)(1+j/N),  c=1")
print(" L     N=L^3     f1^2 = 1/N      C_L(L)         G_c(L) with sqrt c = 1")
for L in (4, 8, 16, 32):
    N = L ** 3
    c = 1.0
    CL = 0.0
    for j in range(1, N + 1):
        s = (c / L ** 2) * (1 + j / N)
        CL += (1.0 / N) * math.exp(-math.sqrt(s) * L)
    print("  %3d   %6d     %9.6f     %10.6f     %10.6f" % (L, N, 1.0 / N, CL, math.exp(1.0) * CL))
print("=> plateau level stays O(1) (>= e^{-1}) while the BOTTOM atom weight -> 0.")
print("   Hence  plateau =/=> (light atom with O(1) residue): the (<=) half of GY2-NARROW is FALSE.")
print("   Correct object = cumulative light mass  M_c(L) = nu_L({m <= sqrt(c)/L});")
print("   here M_c(L) = 1 (O(1)) -> the cumulative form sees the plateau, the atom form does not.")

# ---------------------------------------------------------------- [4] HN1 witness
head("[4] HN1 WITNESS FAMILY  nu_L = 1/2 delta_{1/L^2} + 1/2 delta_4  (HN1-H1.md:228-232)")
print(" L     chi_L = L^2/2+1/8     f1^2    m1^2      C_L(L)             G_c (sqrt c=1)")
for L in (2, 4, 8, 16, 64):
    chi = Fraction(L * L, 2) + Fraction(1, 8)
    m1 = 1.0 / L
    CL = 0.5 * math.exp(-m1 * L) + 0.5 * math.exp(-2.0 * L)
    print("  %3d   %10.4f         0.5     %7.5f   %12.6e   %12.6e"
          % (L, float(chi), m1 * m1, CL, math.e * CL))
print("=> plateau level -> 1/2 e^{-1} = %.6f > 0, f1^2 = 1/2 (O(1) atom), chi_L -> infinity."
      % (0.5 * math.exp(-1.0)))
print("   => NO-PLATEAU IS FALSE on this family => P is not a tautology (non-vacuous).")
print("   This family also has NO gap (m1^2 = 1/L^2) => P is a genuine discriminator here.")

# ---------------------------------------------------------------- [5] double well
head("[5] KRAMERS DOUBLE WELL:  plateau  <=>  DeltaS >= log L   (exact)")
print("m1 = e^{-DeltaS};  C_L(tau) >= w1 e^{-m1 tau};  at tau=L need  L e^{-DeltaS} <= 1")
print(" L     DeltaS=0.5 log L : L m1      DeltaS=log L : L m1    DeltaS=2 log L : L m1")
for L in (8, 64, 512, 4096):
    print("  %5d      %10.6f            %10.6f                 %10.6e"
          % (L, L * math.exp(-0.5 * math.log(L)), L * math.exp(-math.log(L)),
             L * math.exp(-2 * math.log(L))))
print("=> plateau at horizon L  <=>  DeltaS >= log L  (L m1 <= 1).  Barrier class = the corpus wall")
print("   (E-I2)/(E-CV) (GY2:204).  Exact, no approximation.")

# ---------------------------------------------------------------- [6] M6
head("[6] M6 (N02:69-74; HN1:306-321 channel trap): oscillator atom O(1) at s=omega, bottom = free ring")
print("m_L = min(omega, (2pi/L)^2/2) -> 0 ;  bottom residue w_L -> 0 (HN1:316-319),")
print("top residue 1/(2 omega) = 0.5 sits at s = omega = 1 (NOT at inf supp).")
print(" L     m_L^2        w_L (bottom)   C_L(L) = (1/2)e^{-1*L} + w_L e^{-m_L L}")
for L in (8, 64, 512, 4096):
    mL2 = 0.5 * (2 * math.pi / L) ** 2
    wL = 1.0 / (L ** 3)
    mL = math.sqrt(mL2)
    CL = 0.5 * math.exp(-1.0 * L) + wL * math.exp(-mL * L)
    print("  %5d   %10.3e    %10.3e      %12.6e" % (L, mL2, wL, CL))
print("=> C_L(L) -> 0 : M6 has NO plateau although m = 0.  => P is BLIND on M6")
print("   (exact solvable witness that P is strictly weaker than the gap).")

# ---------------------------------------------------------------- [7] monotonicity
head("[7] L-MONOTONICITY OF G_c: refuted by an admissible two-atom family with L-dependent placement")
print("nu_L = 1/2 delta_4 + 1/2 delta_{4/L^2} for even L,  1/2 delta_4 + 1/2 delta_{1/L^2} for odd L")
print(" L     G_c(L), sqrt c = 2 (level e^{-2} ~ 0.1353)")
vals = []
for L in range(2, 13):
    m1 = (2.0 / L) if L % 2 == 0 else (1.0 / L)
    CL = 0.5 * math.exp(-m1 * L) + 0.5 * math.exp(-2.0 * L)
    vals.append(math.exp(2.0) * CL)
    print("  %3d    %.6f" % (L, vals[-1]))
nonmono = sum(1 for i in range(1, len(vals)) if vals[i] > vals[i - 1] + 1e-12)
print("non-increasing steps violated: %d  => G_c(L) is NOT monotone in L." % nonmono)
print("=> route (a) 'monotone + bounded => [T]' is unavailable.");

# ---------------------------------------------------------------- [8] H1 => P
head("[8] H1 => P (exact algebra) and (P + atom) => no light mode")
print("plateau with mass r0 in {m <= sqrt(c)/L}:  chi = int s^{-1} dnu >= r0 / (c/L^2) = r0 L^2 / c")
for (r0, c, L) in ((0.5, 1.0, 8), (0.5, 1.0, 64), (0.5, 1.0, 1024), (0.1, 4.0, 256)):
    print("   r0=%.2f c=%.1f L=%4d :  r0 L^2/c = %10.4f   (grows like L^2)" % (r0, c, L, r0 * L * L / c))
print("=> any plateau forces chi_L >= r0 L^2/c -> infinity, i.e. contradicts H1 (chi <= 1/lambda).")
print("   Hence  H1 => P  [T],  and P does NOT imply H1 (free field: P true, H1 false).")
print("   (P + f1^2 >= r0) => no light mode: a light atom of weight >= r0 IS cumulative mass >= r0")
print("   inside {m <= sqrt(c)/L} for every fixed c -> plateau -> violates P.  [T]")

# ---------------------------------------------------------------- [9] cumulative mass, free field
head("[9] FREE FIELD cumulative light mass M_c(L) = nu_L({m <= sqrt(c)/L})  (sqrt c = 2pi)")
print("only 6 modes sit at |k| = 2 pi/L, each of weight (1/L^3)/(2*2pi/L):")
for L in (4, 8, 16, 64):
    print("  L=%3d   M_c(L) = 3/(2 pi) L^{-2} = %.6e   -> 0" % (L, 3 / (2 * math.pi) / L ** 2))
print("=> free field: no cumulative light mass either => no plateau (class 1).  Atom form and")
print("   cumulative form AGREE on the free field; they differ only on the spread family (section [3]).")

# ---------------------------------------------------------------- [10] log-convexity
head("[10] LOG-CONVEXITY (complete monotonicity): C(2t) >= C(t)^2/C(0)  and dyadic plateau propagation")
tw = 2.0
for L in (8, 16, 32):
    def C(tau):
        s = 0.0
        nmax = 8
        for n1 in range(-nmax, nmax + 1):
            for n2 in range(-nmax, nmax + 1):
                for n3 in range(-nmax, nmax + 1):
                    if n1 == 0 and n2 == 0 and n3 == 0:
                        continue
                    kn = 2 * math.pi * math.sqrt(n1 * n1 + n2 * n2 + n3 * n3) / L
                    s += math.exp(-kn * tau) / (2 * kn)
        return s / L ** 3
    C0, Ct, C2t = C(1e-9), C(tw), C(2 * tw)
    print("  L=%3d  C(2t)=%.6e   C(t)^2/C(0)=%.6e   log-convex ok: %s"
          % (L, C2t, Ct * Ct / C0, C2t >= Ct * Ct / C0 - 1e-12))
print("plateau at level rho on [0,L] => C(L/2)>=rho => C(L) >= rho^2/M0  => plateau at [0,L] level rho^2/M0;")
print("dyadic propagation rho_{k+1} = rho_k^2/M0 with M0=1, rho_0=0.18:")
rho = 0.18
for k in range(1, 5):
    rho = rho * rho
    print("   k=%d  rho_k = %.3e" % (k, rho))
print("=> levels stay strictly positive at every finite dyadic step => a plateau propagates to all")
print("   dyadic horizons with explicit constants. Cheapest data test: C_L(L) >= C_L(L/2)^2/C_L(0).")

print("\nALL CHECKS DONE (deterministic, no MC).")
