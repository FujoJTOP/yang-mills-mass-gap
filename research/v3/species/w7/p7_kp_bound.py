"""p7_kp_bound.py -- P3 total attack: is the KP per-rep bound
    c_r * dim r  <=  beta^{|r|} / (|r|! * e^beta)
a REAL blocker at beta = 2.5 ?

Reads: research/unconditional/P3-FULLREP-KP.md (lines 29,31,49,76,85,90,96,121).
Target file: research/v3/species/P7-KP.md

Exact SU(2) input (N=2): with exponent  (beta/2N)(trU + trU^dag) = beta*cos(theta),
the character coefficient is
    c_n(beta) = I_n(beta) - I_{n+2}(beta),      dim r_n = n+1,   |r| = n
(fundamental n=1, dim 2, c_1 = beta/2 + O(beta^3)).
Derived: e^{z cos th} = I_0 + 2 sum_k I_k cos(k th); sin th sin((m+1)th) = (cos(m th)-cos((m+2)th))/2;
Haar = (2/pi) sin^2 th dth on [0,pi]; int_0^pi cos(a th)cos(b th) dth = (pi/2) delta_ab (a,b>=1).
Verified below against direct numerical projection.
"""
import mpmath as mp

mp.mp.dps = 30

BETAS = [mp.mpf("0.018394"), mp.mpf("0.1"), mp.mpf("1"), mp.mpf("2.5")]
NULLS = list(range(1, 9))


# ---------------------------------------------------------------- exact SU(2)
def c_exact(n, beta):
    """character coefficient of chi_{r_n}, n = |r| = 2j boxes, dim = n+1."""
    z = beta
    return mp.besseli(n, z) - mp.besseli(n + 2, z)


def c_projection(n, beta):
    """independent check: c_n = int dU e^{beta Re tr U} chi_n(U)^*."""
    f = lambda th: (2 / mp.pi) * mp.sin(th) ** 2 * mp.e ** (beta * mp.cos(th)) * (
        mp.sin((n + 1) * th) / mp.sin(th) if n > 0 else mp.mpf(1)
    )
    # n=0: chi_0 = 1 -> above sin((n+1)th)/sin(th) = sin(th)/sin(th) = 1, fine
    return mp.quad(f, [0, mp.pi])


def dim_r(n):
    return n + 1


# ------------------------------------------------------- the two readings
def R_file(n, beta):
    """P3-FULLREP-KP.md reading (line 49): sum_{m>=n} b^m/m! <= b^n/n! * e^{+b}."""
    return beta**n / mp.factorial(n) * mp.e**beta


def R_quoted(n, beta):
    """the task's literal string beta^{|r|}/(|r|! e^beta): e^beta in the DENOMINATOR."""
    return beta**n / (mp.factorial(n) * mp.e**beta)


def main():
    print("=" * 96)
    print("P7 / KP:  c_r dim r <= beta^{|r|}/|r|! e^{+-beta}  --  four betas, |r| = 1..8 (SU(2), N=2)")
    print("=" * 96)

    # --- 0. self-check of the exact c_n
    print("\n[0] exact c_n(beta) = I_n(beta) - I_{n+2}(beta)  vs  direct Haar projection integral")
    for beta in BETAS:
        for n in (0, 1, 4):
            a, b = c_exact(n, beta), c_projection(n, beta)
            assert abs(a - b) < mp.mpf("1e-18"), (beta, n, a, b)
    print("    max |closed-form - projection| < 1e-18 over (beta,n) in {4 betas} x {0,1,4}   [E, self-check OK]")
    for beta in BETAS:
        print(f"    beta={float(beta):<9.6g} c_1 = {float(c_exact(1, beta)):<12.6g} "
              f"(beta/2 = {float(beta/2):<12.6g}; c_1/(beta/2)-1 = {float(c_exact(1, beta)/(beta/2)-1):<11.3e})")

    # --- 1. the two RHS tables
    for tag, fn in (("[1A] RHS as in P3-FULLREP-KP.md line 49:  beta^n/n! * e^{+beta}", R_file),
                    ("[1B] RHS as literally quoted in the task:  beta^n/(n! e^{beta})", R_quoted)):
        print("\n" + tag)
        hdr = "    n=|r| " + "".join(f"{float(b):>14.6g}" for b in BETAS)
        print(hdr)
        for n in NULLS:
            print(f"    {n:5d} " + "".join(f"{float(fn(n, b)):>14.6g}" for b in BETAS))

    # --- 2. exact side: c_r dim r, saturation c_r dim r / RHS
    print("\n[2] exact c_n, dim, c_n*dim  (SU(2); dim = n+1)")
    print("    n=|r|  dim   c_n(beta) for beta = " + ", ".join(str(float(b)) for b in BETAS))
    for n in NULLS:
        row = "".join(f"{float(c_exact(n, b)):>14.5g}" for b in BETAS)
        print(f"    {n:5d} {dim_r(n):4d} {row}")
    print("\n[2b] SATURATION  s_n := c_n dim r_n / RHS_A(n)   (1 = bound tight, 0 = bound useless)")
    print("    n=|r| " + "".join(f"{float(b):>14.6g}" for b in BETAS))
    for n in NULLS:
        print(f"    {n:5d} " + "".join(
            f"{float(c_exact(n, b) * dim_r(n) / R_file(n, b)):>14.6g}" for b in BETAS))

    # --- 3. which n "still satisfy" the bound, at beta = 2.5
    print("\n[3] beta = 2.5: does the bound still HOLD (RHS >= c_r dim r) for each |r|?")
    b = mp.mpf("2.5")
    for n in NULLS:
        lhs = c_exact(n, b) * dim_r(n)
        print(f"    |r|={n}: c_r dim r = {float(lhs):<11.5g}  RHS_A = {float(R_file(n, b)):<11.5g}"
              f"  slack = {float(R_file(n, b) / lhs):>9.4g}x   {'HOLDS' if R_file(n, b) >= lhs else 'FAILS'}")
    lhs0, rA0, rB0 = c_exact(0, b), R_file(0, b), R_quoted(0, b)
    print(f"    |r|=0 (trivial rep): c_0 = {float(lhs0):.6g};  RHS_A(0)=e^beta = {float(rA0):.6g}  -> "
          f"{'HOLDS' if rA0 >= lhs0 else 'FAILS'};  RHS_B(0)=e^-beta = {float(rB0):.6g} -> "
          f"{'HOLDS' if rB0 >= lhs0 else 'FAILS'}")
    print("    -> reading [1B] is falsified already at |r|=0 for every beta>0: c_0 ~ 1 > e^-beta.")

    # --- 4. the real KP sum  Phi(a) = sum_n c_n dim r_n e^{a n}
    print("\n[4] KP per-plaquette rep activity  Phi(a) = sum_{n>=1} c_n dim r_n e^{a n},  a = 1")
    print("    beta      Phi_exact(1)   majorant e^b(e^{b e}-1)   term-peak |r|   frac(|r|<=2)   frac(|r|<=4)")
    for beta in BETAS:
        a = mp.mpf(1)
        terms = {}
        for n in range(1, 60):
            t = c_exact(n, beta) * dim_r(n) * mp.e**(a * n)
            if abs(t) < mp.mpf("1e-40"):
                break
            terms[n] = t
        Phi = sum(terms.values())
        maj = mp.e**beta * (mp.e**(beta * mp.e**a) - 1)
        peak = max(terms, key=lambda k: terms[k])
        f2 = sum(v for k, v in terms.items() if k <= 2) / Phi
        f4 = sum(v for k, v in terms.items() if k <= 4) / Phi
        print(f"    {float(beta):<9.6g} {float(Phi):>14.6g} {float(maj):>23.6g} {peak:>15d} "
              f"{float(f2):>14.6g} {float(f4):>14.6g}")

    # --- 5. the KP hypothesis itself: lambda_a Phi(a) < 1  and  u*K < 1
    print("\n[5] KP HYPOTHESIS (P3-FULLREP-KP.md lines 90,96,121):  lambda_a Phi(a) < 1  [and Thm: u K < 1]")
    lam = 20 * mp.e                       # Delta = 20  (line 116)
    N = 2
    print("    lambda_a = e*Delta = 20e = %.4f ;  C_g = 2N e^b ;  K = 2 lambda C_g e ;  u = c_f" % float(lam))
    print("    beta        Phi(1)      lambda*Phi   u=c_1        K          u*K      verdict(KP hyp)")
    for beta in BETAS:
        a = mp.mpf(1)
        Phi = sum(c_exact(n, beta) * dim_r(n) * mp.e**(a * n) for n in range(1, 40))
        Cg = 2 * N * mp.e**beta
        K = 2 * lam * Cg * mp.e
        u = c_exact(1, beta)
        print(f"    {float(beta):<9.6g} {float(Phi):>10.5g} {float(lam*Phi):>12.5g} {float(u):>12.5g} "
              f"{float(K):>11.5g} {float(u*K):>11.5g}   {'HOLDS' if u*K < 1 else 'FAILS'}")
    # numeric KP radius with these very constants
    f = lambda bb: float(c_exact(1, mp.mpf(bb)) * (2 * lam * (2 * N * mp.e**mp.mpf(bb)) * mp.e)) - 1
    lo, hi = 1e-6, 0.2
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    print(f"    -> with these (crude, unoptimised) constants: u*K<1  <=>  beta < {lo:.5g}"
          f"   (1/48 = {1/48:.5g}; ratio {1/48/lo:.3g}x)")

    # --- 6. does the rep sum converge at beta = 2.5? (the only thing KP's per-rep bound is for)
    print("\n[6] convergence of the full rep sum at beta = 2.5 (the *only* job of c_r dim r <= RHS):")
    beta = mp.mpf("2.5")
    for a in (mp.mpf(1), mp.mpf(2), mp.mpf(4)):
        nmax, tails, peak, peakv = 400, mp.mpf(0), 0, mp.mpf(0)
        S = mp.mpf(0)
        terms = []
        for n in range(1, nmax + 1):
            t = c_exact(n, beta) * dim_r(n) * mp.e**(a * n)
            terms.append((n, t))
            S += t
            if t > peakv:
                peak, peakv = n, t
        last_ratio = terms[-1][1] / terms[-2][1]
        print(f"    a={float(a)}: Phi = {float(S):.6g};  term peak at |r|={peak};  "
              f"t(400)/t(399) = {float(last_ratio):.3g}  (<1 => converged);  "
              f"partials S(|r|<=2) = {float(sum(t for n, t in terms if n <= 2)):.5g}")
    print("    -> sum converges for every a (term_n ~ (beta/2)^n (n+1) e^{a n}/n! : factorial wins).")
    print("    NOTE: RHS_A(n) is INCREASING in beta for every n (see [1A]) -> the bound is TIGHTEST")
    print("          at small beta, LOOSEST at beta=2.5.  'beta small => RHS small => keep small |r|'")
    print("          has the sign of the beta-dependence inverted; and the bound never truncates |r|.")


if __name__ == "__main__":
    main()
