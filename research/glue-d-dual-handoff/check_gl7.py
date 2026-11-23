"""Reproducible numerical checks for GL7 (dual handoff, one RT).

(1) EXACT check (rational arithmetic) of the collapse identity

        (q1m <= q0m  and  q1s <= q0s)  <=>  Rm >= E * max(q1m, c*q1s)

    with q0_i = R_i / E, Rm = c*Rs, E > 0.  Only positivity of E and c is used,
    so random rationals exhaust the algebra (5000 cases, fixed seed).

(2) FRONTIER check: with one-loop t_up = C_RG/lam, both H3s hold iff

        lam >= lam_req := C_RG / ln(min(Rm/q1m, Rs/q1s)).

(3) FALSIFICATION cases: the two mu* = invariant/q1 differ; coverage needs
    mu_up <= min(mu*_m, mu*_s), so the smaller mu* decides.

Inputs consumed by GL7 (listed, not proved here):
    GL7-a  same RT (external obligation)
    CPB-W  weak side incl. a >=2-loop bound on t_up
    CPB-S^m, CPB-S^sigma  strong sides
    C = C_RG = 24*pi^2/11
    hcert_m, hcert_sigma  the two GL2b bridges

Exit code 0 iff all checks pass.  Run with:  python -I check_gl7.py
"""
from fractions import Fraction as F
import math
import random

C_RG = 24 * math.pi ** 2 / 11


def exact_collapse(n=5000, seed=20261122):
    rng = random.Random(seed)

    def pos():
        return F(rng.randint(1, 50), rng.randint(1, 50))

    bad = 0
    for _ in range(n):
        Rs, c, q1m, q1s, E = pos(), pos(), pos(), pos(), pos()
        Rm = c * Rs
        lhs = (q1m <= Rm / E) and (q1s <= Rs / E)
        rhs = Rm >= E * max(q1m, c * q1s)
        bad += (lhs != rhs)
    assert bad == 0, f"collapse identity failed in {bad}/{n} rational cases"
    return n


def frontier_holds(Rm, Rs, q1m, q1s, lam):
    t = C_RG / lam
    return (q1m <= Rm * math.exp(-t)) and (q1s <= Rs * math.exp(-t))


def check_frontier(n=2000, seed=7):
    rng = random.Random(seed)
    bad = 0
    for _ in range(n):
        Rm = 10 ** rng.uniform(-1, 1.3)
        Rs = Rm / rng.uniform(1.5, 5.0)          # Rs = sqrt(sigma)/Lambda
        q1m = Rm * rng.uniform(0.01, 0.2)
        q1s = Rs * rng.uniform(0.01, 0.2)
        lam_req = C_RG / math.log(min(Rm / q1m, Rs / q1s))
        # Away from the exact boundary (where equality holds and rounding
        # decides the float comparison), the equivalence must be exact.
        for lam in (lam_req * 0.98, lam_req * 0.99,
                    lam_req * 1.01, lam_req * 1.02):
            got = frontier_holds(Rm, Rs, q1m, q1s, lam)
            want = lam >= lam_req
            bad += (got != want)
    assert bad == 0, f"frontier equivalence failed in {bad} sampled points"
    return n


def falsification_report():
    """The two mu* differ; the smaller decides coverage."""
    Rm, Rs = 8.0, 2.0                 # m/Lambda, sqrt(sigma)/Lambda
    q1m, q1s = 1.0, 1.0
    mu_m, mu_s = Rm / q1m, Rs / q1s   # in units of Lambda
    cases = [
        # (t_up, expected H3_m, expected H3_s)
        (0.60, True, True),           # mu_up = 1.82 <= min(8,2)
        (1.00, True, False),          # 2.72 in (2,8]: mass yes, string no
        (2.10, False, False),         # 8.17 > 8: neither
    ]
    print(f"  mu*_m/Lambda = {mu_m:.3f},  mu*_sigma/Lambda = {mu_s:.3f}")
    for t, em, es in cases:
        mu_up = math.exp(t)
        hm = mu_up <= mu_m
        hs = mu_up <= mu_s
        print(f"  t_up={t:4.2f}  mu_up/Lambda={mu_up:6.3f}  "
              f"H3_m={hm!s:5}  H3_s={hs!s:5}  dual={hm and hs}")
        assert hm == em and hs == es, "falsification classification wrong"
    assert mu_s < mu_m, "with sqrt(sigma) < m the string mu* must be smaller"

    # strict two-loop (RB2, handoff-proof): t_up^{2L}(6) = 2.639
    t2 = 2.639
    q0m, q0s = Rm * math.exp(-t2), Rs * math.exp(-t2)
    print(f"  strict 2-loop t_up(lam=6)={t2}:  q0^m={q0m:.3f}  q0^sigma={q0s:.3f}")
    print(f"  => thresholds at lam=6: q1^m <= {q0m:.3f}, q1^sigma <= {q0s:.3f}")
    assert abs(q0s - 0.143) < 5e-3, "strict 2-loop sigma threshold should be ~0.143"
    return len(cases)


def main():
    print("GL7 consumed inputs: SameRT; CPB-W (>=2-loop t_up); CPB-S^m; "
          "CPB-S^sigma; C=C_RG=24*pi^2/11; hcert_m; hcert_sigma")
    print(f"C_RG = {C_RG:.6f}")
    n1 = exact_collapse()
    print(f"[1] exact rational collapse: {n1} cases OK")
    n2 = check_frontier()
    print(f"[2] frontier lambda_req equivalence: {n2} sampled points OK")
    print("[3] falsification: two mu* differ (smaller decides)")
    n3 = falsification_report()
    print(f"[3] {n3} classification cases OK")
    print("ALL GL7 NUMERICAL CHECKS PASS")


if __name__ == "__main__":
    main()
