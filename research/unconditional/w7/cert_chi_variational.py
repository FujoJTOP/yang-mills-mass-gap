"""
cert_chi_variational.py -- witness for CERT-CHI-VARIATIONAL (fork directive).

Goal: a genuinely action-specific variational UPPER bound on Var(a_{k->0}) = chi
(the susceptibility half of H1), and test whether it can reach Var <= 1/lambda at
beta ~ 2.5, or reduces to H1.

The object: decompose Var(a_{k->0}) for a = sum_x v_x X_x (X = su(2) log-link,
sum v_x^2 = 1) into
    Var(a) = [ single-link / diagonal part ]  +  [ connected off-diagonal sum ]
The diagonal part is a genuine action-specific computable (single-link marginal
under the Wilson weight). The off-diagonal is the connected two-point sum.

SU(2) conventions: U = cos(alpha) I + i sin(alpha) (n.sigma), class angle
alpha in [0,pi], tr U = 2 cos(alpha), Wilson weight e^{beta Re tr U} = e^{2 beta cos alpha},
Haar density propto sin^2(alpha) on [0,pi].

Run: py -3.12 research/unconditional/w7/cert_chi_variational.py
"""
import numpy as np
from scipy import integrate

np.set_printoptions(precision=6, suppress=True)


def moments_wilson(beta, n=4000):
    """E[alpha^k] under density ∝ sin^2(alpha) e^{2 beta cos alpha} on [0,pi]."""
    a = np.linspace(0.0, np.pi, n)
    w = np.sin(a) ** 2 * np.exp(2.0 * beta * np.cos(a))
    w /= np.trapezoid(w, a)
    m1 = np.trapezoid(a * w, a)
    m2 = np.trapezoid(a * a * w, a)
    return a, w, m1, m2, m2 - m1 ** 2


def main():
    print("=" * 74)
    print("PART 1  single-link (diagonal) variance sigma_link^2(beta), SU(2)")
    print("=" * 74)
    print(" density ~ sin^2(alpha) e^{2 beta cos alpha}; alpha = class angle in [0,pi]")
    print("  (Haar = beta=0; orders toward alpha=0 as beta grows -> variance DROPS)")
    print(f"{'beta':>6} {'E[alpha]':>10} {'E[alpha^2]':>12} {'Var(alpha)':>12}")
    sig2 = {}
    for beta in [0.0, 0.5, 1.0, 1.5, 2.5]:
        a, w, m1, m2, var = moments_wilson(beta)
        sig2[beta] = var
        print(f"{beta:>6.2f} {m1:>10.5f} {m2:>12.5f} {var:>12.6f}")
    print("\n  Closed-form Haar check (beta=0): Var = pi^2/12 - 1/2 =",
          np.pi**2 / 12 - 0.5)
    print("  -> diagonal sigma_link^2 is a genuine, computable, action-specific O(1) number.")

    print("\n" + "=" * 74)
    print("PART 2  the decomposition  Var(a) = diagonal + connected off-diagonal")
    print("=" * 74)
    print("""
 For a = sum_x v_x X_x with |v|=1:
   Var(a) = sum_x v_x^2 * sigma_link^2                     (DIAGONAL: computable)
          + sum_{x!=y} v_x v_y Cov(X_x, X_y)              (OFF-DIAGONAL: connected)

 - beta = 0 (pure Haar, strong coupling): links independent => off-diagonal = 0,
   Var(a) = sigma_link^2 / V_eff (compact confinement: the zero mode is PINNED).
 - beta > 0: the Wilson weight CORRELATES the links; the off-diagonal term grows
   (the modes align toward the weak-coupling / gapless direction).
""")
    print(" Numerically the split is exact; the whole difficulty lives in the")
    print(" off-diagonal (connected two-point) term.")

    print("\n" + "=" * 74)
    print("PART 3  the sharpest bound on the off-diagonal (and why it is H1)")
    print("=" * 74)
    print("""
 KL / spectral representation:  chi = Var(a) = sum_n |<0|X|n>|^2 / (E_n - E_0)
   (over excited states). If a mass gap Delta = E_1 - E_0 > 0 exists,
       Var(a) <= (1/Delta) * sum_n |<0|X|n>|^2 <= Var(X) / Delta .
   This is the sharpest spectral upper bound on the off-diagonal, and it NEEDS
   Delta > 0  ->  the mass gap  =  H1.
 So  chi <= Var(X)/Delta  is CERTIFIED only GIVEN the gap; without a gap, Delta=0,
 and the bound is +infinity (the off-diagonal is unbounded).
""")

    # illustrate the KL bound with a 3-level toy at a few gaps
    print(" Toy KL check (3 levels E=[0,Delta,2Delta], overlaps by |<0|X|n>|^2):")
    ov = np.array([0.0, 0.6, 0.4])  # overlaps; Var(X) = sum |X|^2
    VarX = ov.sum()
    for Delta in [0.1, 0.5, 1.0]:
        E = np.array([0.0, Delta, 2 * Delta])
        chi = np.sum(ov[1:] / E[1:])
        bound = VarX / Delta
        print(f"   Delta={Delta:>4}: chi={chi:>8.4f}  Var(X)/Delta={bound:>8.4f}"
              f"  (chi<=bound: {chi <= bound + 1e-12})  |Delta->0 => bound->inf")

    print("\n" + "=" * 74)
    print("PART 4  action-specific inequalities that COULD bound the off-diagonal")
    print("=" * 74)
    print("""
 - Gaussian domination (Froehlich-Simon-Spencer): bounds E[W(C)W(C')] for WILSON
   LOOPS W (characters), NOT for the log-link X. It does not control Cov(X,X).
 - FKG / positive association: needs ferromagnetic (increasing) coupling. The
   full-representation SU(N) Wilson weight is a character SUM with alternating
   signs -> NO FKG -> no sign control on the off-diagonal.
 - Brascamp-Lieb: needs log-concavity; the Savvidy Hessian is indefinite at
   beta~2.5 (rank-one negative) -> BL holds only on the positive fiber, softened
   by (1 - g|B|alpha)^{-1} (BLAST-OPERATOR-INEQ / BLAST-YM-NONLINEAR).
 => no available action-specific inequality bounds the X off-diagonal below 1/lambda
    at beta=2.5; each reduces to the gap.
""")

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print("""
 The sharpest action-specific variational bound = [diagonal sigma_link^2 (real,
 computable, O(1))] + [off-diagonal <= Var(X)/Delta (KL, needs the gap)].
 The diagonal is genuine but insufficient (O(1)/V_eff, does not survive the
 weak-coupling growth); the off-diagonal bound requires Delta = E_1-E_0 > 0 = the
 gap = H1. So reaching Var <= 1/lambda at beta=2.5 REDUCES TO H1. A genuinely
 action-specific bound exists and is non-vacuous on the diagonal, but the
 off-diagonal (where the crossover lives) is H1. Labels [T] decomposition + KL
 bound; [V-num] sigma_link^2; [O] off-diagonal = H1. No QED; Clay 0%.
""")
    print(" single-link diagonal variance at beta=2.5:", round(sig2[2.5], 6),
          " (Haar:", round(sig2[0.0], 6), ")")


if __name__ == "__main__":
    main()
