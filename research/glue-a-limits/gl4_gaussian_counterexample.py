"""GL4-C sharpness check: weak limits do NOT preserve clustering without a
uniform-in-a bound.

The 4D Euclidean massive scalar propagator is

    G_m(x) = m * K_1(m r) / (4 pi^2 r),     r = |x|,

where K_1 is the modified Bessel function, evaluated below from its exact
integral representation

    K_1(z) = int_0^inf exp(-z cosh t) cosh t dt .

For every fixed r > 0 one has G_m(x) -> 1/(4 pi^2 r^2) as m -> 0 (K_1(z) ~ 1/z).
The Gaussian measures with covariance G_{m_n}, m_n -> 0, are reflection
positive and converge weakly to the massless Gaussian free field, whose
correlations decay like 1/r^2 -- not exponentially.  Hence GL4-C needs the
a-uniform bound (U_a); RP and weak convergence alone do not give it.

Run with -I (isolated): no local imports on purpose. Standard library only.
"""
import math


def k1(z: float) -> float:
    """K_1(z) for z > 0 via Simpson quadrature of the exact integral."""
    if z <= 0.0:
        raise ValueError("z must be positive")
    # integrand decays like exp(-z e^t / 2); t_max = 14 is ample for z >= 1e-3
    t_max = 14.0
    n = 20000
    h = t_max / n

    def f(t: float) -> float:
        return math.exp(-z * math.cosh(t)) * math.cosh(t)

    s = f(0.0) + f(t_max)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(i * h)
    return s * h / 3.0


def G(m: float, r: float) -> float:
    return m * k1(m * r) / (4.0 * math.pi ** 2 * r)


def main() -> None:
    ms = (1.0, 0.25, 0.0625, 0.015625)
    print("4D massive Gaussian propagator  G_m(r) = m K_1(m r) / (4 pi^2 r)")
    print()
    print(f"{'r':>4}  " + "  ".join(f"m={m:<8g}" for m in ms)
          + f"  {'1/(4pi^2 r^2)':>13}  {'ratio(1/64)':>11}")
    for r in (1.0, 2.0, 5.0):
        vals = [G(m, r) for m in ms]
        massless = 1.0 / (4.0 * math.pi ** 2 * r * r)
        print(f"{r:4.1f}  " + "  ".join(f"{v:<10.6f}" for v in vals)
              + f"  {massless:13.6f}  {vals[-1] / massless:11.6f}")
    print()
    ref = 1.0 / (4.0 * math.pi ** 2)
    print(f"sanity: G_(1/64)(r=1) = {G(1/64, 1):.8f},  1/(4 pi^2) = {ref:.8f}")
    print()
    print("Every member is reflection positive (Gaussian), and the family")
    print("converges weakly.  The pointwise limit 1/(4 pi^2 r^2) > 0 does not")
    print("decay exponentially, so clustering is lost in the limit.")
    print("=> GL4-C's hypothesis (U_a) -- a uniform-in-a exponential bound --")
    print("   is necessary, not cosmetic.  Risk R28 in REPORT.md.")


if __name__ == "__main__":
    main()
