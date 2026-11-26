"""
mu4_lattice.py -- contact term on the lattice: is  mu^2 = delta^3_perp(0) = 1/a^2 ?

Object (RV-2 (star), x on C):   oint_C dy_nu delta^4(x-y) |_{x=y}  =  t_nu * delta^3_perp(0)
Reduction constant (RV-2 sec.4): (sigma/2) * delta_C(x)  =  g^2 N * delta^3_perp(0)
                                 => sigma = 2 g^2 N * delta^3_perp(0)/delta_C      (*)

Fixed physical volume:  L*a = 1  =>  a = 1/L,   L = 4, 6, 8 (then 16, 32 for the trend).

Three discretisations, written down BEFORE any number:
  A  delta^4 -> delta_{x,y}        (Kronecker, dimensionless), oint dy -> plain sum  : K_A = 1
  B  delta^4 -> a^-4 delta_{x,y}   , oint dy -> a * sum                            : K_B = a^-3
  C  B divided by the 1D curve delta delta_C = a^-1  (RV-1 sec.8.4 note: the missing L^-1)
                                                                                  : K_C = a^-2

Run: py -3.12 mu4_lattice.py
"""
import numpy as np

N = 3                                    # SU(N); the invariant in (star) is g^2 N
b0 = 11.0 * N / (48 * np.pi ** 2)        # 11N/48pi^2
b1 = 34.0 * N ** 2 / (3 * (16 * np.pi ** 2) ** 2)


def contact(L):
    a = 1.0 / L                          # fixed physical volume L*a = 1
    K_A = 1.0                            # dimensionless Kronecker
    K_B = a * a ** -4                    # one line measure a, four 1/a from delta^4
    dC = 1.0 / a                         # 1D curve delta on the lattice
    return a, K_A, K_B, K_B / dC


def g2_two_loop(L):
    """a*Lambda = (b0 g^2)^(-b1/(2b0^2)) exp(-1/(2 b0 g^2)) = 1/L  =>  g^2(L), x=1/(2b0g^2)."""
    c = b1 / (2 * b0 ** 2)
    t = 1.0 / L

    def f(x):                            # decreasing in x
        return (2 * x) ** c * np.exp(-x) - t

    lo, hi = 0.2, 80.0
    assert f(lo) > 0 > f(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    x = 0.5 * (lo + hi)
    return 1.0 / (2 * b0 * x), x, t


def G0(L):
    """free massless propagator at coincidence (a=1), zero mode dropped -> 0.154933."""
    c = 4 * np.sin(np.pi * np.arange(L) / L) ** 2
    th2 = (c[:, None, None, None] + c[None, :, None, None]
           + c[None, None, :, None] + c[None, None, None, :])
    th2[0, 0, 0, 0] = np.inf
    return float(np.sum(1.0 / th2) / L ** 4)


print("=" * 78)
print("0. lattice normalisations: a = 1/L, fixed physical volume L*a = 1")
print("-" * 78)
for L in (4, 6, 8):
    a = 1.0 / L
    print(f"  L={L}: a={a:.6f}   delta^4(0)=(1/V)sum_k 1 = {L**4:8.1f} = a^-4 = {a**-4:8.1f}"
          f"   G0(a=1) = {G0(L):8.5f}  -> a^-2 * (0.15493 textbook)")

print()
print("=" * 78)
print("1. coincidence contact term  K = oint_C dy_nu delta^4(x-y)|_{x=y}")
print("-" * 78)
print(f"{'L':>3} {'a':>9} {'K_A':>7} {'K_B=a^-3':>11} {'dC=a^-1':>9} {'K_C=a^-2':>10}")
for L in (4, 6, 8):
    a, KA, KB, KC = contact(L)
    print(f"{L:>3} {a:>9.6f} {KA:>7.4f} {KB:>11.4f} {1/a:>9.4f} {KC:>10.4f}")

print()
print("=" * 78)
print("2. a-power of sigma = 2 g^2 N K   (2 g^2 N = 2 fixed: only the a-power is on trial)")
print("-" * 78)
print(f"{'L':>3} {'a':>9} {'sigma_A':>9} {'sigma_B':>10} {'sigma_C':>10}"
      f" {'sig_B a^3':>10} {'sig_C a^2':>10}")
Ls = np.array([4.0, 6.0, 8.0, 16.0, 32.0])
sB, sC, sA = [], [], []
for L in Ls.astype(int):
    a, KA, KB, KC = contact(int(L))
    sA.append(2 * KA); sB.append(2 * KB); sC.append(2 * KC)
    if L <= 8:
        print(f"{int(L):>3} {a:>9.6f} {2*KA:>9.3f} {2*KB:>10.3f} {2*KC:>10.3f}"
              f" {2*KB*a**3:>10.4f} {2*KC*a**2:>10.4f}")
for nm, s in (("A", sA), ("B", sB), ("C", sC)):
    p = np.polyfit(np.log(Ls), np.log(np.array(s)), 1)[0]   # sigma ~ L^p
    print(f"  sigma_{nm} ~ L^{p:+.4f}   =  a^{{-{p:.4f}}}   [dimension L^{{-{p:.4f}}}]")

print()
print("=" * 78)
print("3. two-loop beta, physical volume fixed in units of Lambda:  a*Lambda = 1/L")
print("-" * 78)
rows = []
print(f"{'L':>3} {'a':>9} {'x=1/(2b0g^2)':>13} {'g^2':>8} {'beta=2N/g^2':>12}"
      f" {'sig a^2=2g^2N':>13} {'R=sig a^2/e^-x':>14} {'(aL)^2':>8} {'2g^2N/(aL)^2':>13}")
for L in Ls.astype(int):
    a = 1.0 / L
    g2, x, t = g2_two_loop(L)
    sa2 = 2 * g2 * N
    R = sa2 * np.exp(x)                  # lead's ratio  sigma a^2 / exp(-1/(2 b0 g^2))
    rows.append((L, a, x, g2, 2 * N / g2, sa2, R, t ** 2, sa2 / t ** 2))
    print(f"{L:>3} {a:>9.6f} {x:>13.6f} {g2:>8.4f} {2*N/g2:>12.4f}"
          f" {sa2:>13.6f} {R:>14.3f} {t**2:>8.5f} {sa2/t**2:>13.2f}")

print()
print("-" * 78)
print("4. is R = sigma a^2 / e^(-1/(2b0g^2)) constant (1/Lambda^2) or running (1/a^2)?")
print("-" * 78)
L3 = Ls[:3]
R3 = np.array([r[6] for r in rows[:3]])
p = np.polyfit(np.log(L3), np.log(R3), 1)[0]
print(f"  R at L=4,6,8 : {R3[0]:.3f}, {R3[1]:.3f}, {R3[2]:.3f}"
      f"   max/min = {R3.max()/R3.min():.4f}  (constant => 1.000);  fit R ~ L^{p:+.3f}")
print(f"  sigma a^2 = 2g^2N itself : {rows[0][5]:.3f}, {rows[1][5]:.3f}, {rows[2][5]:.3f}"
      f"   vs  1/L^2 = (aL)^2 : {rows[0][7]:.5f}, {rows[1][7]:.5f}, {rows[2][7]:.5f}")
print(f"  ratio 2g^2N/(aL)^2 : {rows[0][8]:.1f}, {rows[1][8]:.1f}, {rows[2][8]:.1f}"
      f"   (diverges ~L^2  =>  ',1/Lambda^2' is excluded by rv-2's own constant)")

print()
print("=" * 78)
print("5. summary")
print("-" * 78)
print("  contact term: 4 powers of 1/a from delta^4, minus 1 from the line measure dy  =>  a^-3.")
print("  sigma = 2g^2N * delta^3_perp(0)          ~ a^-3 : [L^-3], not a string tension.")
print("  sigma = 2g^2N * delta^3_perp(0)/delta_C  ~ a^-2 : [L^-2] = 2g^2N/a^2  ->  sigma ~ 1/a^2.")
print("  sigma ~ a^0 would need K dimensionless with no leftover a anywhere: that is a=1 by hand,")
print("  i.e. exactly the statement sigma a^2 = 2g^2N written with the a suppressed.")
print(f"  ratio test: R runs as L^{p:.2f} (not constant)  =>  1/a^2, not 1/Lambda^2.")
