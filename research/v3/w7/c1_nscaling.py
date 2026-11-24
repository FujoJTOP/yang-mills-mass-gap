"""C1 witness: N-scaling of the Shen-Zhu-Zhu (CMP 400, 2023) strong-coupling
Log-Sobolev / Poincare constant K_S, and of the mass-gap correlation-decay rate c_N.

Source: arXiv:2204.12737, Assumption 1 (K_S), Lemma 4.1 (|Hess_S|<=8(d-1)N|beta|),
Theorem 4.2 / eq.(c1) Bakry-Emery, and the remark after the mass-gap Corollary:
    c_N ~ K_S / ( d_g * ( a_{e,e} + 6(d-1) a_{e,ebar} ) ),   d=4, SU(N).
with  d_g = N^2-1,  a_{e,ebar} = N|beta| sqrt(d_g),
      a_{e,e} = 2(d-1) N|beta| ( sqrt(d_g) + 6 sqrt(N) )      [gamma=3sqrt2].
All in the paper's own 't Hooft-scaled normalization S = N beta Re sum tr Q_p.
One job, numpy only.
"""
import numpy as np


def K_S(N, beta, d=4):
    # SU(N): (N+2)/2 - 1 - 8(d-1) N |beta|
    return (N + 2) / 2.0 - 1.0 - 8.0 * (d - 1) * N * abs(beta)


def c_N(N, beta, d=4):
    dg = N * N - 1.0
    a_eebar = N * abs(beta) * np.sqrt(dg)
    a_ee = 2.0 * (d - 1) * N * abs(beta) * (np.sqrt(dg) + 6.0 * np.sqrt(N))
    denom = dg * (a_ee + 6.0 * (d - 1) * a_eebar)
    return K_S(N, beta, d) / denom


def loglog_slope(Ns, ys):
    x = np.log(np.asarray(Ns, float))
    y = np.log(np.asarray(ys, float))
    A = np.vstack([x, np.ones_like(x)]).T
    slope, _ = np.linalg.lstsq(A, y, rcond=None)[0]
    return slope


print("=" * 72)
print("SZZ N-scaling. SU(N), d=4.  K_S > 0  <=>  |beta| < 1/48  (N-independent).")
print("=" * 72)

for beta in [1.0 / 96.0, 0.9 / 48.0]:
    print(f"\n-- beta = {beta:.6f}  (1/48 = {1/48:.6f}) --")
    print(f"{'N':>5} {'K_S(N)':>14} {'K_S/N':>12} {'c_N(N)':>14} {'c_N*N^3':>12}")
    Ns = [2, 3, 4, 5, 10, 20, 50, 100]
    kvals, cvals = [], []
    for N in Ns:
        ks = K_S(N, beta)
        c = c_N(N, beta)
        kvals.append(ks); cvals.append(c)
        print(f"{N:>5} {ks:>14.6e} {ks/N:>12.6e} {c:>14.6e} {c*N**3:>12.6e}")
    # exclude N=2,3 from asymptotics for c_N (O(1/N) corrections)
    big = [N for N in Ns if N >= 10]
    cs = [c_N(N, beta) for N in big]
    ks = [K_S(N, beta) for N in big]
    print(f"  slope d log(K_S)/d log N (N>=10) = {loglog_slope(big, ks):+.4f}  (expect +1)")
    print(f"  slope d log(c_N)/d log N (N>=10) = {loglog_slope(big, cs):+.4f}  (expect -3)")
    print(f"  c_N * N^3 -> {cs[-1]*big[-1]**3:.4f}  (bounded => polynomial, NOT exponential)")

print("\n" + "=" * 72)
print("Exponential test: ratio c_N(N+1)/c_N(N) should approach 1 (poly), not e^-C<1.")
print("=" * 72)
beta = 1.0 / 96.0
N = 2
while N <= 12000:
    r = c_N(2 * N, beta) / c_N(N, beta)
    print(f"  c_N({2*N})/c_N({N}) = {r:.4f}   (poly ~ 2^-3 = 0.125; exp would -> 0 fast)")
    N *= 3

print("\nVerdict data: K_S grows ~ N; c_N decays ~ N^-3; both polynomial.")
