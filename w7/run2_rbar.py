"""RUN2 -- how to estimate Rbar = M1/M0 for the residue-level-gap bridge. NO sampling.

Definitions (verbatim corpus, not invented):
  nu on s = m^2 (m = energy);  M_k = int s^k dnu(s)              ch02-squeeze.tex:12-13
  C(t) = int e^{-m t} dnu(m) = sum_n f_n^2 e^{-m_n t}            PROOF-P3:19
  M0 = C(0) = total spectral weight;  M1 = C''(0) = int m^2 dnu  ch02:12-13
  Rbar := M1/M0 = <m^2>_nu  (mass^2)                             ch05-theorems.tex:426-428
  bridge: f1^2 >= M0 (Lstar - Rbar)/(Lstar - m1^2), Lstar = rho^2 m1^2
  =>  w1 = f1^2/M0 >= (rho^2 - rbar)/(rho^2 - 1),   rbar := Rbar/m1^2   P8-TRUTH:36-38

What is EXACT and zero-cost (no ansatz):
  (R0) C(0) - C(1) = M0 <1 - e^{-m}>                              one line, [T]
  (R1) <m> >= -ln<e^{-m}> = -ln(C1/C0)            (Jensen, -ln convex)      [T]
  (R2) Rbar = <m^2> >= <m>^2 >= [-ln(C1/C0)]^2    (Jensen, x^2 convex)      [T]
What is NOT zero-cost: M1 itself.  P8:39 calls it "C(t) 的时间一阶矩/m_eff 加权";
that functional is int C(t) m^2 dt = M0 <m>  ->  <m> (a MASS, not a mass^2).
Using it as Rbar underestimates rbar by ~m1 and makes the bridge falsely strong.
Measured below.

ASCII-only output.  Reads archives only; nothing is sampled.
"""
import numpy as np
from scipy.optimize import nnls

RHO2 = 2.48          # P8-TRUTH:35 anchor (mu2/mu1, quenched SU(3), [E]/[J])
M1_FIX = 1.50        # U14:54 fixed-m used for the corpus w1 table


def conn(S):
    m = S.mean()
    return np.array([(S * np.roll(S, -t, axis=1)).mean() - m * m
                     for t in range(S.shape[1])])


def neff(C, T):
    """log-derivative effective mass, t=1..T/2-1 (periodic-safe)."""
    out = {}
    for t in range(1, T // 2):
        if min(C[t - 1], C[t], C[t + 1]) > 0:
            r = (C[t - 1] + C[t + 1]) / (2.0 * C[t])
            out[t] = float(np.arccosh(r)) if r >= 1.0 else np.nan
    return out


def nnls_rbar(C, T, grid):
    """non-negative spectral reconstruction on t=0..T/2; Rbar = sum f2 m^2 / sum f2."""
    ts = np.arange(0, T // 2 + 1)
    A = np.exp(-np.outer(ts, grid))
    f2, _ = nnls(A, C[ts])
    w = f2.sum()
    return (float((f2 * grid ** 2).sum() / w) if w > 0 else np.nan), f2, grid


def report(tag, S, levels, beta, L, T, N):
    print(f'\n===== {tag}  beta={beta}  L={L}  T={T}  N={N} =====')
    print(f'  {"s":>2s} {"C0":>8s} {"C1/C0":>7s} {"1-C1/C0":>8s} {"<m>_LB":>7s} '
          f'{"rbar_jLB":>9s} {"bridge@jLB":>11s} | {"P8-lin":>7s} {"rbar(lin)":>9s} '
          f'{"w1(lin)":>8s} | {"Rbar_nnls":>9s} {"rbar_nnls":>9s} {"w1_nnls":>8s}')
    for s in levels:
        C = conn(S[s])
        r1 = C[1] / C[0]
        mLB = -np.log(r1)                       # <m> >= -ln<e^-m>   [T]
        rbar_lb = mLB ** 2 / M1_FIX ** 2        # rbar >= <m>^2/m1^2 [T]
        w1_lb = (RHO2 - rbar_lb) / (RHO2 - 1)
        # P8:39 literal functional: int C(t) m_eff(t)^2 dt / M0  ->  <m>   (dim: mass)
        me = neff(C, T)
        ts = np.array(sorted(me))
        lin = float((C[ts] * np.array([me[t] for t in ts]) ** 2).sum() / C[ts].sum()) \
            if len(ts) else np.nan
        rbar_lin = lin / M1_FIX ** 2            # what you get if you FORGET the extra mass
        w1_lin = (RHO2 - rbar_lin) / (RHO2 - 1)
        Rb, f2, g = nnls_rbar(C, T, np.linspace(1.0, 6.0, 51))
        rbar_n = Rb / M1_FIX ** 2
        w1_n = (RHO2 - rbar_n) / (RHO2 - 1)
        print(f'  {s:2d} {C[0]:8.4g} {r1:7.4f} {1-r1:8.4f} {mLB:7.3f} {rbar_lb:9.3f} '
              f'{w1_lb:11.3f} | {lin:7.3f} {rbar_lin:9.3f} {w1_lin:8.3f} | '
              f'{Rb:9.3f} {rbar_n:9.3f} {w1_n:8.3f}')
    # ansatz spread of the NNLS estimate
    print('  NNLS grid spread (s=0): ', end='')
    C = conn(S[levels[0]])
    for g in (np.linspace(1.2, 3.0, 19), np.linspace(1.0, 4.0, 31),
              np.linspace(1.0, 6.0, 51), np.linspace(1.0, 10.0, 91)):
        Rb, _, _ = nnls_rbar(C, T, g)
        print(f'[{g[0]:.1f},{g[-1]:.0f}] Rbar={Rb:6.2f} rbar={Rb/M1_FIX**2:5.2f}  ',
              end='')
    print()
    # two-component closure: does a (m1, m_UV) pair exist with the corpus w1?
    w1c = 0.577                              # U14:59  beta=2.2 L=6 s=0, fixed m=1.50
    res = conn(S[levels[0]])[1] / conn(S[levels[0]])[0]
    x = res - w1c * np.exp(-M1_FIX)
    print(f'  two-component closure with w1={w1c}: 1-C1/C0 budget = {res:.4f}; '
          f'w1*e^-m1 = {w1c*np.exp(-M1_FIX):.4f} -> '
          + (f'm_UV = {-np.log(x/(1-w1c)):.2f} (feasible)' if x > 0 else
             'INFEASIBLE (w1 too large for the measured C1/C0)'))


def main():
    d = np.load('D:/Dev/math/research/unconditional/f1sq_b2.2_L6.npz')
    lv = [int(x) for x in d['levels']]
    report('archive', {s: d[f'S{s}'] for s in lv}, lv,
           float(d['beta']), int(d['L']), d['S0'].shape[1], int(d['N']))

    for L in (4, 8):
        f = f'D:/Dev/math/research/unconditional/w7/u14_b2.2_L{L}_b2.2.npz'
        try:
            e = np.load(f)
            l2 = [int(x) for x in e['levels']]
            report(f'fresh L={L}', {s: e[f'S{s}'] for s in l2}, l2,
                   float(e['beta']), int(e['L']), e[f'S{l2[0]}'].shape[1], int(e['N']))
        except FileNotFoundError:
            print(f'\n[MISSING] {f}')

    p = 'D:/Dev/math/research/lattice-cpb/pilot_b2.5.npz'
    try:
        e = np.load(p)
        l2 = [int(x) for x in e['smears']]
        print('\n(pilot b2.5: N=40 only -> indicative, NOT a measurement)')
        report('pilot', {s: e[f'O_{s}'] for s in l2}, l2,
               2.5, int(e['L']) if 'L' in e.files else 8,
               e[f'O_{l2[0]}'].shape[1], 40)
    except FileNotFoundError:
        print(f'[MISSING] {p}')

    print('\n===== verdict on Rbar =====')
    print('  [T] definition: Rbar = C\'\'(0)/C(0) = <m^2>_nu, mass^2 (ch05:426-428, ch02:12-13)')
    print('  [T] zero-cost exact: <1-e^-m> = 1-C1/C0 ; <m> >= -ln(C1/C0) ; Rbar >= <m>^2')
    print('  [T] DIRECTION: the bridge needs an UPPER bound on Rbar (Lstar > Rbar).  The')
    print('      zero-cost inequalities are LOWER bounds -> they only cap the headroom')
    print('      rho^2 - rbar, they cannot certify w1 > 0.')
    C0 = conn(np.load('D:/Dev/math/research/unconditional/f1sq_b2.2_L6.npz')['S0'])
    print(f'  [E] degeneracy test (Cauchy-Schwarz, exact): (C1/C0)^2 = {(C0[1]/C0[0])**2:.4f}'
          f' <= C2/C0 = {C0[2]/C0[0]:.4f} -> measure is NON-degenerate, so Rbar > <m>^2.')
    print('  [C] Rbar itself needs an ansatz (spectral reconstruction); no zero-cost value;')
    print('      the NNLS estimate grows without bound with the assumed UV end of the grid')
    print('      (Rbar = 6.8/10/20/52 for m_max = 3/4/6/10) -> Rbar is UV-dominated, i.e. a')
    print('      cutoff quantity, not an IR one.  No (L,a)-uniform premise can be built on it.')
    print('  [T] P8:39 ("time first moment / m_eff weighted") evaluates to <m>, a MASS;')
    print('      feeding it as Rbar divides by m1^2 without the extra power -> rbar too')
    print('      small -> w1 bound too large (column "w1(lin)"). It is NOT Rbar.')


if __name__ == '__main__':
    main()

