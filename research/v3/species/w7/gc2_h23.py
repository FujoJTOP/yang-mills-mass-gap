# -*- coding: utf-8 -*-
"""GC2-H23 — gamma_ent from B_n (exact ints), H4 gate recompute, H2 de-degeneration.

Pure arithmetic: integers, ln, n-th roots.  No MC, no randomness, no numpy.
Inputs: B_n, a_n from GV2-KAPPAC-ALGO.md:118-125 (exact integers).
"""
import math

B = [1, 20, 514, 14748, 450255, 14305816]
A = [1, 20, 458, 11132, 280175, 7218768]
LN3 = math.log(3.0)
WINDOW1 = 0.5 + 1.0 + LN3          # 2.5986122886681098


def hr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def solve3(M, v):
    M = [row[:] + [v[i]] for i, row in enumerate(M)]
    for c in range(3):
        p = max(range(c, 3), key=lambda r_: abs(M[r_][c]))
        M[c], M[p] = M[p], M[c]
        for r_ in range(3):
            if r_ != c:
                f = M[r_][c] / M[c][c]
                for k in range(c, 4):
                    M[r_][k] -= f * M[c][k]
    return [M[i][3] / M[i][i] for i in range(3)]


def extrap(seq):
    s = [math.log(x) / (i + 1) for i, x in enumerate(seq)]
    out = {}
    out["M1 endpoint (1/n)ln x_n"] = s[5]
    r = (s[5] - s[4]) / (s[4] - s[3])
    out["M2 geom-tail of increments"] = s[5] + (s[5] - s[4]) * r / (1 - r)
    xs = [1.0 / n for n in (3, 4, 5, 6)]
    ys = [s[n - 1] for n in (3, 4, 5, 6)]
    mx = sum(xs) / 4; my = sum(ys) / 4
    b3 = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    out["M3 LSQ g+b/n (n=3..6)"] = my - b3 * mx
    out["M4 exact g+b/n+c/n^2"] = solve3([[1.0, 1.0 / n, 1.0 / n ** 2] for n in (4, 5, 6)],
                                         [s[n - 1] for n in (4, 5, 6)])[0]
    out["M5 exact g n+th ln n+c"] = solve3([[float(n), math.log(n), 1.0] for n in (4, 5, 6)],
                                           [math.log(seq[n - 1]) for n in (4, 5, 6)])[0]
    y6 = [math.log(seq[n - 1] / seq[n - 2]) for n in (4, 5, 6)]
    x6 = [1.0 / n for n in (4, 5, 6)]
    mx6 = sum(x6) / 3; my6 = sum(y6) / 3
    bb = sum((x - mx6) * (y - my6) for x, y in zip(x6, y6)) / sum((x - mx6) ** 2 for x in x6)
    out["M6 ratio LSQ (n=4..6)"] = my6 - bb * mx6
    out["M7 Aitken(s4,s5,s6)"] = s[3] - (s[4] - s[3]) ** 2 / ((s[5] - s[4]) - (s[4] - s[3]))
    x5_, x6_ = 1 / 5, 1 / 6
    out["M8 2pt 1/n (n=5,6)"] = (s[4] * x6_ - s[5] * x5_) / (x6_ - x5_)
    # M9 s_n = g - a n^{-p}: exact 3-point solve, p by bisection
    fp = lambda p: (4.0 ** (-p) - 5.0 ** (-p)) / (5.0 ** (-p) - 6.0 ** (-p))
    target = (s[3] - s[4]) / (s[4] - s[5])
    lo, hi = 0.05, 4.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if fp(mid) < target:
            lo = mid
        else:
            hi = mid
    p_ = 0.5 * (lo + hi)
    a_ = -(s[3] - s[4]) / (4.0 ** (-p_) - 5.0 ** (-p_))
    out["M9 exact g-a/n^p"] = s[5] + a_ * 6.0 ** (-p_)
    return s, out, p_


hr("S1  sequences, roots, ratios  (exact ints -> float log only)")
print(f"{'n':>2} {'B_n':>10} {'sB=(1/n)lnB':>13} {'rB':>10} "
      f"{'a_n':>9} {'sA=(1/n)lna':>13} {'rA':>10} {'B/a':>8}")
for i in range(6):
    n = i + 1
    print(f"{n:>2} {B[i]:>10} {math.log(B[i])/n:>13.9f} "
          f"{(B[i]/B[i-1] if i else float('nan')):>10.6f} {A[i]:>9} "
          f"{math.log(A[i])/n:>13.9f} {(A[i]/A[i-1] if i else float('nan')):>10.6f} "
          f"{B[i]/A[i]:>8.5f}")

hr("S2  extrapolations  [E/M, not T; no existence proof for the limit]")
sB, mB, pB = extrap(B)
sA, mA, pA = extrap(A)
print(f"{'method':<28} {'gamma_B (from B_n)':>19} {'gamma_a (from a_n)':>19}")
for k in mB:
    print(f"{k:<28} {mB[k]:>19.9f} {mA[k]:>19.9f}")
print(f"M9 correction exponent p: B {pB:.4f}, a {pA:.4f}")
vb = sorted(mB.values()); va = sorted(mA.values())
print(f"spread B: [{vb[0]:.6f}, {vb[-1]:.6f}]  median {vb[len(vb)//2]:.6f}")
print(f"spread a: [{va[0]:.6f}, {va[-1]:.6f}]  median {va[len(va)//2]:.6f}")
print(f"per-site normalisation: gamma_a - ln6/6 = {(math.log(A[5])/6 - math.log(6)/6):.6f}")
print(f"                        gamma_B - ln6/6 = {(math.log(B[5])/6 - math.log(6)/6):.6f}")

hr("S3  H4 gate:  rho = exp(gamma - c0 q) <= 1/3  <=>  c0 q >= gamma + ln3")
print(f"ln3 = {LN3:.12f};  quoted CPB-S window 0.5+1+ln3 = {WINDOW1:.12f}")
print(f"{'gamma':>9} {'gate g+L3':>10} {'window .5+g+L3':>14} {'win-2.5986':>10} "
      f"{'qmin':>8} {'q @c0=2.5986':>13} {'rho(c0q=1)':>11}")
for g in [1.0, 2.333406, 2.632032, sB[5], 3.2, vb[-1]]:
    gate = g + LN3; win = 0.5 + g + LN3
    print(f"{g:>9.6f} {gate:>10.6f} {win:>14.6f} {win-WINDOW1:>10.6f} "
          f"{gate/win:>8.6f} {gate/WINDOW1:>13.6f} {math.exp(g-1):>11.5f}")
print("identity: (0.5+g+ln3) - (0.5+1+ln3) = g - 1  (deficit of the quoted window)")
print("structural margin at the *true* window is exactly 0.5 for every g")

hr("S4  H2 de-degenerated: c1 exp(-c0q|C|) <= |w| <= c2 exp(-c0q|C|),  R=c2/c1")
print("AS with c1,c2:  |sum w| >= exp(-c0q|C0|)[c1 - c2 rho/(1-rho)]")
print("need >= (c1/2) exp(...)  <=>  rho <= 1/(1+2R)  <=>  c0 q >= gamma + ln(1+2R)")
print(f"{'R=c2/c1':>8} {'ln(1+2R)':>10} {'gate@g=2.746':>13} {'window@g=2.746':>15} "
      f"{'excess over ln3':>15}")
for R in (1.0, 1.25, 1.5, 2.0, math.e, 3.0, 5.0, 10.0):
    c = math.log(1 + 2 * R)
    print(f"{R:>8.3f} {c:>10.6f} {sB[5]+c:>13.6f} {0.5+sB[5]+c:>15.6f} {c-LN3:>15.6f}")
print(f"max R admissible at c0q = 2.5986, gamma = {sB[5]:.6f}: "
      f"R <= {(math.exp(2.5986-sB[5])-1)/2:.6f}  (needs R<1: impossible)")

hr("S4b unified constant: H3 prefactor P and H2 amplitude ratio R enter as ln(1+2PR)")
print("A_j <= P e^{gamma j} c2 e^{-c0q(|C0|+j)}  =>  gate c0 q >= gamma + ln(1+2PR)")
print(f"{'P':>4} {'R':>5} {'ln(1+2PR)':>10} {'gate@g=2.746':>13} {'window@g=2.746':>15}")
for P in (1.0, 1.5, 2.0):
    for R in (1.0, 2.0, math.e):
        c = math.log(1 + 2 * P * R)
        print(f"{P:>4.1f} {R:>5.2f} {c:>10.6f} {sB[5]+c:>13.6f} {0.5+sB[5]+c:>15.6f}")

hr("S5  2D analogue: fixed-root connected cell sets in Z^2 (cross-check GO1:252)")
level = {frozenset([(0, 0)])}
c2 = [1]
for _ in range(2, 9):
    nxt = set()
    for S in level:
        for (x, y) in S:
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                c = (x + dx, y + dy)
                if c not in S:
                    nxt.add(S | {c})
    level = nxt
    c2.append(len(level))
print("counts n=1..8 :", c2, " (GO1:252 lists 1,4,18,76,315,1296)")
print("roots         :", " ".join(f"{math.log(c2[i])/(i+1):.6f}" for i in range(len(c2))))
print("ratios        :", " ".join(f"{c2[i]/c2[i-1]:.4f}" for i in range(1, len(c2))))
print(f"gamma_2D floor (n=8) = {math.log(c2[7])/8:.6f}")

hr("S6  fiducial rho (why the gate bites)")
for c0, q, g in ((1.0, 1.0, 1.0), (1.0, 1.0, 2.746), (2.5986, 1.0, 1.0),
                 (2.5986, 1.0, 2.746), (4.3446, 1.0, 2.746)):
    print(f"  c0={c0:>7.4f} q={q:.2f} gamma={g:>7.4f}  rho={math.exp(g-c0*q):>12.6f}")
print("\ndone.")
