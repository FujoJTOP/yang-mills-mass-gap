"""ur4_analyze.py -- read ur4_enum_n6.json, build the m_n table, test convergence of m_n^(1/n),
and cross-check the C03 tree-sum 0.013635 (Bethe vs true lattice trees)."""
import json, sys, math

fn = sys.argv[1] if len(sys.argv) > 1 else "ur4_enum_n6.json"
d = json.load(open(fn))
res = {int(k): v for k, v in d["res"].items()}
N = max(res)

print(f"file={fn}  pool={d['pool']}  R={d['R']}  elapsed={d['elapsed']:.1f}s\n")

hdr = (f"{'n':>2} {'#connected sets':>16} {'trees(|phi|=1)':>15} {'m_n':>14} "
       f"{'m_n^(1/n)':>10} {'m_n/m_(n-1)':>12} {'max|phi|':>8} {'max^(1/n)':>10} {'#classes':>8}")
print(hdr)
print("-" * len(hdr))
prev = None
rows = []
for n in range(1, N + 1):
    r = res[n]
    t = r["dist"].get("1", 0)
    ratio = r["m_n"] / prev if prev else float("nan")
    rows.append((n, r["sets"], t, r["m_n"], r["mroot"], ratio, r["maxabs"], r["maxroot"], r["nclasses"]))
    print(f"{n:>2} {r['sets']:>16} {t:>15} {r['m_n']:>14} {r['mroot']:>10.6f} {ratio:>12.4f} "
          f"{r['maxabs']:>8} {r['maxroot']:>10.6f} {r['nclasses']:>8}")
    prev = r["m_n"]

print("\n-- increments of  m_n^(1/n) --")
print("  " + "  ".join(f"{rows[i][4]-rows[i-1][4]:.4f}" for i in range(1, len(rows))))

print("\n-- animal-count rate  a_n^(1/n)  vs  mass rate  m_n^(1/n) --")
for n, s, t, m, mr, *_ in rows:
    print(f"n={n}: a_n^(1/n)={s**(1/n):.6f}   t_n^(1/n)={t**(1/n):.6f}   "
          f"m_n^(1/n)={mr:.6f}   mean|phi|={m/s:.6f}")

# --- extrapolations of the limit of m_n^(1/n) ---
print("\n-- extrapolation of  lam = lim m_n^(1/n) --")
vals = [(n, mr) for n, s, t, m, mr, *_ in rows]
if len(vals) >= 3:
    # (i) ratio test: m_n/m_(n-1) for the last step
    lam_ratio = rows[-1][5]
    print(f"(i)  ratio test  m_{rows[-1][0]}/m_{rows[-1][0]-1} = {lam_ratio:.4f}")
    # (ii) linear fit v_n = lam - b/n  on the last three points
    (n1, v1), (n2, v2), (n3, v3) = vals[-3], vals[-2], vals[-1]
    # solve for lam,b by least squares on 1/n
    xs = [1.0 / a for a, _ in vals[-3:]]
    ys = [b for _, b in vals[-3:]]
    nn = len(xs)
    sx = sum(xs); sy = sum(ys); sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    b = (nn * sxy - sx * sy) / (nn * sxx - sx * sx)
    a = (sy - b * sx) / nn
    print(f"(ii) affine fit  m_n^(1/n) = lam + b/n  on n={n1},{n2},{n3}:  lam = {a:.3f} , b = {b:.2f}")
    # (iii) same fit on ratio sequence
    rr = [x[5] for x in rows[-3:]]
    xs2 = [1.0 / x[0] for x in rows[-3:]]
    nn = 3; sx = sum(xs2); sy = sum(rr); sxx = sum(x * x for x in xs2); sxy = sum(x * y for x, y in zip(xs2, rr))
    b2 = (nn * sxy - sx * sy) / (nn * sxx - sx * sx)
    a2 = (sy - b2 * sx) / nn
    print(f"(iii) affine fit on ratios m_n/m_(n-1) = lam + b/n :  lam = {a2:.3f} , b = {b2:.2f}")

# --- C03 tree-sum cross-check: Bethe  R = u(1+R)^20  vs true lattice tree sum ---
print("\n-- C03 cross-check: tree sums at u = 0.0104 (beta_W = 1/12  ->  u = beta/(2N^2), N=2) --")
u = 0.0104
tt = [res[n]["dist"].get("1", 0) for n in range(1, N + 1)]


def bethe(u, delta=20, it=200):
    R = u
    for _ in range(it):
        R = u * (1 + R) ** delta
    return R


R_bethe = bethe(u)
part = sum(tt[n - 1] * u ** n for n in range(1, N + 1))
print(f"Bethe fixed point  u(1+R)^20 = {R_bethe:.8f}   (C03/tt1 record 0.01363543)")
print(f"true lattice tree sum  sum_{{n<=6}} t_n u^n = {part:.8f}")
print(f"missing tail estimate: next term  t_7 u^7  with t_7 ~ {tt[-1]*20:.0f} -> {tt[-1]*20*u**7:.2e}")
print(f"ratio true/Bethe = {part/R_bethe:.6f}")

# --- decision numbers for TT2 (4.3): q = lam_Ur * u  vs 1 ---
print("\n-- TT2 (4.3) decision  q = lam_Ur * u  at u = 0.625 (beta = 2.5, N = 2) --")
for n, s, t, m, mr, *_ in rows:
    print(f"n={n}:  m_n^(1/n) * 0.625 = {mr*0.625:.4f}")
print(f"1/u = {1/0.625:.4f}")
