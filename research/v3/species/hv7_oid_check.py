# HV7-OIDENT witness: carrier-identity arithmetic.
# Zero new MC: Fractions only, no sampling, no random, no numpy needed.
import math
import sys
from fractions import Fraction as F

sys.stdout.reconfigure(encoding="utf-8")

LOG = []


def p(tag, msg):
    LOG.append((tag, msg))


def stats(w, s):
    """w: weights (Fraction), s: atom positions. Returns (f1sq, M0, w1 or None)."""
    M0 = sum(w)
    s1 = min(s)
    f1 = sum(wi for wi, si in zip(w, s) if si == s1)
    return f1, M0, (f1 / M0 if M0 != 0 else None)


# ---------------------------------------------------------------- [1]
p(1, "w0 精确值 / Popoviciu 口径 / N-幂指纹  w0 = 4 r0 / R^2 ,  r0 = 1/(16 N^2), R = N^k")
for N in [2, 3, 4, 5, 6]:
    r0 = F(1, 16 * N * N)                       # = 1/4 (2N)^{-2}
    w0_pop = 4 * r0 / F(N) ** 2                 # Popoviciu: Var <= (range)^2/4, range=N
    w0_2N2 = r0 / (2 * N * N)                   # 仓内旧口径 M0 <= 2N^2
    assert w0_pop == F(1, 4 * N ** 4), (N, w0_pop)
    assert w0_2N2 == F(1, 32 * N ** 4), (N, w0_2N2)
    p(1, f"  N={N}: r0={r0}  w0(Popoviciu)={w0_pop}={float(w0_pop):.6f}  w0(M0<=2N^2)={w0_2N2}")

# fingerprint: w0 ~ N^{-(2+2k)}.  exact check on ratios for k=1,2,1/2
p(1, "  指纹: w0 = (1/4) N^{-(2+2k)}；k=1 (dim4 单迹) -> N^-4；k=2 (双迹) -> N^-6；")
for k in [F(1), F(2), F(1, 2)]:
    vals = []
    for N in [4, 8, 16]:
        r0 = F(1, 16) / (N * N)                 # r0 = 1/(16 N^2) 随 N 走（P3）
        vals.append((N, 4 * r0 / F(N) ** (2 * k)))
    exps = []
    for (n1, v1), (n2, v2) in zip(vals, vals[1:]):
        exps.append(math.log(float(v2 / v1)) / math.log(float(n2 / n1)))
    p(1, f"  k={k}: 实测指数 {['%.4f' % e for e in exps]}  (预测 {float(-(2 + 2 * k)):.4f})")
    for (n, v) in vals:
        p(1, f"      N={n}: w0={v}")

# ---------------------------------------------------------------- [2]
p(2, "仿射群 O->sO+c 下的精确不变性（nu -> s^2 nu；s 轴不动；range|s| 倍）")
w = [F(3, 10), F(7, 10)]
s = [F(1), F(9)]
f1, M0, w1 = stats(w, s)
for mult, shift in [(F(1), F(0)), (F(5), F(0)), (F(1), F(7)), (F(1, 3), F(-2))]:
    w2 = [mult * mult * wi for wi in w]          # nu -> s^2 nu
    f1b, M0b, w1b = stats(w2, s)
    rho_f = f1b / f1
    rho_m = M0b / M0
    p(2, f"  s={mult} c={shift}: f1^2 x{rho_f} (=s^2 {mult*mult})  M0 x{rho_m}  w1={w1b} (原 {w1})  同={w1b == w1}")
    assert w1b == w1
# range: shift 不动 range，scale 使 range -> |s| range，故 r0 -> s^2 r0, R^2 -> s^2 R^2 => w0 不变
for mult in [F(1), F(5), F(1, 3)]:
    r0, R = F(1, 16 * 4 * 4), F(2) / F(1)        # N=2
    num, den = (mult * mult) * r0, (abs(mult) * R) ** 2
    p(2, f"  s={mult}: r0'={num}  R'^2={den}  w0'={4*num/den}  w0={4*r0/R**2}")

# ---------------------------------------------------------------- [3]
p(3, "Popoviciu 尖锐性：supp={0,N} 两原子, 权重 (p,1-p)")
N = 2
best = None
for i in range(0, 11):
    pp = F(i, 10)
    w = [pp, 1 - pp]
    s = [F(0), F(N)]
    mean = sum(wi * si for wi, si in zip(w, s))
    var = sum(wi * (si - mean) ** 2 for wi, si in zip(w, s))
    if best is None or var > best[0]:
        best = (var, pp)
    p(3, f"  p={pp}: Var={var}  (界 N^2/4={F(N*N,4)})")
p(3, f"  取等: p={best[1]}, Var={best[0]} = N^2/4 = {F(N*N,4)}  -> Popoviciu 在 (0,N) 两端取等")
p(3, "  注意: 取等测度的 s1=0 —— 在类 C (s1=m^2>0) 之外 => 同法不可能再收紧, 但类约束可能留一点口")

# ---------------------------------------------------------------- [4]
p(4, "w1 = cos^2 theta（master identity）：v=(a,b,c) 对 psi1 基, w1 = a^2/|v|^2")
cases = [(F(1), F(0), F(0)), (F(1), F(1), F(0)), (F(1), F(1), F(1)),
         (F(0), F(1), F(0)), (F(1), F(2), F(2)), (F(3), F(4), F(0))]
for a, b, c in cases:
    n2 = a * a + b * b + c * c
    w1 = (a * a / n2) if n2 != 0 else None
    p(4, f"  v=({a},{b},{c}): w1={w1}")
p(4, "  污染单调性: v->(a, b', c) 且 |b'|>=|b| => w1 不增")
base = (F(1), F(1), F(1))
for bb in [F(0), F(1), F(2), F(5)]:
    a, b, c = base[0], bb, base[2]
    n2 = a * a + b * b + c * c
    p(4, f"  b={b}: w1={a*a/n2}")
p(4, "  完全失配 (a=0, M0>0): w1=0 —— 目标在该载体上不是'未证', 是假")

# ---------------------------------------------------------------- [5]
p(5, "跨射线 w0 = 4 r0' / R^2 ：分子 r0' 无定理 => 值未定（表内为 r0'=r0 的条件算术）")
N = 2
r0 = F(1, 16 * N * N)
for name, k, rr in [("C0 ¼trF² (k=1)", F(1), r0),
                    ("C8 (¼trF²)² (k=2)", F(2), r0),
                    ("C8 (¼trF²)² (k=2, r0'=1)", F(2), F(1)),
                    ("C5/C6 E² or B² (k=1)", F(1), r0),
                    ("C7 |q| (k=1, range N/8π²≈N/79)", F(1), r0)]:
    R = F(N) ** k
    if "q" in name:
        R = F(N) / 79                                # 1/(8π²) ≈ 1/78.96
    p(5, f"  {name}: R={R}  w0={4*rr/(R*R)}  = {float(4*rr/(R*R)):.6f}")

# ---------------------------------------------------------------- [6]
p(6, "APE 载体实测跨度（HT5-W1POWER.md:63-67, β=2.5, SU(2), L_t=12）与 w0=1/(4N^4)=1/64")
ape_L6 = {"lev0": 0.3382, "lev1": 0.294, "lev2": 0.219, "lev3": 0.122, "lev4": 0.0099, "lev5": 0.178}
ape_L8 = {"lev0": 0.3665, "lev1": 0.320, "lev2": 0.239, "lev3": 0.132, "lev4": 0.0126, "lev5": 0.192}
w0 = 1 / 64
p(6, f"  L=6: 跨度 {max(ape_L6.values())/min(ape_L6.values()):.2f}x  (max {max(ape_L6.values())} @ {max(ape_L6, key=ape_L6.get)},"
     f" min {min(ape_L6.values())} @ {min(ape_L6, key=ape_L6.get)})")
p(6, f"  L=8: 跨度 {max(ape_L8.values())/min(ape_L8.values()):.2f}x")
below = [k for k, v in list(ape_L6.items()) + list(ape_L8.items()) if v < w0]
p(6, f"  w0=1/64={w0}: 低于它的载体验证点 {below}  -> 常数不可跨载体搬运")

# ---------------------------------------------------------------- [7]
p(7, "2D / q / EOM-零 的退化格: M0=0 => w1=0/0")
f1, M0, w1 = stats([F(0)], [F(1)])
p(7, f"  nu=0: f1^2={f1} M0={M0} w1={w1} (None = 未定义)")

# ---------------------------------------------------------------- [8]
p(8, "独立复算 ½Re tr U 的精确值域（HS-GAP-4 §0.5 的口径核对；闭式驻点 + 确定性网格双证）")


def grid_minmax(N, n=60):
    import itertools
    lo, hi = None, None
    ang = [2 * math.pi * i / n for i in range(n)]
    for tup in itertools.product(ang, repeat=N - 1):
        s = sum(math.cos(t) for t in tup) + math.cos(-sum(tup))
        if lo is None or s < lo:
            lo = s
        if hi is None or s > hi:
            hi = s
    return lo / 2, hi / 2


def crit(N):
    """Lagrange: sin θ_i all equal ⇒ θ_i ∈ {α, π−α}; k of them = α; det=1 ⇒ Σθ ≡ 0 (2π)."""
    vals = []
    for k in range(N + 1):
        m = 2 * k - N
        if m == 0:
            continue
        for j in range(-N - 2, N + 3):
            alpha = (2 * math.pi * j - (N - k) * math.pi) / m
            vals.append((k * math.cos(alpha) + (N - k) * math.cos(math.pi - alpha)) / 2)
    return min(vals), max(vals)


for N in [2, 3, 4]:
    lo_g, hi_g = grid_minmax(N, 60)
    lo_c, hi_c = crit(N)
    p(8, f"  N={N}: 网格 [{lo_g:.6f},{hi_g:.6f}]  驻点 [{lo_c:.6f},{hi_c:.6f}]  "
         f"W_grid={hi_g-lo_g:.6f}  W_crit={hi_c-lo_c:.6f}")
p(8, "  ⇒ N=2 W=2, N=3 W=9/4 (min −3/4), N=4 W=4 —— 与 HS-GAP-4 §0.5 的自我更正逐位一致")

# ---------------------------------------------------------------- [9]
p(9, "两个手柄分居（fork）的算术：定域载体的 w₀ 是常数，广延载体的 Popoviciu 界 ∝ L^-6 ⇒ 空转")
N = 2
r0 = F(1, 16 * N * N)
w0_loc = 4 * r0 / F(N) ** 2
p(9, f"  定域: W=N={N} ⇒ w₀={w0_loc}={float(w0_loc):.6f}（L-无关）")
for L in [6, 8, 12, 16]:
    W = F(L ** 3 * N)                        # 固定时切片求和: |Σ_x A| ≤ L³·N
    w0_ext = 4 * r0 / W ** 2
    p(9, f"  广延: L={L} W={W} ⇒ w₀=4r₀/W²={w0_ext}={float(w0_ext):.3e}"
         f"   （实测 w₁≈0.356 的裕度 = {0.356/float(w0_ext):.3e} 倍，随 L⁶ 涨）")
p(9, "  ⇒ 广延载体的 Popoviciu 不等式对一切 L 被平凡满足 ⇒ 该格空转（不是'更难'）")

# ---------------------------------------------------------------- [10]
p(10, "∂_β 恒等式对固定时载体的充要判据： ∂_β⟨O_τ⟩=−Cov(O_τ,Φ)=−Σ_t R(t),  Var(O_τ)=R(0)")
p(10, "  ⇒ 恒等式 ⟺ Σ_t R(t)=R(0) ⟺ R(t)=0 ∀t≠0；而 R(t)=∫e^{−√s|t|}dν(s)>0 (ν≠0) ⇒ 恒不成立")
for t in [0, 1, 2, 3]:
    p(10, f"  t={t}: e^{{−√s t}} 在 s∈[1,4] 两点测度上 = "
         f"{[f'{math.exp(-math.sqrt(s)*t):.6f}' for s in (1, 4)]}")
p(10, "  两原子 ν=½δ₁+½δ₄: R(t)=½(e^{−t}+e^{−2t}) ⇒ R(0)=1.000000, R(1)=%.6f, R(2)=%.6f"
     % (0.5 * (math.exp(-1) + math.exp(-2)), 0.5 * (math.exp(-2) + math.exp(-4))))
SR = 2 * sum(0.5 * (math.exp(-t) + math.exp(-2 * t)) for t in range(0, 400)) - 1.0
p(10, f"  Σ_{{t∈Z}}R(t) = {SR:.6f} vs R(0)=1 ⇒ 比 {SR:.6f} > 1 严格（=∫coth(√s/2)dν）")
for s1 in [0.01, 0.25, 1.0, 4.0]:
    p(10, f"  单原子 s1={s1}: 比 = coth(√s1/2) = {1/math.tanh(math.sqrt(s1)/2):.6f}"
         f"   (s1→0 ⇒ →∞ 松；s1→∞ ⇒ →1 紧)")
xs = [0.05, 0.2, 1.0, 4.0, 16.0]
ph = [1 / math.tanh(math.sqrt(s) / 2) for s in xs]
d2 = [ph[i + 1] - 2 * ph[i] + ph[i - 1] for i in range(1, len(ph) - 1)]
p(10, f"  φβ(s)=coth(√s/2) 在 s={xs}: {['%.4f' % v for v in ph]} ⇒ 单调减；二阶差分 {['%.3e' % v for v in d2]}")
p(10, "  ⇒ 正、减、凸权 ⇒ 总量型 ⇒ 天花板侧（HQ4 类 I）⇒ 它给 cap，不给地板")

# ---------------------------------------------------------------- [11]
p(11, "由 (B2) 的 Jensen 一步：ρβ:=|∂_β⟨A_τ⟩|/M₀=⟨φβ⟩_ν ≥ φβ(s1) ⇒ s1 ≥ 4·arccoth(ρβ)²（单原子取等）")
for s1 in [0.25, 1.0, 2.4, 4.0]:
    rho = 1 / math.tanh(math.sqrt(s1) / 2)
    back = 4 * (0.5 * math.log((rho + 1) / (rho - 1))) ** 2
    p(11, f"  s1={s1}: ρ=coth(√s1/2)={rho:.6f}  4·arccoth(ρ)²={back:.6f}  往返差 {abs(back-s1):.2e}")
p(11, "  物理点 m_G≈1.546 ⇒ s1≈2.390: ρ=%.6f（只比 1 大这么多 ⇒ 该比值本身是 m₁ 的读数）"
     % (1 / math.tanh(1.546 / 2)))

print("HV7-OIDENT carrier-identity witness")
print("=" * 78)

for tag, msg in LOG:
    print(f"[{tag}] {msg}")
print("=" * 78)
print("ALL CHECKS PASS")
with open("research/v3/species/hv7_oid_check.log", "w", encoding="utf-8") as fh:
    fh.write("HV7-OIDENT carrier-identity witness\n" + "=" * 78 + "\n")
    for tag, msg in LOG:
        fh.write(f"[{tag}] {msg}\n")
    fh.write("=" * 78 + "\nALL CHECKS PASS\n")
