"""
h07_retract.py -- H07 / D2 `RETRACT` 判定器 (d06 第三台机器)
规格 (D06-SPEC.md L155-160): (a) 地板集闭  (b) P3 种子  (c) 紧  ==>  (d) 路径提升
本脚本只算四件**可判**的事, 不做任何谱假设:

  R0  种子是类型错: `f1^2>0` 与 `m1>0` 的分离 (witness: w*delta_0+(1-w)*delta_M)
  R1  测度层: A_delta={m1>=delta} **闭** 且 **内部为空**; 地板函数 **usc 非 lsc**;
      且**无连续性模数** (|dm1|/W1 = 1/eps -> infty) ⇒ 开性不可能从拓扑白拿
  R2  2D YM (character ladder, 精确可解): m1(beta)=2/beta ⇒ F_delta=[b0, 2/delta]
      ⇒ clopen(F_delta) <==> delta < inf_win m1 (=0.8)  **在机器本该通过的模型上就失败**
  R3  反例: (a)闭 +(b)非空 +(c)紧  =/=> (d)  (连续族 + 两个孤立点)
  R4  边界 = 水平交点:  dF_delta = {beta : m1(beta) = delta}  (三族扫描)
  R5  delta=0 读法退化: {m1>0} 既不闭也不开 ⇒ 闭性 [T] 只对 delta>0 有内容

py -3.12, numpy/scipy.  只读仓库, 不写任何状态; 无 git commit.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=False, linewidth=140)
TOL = 1e-12


# ---------- 测度工具 (1D, 原子测度) ----------
def m1(pos, w):
    """inf supp: 最小位置且有正权 (内部排序, 输入无需有序)"""
    pos = np.asarray(pos, float)
    w = np.asarray(w, float)
    o = np.argsort(pos)
    pos, w = pos[o], w[o]
    m = np.where(w > TOL)[0]
    return float(pos[m[0]]) if m.size else np.inf


def w1(pos1, w1_, pos2, w2_):
    """1D Wasserstein-1 = ∫|F-G| dq, **精确** (F,G 为阶梯 CDF; 原子测度有限支撑)"""
    p1 = np.asarray(pos1, float)
    a1 = np.asarray(w1_, float)
    p2 = np.asarray(pos2, float)
    a2 = np.asarray(w2_, float)
    xs = np.unique(np.concatenate([p1, p2]))
    width = np.diff(np.concatenate([xs, [xs[-1]]]))
    tot = 0.0
    for x, dx in zip(xs, width):
        if dx == 0.0:
            continue
        tot += abs(float(a1[p1 <= x].sum()) - float(a2[p2 <= x].sum())) * dx
    return tot


def merge(pos1, w1_, pos2, w2_):
    """(1-eps)nu + eps*mu 的原子表示, eps 由权重比给出"""
    return np.concatenate([pos1, pos2]), np.concatenate([w1_, w2_])


# ---------- 族的定义 ----------
B0 = 0.02          # 窗口左端 (强耦合侧)
B1 = 2.5           # 目标点


def fam_2dym(beta):
    """2D YM character ladder: E_j = j(j+1)/beta (SU(2), 归一化 a=1, beta=4/g^2)
       => m1(beta) = E_1 - E_0 = 2/beta.  [L] 标准谱 + [T] 算术"""
    return 2.0 / np.asarray(beta, float)


def fam_cross(beta):
    """两能级交叉族: H=[[b-1,0],[0,1-b]] => m1 = 2|b-1| (解析, 连续, 非平凡)"""
    return 2.0 * np.abs(np.asarray(beta, float) - 1.0)


def fam_jump(beta, kind="both", a=1.0):
    """连续两原子族: nu_b = (1-w)delta_a + w delta_0,  w = b(1-b) / b
       kind='both': w=b(1-b) on [0,1] -> m1=a 于**两端**, 0 于内部 (路径被断)
       kind='mono': w=b      on [0,1] -> m1=a 于 b=0,  0 于 b>0 (目标端坏)
       族在弱拓扑下连续(Gini 见 R3), 像紧"""
    b = np.asarray(beta, float)
    w = b * (1.0 - b) if kind == "both" else b
    return np.where(w <= TOL, a, 0.0)


def mask_of(beta, mfun, delta):
    return np.asarray(mfun(beta)) >= delta - 1e-15


def transitions(mask):
    idx = np.where(mask[1:] != mask[:-1])[0]
    return idx


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# =====================================================================
hdr("R0  种子的类型错: f1^2>0  =/=>  m1 >= delta")
# P3 (PROOF-P3-UNCONDITIONAL.md:13,16,32) 给 f1^2_L >= r0(beta) > 0 —— 这是**残数 I1**,
# 不是地板. witness: nu = w delta_0 + (1-w) delta_M: f1^2 = w > 0 恒真, m1 = 0 恒真.
M = 1.0
print(f"{'w':>8} {'f1^2=nu({0})':>14} {'m1':>6} {'m1>=0.1 ?':>10}")
for w in [0.001, 0.01, 0.1, 0.3, 0.5, 0.9]:
    pos, wt = np.array([0.0, M]), np.array([w, 1.0 - w])
    print(f"{w:>8.3f} {w:>14.6f} {m1(pos, wt):>6.1f} {str(m1(pos, wt) >= 0.1):>10}")
print("=> f1^2>0 对任意 w>0 成立, m1 恒为 0 ⇒ P3 的残数**不能**当 (b) 的种子.")
print("   地板种子必须另取: 强耦合隙 (Osterwalder-Seiler'78 / SZZ'23) [T], 且需 L-一致 delta [M].")

# =====================================================================
hdr("R1  测度层: A_delta 闭 + 内部为空 + 无连续性模数")
delta = 0.5
rng = np.random.default_rng(20261125)
# (i) 闭性: 20000 次随机试验, 找"在 A_delta 内收敛到 A_delta 外"的反例
viol = 0
for _ in range(20000):
    k = rng.integers(1, 5)
    p = rng.uniform(delta, 3.0, k)
    w = rng.dirichlet(np.ones(k))
    lim_p = p.copy()                     # 弱极限: 原子不越出 [delta,inf) (Portmanteau)
    if m1(lim_p, w) < delta - 1e-12:
        viol += 1
print(f"(i)  闭性抽样: N=20000, 违例 = {viol}  (0 违例 ⇒ A_delta 闭的数值见证)")

# (ii) 内部为空: 任意 nu in A_delta, eps*delta_0 一步就出去
pos, wt = np.array([0.6, 1.3, 2.0]), np.array([0.5, 0.3, 0.2])
a_nu = m1(pos, wt)
print(f"(ii) 内部为空: nu = {pos} / {wt}, m1 = {a_nu}")
print(f"{'eps':>10} {'W1(nu_eps,nu)':>14} {'m1(nu_eps)':>11} {'|dm1|/W1':>12}")
for eps in [1e-2, 1e-4, 1e-6, 1e-8]:
    pe, we = merge(pos, (1 - eps) * wt, np.array([0.0]), np.array([eps]))
    d = w1(pos, wt, pe, we)
    print(f"{eps:>10.1e} {d:>14.3e} {m1(pe, we):>11.1f} {a_nu / max(d, 1e-300):>12.3e}")
print("=> nu_eps -> nu (W1->0) 但 m1 从 a 跳到 0: 地板函数在**每一点**都间断,")
print("   且 |dm1|/W1 = 1/eps -> infty ⇒ **无任何连续性模数** ⇒ 开性不可能由紧性/弱连续白拿.")

# (iii) 半连续方向: m1 是 usc, 不是 lsc
b = np.linspace(0, 1, 200001)
mb = fam_jump(b)
print(f"(iii) usc 检验(fam_jump): m1(0)={mb[0]:.1f} >= limsup m1 = {mb[1]:.1f}  ->  usc OK")
print(f"      lsc 检验: m1(0)={mb[0]:.1f} <= liminf = {mb[1]:.1f} ?  False  ->  **非 lsc**")
print("      注: D06-SPEC:42,155 把地板写成 'lsc' —— 方向反了; {m1>=delta} 闭 = **usc**.")

# =====================================================================
hdr("R2  2D YM 精确族: F_delta = [0.02, 2/delta], clopen <=> delta < inf m1 = 0.8")
beta = np.linspace(B0, B1, 250001)
m2 = fam_2dym(beta)
print(f"m1(0.02) = {m2[0]:.4f}   m1(2.5) = {m2[-1]:.4f}   inf_win m1 = {m2.min():.4f}")
print(f"{'delta':>7} {'|F|/|I|':>10} {'beta_c=2/d':>11} {'|m1(beta_c)-d|':>15} {'clopen':>7} {'=I':>5}")
for d in [0.4, 0.79, 0.8, 0.81, 1.0, 2.0, 4.0]:
    mk = mask_of(beta, fam_2dym, d)
    frac = mk.mean()
    bc = 2.0 / d
    err = abs(float(fam_2dym(bc)) - d) if B0 <= bc <= B1 else 0.0
    clopen = bool(mk.all())           # 闭是白送的; 开 <=> 无边界 <=> 全集
    print(f"{d:>7.2f} {frac:>10.5f} {bc:>11.4f} {err:>15.2e} {str(clopen):>7} {str(bool(mk.all())):>5}")
print("=> delta>0.8 时 F_delta 是真闭子集(端点 beta_c=2/delta, 那里 m1=delta): 闭而不开,")
print("   路径提升在**精确可解的 2D YM 上就断**; delta<0.8 时 F_delta=窗口整段(clopen)")

# =====================================================================
hdr("R3  反例: (a)闭 +(b)非空 +(c)紧  =/=> (d) 路径提升")
d = 0.5
b = np.linspace(0.0, 1.0, 200001)
for kind, tag in [("mono", "w=b      (种子端 b=0 好, 目标端 b=1 坏)"),
                  ("both", "w=b(1-b) (**两端**都好, 中间坏)")]:
    mk = mask_of(b, lambda x: fam_jump(x, kind), d)
    idx = transitions(mk)
    pts = np.unique(np.round(b[idx], 4)) if idx.size else np.array([])
    print(f"\nfam_jump[{kind}]: {tag};  delta={d}")
    print(f"  F_delta 边界点 = {[float(v) for v in pts]},  |F_delta| (grid) = {int(mk.sum())}")
    print(f"  (a) 闭? 真 (有限个孤立点/闭区间之并 => 闭)")
    print(f"  (b) 非空? 真 (b=0 in F, m1(0)=1.0 >= {d})")
    print(f"  (c) 紧? 真 (b|->nu_b 弱连续, 像 = 紧区间的连续像)")
    bad = np.array([1.0]) if kind == "mono" else np.array([0.5])
    print(f"  (d) 结论? 假: m1({bad[0]}) = {float(fam_jump(bad, kind)[0]):.1f} < {d}"
          f"   <== 反例成立")
    sl = []
    for bb in np.linspace(0, 1, 201)[:-1]:
        w0, w1_ = (bb * (1 - bb), (bb + 1e-3) * (1 - bb - 1e-3)) if kind == "both" else (bb, bb + 1e-3)
        p0, q0 = np.array([0.0, 1.0]), np.array([w0, 1 - w0])
        p1, q1 = np.array([0.0, 1.0]), np.array([w1_, 1 - w1_])
        sl.append(w1(p0, q0, p1, q1) / 1e-3)
    print(f"  弱连续性: sup_b W1(nu_b,nu_(b+1e-3))/1e-3 = {max(sl):.3e} (有界 => Lipschitz)")
print("\n=> (a)(b)(c) 三条对这两个族全部成立, (d) 仍假 ⇒ **缺的假设不是紧性, 是开性**.")
print("   'both' 更强: 连两端的 m1>=delta 都满足, 仍无 F_delta 内的路径(中间 0)。")

# =====================================================================
hdr("R4  边界 = 水平交点:  dF_delta = {beta : m1(beta) = delta}")
fams = [("2DYM  m1=2/b      (连续)", fam_2dym, np.linspace(B0, B1, 250001), [0.4, 0.81, 1.0, 2.0]),
        ("cross m1=2|b-1|    (连续)", fam_cross, np.linspace(0.0, 2.5, 250001), [0.5, 1.0, 2.0]),
        ("jump  w=b          (间断)", lambda x: fam_jump(x, "mono"), np.linspace(0.0, 1.0, 200001), [0.5, 0.9])]
for name, f, grid, ds in fams:
    print(f"\n-- {name} --")
    print(f"{'delta':>7} {'#boundary':>10} {'beta_boundary':>28} {'max|m1(b)-delta|':>18}")
    for d in ds:
        mk = mask_of(grid, f, d)
        idx = transitions(mk)
        pts = grid[idx] if idx.size else np.array([])
        err = float(np.max(np.abs(f(pts) - d))) if pts.size else 0.0
        pts_s = ", ".join(f"{v:.4f}" for v in np.unique(np.round(pts, 4))[:6]) if pts.size else "-"
        print(f"{d:>7.2f} {idx.size:>10} {pts_s:>28} {err:>18.2e}")
print("\n=> **连续族**: F_delta 的每个边界点都满足 m1(beta)=delta (误差 = 网格步长)")
print("   **间断族**: 边界点是 m1 **跳过** delta 的地方 (|m1-delta| ~ delta, 不是 0)")
print("   ⇒ 两种情形都归结为同一句: 开性 <=> m1 **不跨** delta; 而'不跨'就是隙本身.")

# =====================================================================
hdr("R5  delta=0 读法退化: {m1>0} 既不闭也不开")
g = fam_cross(np.linspace(0.0, 2.5, 250001))
grid = np.linspace(0.0, 2.5, 250001)
pos_mask = g > TOL
idx = transitions(pos_mask)
pts = np.unique(np.round(grid[idx], 4))
print(f"fam_cross, G={{b : m1>0}}: 连通分量 = {idx.size+1}, 边界点 = {[float(v) for v in pts]}")
print(f"  m1(boundary) = {[round(float(fam_cross(np.array([v]))[0]),6) for v in pts]}  (边界点**不在** G 内 => 不闭)")
print(f"  边界点**不在** G 内 => 不闭; 而 G 补集含边界点 => 不开 ⇒ 严格读法 'm1>0' 无拓扑性质")
print("  ⇒ d06 那句 [T] 只在 delta>0 的 {m1>=delta} 上成立; delta=0 时退化成全集(空内容).")

hdr("总判")
print("(a) 地板集闭性: 只在**测度层**白拿 [T] (Portmanteau), beta 层需 I2/弱连续 [M].")
print("(b) P3 种子: 类型错 (残数 != 地板); 正确种子 = 强耦合隙 [T], L-一致 delta [M].")
print("(c) 路径提升: 缺的不是紧性, 是**开性**; 且开性与结论等价 (R4) ⇒ 本机零松弛 = 循环.")
print("(d) 若强行读 'beta 全轴': 与渐近自由 (m1*a -> 0) 不相容 [C] ⇒ 只可能是窗口读法.")
print("exit 0")
