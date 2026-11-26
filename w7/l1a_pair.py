"""l1a_pair.py -- L1 第一步: 在 L=2/SU(2)/beta=2.5 上测 (||Pi_perp H phi_1||^2 , 1-f_1^2) 这一对.

[A] 精确恒等式 (phi_1 与 psi 的两种读法) + 两条推论.  [T] + [E]
[B] 2/3 层精确表: LHS/(1-f_1^2) 到底是什么量 (sigma_bar^2 + f_1^2(lam_bar-lam_1)^2).
[C] L=2 可测性: L_t=2 只剩 1 个时间间隔; f_1^2 是 residue, 需要 >=3.
[D] 自由场在有限 L 上的"原子分数" f_1^2_free(L) -- 酸测: 必须 ->0.
[E] N11 §4 自带的死亡判据"散点无界 => (II) 假": 构造 f_1^2 固定而 LHS/(1-f_1^2) 无界的族.
[F] r04 维数账机检 (192 / 144 / 143 / 1 orbit).

py -3.12 w7/l1a_pair.py > w7/l1a_pair.log
标签 [T]/[E]/[V]. 不写 docs/, 不 git commit.
"""
import itertools
import numpy as np

P = print
P("=" * 96)
P("L1A-PAIR -- 测 ||Pi_perp H phi_1||^2  vs  (1-f_1^2)   [py -3.12, numpy]")

# ---------------------------------------------------------------- [A] identity
P("=" * 96)
P("[A] 恒等式 (两种读法都机检)")
P(" 读法 R-eig : phi_1 = H 在 0^{++} 道的**精确**本征矢,  H phi_1 = lam_1 phi_1")
P("              => Pi_perp H phi_1 = 0 **恒等**,  ||.||^2 = 0.0")
P(" 读法 R-trial: phi_1 = psi := (O-<O>)Omega/||.||,  f_1^2 = |<phi_1|psi>|^2 = p_1")
P("              => 1-f_1^2 = ||Pi_perp psi||^2 = q")
P(" 恒等式(R2): ||Pi_perp H phi_1||^2 = Var_phi1(H) = q * [sigma_bar^2 + f_1^2 (lam_bar - lam_1)^2]")
P(" 注意: Pi_perp 必须相对 phi_1 **自己**; 若 Pi_perp=1-|G><G| 而 phi_1=psi, 得第三种 (R3).")

rng = np.random.default_rng(7)
w1 = w2 = w3 = 0.0
for trial in range(6):
    n = 5 + trial
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    H = (M + M.conj().T) / 2
    w, V = np.linalg.eigh(H)                      # ascending
    G = V[:, 0]
    psi = rng.normal(size=n) + 1j * rng.normal(size=n); psi /= np.linalg.norm(psi)
    p = np.abs(V.conj().T @ psi) ** 2; p /= p.sum()
    f1sq, q = p[0], 1 - p[0]
    lam_bar = (p[1:] @ w[1:]) / q
    sig_bar = (p[1:] @ (w[1:] - lam_bar) ** 2) / q
    # (R1) phi_1 = G, Pi_perp = 1-|G><G|  => 0 identically
    P1 = np.eye(n) - np.outer(G, G.conj())
    r1 = np.linalg.norm(P1 @ (H @ G)) ** 2
    # (R2) phi_1 = psi, Pi_perp = 1-|psi><psi|  => Var_psi(H)
    P2 = np.eye(n) - np.outer(psi, psi.conj())
    r2 = np.linalg.norm(P2 @ (H @ psi)) ** 2
    # (R3) phi_1 = psi, Pi_perp = 1-|G><G|
    r3 = np.linalg.norm(P1 @ (H @ psi)) ** 2
    w1 = max(w1, abs(r1))
    w2 = max(w2, abs(r2 - q * (sig_bar + f1sq * (lam_bar - w[0]) ** 2)))
    w3 = max(w3, abs(r3 - q * (sig_bar + lam_bar ** 2)))
P(f"  [E] 6 组随机 Hermite H (5..10 维), 3 种读法:")
P(f"     (R1) ||(1-|G><G|)H G||^2                 = 0      : max dev = {w1:.2e}")
P(f"     (R2) ||(1-|psi><psi|)H psi||^2 = q[...]           : max dev = {w2:.2e}")
P(f"     (R3) ||(1-|G><G|)H psi||^2     = q(sb+lbar^2)     : max dev = {w3:.2e}")
P("  [T] 证明(R2): Var = sum p_k lam_k^2 - (sum p_k lam_k)^2, p_1=f_1^2, q=1-f_1^2,")
P("      按 (lam_bar, sigma_bar^2) 重组 k>=2 => 交叉项恰为 p_1 q (lam_bar-lam_1)^2.")
P("  [T] 证明(R3): sum_k p_k lam_k^2 - f_1^2 lam_1^2 = q(lam_bar^2 + sigma_bar^2).")

# ---------------------------------------------------------------- [A2] corollaries
P("-" * 96)
P("[A2] 三条推论 (方向)")
P(" (C0) (R1) 读法下 LHS = 0 **精确** => (II) 就是 '0 <~ 1-f_1^2', 恒真、零内容")
P("      = N11-FALLBACK §6 自否证条件 1 的字面触发 ('不等式在 f_1^2->0 时恒真 => 空洞, L1 立死')")
P(" (C1) (R2) 读法: LHS/(1-f_1^2) = sigma_bar^2 + f_1^2 (lam_bar - lam_1)^2 >= 0")
P("      => (II) <=> '激发层的 (2 阶矩 + 均值偏移) <~ 1' (H 以 m_1 为单位)")
P("      <=> 层距的**上**界 (Lambda - m_1 <~ m_1/sqrt(f_1^2)).  I2 要的是**下**界 => 反号.")
P("      (R3) 读法: LHS/(1-f_1^2) = sigma_bar^2 + lam_bar^2 >= Lambda_star^2  => 同样只给上界.")
P(" (C2) 一般原理 [T]: 方差/二阶矩只给'谱支撑位置'的**上**界 (Lambda_star <= sqrt(LHS/q));")
P("      而'边离开底' (I2) 是下界句.  => 任何方差型耦合量都不可能证 I2 (与 N11 §0.1 (S1) 同源).")
P(" (C3) f_1^2 端点: q->0 时 LHS = q[...] -> 0, (II) **最松**(趋于恒真); q->1 时 (II) 最紧.")
P("      L1 要用的恰是 f_1^2->1 那一端 => 不等式在**最需要它的地方最空洞**.")
P(" [E] 例(R2): f_1^2=0.6, Lambda_star=3, m_1=1 => f_1^2(Lambda_star-m_1)^2 = 0.960000")

# ---------------------------------------------------------------- [B] level tables
P("=" * 96)
P("[B] LHS/(1-f_1^2) 的精确值: 2 层 / 3 层 (H/m_1 单位, m_1=1)")
P(f"  {'levels':26s} {'f_1^2':>7s} {'lam/m_1':>22s} {'LHS/(1-f_1^2)':>15s} {'(II) <~1?':>10s}")
rows = [
    ("2 层, Delta/m_1=0.5", [1.0, 1.5], 0.90),
    ("2 层, Delta/m_1=0.5", [1.0, 1.5], 0.50),
    ("2 层, Delta/m_1=0.5", [1.0, 1.5], 0.05),
    ("2 层, Delta/m_1=2.0", [1.0, 3.0], 0.90),
    ("2 层, Delta/m_1=2.0", [1.0, 3.0], 0.50),
    ("2 层, Delta/m_1=2.0", [1.0, 3.0], 0.05),
    ("2 层, Delta/m_1=10.0", [1.0, 11.0], 0.50),
    ("3 层, 1,2,3", [1.0, 2.0, 3.0], 0.50),
    ("3 层, 1,2,3", [1.0, 2.0, 3.0], 0.05),
]
for name, lam, f1 in rows:
    lam = np.array(lam, float)
    n = len(lam)
    if n == 2:
        p = np.array([f1, 1 - f1])
    else:
        p = np.array([f1, (1 - f1) / 2, (1 - f1) / 2])
    q = 1 - p[0]
    lb = (p[1:] @ lam[1:]) / q
    sb = (p[1:] @ (lam[1:] - lb) ** 2) / q
    ratio = sb + p[0] * (lb - lam[0]) ** 2
    ok = "YES" if ratio <= 1 else "NO"
    P(f"  {name:26s} {f1:7.2f} {str(np.round(lam,3)):>22s} {ratio:15.6f} {ok:>10s}")
P(" [T] 2 层闭式: ratio = f_1^2 (Delta/m_1)^2  => (II) <=> Delta <= m_1/sqrt(f_1^2):")
P("     f_1^2->1 时**收紧**, f_1^2->0 时**放松**.  与 'f_1^2->1 时隙开' 的要求反号.")

# ------------------------------------------------- [C] L=2 testability (residue)
P("=" * 96)
P("[C] L=2 能不能测 f_1^2 (residue)?")
P("  f_1^2 = nu({m_1^2}) 是 KL 测度的**原子权重**, 提取需要: C(t) 的 >=3 个时间点 + 多态拟合/GEVP 平台.")
Lt = 2
P(f"  L=2 的转移矩阵时间方向: L_t = {Lt} 个切片 => 可用时间间隔 = L_t-1 = {Lt-1} 个 (t=0,1)")
P("  => 只有 C(0), C(1): 一个比值, 0 个平台 => m_1 与 f_1^2 **不可分离** (自由度不足)")
P("  语料同判: R3-ATOM-FUNCTIONAL:256 / F06-KLATOM:57 => 有限 L 的 'f_1^2>0' 按定义平凡, vacuous")
P("  语料自己给的提取门槛: R4-RESIDUE-PROTOCOL 计划 L=8,12, N=2000, 6--12 CPU-h (L<=8 未解出, M1-4D-RUN)")

# ---------------------------------------- [D] free-field acid test: f_1^2_free(L)
P("=" * 96)
P("[D] 酸测: 自由场(有质量,m=1/2)在 L^3 盒上的'原子分数' f_1^2_free(L)")
P("    O_0 = sum_x :phi_x^2:,  Z(t)=<O_0(t)O_0(0)>_c = sum_k (1/(2 w_k^2)) e^{-2 w_k t}")
P("    f_1^2_free(L) := [k=0 层权重] / [总权重].  自由场必须 -> 0; 在 L=2 上它必须是 O(1) 才叫 vacuous.")
P(f"  {'L':>3s} {'d':>2s} {'m':>4s} {'#levels':>8s} {'E_min/2w':>9s} {'f_1^2_free':>11s}")
for m in (1.0, 0.5):
    for L in (2, 4, 6, 8, 12):
        ks = np.array(list(itertools.product(*[2 * np.pi * np.arange(L) / L] * 3)))
        w2 = m * m + 4 * np.sum(np.sin(ks / 2) ** 2, axis=1)
        w = np.sqrt(w2)
        wt = 1.0 / (2.0 * w2)
        lev = {}
        for E, a in zip(2 * w, wt):
            lev[round(E, 10)] = lev.get(round(E, 10), 0.0) + a
        Es = np.array(sorted(lev)); Ws = np.array([lev[e] for e in Es])
        f1 = Ws[0] / Ws.sum()
        P(f"  {L:3d} {3:2d} {m:4.2f} {len(Es):8d} {Es[0]/2:9.4f} {f1:11.6f}")
P("  [E] L=2: f_1^2_free ~ 0.5 (≈ 总权重的一半)  => 自由场在 L=2 上报出 O(1) 的'原子'")
P("      => L=2 的 f_1^2 测量 **分不开 YM 与自由场** = N11 §6 自否证条件 1 / CLAUDE.md §1 酸测 FAIL")
P("      需要的 L: 让 f_1^2_free 落到靶以下; 表里给出衰减, 语料的提取门槛是 L=8,12.")

# ------------------------------- [E] N11 §4 kill criterion: unbounded scatter
P("=" * 96)
P("[E] N11 §4 自带判据 '散点无界 => (II) 假': 固定 f_1^2, 拉大激发层")
P("    族: nu = f_1^2 delta_{m_1^2} + (1-f_1^2) delta_{Lambda^2}  (两层, 精确可算)")
P(f"  {'f_1^2':>6s} {'Lambda/m_1':>11s} {'LHS/m_1^2':>11s} {'(1-f_1^2)':>10s} {'ratio=LHS/m_1^2/q':>18s}")
for f1 in (0.5, 0.9):
    for Lam in (1.5, 3.0, 10.0, 100.0):
        q = 1 - f1
        lhs = f1 * q * (Lam - 1.0) ** 2        # m_1=1, 2 层: Var = q[f_1^2 (lam_2-m_1)^2]
        P(f"  {f1:6.2f} {Lam:11.1f} {lhs:11.4f} {q:10.2f} {lhs/q:18.4f}")
P("  [T]+[E] ratio = f_1^2 (Lambda-m_1)^2/m_1^2: 对 Lambda **无界**  => 不存在与理论无关的常数 C")
P("      => (II) 作为普适不等式 **假** (N11 §4 判据触发); 作为单个 (L,beta) 的数 = 无内容.")

# ------------------------------------------------------------ [F] r04 bookkeeping
P("=" * 96)
P("[F] |V| r04 维数账机检 (L=2, 4D, SU(2))")
L, D, N = 2, 4, 3
ns, nl = L ** D, D * L ** D
npq = len([(x, mu, nu) for x in itertools.product(range(L), repeat=D)
           for mu in range(D) for nu in range(mu + 1, D)])
P(f"  sites={ns}  links={nl}  link-components={N*nl}  plaquettes={npq}  gauge-gen={N*ns}")
P(f"  物理维数 = {N*nl} - {N*ns} = {N*nl-N*ns} ;  条件化 (mod <v,X>=0) = {N*nl-N*ns-1}")
P(f"  => '143 维' 是**全环面规范商位形空间**的维数, 与 0^{{++}} 道的谱无关")
P(f"  L_t = 2 => 时间间隔 1 个 => 不能测 residue (见 [C])")
P("  LOOP-RUNBOOK:4594 `L=2` 的 [T] 由 r04 判死(143 维, N~3e8, ~35 d 单核) -- 那是 G'' 证书, 不是 f_1^2")
P("=" * 96)
P("[VERDICT] 见 research/v3/species/L1A-PAIR.md")
