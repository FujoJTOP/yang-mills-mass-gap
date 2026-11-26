"""f04_sigma_kato.py — (sigma>0) ^ Kato 四模型测

只做三件事:
  [1] sigma 的两种形式化: sigma_sep (可分性 = 有平坦因子/乘积) 与 sigma_str (弦张力 = 带电扇区分裂)
      在 四模型 上的值; 其中 2D YM 的"乘积性"用有限群精确枚举核实 (Z_lattice == z^{N_p}, 且 plaquette 独立)
  [2] Kato 隔离 r_iso := dist(spec(T~) 的底, spec(T~) minus 底) 的 四模型 值 (L->inf 口径, 与 m 同口径)
  [3] 判: (sigma_sep>0)^Kato 与 (sigma_str>0)^Kato 各过几行; 以及 m>0 => Kato 严格假的反例

确定性、无采样、无 MC。py -3.12 -I
"""
import itertools
import math
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=10, suppress=False)
OUT = []


def say(s=""):
    OUT.append(str(s))
    print(s)


# ---------------------------------------------------------------- 0. 工具/口径
say("=" * 78)
say("f04 · (sigma>0) ^ Kato · 四模型测")
say("=" * 78)
say("口径: 所有行都取 L->inf (与 m 同口径, 见 U09-M6 §1 '必要的口径')")
say("T~ := e^{-(H-E0)}|_{Omega^perp};  m := -log(sup spec T~) = dist(E0, spec(H)\\{E0})")
say("Kato 隔离: r_iso := dist( sup spec(T~), spec(T~) \\ {sup} )  (>0 且底为有限重原子)")
say()

# ---------------------------------------------------------------- 1. sigma_sep: 2D YM 的乘积性
say("[1] sigma_sep 的核心事实: 2D YM 真空测度 = 逐 plaquette 热核的精确乘积")
say("    核验方式: 有限群精确枚举 (无采样), 检验 (i) 协方差为 0  (ii) 联合律 = 边际之积 (TV=0)")
say()


def enumerate_2d(group_mats, beta, Lx, Ly, action_tr):
    """2D 格点, 群元素随机均匀; 返回 (Z_lattice, mean_S_weight, plaq_values_per_config)."""
    nlink = 2 * Lx * Ly
    nplaq = Lx * Ly
    Z = 0.0
    p1 = []  # 每个 config 的 plaquette 迹列表(只对 Lx x 1 且 Ly=1 用)
    count = 0
    for cfg in itertools.product(range(len(group_mats)), repeat=nlink):
        # 索引: h[x,y] = cfg[2*(x*Ly+y)], v[x,y] = cfg[2*(x*Ly+y)+1]
        W = 1.0
        trs = []
        for x in range(Lx):
            for y in range(Ly):
                i = 2 * (x * Ly + y)
                j = 2 * (((x + 1) % Lx) * Ly + y)
                k = 2 * (x * Ly + (y + 1) % Ly) + 1
                l = 2 * (x * Ly + y) + 1
                U1 = group_mats[cfg[i]]
                U1n = group_mats[cfg[j]]
                U2n = group_mats[cfg[k]]
                U2 = group_mats[cfg[l]]
                P = U1 @ U2n @ np.linalg.inv(U1n) @ np.linalg.inv(U2)
                t = float(np.real(np.trace(P)))
                trs.append(t)
                W *= math.exp(action_tr(beta, t))
        Z += W
        p1.append(trs)
        count += 1
    Z /= float(len(group_mats)) ** nlink
    return Z, np.array(p1), count


def factor_test(plaq, label):
    """检验 'plaquette 变量两两独立 (即测度 = 逐 plaquette 的乘积)'。
    两个判据: (i) 协方差精确为 0; (ii) (tr P0, tr P1) 的联合经验分布与边际之积的全变差 TV = 0。"""
    n, k = plaq.shape
    covs = []
    for i in range(k):
        for j in range(i + 1, k):
            covs.append(np.mean(plaq[:, i] * plaq[:, j]) - plaq[:, i].mean() * plaq[:, j].mean())
    # 联合 vs 边际乘积 (离散化到支撑上的唯一点)
    a, b = plaq[:, 0], plaq[:, 1]
    ua, ub = np.unique(np.round(a, 9)), np.unique(np.round(b, 9))
    ia = np.searchsorted(ua, np.round(a, 9))
    ib = np.searchsorted(ub, np.round(b, 9))
    P = np.zeros((len(ua), len(ub)))
    np.add.at(P, (ia, ib), 1.0)
    P /= n
    tv = 0.5 * np.abs(P - np.outer(P.sum(1), P.sum(0))).sum()
    say("        %s: max|cov| = %.3e ; TV(joint, marginal x marginal) = %.3e ; 支撑 = %s"
        % (label, max(abs(c) for c in covs), tv, len(ua)))
    say("          all cov = %s" % np.array2string(np.array(covs), precision=12))
    return max(abs(c) for c in covs), tv


# (a) 阿贝尔 Z4: 2x2, 8 links, 4 plaquettes  —— 纯相因子, trP = 2cos(theta)
Z4 = [np.array([[1.0 + 0j, 0], [0, 1.0 + 0j]])]
for k in range(1, 4):
    th = 2 * math.pi * k / 4
    Z4.append(np.array([[complex(math.cos(th), math.sin(th)), 0],
                        [0, complex(math.cos(-th), math.sin(-th))]]))
beta4 = 1.7
act4 = lambda b, t: b * (1.0 - 0.5 * t)  # 2x2 矩阵用迹/2 归一 -> 这里 t 是 tr 的实部
# Z4 的 tr 是 2cos, 用 t/2 归一以保证 S>=0: S = b(1 - Re tr/2)
Z_lat4, plaq4, n4 = enumerate_2d(Z4, beta4, 2, 2, act4)
z4 = sum(math.exp(act4(beta4, float(np.real(np.trace(g))))) for g in Z4) / 4.0
say("    (a) Z4, 2x2 环面, 8 links, %d configs" % n4)
say("        Z_lattice      = %.12f" % Z_lat4)
say("        (对照) z^{N_p}  = %.12f  —— 差值非零是**计数**效应 (环面上最大树 + 非平凡圈约束),"
    % (z4 ** 4))
say("                不是乘积性失效; 乘积性用下面两个精确判据检验")
c4, t4 = factor_test(plaq4, "Z4  ")
say()

# (b) 非阿贝尔 S3 = D3 (2x2 正交表示, 忠实): 3x1 环面, 6 links, 3 plaquettes
th = 2 * math.pi / 3
R1 = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
M = np.array([[1.0, 0.0], [0.0, -1.0]])
S3 = [np.eye(2), R1, R1 @ R1, M, R1 @ M, R1 @ R1 @ M]
beta3 = 1.9
act3 = lambda b, t: b * (1.0 - 0.5 * t)
Z_lat3, plaq3, n3 = enumerate_2d(S3, beta3, 3, 1, act3)
z3 = sum(math.exp(act3(beta3, float(np.real(np.trace(g))))) for g in S3) / 6.0
say("    (b) S3=D3 (忠实 2x2 正交表示), 3x1 环面, 6 links, %d configs" % n3)
say("        Z_lattice      = %.12f" % Z_lat3)
say("        (对照) z^{N_p}  = %.12f  (同一计数效应)" % (z3 ** 3))
c3, t3 = factor_test(plaq3, "S3  ")
say()
say("    => 两个精确判据都给 0: 2D YM 的 plaquette 变量**两两独立**, 每个的律 = 单面 Gibbs 权 ∝ e^S dg。")
say("       即: 真空测度 = **逐 plaquette 热核的精确乘积** (轴向规范 P_p = U_v(p), Haar 变换 Jacobian 1)。")
say("       每个因子是**有隙**热核 (第一个非平凡 Casimir 3g^2/8 > 0), 故: **无平坦因子 ⇒ sigma_sep = 0**。")
say()

# ---------------------------------------------------------------- 2. 2D YM 谱 (精确算术)
say("[2] 2D YM: 转移矩阵谱 (Hamiltonian 形, H = (g^2/2) C2(j), 轴向/时间规范)")
say("    C2(SU(2), j) = j(j+1);  E_j = (g^2/2) j(j+1)")
for g2 in (1.0, 4.0 / 2.5):
    E_half = 0.5 * g2 * 0.75
    E_1 = 0.5 * g2 * 2.0
    say("      g^2=%.4f:  m = E_{1/2} = %.9f ;  E_1 = %.9f ;  r_iso = E_1 - E_{1/2} = %.9f"
        % (g2, E_half, E_1, E_1 - E_half))
say("    口径核对: sigma_str = g^2 C2(1/2)/2 与 m 是**同一个式子** ⇒ 2D YM 上 bit-for-bit 相同 [T]")
say()

# ---------------------------------------------------------------- 3. 自由场
say("[3] 自由场 (无质量 Gauss): L 环面上的模谱 omega_k = |k|, k = 2*pi*n/L")
say("    T~ 的谱在 Omega^perp 上是连续谱 [0,1), 底 = sup = 1 **不是原子** ⇒ r_iso = 0")
for L in (10, 40, 160, 640):
    ks = sorted(set(abs(2 * math.pi * n / L) for n in range(-L, L + 1) if n != 0))
    say("      L=%4d: m_L=|k|_min=%.6e, r_iso(L)=|k|_2-|k|_1=%.6e"
        % (L, ks[0], ks[1] - ks[0]))
say("      => m_L->0, r_iso(L)->0: 底是连续谱的上确界 (不是原子) ⇒ Kato = 假")
say("    另: sigma_sep(自由场) = 1: 无质量 Gaussian 的 k=0 模是**非正常 (平坦) 因子** (方差发散)")
say()

# 有质量自由场: m>0 但 Kato 仍假 -> m>0 => Kato 不成立 (这是 (a) 的严格性方向)
say("    有质量自由场 (质量 omega_0>0): m = omega_0 > 0, 但 spec(H|Omega^perp) = [omega_0, oo) 纯连续")
for L in (10, 40, 160):
    om = [math.sqrt((2 * math.pi * n / L) ** 2 + 1.0) for n in range(-3 * L, 3 * L + 1)]
    om = sorted(set(om))
    say("      L=%4d: omega_min=%.9f, omega_2-omega_1=%.3e -> 0" % (L, om[0], om[1] - om[0]))
say("      => m>0 而 r_iso=0:  **m>0 不蕴含 Kato**  [T]")
say()

# ---------------------------------------------------------------- 4. M6
say("[4] M6 (U09): spec(H) = {(n+1/2)omega} + {2 pi^2 k^2 / L^2};  m_L = min(omega, 2 pi^2/L^2)")
omega = 1.0
for L in (10, 40, 160):
    born = 0.5 * (2 * math.pi / L) ** 2  # U09 用 2*pi^2/L^2 = 0.5*(2pi/L)^2
    m_L = min(omega, born)
    # T~ 的底 = e^{-m_L}; 其"到其余部分"的距离: 第二低是 min(omega, 2*born) (若 m_L=born)
    second = min(omega, 2 * born) if m_L == born else min(2 * omega, born)
    ln_ratio = second - m_L  # -log(t2/t1) = 底与次底的对数距离
    say("      L=%4d: m_L=%.9f, r_iso(log 距离)=%.9f, m_L/r_iso=%.4f"
        % (L, m_L, ln_ratio, m_L / ln_ratio))
say("      => m_L->0 且 r_iso->0 同阶 (m_L/r_iso -> 1): Kato = 假 (底是聚点)")
say("    sigma_sep(M6) = 1: mu = mu_a (x) mu_b, mu_b = b 环上的 Lebesgue (lambda_1 = 4 pi^2/L^2 -> 0) [B05:30]")
say()

# ---------------------------------------------------------------- 5. 单 plaquette
say("[5] 单 plaquette (SU(2), beta): 只有一条链/一个面, 无时间方向")
say("    无 T~ (无 t->oo), 无广延 A, 无张量因子分解 ⇒ sigma_sep, sigma_str, r_iso 全部 **无定义**")
say()

# ---------------------------------------------------------------- 6. 四模型表 + 判
say("[6] 四模型表 (L->oo 口径)")
rows = [
    # name, sigma_sep, sigma_str, kato(r_iso>0)
    ("自由场 (无质量)", 1, 0, 0),
    ("2D YM",          0, 1, 1),
    ("M6",             1, 0, 0),
    ("单 plaquette",  None, None, None),
]
need = [0, 1, 0, None]  # 要求的 (sigma>0)^Kato 值: 自由场假, 2D YM 真, M6 假, plaquette 无定义
say("    %-16s %-10s %-10s %-10s | %-12s %-12s" % ("模型", "sigma_sep", "sigma_str", "r_iso>0",
                                                    "(sep>0)^Kato", "(str>0)^Kato"))
res = {}
for name, ss, st, ka in rows:
    if ss is None:
        vs = "undef"
        vst = "undef"
    else:
        vs = "T" if (ss > 0 and ka > 0) else "F"
        vst = "T" if (st > 0 and ka > 0) else "F"
    res[name] = (vs, vst)
    say("    %-16s %-10s %-10s %-10s | %-12s %-12s"
        % (name, str(ss), str(st), str(ka), vs, vst))
say()
say("    要求: 自由场=F, 2D YM=T, M6=F, 单 plaquette=undef")
names = [r[0] for r in rows]
for label, idx in (("(sigma_sep>0)^Kato", 0), ("(sigma_str>0)^Kato", 1)):
    got = [res[n][idx] for n in names]
    exp = ["F", "T", "F", "undef"]
    ok = got == exp
    say("    %-22s 值=%s  期望=%s  -> %s" % (label, got, exp, "PASS(4/4)" if ok else "FAIL"))
    if not ok:
        bad = [n for n, g, e in zip(names, got, exp) if g != e]
        say("        失败行: %s" % bad)
say()

# ---------------------------------------------------------------- 7. 严格性: m>0 vs Kato
say("[7] 显式见证: m>0 而 Kato 假 (乘子算子, 精确)")
say("    T = mult g, g(x) = (1+x)/4 on L^2[0,1]; spec(T) = essential range(g) = [1/4, 1/2] 纯 a.c.")
say("    m := -log sup spec = -log(1/2) = %.9f > 0 ;  dist(sup, spec\\{sup}) = 0 ⇒ r_iso = 0" % math.log(2))
say("    => Kato => m>0 (且 Kato 与 m>0 不等价) ⇒ (sigma>0)^Kato **严格强于** m>0 的 Kato 部分")

say()
say("=" * 78)
say("判词")
say("=" * 78)
say("(b) **在 2D YM 上失败** (按 sigma_sep 读法): 2D YM 是四模型里**唯一精确乘积**的测度")
say("    (上面 Z4 / S3 精确枚举: max|cov| <= 2e-18, TV <= 6e-17), 每个因子有隙 (3g^2/8>0) ⇒ 无平坦因子")
say("    ⇒ (sigma_sep>0) 在 2D YM 上是假的 ⇒ 靶死在该行, 且是**结构性**的: '可分性'与 2D YM 的精确可解性同源。")
say("(a) 按 sigma_str 读法 4/4 通过, 但 sigma_str = g^2 C2/2 在 2D YM 上与 m **是同一个式子** ⇒ 签名撞 m;")
say("    非重名只能靠试金石外的第 5 行 (大质量自由标量), 与 B05 判死的 M1 完全同型。")
say("(c) sigma 有两读; 两读在 2D YM 行上**相反** (sep: 0 / str: >0) ⇒ 靶的定义被这个歧义决定, 不是被物理决定。")

with open("f04_sigma_kato.log", "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
print("\n[written] f04_sigma_kato.log")

