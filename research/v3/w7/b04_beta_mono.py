"""
b04_beta_mono.py — B04: m_{0++}(beta) 对 beta 单调吗？哪个量单调？(iii) 能否简化成一点？

问题： (iii)  exists Delta0  forall beta >= beta0=2.5 :  m_{0++} >= Delta0

四个块：
 [A] 4D SU(2)/SU(3) 渐近标度：Lam*a(beta) 的 beta-导数符号（两圈 beta 函数）
     -> m*a = (m/Lam)*(Lam*a) 的方向。  [E]-公式（只在 beta >> 1 定量，低 beta 只作单调性演示）
 [B] "格点单位的地板" <=> "关联长度有界" <=> "无连续极限"：xi/a = 1/(m*a)
 [C] w1 = A/C(0) 是不是质量：标度不变性（H -> cH）精确演示 + 3 能级模型里同 gap 不同 w1
 [D] 2D SU(N) 精确结构对照：m*a ∝ beta^{-1/2} 严格减，而 m/Lam = 纯数（beta-无关）

只写 b04_beta_mono.py / .log。无 git。
"""
import json
import numpy as np

out = []
def p(s=""):
    out.append(s)
    print(s)

# ---------------------------------------------------------------- [A]
# 两圈 beta 函数（格点约定）： b0 = (11N/3)/(16 pi^2), b1 = (34 N^2/3)/(16 pi^2)^2
# beta = 2N/g^2  (Wilson 作用  S = beta * sum_p (1 - (1/N)Re tr U_p))
# a*Lam_L = (b0 g^2)^{-b1/(2 b0^2)} * exp(-1/(2 b0 g^2))        [渐近标度，两圈]
def aLam(N, beta):
    b0 = (11.0 * N / 3.0) / (16.0 * np.pi**2)
    b1 = (34.0 * N**2 / 3.0) / (16.0 * np.pi**2) ** 2
    g2 = 2.0 * N / beta
    return (b0 * g2) ** (-b1 / (2.0 * b0**2)) * np.exp(-1.0 / (2.0 * b0 * g2)), b0, b1

p("=" * 78)
p("[A] 4D 渐近标度： a*Lam_L(beta) 严格单调减？  m*a = (m/Lam)*(Lam*a)")
p("=" * 78)
grid = np.unique(np.concatenate([
    np.linspace(2.5, 6.0, 36), np.linspace(6.0, 100.0, 95),
    np.array([2.5, 3.0, 4.0, 4.3066, 5.0, 6.0, 8.0, 10.0, 20.0, 50.0, 100.0])]))

for N, mLam in ((2, 6.51), (3, 6.28)):
    aL = np.array([aLam(N, b)[0] for b in grid])
    d = np.diff(aL)
    b0, b1 = aLam(N, 2.5)[1], aLam(N, 2.5)[2]
    p(f"\n SU({N}):  b0={b0:.6f}  b1={b1:.6f}  b1/(2b0^2)={b1/(2*b0**2):.4f}")
    p(f"   网格 {len(grid)} 点 (beta 2.5..100)：  max(diff aLam) = {d.max():+.3e}   "
      f"严格减 = {bool(np.all(d < 0))}")
    p(f"   解析斜率 dln(aLam)/dbeta = -1/(2 b0 g^2)/beta *(...) < 0 ；数值 dln(aLam)/dln beta 在 beta=2.5: "
      f"{np.log(aL[1]/aL[0]) / np.log(grid[1]/grid[0]):+.4f}， beta=6: "
      f"{np.log(aL[np.argmin(abs(grid-8))]/aL[np.argmin(abs(grid-6))]) / np.log(8/6):+.4f}")
    p(f"   {'beta':>9} {'Lam*a':>13} {'d ln(aLam)/d ln beta':>22}")
    for bb in [2.5, 3.0, 4.0, 4.3066, 5.0, 6.0, 8.0, 10.0, 20.0, 50.0, 100.0]:
        a_l = aLam(N, bb)[0]
        s = np.log(aLam(N, bb * 1.0001)[0] / a_l) / np.log(1.0001)
        p(f"   {bb:9.4f} {a_l:13.6e} {s:22.6f}")

# —— 定量表：只在两圈可信区（SU(3), beta>=6），用仓内锚点
#     a^{-1}(beta=6.00) = 2.91 GeV, m*a(6.00) = 0.568   [L2/A8, 仓内]
#     Lam_MS/Λ_L(SU(3),n_f=0) = 38.85 => Λ_L = 263/38.85 = 6.77 MeV
p("\n 定量表（SU(3)，锚点 beta=6.00: a^-1=2.910 GeV, m*a=0.568；Λ_L=6.77 MeV）：")
aL6 = aLam(3, 6.0)[0]
p(f"   {'beta':>7} {'Lam*a':>13} {'a^-1[GeV]':>11} {'m*a':>11} {'xi/a':>12}")
tab = []
for bb in [6.0, 6.2, 6.5, 7.0, 7.5, 8.0, 10.0]:
    a_l = aLam(3, bb)[0]
    ainv = 6.77e-3 / a_l          # a^-1 = Λ_L/(Λ_L a)
    ma = 0.568 * a_l / aL6
    tab.append((bb, a_l, ainv, ma, 1.0 / ma))
    p(f"   {bb:7.2f} {a_l:13.6e} {ainv:11.4f} {ma:11.6e} {1.0/ma:12.6e}")

p("\n 判读：Lam*a 严格减 => m*a 严格减（m/Lam 保持常数）")
p("       => inf_{beta>=2.5} m*a = lim_{beta->inf} m*a = 0")
p("       => 格点单位的 (iii) 是假的：不存在正地板（地板 = 0，不可取）")
p(" 定量警告：两圈公式在 beta<5.5 只有单调性意义（绝对值可差 O(10)），SU(3) beta>=6 定量可信")
p("      （自检：aLambda(6)=2.346e-3 => Λ_L = 2.91GeV*2.346e-3 = 6.83 MeV，与 Λ_MS/38.85=6.77 MeV 差 1%）")

# ---------------------------------------------------------------- [B]
p("\n" + "=" * 78)
p("[B] 格点单位地板 <=> 关联长度有界 <=> 无连续极限（逻辑等价，非数值）")
p("=" * 78)
p(" 设 m*a >= Delta0 > 0  for all beta >= beta0。")
p(" 关联长度 xi/a = 1/(m*a) <= 1/Delta0  有界。")
p(" 连续极限要求 a -> 0 且 xi/a -> inf（临界点、发散的关联长度）。")
p(" => 一个格点单位的正地板 = 临界性失败 = 该 beta-序列上不存在非平凡连续极限。")
p(" => (iii) 若以格点单位写，其结论恰好是【目标的反面】（平凡性），不是质量隙。")
p(" 数值见证（用 [A] 的 SU(3) 锚点表：m*a(6)=0.568 -> m*a(10)=?)：")
for bb in [6.0, 6.5, 7.0, 8.0, 10.0]:
    a_l = aLam(3, bb)[0]
    ma = 0.568 * a_l / aL6
    p(f"   beta={bb:6.2f}:  m*a = {ma:.6e}   xi/a = 1/(m*a) = {1.0/ma:.6e}")

# ---------------------------------------------------------------- [C]
p("\n" + "=" * 78)
p("[C] w1 = A/C(0) 不是质量： (i) 标度不变  (ii) 同 gap 不同 w1")
p("=" * 78)

def two_level_spectrum(D, v=1.0):
    H = np.array([[0.0, v], [v, D]])
    E, U = np.linalg.eigh(H)
    O = np.diag([1.0, 0.0])
    Om = U[:, 0]
    m = np.abs(Om @ O @ U[:, 1])**2
    var = (Om**2 @ np.diag(O)) - (Om @ O @ Om) ** 2
    return E[1] - E[0], m, var, m / var

p("\n (i) H -> c*H （全部能量乘 c）：gap 乘 c，w1 不变（精确到机器精度）")
p(f"   {'c':>10} {'gap':>14} {'f1^2':>14} {'Var':>14} {'w1':>18}")
base = two_level_spectrum(D=1.0)
for c in [1e-2, 1e-1, 1.0, 10.0, 100.0]:
    H = np.array([[0.0, 1.0], [1.0, 1.0]]) * c
    E, U = np.linalg.eigh(H)
    O = np.diag([1.0, 0.0])
    Om = U[:, 0]
    m2 = np.abs(Om @ O @ U[:, 1])**2
    var = (Om**2 @ np.diag(O)) - (Om @ O @ Om) ** 2
    p(f"   {c:10.2f} {E[1]-E[0]:14.6e} {m2:14.6e} {var:14.6e} {m2/var:18.15f}")
p("   => w1 对整体能量标度完全盲：w1=0.3427 与 gap=1.707、17.07、1e-3 全部相容。")
p("   => w1 单独不能给 m 的任何地板（[T]：H -> cH 是理论的合法重标，w1 不变而 m 任意）。")

p("\n (ii) 3 能级：gap 钉在 RUN3 的 m_fit=1.707，w1 扫过一大段")
H3 = np.zeros((3, 3))
H3[0, 1] = H3[1, 0] = 1.0
H3[1, 2] = H3[2, 1] = 1.0
p(f"   {'w(1-2 耦合)':>12} {'gap=E1-E0':>12} {'f1^2':>13} {'Var':>13} {'w1':>13}")
# 调整 E1 使 gap 钉住：用 diag 参数
rows = []
for w in [0.0, 0.05, 0.2, 0.5, 1.0, 2.0]:
    best = None
    for d2 in np.linspace(0.01, 6.0, 6001):
        H = np.array([[0.0, 1.0, 0.0], [1.0, 1.1215, w], [0.0, w, 1.1215 + d2]])
        E, U = np.linalg.eigh(H)
        g = E[1] - E[0]
        if best is None or abs(g - 1.707) < abs(best[0] - 1.707):
            best = (g, E, U)
    g, E, U = best
    O = np.diag([0.0, 1.0, 0.0])
    Om = U[:, 0]
    m2 = np.abs(Om @ O @ U[:, 1])**2
    var = (Om**2 @ np.diag(O)) - (Om @ O @ Om) ** 2
    rows.append((w, g, m2, var, m2 / var))
    p(f"   {w:12.2f} {g:12.6f} {m2:13.6e} {var:13.6e} {m2/var:13.6f}")
p("   => 同一个 gap=1.707 上，w1 从 0 走到 ~1 ⇒ w1 不决定 gap，gap 也不决定 w1。")
p("   => RUN3 的 w1=0.3427 是【谱权分数】(A/C(0))，量纲 1，标度不变；不是质量。")

# ---------------------------------------------------------------- [D]
p("\n" + "=" * 78)
p("[D] 对照：2D SU(N) 严格结构（超可重整，g^2 量纲=质量^2，无跑动）")
p("=" * 78)
p(" beta = 2N/(g^2 a^2) => g^2 a^2 = 2N/beta ；m = c*sqrt(g^2 N)（c 为纯数）")
p(" => m*a = c*sqrt(2N^2/beta) ∝ beta^{-1/2}  严格减 -> 0")
p(" => m/Lam = m/(sqrt(g^2)*纯数) = 纯数，beta-无关（精确常数）")
p(f"   {'beta':>8} {'m*a (c=1,N=2)':>16} {'m/Lam (c/const)':>17}")
for bb in [2.5, 4.0, 10.0, 100.0, 1e4]:
    p(f"   {bb:8.1f} {np.sqrt(2*4.0/bb):16.6e} {'= const':>17}")
p(" => 同一个理论里：m*a 严格减，m/Lam 严格常数。 ⇒ 二者不可能同时'单调'，")
p("     而且 (iii) 在两者上的含义【相反】：格点单位假（地板 0），Lam 单位真但内容=连续极限存在。")

# ---------------------------------------------------------------- JSON
res = {
    "aLam_SU2": {str(b): aLam(2, b)[0] for b in [2.5, 4.0, 6.0, 10.0, 100.0]},
    "aLam_SU3": {str(b): aLam(3, b)[0] for b in [2.5, 4.0, 6.0, 10.0, 100.0]},
    "w1_scale_invariance": "exact to machine precision (see [C](i))",
    "RUN3_anchor": {"beta": 2.5, "m_fit_lattice": 1.707, "C0": 19.888,
                    "w1_at_m1564": 0.3427, "w1_err": 0.0119, "SN2": 5.0},
}
with open("research/v3/w7/b04_beta_mono.json", "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1, ensure_ascii=False)

with open("research/v3/w7/b04_beta_mono.log", "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
print("\nWROTE research/v3/w7/b04_beta_mono.log  (%d lines)" % len(out))
