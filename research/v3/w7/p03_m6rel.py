"""P03-M6REL — M6 上 `G'' = Delta_dyn + W_full` 这条分解式成立吗？(W_full = 0 是关键)

对象（逐字，见 species/P03-M6REL.md）:
  T01:35 / ch05 eq:split :  G''(0) = E[H_vv + Phi''|xi=0] - Var(G_v + Phi'|xi=0) = Delta_dyn + W_geo + X
  app:792  :  [G''] = the coarse curvature V_eff''(k->0)          <- G'' 的定义（不是别的量）
  app:795  :  W_geo = E[Phi''] - Var(Phi'),  Phi = -log J, dU = J dX
  app:796  :  W_full = W_geo + X,  X = -2 Cov(G, Phi')
  U09:50-59:  G''(q) = -d_q^2 log rho(q) = -2 d_q^2 log|psi_0| = 4(E_0 - V(q)) + 2(psi_0'/psi_0)^2
  U09:68   :  带权测度 rho = J|psi_0|^2  =>  G'' = 4(E_0-V) + 2(psi_0'/psi_0)^2 - (log J)''
  U09:23-28:  M6: H = (p_a^2+p_b^2)/2 + (1/2)w^2 a^2,  rho_a = |psi_0^HO|^2 => G''(0) = 2w  [T]

本脚本核 5 件事（全部确定性，无采样）:
 [1] HO 格点基态: G''(0) = 2w  (数值), 且 = 4(E_0 - V(0)) + 2(psi_0'/psi_0)^2
 [2] 同一 M6 上 "Delta_dyn = E[H|0] - Var(G|0)" 的三种读法给出 2w / w^2 / 0 —— 数值对照
 [3] Phi == 0 (M6 无群无 Haar, J=1) => W_geo = 0, X = 0 精确; 自由场同样 W_full = 0
 [4] 恒等式 (i) 的精确性: 对一般 Gibbs 测度 e^{-V} 验证 V_eff'' = <V_aa> - Var(V_a) (非 Gauss)
 [5] 反例: V = a^4/4 (非 Gauss), G''(0) = 1.683 但 E[V''|0] - Var(V'|0) = 0
     => (i) 的 "V" 必须是测度的对数密度, 不能用 Hamiltonian 势 —— 这就是 M6 的口径裂缝

用法: py -3.12 research/v3/w7/p03_m6rel.py  (log 落盘 p03_m6rel.log)
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

NA = 401
AMAX = 8.0
LINES = []


def P(s=""):
    print(s)
    LINES.append(s)


def ground_state(Vfun, n=NA, amax=AMAX):
    """1D 三对角精确对角化: H = -1/2 d^2/da^2 + V(a)。返回 (a, h, E0, psi0)。"""
    h = 2.0 * amax / (n - 1)
    a = -amax + h * np.arange(n)
    diag = 1.0 / h**2 + Vfun(a)
    off = -0.5 / h**2 * np.ones(n - 1)
    w, v = eigh_tridiagonal(diag, off, select="i", select_range=(0, 1))
    psi = v[:, 0] / np.sqrt(h)          # L2 归一
    if psi[n // 2] < 0:
        psi = -psi
    return a, h, w[0], psi


def gpp_fd(a, h, psi, q=0.0):
    """G''(q) = -d_q^2 log|psi|^2, 中心差分; q 必须在格点上。"""
    i = int(round((q + a[0]) / h))
    f = np.log(psi**2)
    return -(f[i + 1] - 2.0 * f[i] + f[i - 1]) / h**2, i


def schrod(a, h, psi, E0, Vfun, q=0.0):
    """G''(q) = 4(E_0 - V(q)) + 2(psi'/psi)^2(q), psi' 中心差分。"""
    i = int(round((q + a[0]) / h))
    r = (psi[i + 1] - psi[i - 1]) / (2.0 * h) / psi[i]
    return 4.0 * (E0 - Vfun(np.array([a[i]]))[0]) + 2.0 * r**2, r


def head(t):
    P("=" * 78)
    P(t)
    P("=" * 78)


# ---------------------------------------------------------------- [1] M6: G''=2w
head("[1] M6 = osc(w) (x) free(b on circle L).  rho_a = |psi_0^HO|^2  ==>  G''(0) = 2w ?")
P(f"grid: N={NA}, a in [-{AMAX},{AMAX}], h={2*AMAX/(NA-1):.4f}, w varies; 精确 E0 = w/2")
P(f"{'w':>6} {'E0(num)':>10} {'w/2':>8} {'G\"(FD)':>10} {'2w':>8} "
  f"{'4(E0-V(0))':>11} {'+2(psi0/psi0)^2':>16} {'|diff|':>9}")
m6 = {}
for wv in (0.5, 1.0, 2.0, 3.0):
    Vh = lambda x, wv=wv: 0.5 * wv**2 * x**2
    a, h, E0, psi = ground_state(Vh)
    g, _ = gpp_fd(a, h, psi)
    s0, r0 = schrod(a, h, psi, E0, Vh)
    m6[wv] = dict(E0=E0, g=g, s=s0)
    P(f"{wv:6.2f} {E0:10.6f} {wv/2:8.4f} {g:10.6f} {2*wv:8.4f} "
      f"{4*(E0-0.0):11.6f} {s0:16.6f} {abs(g-2*wv):9.2e}")
P("")
P("=> [T] G''(0) = 2w 精确（解析: rho_a 是精确 Gaussian）; 数值与 4(E0-V(0)) 同值（psi_0'(0)=0）.")
# 非对称点: (psi_0'/psi_0)^2 项必须出现
wv = 1.0
Vh = lambda x: 0.5 * 1.0**2 * x**2
a, h, E0, psi = ground_state(Vh)
q = 13 * h
gq, _ = gpp_fd(a, h, psi, q)
sq, rq = schrod(a, h, psi, E0, Vh, q)
P(f"   非原点 q={q:.4f}: G''={gq:.6f}, 4(E0-V(q))={4*(E0-0.5*q**2):.6f}, "
  f"(psi'/psi)^2={rq**2:.6f}, 公式和={sq:.6f}  (|diff|={abs(gq-sq):.2e})")
P("   => (star) 式在 q!=0 也需要 2(psi_0'/psi_0)^2; 只在 psi_0'(0)=0 处退化成 4(E0-V(0)).")

# ------------------------------------------- [2] Delta_dyn 的三种读法（同一模型!）
head("[2] 同一 M6 上 Delta_dyn = E[H|0] - Var(G|0) 的读法（H=Hess V, G=d_a V, 条件化在 a=0）")
P("(R1) 有效作用量  V_eff = -log rho_a = w a^2 + const  ==>  H_eff = 2w, G_eff = 2w a")
P("     Var(G_eff|0) = 0 (a=0 纤维上 G_eff == 0)  ==>  Delta_dyn = 2w        [= G'', 定义式重复]")
P("(R2) Hamiltonian 势  V_osc = (1/2)w^2 a^2            ==>  H_osc = w^2, G_osc = w^2 a")
P("     Var(G_osc|0) = 0  ==>  Delta_dyn = w^2   (w=1: 1 != 2 = G'')")
P("(R3) 经典/静态 Gibbs 测度 rho ~ exp(-V_osc)  ==>  Var(a)=1/w^2, Stein 等分 ==> Delta_dyn = 0")
P("")
P(f"{'w':>6} {'G\"=2w':>8} {'R1: 2w':>8} {'R2: w^2':>9} {'R3: 0':>7} "
  f"{'2w-w^2 (R2 残差)':>17} {'w=2 时 R1==R2?':>16}")
for wv, d in m6.items():
    P(f"{wv:6.2f} {d['g']:8.5f} {2*wv:8.4f} {wv**2:9.4f} {0.0:7.4f} "
      f"{2*wv-wv**2:17.4f} {abs(2*wv-wv**2)<1e-12:>16}")
P("")
P("=> 只有 w=2 时 R1 与 R2 相撞。w=1: G''=2.000000, E[H_osc|0]-Var(G_osc|0)=1.000000")
P("   ==> 'G'' = Delta_dyn + W_full' 若用模型自己的 Hamiltonian 势, 残差 = 2w-w^2 != 0,")
P("       而 Phi==0 又强制 W_geo = X = 0. 两个陈述不能同真（除非把 V 换成 V_eff）。")

# ------------------------- [3] Phi == 0 on M6 ==> W_geo = X = 0; 自由场也一样
head("[3] M6 的 measure 是 Lebesgue: J = 1, Phi = -log J = 0  ==> W_geo = X = 0 精确")
rng = np.random.default_rng(0)
n = 200000
Phi = np.zeros(n)          # Phi == 0 identically
Phi_p = np.zeros(n)        # Phi' == 0
Phi_pp = np.zeros(n)       # Phi'' == 0
G = rng.normal(size=n)
W_geo = Phi_pp.mean() - Phi_p.var()
X = -2.0 * np.mean((G - G.mean()) * (Phi_p - Phi_p.mean()))
P(f"Phi==0 (M6, 与自由场同):  W_geo = E[Phi''] - Var(Phi') = {W_geo:+.1e}")
P(f"                          X     = -2 Cov(G,Phi')      = {X:+.1e}")
P(f"                          W_full = W_geo + X          = {W_geo + X:+.1e}")
P("")
P("自由场酸测: 平直 Gaussian 边缘 (w -> 0) ==> G'' -> 0, 且同样 Phi==0 => W_full = 0:")
for wv in (1e-1, 1e-2, 1e-3):
    Vh = lambda x, wv=wv: 0.5 * wv**2 * x**2
    a, h, E0, psi = ground_state(Vh, n=2001, amax=60.0)
    g, _ = gpp_fd(a, h, psi)
    P(f"   w={wv:8.1e}:  G''(0) = {g:.6e}   (2w = {2*wv:.1e})   W_full = 0 (Phi==0)")
P("")
P("=> W_full = 0 在 M6 与自由场上【同时】成立 ==> 该条件对 M6 vs 自由场零分辨力;")
P("   两者的差别全部落在 Delta_dyn 上 (自由场 0 / M6 2w)。")

# --------------------- [4] 恒等式 (i) 对一般 Gibbs 测度的精确性（非 Gauss 检验）
head("[4] 恒等式 (i): 对任意 Gibbs 测度 e^{-V(a,y)} 验证  V_eff''(0) = <V_aa>(0) - Var(V_a)(0)")
lam, gcoupl, clin = 0.5, 0.3, 0.7


def V(a, y):
    return 0.5 * a**2 + 0.5 * y**2 + lam * y**4 + gcoupl * a**2 * y**2 + clin * a * y


yg = np.linspace(-6.0, 6.0, 400001)
dy = yg[1] - yg[0]


def logZ(a):
    w = np.exp(-V(a, yg))
    return np.log(np.trapezoid(w, dx=dy)) if hasattr(np, "trapezoid") else np.log(np.trapz(w, dx=dy))


eps = 1e-3
Veff2 = -(logZ(eps) - 2.0 * logZ(0.0) + logZ(-eps)) / eps**2
w0 = np.exp(-V(0.0, yg))
w0 = w0 / (np.trapezoid(w0, dx=dy) if hasattr(np, "trapezoid") else np.trapz(w0, dx=dy))
Vaa0 = 1.0 + 2.0 * gcoupl * yg**2          # d_a^2 V at a=0
Va0 = clin * yg                            # d_a V = a(1+2 g y^2) + c y -> c y at a=0
EVaa = float((w0 * Vaa0).sum() * dy)
mV = float((w0 * Va0).sum() * dy)
VarVa = float((w0 * (Va0 - mV) ** 2).sum() * dy)
P(f"V = a^2/2 + y^2/2 + {lam} y^4 + {gcoupl} a^2 y^2 + {clin} a y  (非 Gauss, 二维, Var(V_a|0)!=0)")
P(f"   -d_a^2 log Z(a)|_0              = {Veff2:+.8f}")
P(f"   <V_aa|0> - Var(V_a|0)           = {EVaa:+.8f} - {VarVa:.8f} = {EVaa-VarVa:+.8f}")
P(f"   |diff| = {abs(Veff2-(EVaa-VarVa)):.2e}   => [T] 恒等式在 V = 测度对数密度时精确")
P("")

# ------------------- [5] 但 V 换成 Hamiltonian 势, 恒等式立刻不成立（非 Gauss 见证）
head("[5] 反例: V = a^4/4 (单坐标, 测度 rho = |psi_0|^2)。G''(0) = 1.683, 而 E[V''|0]-Var(V'|0) = 0")
Va4 = lambda x: 0.25 * x**4
a, h, E0, psi = ground_state(Va4, n=2001, amax=10.0)
g0, _ = gpp_fd(a, h, psi)
s0, r0 = schrod(a, h, psi, E0, Va4)
EVpp = 0.0        # V''(0) = 3a^2 = 0
VarVp = 0.0       # a=0 纤维是单点 ==> Var = 0
P(f"   G''(0) num         = {g0:+.6f}   (U09:106 报 1.683)")
P(f"   4(E0 - V(0))       = {4*(E0-0.0):+.6f}   (E0 = {E0:.6f})")
P(f"   公式 (star) 和      = {s0:+.6f}")
P(f"   E[V''|0]-Var(V'|0) = {EVpp-VarVp:+.6f}   <== Hamiltonian 势给 0, 与 G'' 差 {abs(g0):.4f}")
P("")
P("=> (i) 里的 V 必须是【测度的对数密度】(Gibbs 作用量), 不是 Hamiltonian 势.")
P("   M6 的对数密度是 V_eff = w a^2, 而模型写下的势是 (1/2)w^2 a^2 —— 二者不同 (w!=2 时).")
P("   这就是 'Delta_dyn = 2w' 与 'Delta_dyn = w^2' 的来源.")
P("")
head("SUMMARY")
P("M6: Phi == 0 (J=1)  =>  W_geo = 0, X = 0, W_full = 0            [T]")
P("M6: G''(0) = 2w = 4(E0 - V(0)) + 2(psi_0'/psi_0)^2|_{q=0}       [T]")
P("M6: Delta_dyn = G'' - W_full = 2w                                [T, 但这是定义式重复]")
P("M6: E[H_osc|0] - Var(G_osc|0) = w^2 != 2w (w!=2)                 [T]")
P("     ==> 'Delta_dyn = 2w' 只在把 V 取成 V_eff = -log rho 时成立;")
P("         用模型自己的 Hamiltonian 势则残差 = 2w - w^2, 而 W_geo=X=0.")
P("     ==> M6 上 Delta_dyn > 0 在两种读法下都成立 (2w 或 w^2), 靶对 M6 失明;")
P("         但 '用 M6 测 Delta_dyn' 是空的: 平直测度下 Delta_dyn 恒等于 G''.")

open("research/v3/w7/p03_m6rel.log", "w", encoding="utf-8").write("\n".join(LINES) + "\n")
print("\n[log] research/v3/w7/p03_m6rel.log")
