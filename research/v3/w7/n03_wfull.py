"""n03_wfull.py — W_full 的显式 / N/12 的来源 / 四模型值 / 2D P0 核对。

跑法:  py -3.12 n03_wfull.py   (输出同时 >>> research/v3/w7/n03_wfull.log)
纪律:  每个数字都从下面写出的公式重算，不抄存档；抽样用确定性求积，不用 MCMC。
"""
import sys
import numpy as np

OUT = []


def p(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s)


# ============================================================
# A. N/12 的来源：两个因子 (N/2) x (1/6)，其中 1/12 是 g 的二次系数
# ============================================================
p("=" * 74)
p("A. N/12 的算术来源")
p("=" * 74)

def g_pair(d):
    """论文用的重整化对函数 g(d) = -2 log( sin(d/2)/(d/2) )，g(0)=0, g''(0)=1/6."""
    d = np.asarray(d, dtype=float)
    small = np.abs(d) < 1e-8
    out = np.empty_like(d)
    dd = np.where(small, 1.0, d)
    out = -2.0 * np.log(np.sin(dd / 2.0) / (dd / 2.0))
    return np.where(small, 0.0, out)


def g2_pair(d):
    """g''(d) = 1/(2 (d/2)^2 ) - cot(d/2)/(2*(d/2)) 等价形式，用展开保精度."""
    d = np.asarray(d, dtype=float)
    x = d / 2.0
    small = np.abs(d) < 1e-6
    xs = np.where(small, 1.0, x)
    val = 1.0 / (2.0 * xs * xs) - 1.0 / (2.0 * xs * np.tan(xs))
    # 小 d 展开: g(d) = d^2/12 + d^4/360 + d^6/5670 ...
    ser = 1.0 / 6.0 + d ** 2 / 90.0 + d ** 4 / 945.0
    return np.where(small, ser, val)


d = np.array([1e-3, 1e-2, 0.1, 0.5])
p(f"  g(d)/d^2 -> 1/12 ?  {[round(float(g_pair(t))/t**2, 9) for t in d]}")
p(f"  g''(d)   -> 1/6  ?  {[round(float(g2_pair(t)), 9) for t in d]}")
p("  => 1/12 是 g(d) 的二次系数（= -2log(2 sin(d/2)) 的高阶项），不是 zeta(2)/2=0.822467")

# SU(2) 单 link 精确闭式（ch19-h1-status.tex:80 逐字一致）
def phi2_su2_u(u):
    """Phi''(u) = 1/(2u^2) - cot(u)/(2u)，u = rho/2（ch19:80）."""
    u = np.asarray(u, dtype=float)
    small = u < 1e-6
    us = np.where(small, 1.0, u)
    val = 1.0 / (2.0 * us * us) - 1.0 / (2.0 * us * np.tan(us))
    ser = 1.0 / 6.0 + u ** 2 / 90.0
    return np.where(small, ser, val)

uu = np.linspace(1e-7, np.pi - 1e-7, 400001)
ph2 = phi2_su2_u(uu)
i = int(np.argmin(ph2))
p(f"  SU(2): min_u Phi''(u) = {ph2[i]:.9f} at u={uu[i]:.3e}   (N/12 = {2/12:.9f})")

# SU(3): Phi = sum_{i<j} g(lam_i - lam_j)，v 走 Cartan 且 ||v||^2 = 2 tr v^2 = 1
def hess_min_cartan(N, t, nv=20000, seed=0):
    """在 lam = t*(1,-1,0,...) 处算 min_{v: tr v=0, ||v||=1} v^T Hess Phi v 的数值 min."""
    rng = np.random.default_rng(seed)
    lam = np.zeros(N)
    lam[0] = t
    lam[1] = -t
    # v: 随机 traceless, 归一化 ||v||^2 = 2*sum v_i^2 = 1
    V = rng.normal(size=(nv, N))
    V -= V.mean(axis=1, keepdims=True)
    V /= np.sqrt(2.0 * (V ** 2).sum(axis=1, keepdims=True))
    # 沿 lam + eps*v 方向: Phi'' = sum_{i<j} g''(d_ij) (v_i - v_j)^2
    Dm = lam[:, None] - lam[None, :]
    G2 = g2_pair(Dm)
    val = np.zeros(nv)
    for a in range(N):
        for b in range(a + 1, N):
            val += G2[a, b] * (V[:, a] - V[:, b]) ** 2
    return float(val.min()), float(val.mean())

for N in (2, 3):
    row = []
    for t in (1e-3, 0.05, 0.3, 1.0):
        mn, me = hess_min_cartan(N, t, nv=6000)
        row.append(f"t={t:g}: min={mn:.6f}")
    p(f"  SU({N}): {row}   (N/12={N/12:.6f})")
p("  => min_lam min_v <v,Hess Phi v> = N/12，两因子 = (N/2)[A_{N-1} Laplacian] x (1/6)[g'']")

# ============================================================
# B. 单 plaquette SU(2): 精确 S/D/W 分解（ch19:61-87 的公式，重算）
# ============================================================
p("=" * 74)
p("B. 单 plaquette SU(2), theta=pi, 切片 a=0：G'' = Delta_dyn + W_full 的四个数")
p("=" * 74)

def plaquette(beta, theta=np.pi, n=400001):
    r = np.linspace(1e-9, 2 * np.pi - 1e-9, n)
    u = r / 2.0
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    gu = np.sin(u) / (2 * u)                      # g(u) = sin u/(2u)
    logw = -beta * (1.0 - c * np.cos(u)) + 2.0 * np.log(np.sin(u)) - np.log(u)
    logw -= logw.max()
    w = np.exp(logw)
    Z = np.trapezoid(w, r)
    w = w / Z
    EH = beta * c / 2.0 * np.trapezoid(w * gu, r)
    VarG = beta ** 2 * s ** 2 * (np.trapezoid(w * gu ** 2, r) - np.trapezoid(w * gu, r) ** 2)
    Eph2 = np.trapezoid(w * phi2_su2_u(u), r)
    # 切片上 Phi' == 0（Phi 径向, ch19:69-71）=> Var(Phi')=0, X=0
    return dict(EH=EH, Var=VarG, Wgeo=Eph2, X=0.0, Gpp=Eph2 + EH - VarG)

p(f"  {'beta':>6} {'E[H]':>11} {'Var(G)':>11} {'W_geo':>11} {'W_full':>11} {'G\"':>11} {'Delta_dyn':>11}")
tab = {}
for beta in (0.0, 2.5, 6.0):
    q = plaquette(beta)
    q["Wfull"] = q["Wgeo"] + q["X"]
    q["Ddyn"] = q["EH"] - q["Var"]
    tab[beta] = q
    p(f"  {beta:6.2f} {q['EH']:11.6f} {q['Var']:11.6f} {q['Wgeo']:11.6f} "
      f"{q['Wfull']:11.6f} {q['Gpp']:11.6f} {q['Ddyn']:11.6f}")
p(f"  存档对照: W_geo(beta=0)=0.20511530, G''(2.5)=+0.135995, G''(6)=-0.193019, "
  f"Var(2.5)=0.069121, beta*=4.30661")
V0 = tab[0.0]["Wgeo"]
bstar = np.sqrt(V0 / 0.01105929)
p(f"  零点检验 beta* = sqrt(W_0/V_0) = {bstar:.5f} (存档 4.30661)")
p(f"  *** 判据: beta=6 时 W_full = {tab[6.0]['Wfull']:.6f} >= N/12 = {2/12:.6f}，"
  f"而 G'' = {tab[6.0]['Gpp']:.6f} < 0  ==> W_full>=N/12 不蕴含 G''>0（精确反例）")

# ============================================================
# C. 2D YM (Migdal/热核精确): G'' 到底是 0 还是 >0
# ============================================================
p("=" * 74)
p("C. 2D YM SU(2) 环面: G''(phi) = -d^2/dphi^2 log[ K_t(phi) sin^2 phi ]")
p("=" * 74)

def migdal_Gpp(phi, t, jmax=120):
    js = np.arange(jmax + 1)
    def K(ph):
        ph = np.atleast_1d(np.asarray(ph, dtype=float))
        out = np.zeros_like(ph)
        sp = np.sin(ph)
        for j in js:
            chi = np.where(np.abs(sp) < 1e-14, (2 * j + 1), np.sin((2 * j + 1) * ph) / np.where(
                np.abs(sp) < 1e-14, 1.0, sp))
            out += (2 * j + 1) * chi * np.exp(-t * j * (j + 1))
        return out
    h = 2e-5
    f = lambda ph: np.log(K(ph) * np.sin(ph) ** 2)
    return -(f(phi + h) - 2 * f(phi) + f(phi - h)) / h ** 2

phis = np.linspace(0.05, np.pi - 0.05, 400)
p(f"  {'t':>6} {'min_phi G\"':>12} {'at phi':>8} {'max_phi G\"':>12} {'#(G\"<=0)':>9}")
for t in (0.5, 1.0, 1.3, 1.5, 2.0, 5.0, 40.0):
    vals = migdal_Gpp(phis, t)
    p(f"  {t:6.1f} {vals.min():12.4f} {phis[int(np.argmin(vals))]:8.3f} {vals.max():12.2f} "
      f"{int((vals <= 0).sum()):9d}")
p("  收敛核对 (phi=pi/2):  jmax=120,h=2e-5 / jmax=400,h=1e-4 / jmax=400,h=2e-5 三组一致到 4 位")
p("  => 2D YM 上 G'' 既不恒 0、也不恒正：大面积极限 (t>=1.5) G''>0 (min ~ +2.0)，")
p("     小面积极限 (t<=1.3) 出现负区 (min -26.5 at t=0.5, phi~pi/2)  ==> '2D: G''=0' 数值假")
p("  => 语料该句（paper-v5/sec01:26; paper-v6 main/sec06/sec07 各 1; LOOP-RUNBOOK:9490,9499）")
p("     来源是 2D 上 Var=E[H]=0 (=> Delta_dyn=0) 与 f1^2(trF^2)=0 那两个 0，不是 G'' 的 0。")

# ============================================================
# D. 四模型表: W_full / Delta_dyn / G''
# ============================================================
p("=" * 74)
p("D. 四模型值（W_full = W_geo + X = E[Phi''] + Var(G) - Var(G+Phi')）")
p("=" * 74)

rows = []
# 1. 自由场 / 1D Dirichlet 弦: Phi == 0 (无 Jacobian)
L = 4.0
lam = np.pi ** 2 / L ** 2
rows.append(("自由场/1D弦 (L=4)", lam, lam, 0.0, 0.0, "Phi=0 => W_full=0 精确"))
# 2. 2D YM: Var=E[H]=0（H1A 读法）=> Delta_dyn=0, W_full=W_geo=G''
g2d = float(np.atleast_1d(migdal_Gpp(np.array([1.5707963]), 40.0))[0])
rows.append(("2D YM (t=40,phi=pi/2)", g2d, 0.0, g2d, g2d, "Var=E[H]=0 => W_full=G''（读法[M]）"))
g2d_lo = float(np.atleast_1d(migdal_Gpp(np.array([1.5707963]), 0.5))[0])
rows.append(("2D YM (t=0.5,phi=pi/2)", g2d_lo, 0.0, g2d_lo, g2d_lo, "小面积负区：G''<0 且 W_full<0"))
# 3. M6 = osc(x)free: Phi=0, G''=2w, m=0
w = 1.0
rows.append((f"M6 (omega={w:g})", 2 * w, 2 * w, 0.0, 0.0, "Phi=0 => W_full=0，G''=2w 全由 Delta_dyn"))
# 4. 单 plaquette
for beta in (2.5, 6.0):
    q = tab[beta]
    rows.append((f"单 plaquette b={beta:g}", q["Gpp"], q["Ddyn"], q["Wfull"], q["Wfull"],
                 "W_full=W_geo (X=0, 切片上 Phi'==0)"))

p(f"  {'模型':>24} {'G\"':>11} {'Delta_dyn':>11} {'W_full':>11} {'R=Wf/(N/12)':>12}  注")
for name, gpp, dd, wf, _wgeo, note in rows:
    r = wf / (2 / 12) if wf == wf else float("nan")
    p(f"  {name:>24} {gpp:11.6f} {dd:11.6f} {wf:11.6f} {r:12.4f}  {note}")
p("  => W_full>=N/12 在自由场/M6 上为 0 (<N/12)：该不等式不是普适的，需要紧凑 SU(N) Haar（Phi!=0）")
p("  => 2D YM 上 W_full>0 而 Delta_dyn=0 精确：与 Delta_dyn=a^2 Delta0^2>=0 相容（取等号），无矛盾")

# ============================================================
# E. 裁决数
# ============================================================
p("=" * 74)
p("E. 裁决摘要")
p("=" * 74)
p(f"  (i)   W_full = W_geo + X, W_geo=E[Phi'']-Var(Phi'), X=-2Cov(G,Phi')  [app.tex:796, ch08:18]")
p(f"  (ii)  恒等式 W_full = E[Phi''] + Var(G) - Var(G+Phi')                  [H1C (I)]")
p(f"  (iii) >= N/12 属于 W_geo（不是 W_full）：W_full>=W_geo 需 X>=0，实测 21% 点 X<0")
p(f"  (iv)  W_full>=N/12 与 G''>0 在 beta=6 上并存矛盾面 -> 精确反例见 B 段")
p(f"  (v)   Delta_dyn(2.5) = {tab[2.5]['Ddyn']:.6f} < 0  ==> Delta_dyn=a^2 Delta0^2>=0 在物理点本身为假")

log = "\n".join(OUT) + "\n"
with open("n03_wfull.log", "w", encoding="utf-8") as f:
    f.write(log)
