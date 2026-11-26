# u09_m6.py -- U09: 解剖 M6 (t03 的 (a) 型反例)
#   1) 显式实现 M6 = 谐振子(a) (x) 自由(b, 周长 L 的环): 核 G''=2w, m=0
#   2) 核 G''(0) = 4(E0 - V(0)); 并核一般 q 处的  G''(q) = 4(E0-V(q)) + 2(psi0'/psi0)^2
#   3) 判"什么量同时满足 (i) M6=0 (ii) 有隙>0 (iii) 自由场=0"
#      -> 测候选:  G''(0) / Gamma0 = lam_min(-D^2 log rho)(0) / Lambda = Poincare 常数
# 归一化: H = -1/2 d^2 + V  (hbar = m0 = 1) => G'' = 4(E0-V(0)); 若改归一化, 因子随 m0/hbar^2 变
# numpy/scipy only, 确定性 (三对角精确对角化 + 稀疏 2D 对角化 + 广义特征问题)

import math
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


def H1d(Vfun, amin, amax, N):
    """H = -1/2 d^2 + V, Dirichlet 盒; 返回 (a, E[:4], U)"""
    a = np.linspace(amin, amax, N)
    h = a[1] - a[0]
    d = 1.0 / h**2 + Vfun(a)
    e = np.full(N - 1, -0.5 / h**2)
    E, U = sla.eigh_tridiagonal(d, e, select="i", select_range=(0, 3))
    return a, E, U


def gpp0(a, psi):
    """G''(0) = -d^2 log(psi^2) at a=0 (3 点差分, 网格含 0)"""
    i0 = int(np.argmin(np.abs(a)))
    assert abs(a[i0]) < 1e-12
    h = a[1] - a[0]
    f = np.log(psi**2)
    return -(f[i0 + 1] - 2 * f[i0] + f[i0 - 1]) / h**2


def gpp_at(a, psi, q):
    """-d^2 log(psi^2)(q), 5 点差分; 同时返回 4(E0-V)+2(psi'/psi)^2 需要的量"""
    h = a[1] - a[0]
    i = int(np.argmin(np.abs(a - q)))
    f = np.log(psi**2)
    d2 = (-f[i - 2] + 16 * f[i - 1] - 30 * f[i] + 16 * f[i + 1] - f[i + 2]) / (12 * h**2)
    return -d2


def poincare_const(a, psi, h):
    """Lambda = inf_f  (1/2) int |f'|^2 d rho / int f^2 d rho ,  rho = |psi|^2
       离散步: K_ij from ½ Σ_edges (f_{i+1}-f_i)^2/h^2 * w_edge ,  M = diag(rho h)
       返回 (Lambda, 0th=0 核验)"""
    rho = psi**2
    rho = rho / (rho.sum() * h)
    w = 0.5 * (rho[1:] + rho[:-1]) * h          # 边的 rho-质量
    c = 0.5 * w / h**2                          # ½ 因子进刚度
    N = len(rho)
    diag = np.zeros(N)
    diag[:-1] += c
    diag[1:] += c
    K = sp.diags([-c, diag, -c], offsets=[-1, 0, 1], format="csc")
    M = sp.diags([rho * h], offsets=[0], format="csc")
    vals = sla.eigh(K.toarray(), M.toarray(), eigvals_only=True,
                    subset_by_index=[0, 2])
    return vals[1], vals  # 最小 = 0 (常数), 第二 = Lambda


def poincare_circle(L, Nb):
    """rho = 1/L 均匀 -> Lambda = ½(2π/L)^2 (Fourier k=1); 核验离散值"""
    b = np.linspace(0, L, Nb, endpoint=False)
    h = b[1] - b[0]
    rho = np.full(Nb, 1.0 / L)
    w = 0.5 * (rho + np.roll(rho, -1)) * h
    c = 0.5 * w / h**2
    K = np.zeros((Nb, Nb))
    for i in range(Nb):
        j = (i + 1) % Nb
        K[i, i] += c[i]
        K[j, j] += c[i]
        K[i, j] -= c[i]
        K[j, i] -= c[i]
    M = np.diag(rho * h)
    vals = sla.eigh(K, M, eigvals_only=True, subset_by_index=[0, 2])
    return vals[1], 0.5 * (2 * math.pi / L) ** 2


def H2d_osc_free(om, L, na, nb, aa=8.0):
    """H = (p_a^2+p_b^2)/2 + om^2 a^2/2 ; a: Dirichlet [-aa,aa]; b: 周期 [0,L)"""
    a = np.linspace(-aa, aa, na)
    ha = a[1] - a[0]
    b = np.linspace(0, L, nb, endpoint=False)
    hb = b[1] - b[0]
    Va = 0.5 * om**2 * a**2
    da = 1.0 / ha**2 + Va
    ea = np.full(na - 1, -0.5 / ha**2)
    Ta = sp.diags([ea, da, ea], offsets=[-1, 0, 1], format="csr")
    db = np.full(nb, 1.0 / hb**2)
    eb = np.full(nb - 1, -0.5 / hb**2)
    Tb = sp.diags([eb, db, eb], offsets=[-1, 0, 1], format="lil")
    Tb[0, nb - 1] = -0.5 / hb**2      # 周期
    Tb[nb - 1, 0] = -0.5 / hb**2
    Tb = Tb.tocsr()
    H = sp.kron(Ta, sp.identity(nb), format="csr") + sp.kron(sp.identity(na), Tb, format="csr")
    E, U = spla.eigsh(H, k=4, which="SA", tol=1e-10)
    idx = np.argsort(E)
    E, U = E[idx], U[:, idx]
    return a, b, ha, hb, E, U[:, 0].reshape(na, nb)


# =====================================================================
P("=" * 80)
P("U09 / M6 解剖    H = -1/2 d^2 + V   (hbar = m0 = 1)   G''(0) = -d^2 log rho |_0")
P("=" * 80)

# ------------------------------------------------------------------ 1) M6 显式实现
P()
P("=" * 80)
P("[1] M6 显式实现:  V(a) = w^2 a^2/2 (谐振子) (+) b 自由 (周期环, 周长 L)")
P("=" * 80)
om = 1.0
P(f"omega = {om}")
P("解析: psi0(a,b) = psi0_HO(a) * L^{-1/2}  (精确可分, b 方向常数)")
P("      rho_a(a) = |psi0_HO(a)|^2 = sqrt(w/pi) exp(-w a^2)   =>  G''(0) = 2w    [T]")
P("      spec(H) = {(n+1/2)w} + {2 pi^2 k^2 / L^2 /(2)}  ->  L->inf: [w/2, inf) 连续")
P("      m_L = min(w, (2 pi/L)^2/2) -> 0        (L->inf 才是 m=0)                [T]")
P()
P(f"{'L':>6} {'Nb':>4} {'na':>5} {'E0(num)':>12} {'E1(num)':>12} {'m_L=E1-E0':>12} "
  f"{'½(2pi/L)²':>12} {'min(w,·)':>10} {'G0(num)':>12} {'G0-2w':>9}")
rows_m6 = []
for (L, nb, na) in ((2.0, 32, 401), (10.0, 32, 401), (40.0, 32, 401)):
    a, b, ha, hb, E, psi = H2d_osc_free(om, L, na, nb)
    m_num = E[1] - E[0]
    bmode = 0.5 * (2 * math.pi / L) ** 2
    rho_a = (psi**2).sum(axis=1) * hb
    rho_a /= rho_a.sum() * ha
    i0 = int(np.argmin(np.abs(a)))
    f = np.log(rho_a)
    G0 = -(f[i0 + 1] - 2 * f[i0] + f[i0 - 1]) / ha**2
    P(f"{L:6.1f} {nb:4d} {na:5d} {E[0]:12.6f} {E[1]:12.6f} {m_num:12.6f} "
      f"{bmode:12.6f} {min(om, bmode):10.6f} {G0:12.6f} {G0-2*om:9.1e}")
    rows_m6.append((L, m_num, min(om, bmode), G0, bmode))
P()
P("读法: G0 对 L 的依赖 = 0 (三个 L 上同一数, 差 <网格误差);  m_L 随 L^{-2} 掉。")
P("      => G'' 是 IR-调节子无关的局部量, m 是 IR-敏感的整体量。这是'无关'的最硬一条。")
# a-边缘 vs 有隙谐振子: 逐点同一函数
a1, E1s, U1 = H1d(lambda x: 0.5 * om**2 * x**2, -8, 8, 4001)
psi_ho = U1[:, 0]
psi_ho = psi_ho / math.sqrt(np.sum(psi_ho**2) * (a1[1] - a1[0]))
rho_ho = psi_ho**2
a2, b2, ha2, hb2, E2d, psi2d = H2d_osc_free(om, 10.0, 401, 32)
rho_a2 = (psi2d**2).sum(axis=1) * hb2
rho_a2 /= rho_a2.sum() * ha2
a3, E3s, U3 = H1d(lambda x: 0.5 * om**2 * x**2, -8, 8, 401)   # 同网格
psi_ho3 = U3[:, 0]
psi_ho3 = psi_ho3 / math.sqrt(np.sum(psi_ho3**2) * (a3[1] - a3[0]))
rho_ho3 = psi_ho3**2
dmax = np.max(np.abs(rho_a2 - rho_ho3))
P(f"[边缘同一性] 同网格 sup|rho_a^M6 - rho^HO(omega=1)| = {dmax:.3e}   (同一函数 => 同一 G'')")
P(f"[同一 G'']  G''_M6 = {gpp0(a2, np.sqrt(np.maximum(rho_a2,1e-300))):.6f},"
  f"  G''_HO = {gpp0(a3, psi_ho3):.6f}  而 m: M6 -> 0, HO = {E3s[1]-E3s[0]:.6f}")

# ------------------------------------------------------------------ 2) 恒等式
P()
P("=" * 80)
P("[2] G''(0) = 4(E0 - V(0)) 的核验 + 一般 q 的恒等式")
P("=" * 80)
P("推导: -1/2 psi'' + V psi = E0 psi  =>  psi''/psi = 2(V-E0)")
P("      G''(q) = -(log psi^2)'' = -2(psi''/psi - (psi'/psi)^2)")
P("             = 4(E0 - V(q)) + 2 (psi'/psi)^2        [恒等, 一般 q]")
P("      对称点 psi'(0)=0  =>  G''(0) = 4(E0 - V(0))    [t03 的式子]")
P("      (若测度带 Jacobian J: rho = J psi^2 => G'' = 4(E0-V) + 2(psi'/psi)^2 - (logJ)'' )")
P()
models = [
    ("HO  w=1",           lambda x: 0.5 * x**2,                    -8, 8, 4001, 0.0),
    ("HO  w=1 + c=3",     lambda x: 0.5 * x**2 + 3.0,              -8, 8, 4001, 3.0),
    ("quartic a^4/4",     lambda x: x**4 / 4,                      -8, 8, 4001, 0.0),
    ("DW a^4/4-a^2/2",    lambda x: x**4 / 4 - x**2 / 2,           -6, 6, 4001, 0.0),
    ("DW 10(a^2-1)^2",    lambda x: 10.0 * (x**2 - 1) ** 2,        -4, 4, 4001, 10.0),
]
P(f"{'model':20} {'V(0)':>8} {'E0':>12} {'4(E0-V0)':>12} {'G0(num)':>12} {'diff':>9} "
  f"{'m':>10} {'0 允许?':>8}")
tab = []
for (nm, Vf, lo, hi, N, v0) in models:
    a, E, U = H1d(Vf, lo, hi, N)
    psi = U[:, 0] * np.sign(U[int(np.argmin(np.abs(a))), 0])
    G0 = gpp0(a, psi)
    pred = 4 * (E[0] - v0)
    P(f"{nm:20} {v0:8.3f} {E[0]:12.6f} {pred:12.6f} {G0:12.6f} {G0-pred:9.1e} "
      f"{E[1]-E[0]:10.6f} {'YES' if E[0] > v0 else 'NO':>8}")
    tab.append((nm, v0, E[0], pred, G0, E[1] - E[0]))
P()
P("关键: quartic a^4/4 的 V 在 0 处是【极小】, 0 在经典【允许】区, 但 G''=4E0>0。")
P("      => G''>0 只等价于 E0>V(0)（原点低于基态能量=允许区/密度峰），与势垒无关。")
P("      DW 10(a^2-1)^2: V(0)=10 > E0 => G''<0（原点在禁区, rho 在 0 处是两峰间的凹）。")
P()
P("[一般 q 的恒等式核验]  DW a^4/4-a^2/2:")
a, E, U = H1d(lambda x: x**4 / 4 - x**2 / 2, -6, 6, 4001)
psi = U[:, 0]
psi *= np.sign(psi[int(np.argmin(np.abs(a)))])
h = a[1] - a[0]
P(f"  {'q':>6} {'G''(q) num':>13} {'4(E0-V(q))':>13} {'2(psi''/psi)^2':>15} {'sum':>13} {'diff':>9}")
for q in (0.0, 0.3, 0.7, 1.0):
    i = int(np.argmin(np.abs(a - q)))
    qg = a[i]                      # 实际网格点 (5 点差分在此点展开)
    ln = np.log(psi)
    d1 = (ln[i + 1] - ln[i - 1]) / (2 * h)
    Gn = gpp_at(a, psi, qg)
    Vq = qg**4 / 4 - qg**2 / 2
    four = 4 * (E[0] - Vq)
    two = 2 * d1**2
    P(f"  {qg:6.3f} {Gn:13.6f} {four:13.6f} {two:15.6f} {four+two:13.6f} {Gn-(four+two):9.1e}")

# ------------------------------------------------------------------ 3) 候选量
P()
P("=" * 80)
P("[3] 三个候选: G''(0) / Gamma0 = lam_min(-D^2 log rho)(0) / Lambda = Poincare 常数")
P("=" * 80)
P("Lambda := inf_f (1/2) int |grad f|^2 d rho / Var_rho(f)   (rho = |psi0|^2)")
P("基态变换: H = psi0 (-L) psi0^{-1} + E0, L = 1/2 Delta - grad(phi).grad  => spec(H)-E0 = -spec(L)")
P("=> Lambda == m  (不是不等式, 是恒等)   [T]")
P()
P("[3a] 1D 有隙模型的 Lambda vs m:")
P(f"  {'model':22} {'m = E1-E0':>12} {'Lambda(FD)':>12} {'diff':>10}")
for (nm, Vf, lo, hi, N) in (("HO w=1", lambda x: 0.5 * x**2, -8, 8, 2001),
                            ("DW a^4/4-a^2/2", lambda x: x**4 / 4 - x**2 / 2, -6, 6, 2001)):
    a, E, U = H1d(Vf, lo, hi, N)
    psi = U[:, 0]
    hh = a[1] - a[0]
    Lam, vals = poincare_const(a, psi, hh)
    P(f"  {nm:22} {E[1]-E[0]:12.6f} {Lam:12.6f} {Lam-(E[1]-E[0]):10.1e}")
P()
P("[3b] 平坦方向 (环): Lambda_b = ½(2π/L)^2, 与 m 的 b-模一致:")
for L in (2.0, 10.0, 40.0):
    Lam_b, exact = poincare_circle(L, 32)
    P(f"  L={L:5.1f}:  Lambda_b(FD) = {Lam_b:.8f}   exact ½(2π/L)^2 = {exact:.8f}"
      f"   rel.err = {abs(Lam_b-exact)/exact:.1e}")
P()
P("[3c] 乘积测度: Lambda = min(Lambda_a, Lambda_b); M6 的 Lambda 与 2D 谱隙对照:")
P(f"  {'L':>6} {'m_L(2D num)':>13} {'min(w,Lambda_b)':>16} {'Lambda_a(HO)':>13} {'diff':>9}")
a1, E1s, U1 = H1d(lambda x: 0.5 * om**2 * x**2, -8, 8, 2001)
La, _ = poincare_const(a1, U1[:, 0], a1[1] - a1[0])
for (L, m_num, _m_an, _G0, bmode) in rows_m6:
    P(f"  {L:6.1f} {m_num:13.6f} {min(om, bmode):16.6f} {La:13.6f} "
      f"{m_num-min(om,bmode):9.1e}")
P(f"  => Lambda_a = {La:.6f} (= w = {om}); M6 的 Lambda = min(w, ½(2π/L)^2) -> 0 (L->inf)")
P()
P("[3d] Gamma0 = lam_min(-D^2 log rho)(0)  (最坏方向局部曲率):")
P("  M6:  a-方向 2w = 2,  b-方向恒 0  =>  Gamma0 = min(2w, 0) = 0      (满足 (i))")
P("  HO:  2w = 2                                                      (>0, 满足 (ii) 该例)")
P("  free: 0                                                          (满足 (iii))")
P("  2D YM (t03 M2): rho = Haar 平坦 => Gamma0 = 0, 但 m = 3g^2/8 > 0  => 不满足 (ii)  [T]*")
P("  => Gamma0 不是靶 (被 M2 一票否决)。")

# ------------------------------------------------------------------ 4) 判决
P()
P("=" * 80)
P("[4] 判决")
P("=" * 80)
P("(i) M6=0 / (ii) 有隙>0 / (iii) 自由场=0 三者同时满足的量: Lambda (= m, 基态变换恒等)。")
P("    但 Lambda ≡ m, 不是新的充分判据; 它是'把证隙改写成 Poincare 不等式'的形式化。")
P("    单点 jet 类 (G'', Gamma0, 任何 rho_a 在 a=0 的有限阶导数) 被'边缘同一性'判死:")
P("    M6 的 a-边缘与有隙谐振子的 a-边缘是同一函数 => 任何 rho_a-泛函给同一值,")
P("    而 (i)(ii) 要求不同值 (0 vs 2w)。[T]")
P("(ii) 局部最坏方向 Gamma0 只修了'平坦方向'的盲点 (M6 -> 0), 但 M2 上 =0 而 m>0。")
P("(iii) G'' 的机制: 它测的是 E0 - V(0), 即'原点高度', 不是谱的曲率。")
P()

with open("u09_m6.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
