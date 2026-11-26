# -*- coding: utf-8 -*-
"""
D03-ML : m(L) 的 L-依赖 —— 精确模型测试（四模型 + 2 个补强模型）
=================================================================
定义（逐字，全文件统一）:
    T_L = 转移矩阵（时间步长 a=1）;  lam_1 > lam_2 > ...  其特征值（按 |.| 排序）
    m(L) := -ln( lam_2 / lam_1 )            [本文件唯一使用的 m 定义]

纪律（T03 §1）:
    0 维模型（单 plaquette）无 transfer 方向 => m 无定义（不是 0，不是 >0）。
    m = 0 一律按 L -> inf 读。

PART A  2D YM (Migdal/精确): 显式搭 L 条 link 环上的 T => m(L) = L * sigma ?
PART B  M6 (谐振 ⊗ 自由): m_L = min(omega, 2 pi^2 / L^2)
PART C  横场 Ising 链 (PBC): (定义谱隙) vs (奇宇称扇区谱隙) —— 两个数
PART D  自由有质量标量场: 周期 vs 扭曲边界
PART E  单 plaquette: m 无定义
PART F  汇总

确定性: 无采样。numpy 2.4.6 / scipy 1.18.1。
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.integrate import quad

LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


# =====================================================================
# PART A — 2D YM
# =====================================================================
# 结构: 单时间片 = 空间环上的 L 条 link; 物理(规范不变)态 = 和乐 P=U_1...U_L 的函数。
# 相邻两个时间片之间有 L 个 plaquette, 故
#       T[U_{t+1}, U_t] = prod_{l=1..L} w(U_{l,t+1}^dag U_{l,t})
# 规范不变子空间上 => 卷积 w^{*L} => 特征值 = hat_w(R)^L
#       m(L) = -ln(lam_2/lam_1) = L * ( -ln( hat_w(R1)/hat_w(R0) ) ) =: L * sigma
# 即 2D YM 的"隙"是**绕环一圈的通量(torelon)**, 能量 = sigma * L, 严�格递增、且 -> inf。

def zN_weight(N, beta, k):
    return np.exp(beta * np.cos(2 * np.pi * k / N))


def part_A():
    P("=" * 74)
    P("PART A — 2D YM: 环上 L 条 link 显式 T; m(L) 是线性递增的 torelon")
    P("=" * 74)
    for N, beta in [(2, 0.8), (3, 1.3)]:
        w = np.array([zN_weight(N, beta, k) for k in range(N)])
        what = np.real(np.array([np.sum(w * np.exp(-2j * np.pi * r * np.arange(N) / N))
                                 for r in range(N)]))
        o = np.argsort(-what)
        lam1, lam2 = what[o[0]], what[o[1]]
        sigma = -np.log(lam2 / lam1)
        P(f"\n[Z_{N}] beta={beta}  单 link 核特征值 hat_w = {np.round(what, 10)}")
        P(f"        sigma := -ln(hat_w_1/hat_w_0) = {sigma:.12f}   (单位长度的弦张力)")
        P(f"        {'L':>3} {'dim':>6} {'限制维':>8} {'lam1':>18} {'lam2':>18} {'m(L)':>14} {'m(L)/L':>14}")
        for L in (2, 3, 4, 5):
            T = np.array([[1.0]])
            for _ in range(L):
                T = np.kron(T, w.reshape(-1, 1) if False else
                            np.array([[w[(b - a) % N] for b in range(N)] for a in range(N)]))
            dim = N ** L
            basis = []
            for Pv in range(N):
                v = np.zeros(dim)
                for idx in range(dim):
                    s, t = 0, idx
                    for _ in range(L):
                        s += t % N
                        t //= N
                    if s % N == Pv:
                        v[idx] = 1.0
                basis.append(v / np.linalg.norm(v))
            B = np.array(basis).T
            ev = np.sort(np.real(np.linalg.eigvals(B.T @ T @ B)))[::-1]
            mL = -np.log(ev[1] / ev[0])
            P(f"        {L:3d} {dim:6d} {B.shape[1]:8d} {ev[0]:18.8f} {ev[1]:18.8f} "
              f"{mL:14.10f} {mL/L:14.10f}")
        P(f"        => m(L)/L = sigma 逐位成立 => **m(L) = sigma*L（严格递增、线性、-> infinity）**")

    # SU(2) Wilson: 字符系数用 Haar 积分, Peter-Weyl 本征值 = c_R/d_R
    P("\n[SU(2) Wilson] w(U)=exp(beta cos t);  c_jj = (2/pi)int sin^2 t [sin(jj t)/sin t] e^{beta cos t} dt")
    P("              重构自校验: sum_jj jj*c_jj = e^beta  (w(theta=0)=sum_R d_R c_R)")
    P("              传递本征值 hat_w(R) = c_jj/d_R, d_R=jj  => lam2/lam1 = (c_2/2)/c_1")
    for beta in (0.5, 1.0, 2.5, 4.0, 6.0):
        c = []
        for jj in range(1, 41):
            f = lambda t, jj=jj: (2.0 / np.pi) * np.sin(t) ** 2 * (np.sin(jj * t) / np.sin(t)) * np.exp(beta * np.cos(t))
            c.append(quad(f, 0, np.pi, limit=500)[0])
        recon = sum(jj * c[jj - 1] for jj in range(1, 41))
        P(f"   beta={beta}: c_1={c[0]:.10f} c_2={c[1]:.10f} c_3={c[2]:.10f} c_4={c[3]:.10f}")
        P(f"      rec sum= {recon:.8f}  vs e^beta= {np.exp(beta):.8f}  rel.err={abs(recon-np.exp(beta))/np.exp(beta):.1e}"
          f"   PF(c_2/2<=c_1): {c[1]/2 <= c[0]+1e-12}")
        P(f"      m(L) = -ln((c_2/2)/c_1) = {-np.log((c[1]/2)/c[0]):.10f}   [对任意 L 相同: 见 Z_N]")
    P("\n[heat-kernel] w_j=exp(-lam*C2(j)), C2=j(j+1) => hat_w_j = w_j => m = lam*C2(1/2) = 3lam/4 (任意 L)")
    for lam in (0.1, 0.25, 0.5):
        P(f"   lam={lam}: sigma = 3lam/4 = {0.75*lam:.10f}   (L 不出现; 对应 B01/T03 的 m=3g^2/8)")
    P("\n=> 2D YM: **m(L) = sigma*L 严格递增**（唯一递增的精确模型）但 **极限是 infinity**，")
    P("   因为第二态是绕环通量 torelon，不是体态。sigma 才是 L-无关的物理量。")


# =====================================================================
# PART B — M6
# =====================================================================
def part_B(omega=1.0):
    P("\n" + "=" * 74)
    P("PART B — M6 (谐振 ⊗ 自由): 精确谱 => m_L = min(omega, 2 pi^2/L^2)")
    P("=" * 74)
    na, amax = 601, 12.0
    a = np.linspace(-amax, amax, na)
    h = a[1] - a[0]
    main = 1.0 / h ** 2 + 0.5 * omega ** 2 * a ** 2
    off = -0.5 / h ** 2 * np.ones(na - 1)
    Ha = sparse.diags([off, main, off], [-1, 0, 1])
    ev_a = np.sort(eigsh(Ha, k=3, which='SA', maxiter=100000)[0])
    P(f"[a 方向] na={na}, h={h:.5f}: E0={ev_a[0]:.10f} (精确 omega/2={omega/2:.10f}); "
      f"第一激发={ev_a[1]-ev_a[0]:.10f} (精确 omega={omega:.10f})")
    P(f"\n{'L':>5} {'2pi^2/L^2':>14} {'m_L=min':>14} {'m_L 数值':>16} {'m_L/L^0?':>10}")
    for L in (2, 3, 4, 6, 8, 10, 16, 20, 40, 80):
        lad = 2 * np.pi ** 2 / L ** 2
        grid = [n * (ev_a[1] - ev_a[0]) + lad * k * k
                for n in range(3) for k in range(-3, 4) if not (n == 0 and k == 0)]
        mL_num = min(grid)
        P(f"{L:5d} {lad:14.8f} {min(omega, lad):14.8f} {mL_num:16.10f}")
    P("=> 前段 = omega（常），后段 = 2pi^2/L^2（严格递减 -> 0）; **从不递增**")


# =====================================================================
# PART C — 横场 Ising 链 (PBC)
# =====================================================================
def tfim_H(L, g):
    dim = 1 << L
    s = np.arange(dim, dtype=np.int64)
    diag = np.zeros(dim)
    for i in range(L):
        szi = 1 - 2 * ((s >> i) & 1)
        szj = 1 - 2 * ((s >> ((i + 1) % L)) & 1)
        diag -= szi * szj
    H = sparse.diags(diag)
    rows = np.arange(dim)
    for i in range(L):
        H = H + sparse.csr_matrix((-g * np.ones(dim), (rows, s ^ (1 << i))), shape=(dim, dim))
    return H.tocsr()


def tfim_gaps(L, g, nev=3):
    H = tfim_H(L, g)
    ev = np.sort(eigsh(H, k=nev, which='SA', maxiter=300000, tol=1e-13)[0])
    # 奇宇称子空间 (P = prod sx = 全比特翻转): 反对称组合
    dim = 1 << L
    mask = dim - 1
    s = np.arange(dim)
    rep = np.minimum(s, s ^ mask)
    keep = s < (s ^ mask)              # 每个配对取一个代表
    reps = s[keep]
    nrep = len(reps)
    idx = {int(r): i for i, r in enumerate(reps)}
    out = {}
    for sign, tag in ((-1.0, "odd"), (+1.0, "even")):
        rowsP, colsP, valsP = [], [], []
        for i, r in enumerate(reps):
            rowsP += [i, i]
            colsP += [int(r), int(r ^ mask)]
            valsP += [1 / np.sqrt(2), sign / np.sqrt(2)]
        Pm = sparse.csr_matrix((valsP, (rowsP, colsP)), shape=(nrep, dim))
        A = (Pm @ H @ Pm.T).tocsr()
        e = np.sort(eigsh(A, k=3, which='SA', maxiter=300000, tol=1e-13)[0])
        out[tag] = e
    return ev, out


def part_C():
    P("\n" + "=" * 74)
    P("PART C — 横场 Ising 链 (PBC, J=1): 定义谱隙 vs 奇宇称扇区谱隙")
    P("=" * 74)
    for g in (0.5, 1.5):
        P(f"\n[g={g}] 精确 m_inf = 2|1-g| = {2*abs(1-g):.10f}   (kink/单粒子质量)")
        P(f"{'L':>4} {'E0':>15} {'E1(full)':>15} {'m_def':>14} {'E1(even)-E0':>15} {'E1(odd)-E0':>15}")
        rows = []
        for L in (4, 6, 8, 10, 12):
            ev, sec = tfim_gaps(L, g)
            mdef = ev[1] - ev[0]
            mev = sec["even"][1] - ev[0]
            mod = sec["odd"][0] - ev[0]
            rows.append((L, mdef, mev, mod))
            P(f"{L:4d} {ev[0]:15.9f} {ev[1]:15.9f} {mdef:14.8f} {mev:15.8f} {mod:15.8f}")
        dm = [rows[i + 1][1] - rows[i][1] for i in range(len(rows) - 1)]
        de = [rows[i + 1][2] - rows[i][2] for i in range(len(rows) - 1)]
        P(f"  m_def 差分: {np.round(dm, 8)}   全<0: {all(d < 0 for d in dm)}")
        P(f"  偶扇区隙差分: {np.round(de, 8)}   全<0: {all(d < 0 for d in de)}")
        P(f"  L=4->12: m_def {rows[0][1]:.3e} -> {rows[-1][1]:.3e} ; 偶扇区隙 {rows[0][2]:.6f} -> {rows[-1][2]:.6f}")
    P("\n=> g<1 (有序相, 真空二重简并): m_def -> 0 (隧穿劈裂), 真 m=1; 偶扇区隙 -> 2 (=2m, 二 kink 阈)")
    P("   => 环上 PBC 的谱里**根本不存在质量为 m 的单粒子态**(单 kink 不能在环上存活)")
    P("=> g>1 (唯一真空): m_def = 奇扇区隙, 从**上**方严格递减到 m_inf")
    P("=> 两个相都**没有递增支**")


# =====================================================================
# PART D — 自由有质量标量场
# =====================================================================
def part_D(m=1.0):
    P("\n" + "=" * 74)
    P("PART D — 自由有质量标量场 on ring, 精确")
    P("=" * 74)
    P(f"{'L':>5} {'m_L 周期':>14} {'m_L 扭曲':>14}")
    for L in (2, 4, 8, 16, 32, 64):
        P(f"{L:5d} {m:14.10f} {np.sqrt(m**2+(np.pi/L)**2):14.10f}")
    P("=> 周期: m(L) ≡ m 恒等（k=0 模在, L 完全不出现）; 扭曲: 严格递减 -> m")


def part_E():
    P("\n" + "=" * 74)
    P("PART E — 单 plaquette (0 维): m(L) 无定义")
    P("=" * 74)
    P("无 transfer 方向 => 无谱 => m 无定义。同 T03 §1 纪律。")


def part_F():
    P("\n" + "=" * 74)
    P("PART F — 汇总")
    P("=" * 74)
    rows = [
        ("M2 2D YM (Migdal/Z_N 显式)", "sigma*L", "严格递增(线性)", "-> infinity (torelon)"),
        ("M6 谐振⊗自由", "min(omega,2pi^2/L^2)", "不增", "-> 0"),
        ("M5 自由标量(周期)", "m", "常", "= m"),
        ("自由标量(扭曲)", "sqrt(m^2+pi^2/L^2)", "严格递减", "-> m (上方)"),
        ("TFIM g<1 (有序相)", "~A e^{-kL}", "严格递减", "-> 0 (而真隙=1)"),
        ("TFIM g>1 (唯一真空)", "m_inf + 正修正", "严格递减", "-> m_inf (上方)"),
        ("M1 单 plaquette", "无定义", "n/a", "n/a"),
    ]
    P(f"{'模型':<30}{'m(L)':<24}{'单调方向':<16}{'极限':<22}")
    P("-" * 92)
    for r in rows:
        P(f"{r[0]:<30}{r[1]:<24}{r[2]:<16}{r[3]:<22}")
    P("\n=> 递增只出现在 2D YM，而出处是 torelon（∝L），极限 = infinity，")
    P("   所以 'm(L) 递增 => m(inf) >= m(2) > 0' 在唯一递增的实例上**是空话**（m(inf)=inf）。")
    P("=> 有有限极限的实例里，m(L) 全部**不增**（常或严格递减）。")
    P("=> 没有任何精确模型是 '递增 + 有限正极限' = 路线 (a) 的假设。")


def main():
    P("D03-ML: m(L) 的 L-依赖, 精确模型测试")
    P(f"numpy {np.__version__} / scipy {__import__('scipy').__version__}")
    part_A()
    part_B()
    part_C()
    part_D()
    part_E()
    part_F()
    P("\n[DONE]")


if __name__ == "__main__":
    main()
    with open("w7/d03_ml_models.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")

