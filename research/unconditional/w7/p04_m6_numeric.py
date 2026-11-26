"""p04_m6_numeric.py — M6 (谐振子(a) ⊗ 自由粒子(b)) 的独立数值：
    G″(0), E[H], Var(G), Δ_dyn = E[H] − Var(G)
不引用闭式判据（T03 的 4(E0−V0) 只当旁证）。numpy/scipy，确定性。
"""
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags, kron, identity
from scipy.sparse.linalg import eigsh
from scipy.optimize import brentq


def ho1d(omega, A, N=4001):
    """Dirichlet 盒 [-A,A] 上的谐振子：三对角精确对角化 → (x, E0, psi0)"""
    x = np.linspace(-A, A, N)
    d = 0.5 / ((x[1] - x[0]) ** 2) * 2 + 0.5 * omega**2 * x**2
    e = -0.5 / (x[1] - x[0]) ** 2 * np.ones(N - 1)
    w, v = eigh_tridiagonal(d, e)
    psi = v[:, 0].copy()
    psi *= np.sign(psi[N // 2])
    psi /= np.sqrt(x[1] - x[0])
    return x, w[0], psi


def d2_at(x, y, idx, m):
    h = x[1] - x[0]
    return (y[idx + m] - 2.0 * y[idx] + y[idx - m]) / (m * h) ** 2


def dd_rich(x, y, idx=None):
    """−∂²y 的 Richardson 二阶差分（idx 缺省 = 0 点）"""
    idx = len(x) // 2 if idx is None else idx
    return -(4.0 * d2_at(x, y, idx, 1) - d2_at(x, y, idx, 2)) / 3.0


def m6_2d(omega, L, A=7.0, Na=281, Nb=241):
    """a ∈ [-A,A] × b ∈ [-L/2,L/2]，Dirichlet 双方向；返回 (a, E0, rho_a边际)"""
    a = np.linspace(-A, A, Na)
    ha, hb = a[1] - a[0], L / (Nb - 1)
    Ia, Ib = identity(Na, format="csr"), identity(Nb, format="csr")

    def Lap(N, h):
        return diags([1, -2, 1], [-1, 0, 1], shape=(N, N), format="csr") / h**2

    H = (-0.5 * kron(Lap(Na, ha), Ib) - 0.5 * kron(Ia, Lap(Nb, hb))
         + kron(diags(0.5 * omega**2 * a**2, 0, format="csr"), Ib)).tocsr()
    w, v = eigsh(H, k=1, which="SA", tol=0, maxiter=50000)
    psi = v[:, 0].reshape(Na, Nb)
    psi *= np.sign(psi[Na // 2, Nb // 2])
    rho = (psi**2).sum(axis=1)
    return a, float(w[0]), rho / (rho.sum() * ha)


def b_spec(L, N=3001, bc="periodic"):
    x = np.linspace(-L / 2, L / 2, N, endpoint=(bc != "periodic"))
    h = x[1] - x[0]
    if bc == "dirichlet":
        ve = eigh_tridiagonal(2.0 / h**2 * np.ones(len(x)), -1.0 / h**2 * np.ones(len(x) - 1),
                             eigvals_only=True, select="i", select_range=(0, 5))
        return 0.5 * ve                      # 动能 = −½ d²/dx²
    k = 2 * np.pi * np.fft.fftfreq(len(x), d=L / len(x))
    return np.sort(0.5 * k**2)


def main():
    omega, L = 1.0, 10.0
    o = []
    p = o.append
    p("=" * 80)
    p("p04_m6_numeric — M6: H = ½(p_a²+p_b²) + ½ω²a²   (b 自由，a|b 可分)")
    p(f"numpy {np.__version__} | scipy {scipy.__version__} | 锚点 ω={omega} L={L}")
    p("=" * 80)

    # [1] G″(0)
    x, E0a, psi = ho1d(omega, 9.0)
    h = x[1] - x[0]
    i0 = len(x) // 2
    rho = psi**2
    gpp_log = dd_rich(x, np.log(rho))
    gpp_C = dd_rich(x, rho / rho[i0])        # −ρ″/ρ (归一化后同式)
    j1 = np.argmin(np.abs(x - 1.0))
    gpp_log_1 = dd_rich(x, np.log(rho), j1)
    gpp_C_1 = dd_rich(x, rho / rho[j1], j1)
    dpsi = (psi[i0 + 1] - psi[i0 - 1]) / (2 * h) / psi[i0]
    gpp_id = 4 * (E0a - 0.0) + 2 * dpsi**2

    p("\n[1] G″(0)   （1D 谐振子网格 A=9 N=4001，h=%.3e）" % h)
    p(f"  E0(num) = {E0a:.15f}     ω/2 = {omega/2:.15f}")
    p(f"  G″=−∂²log ρ|_0   = {gpp_log:.12f}   2ω={2*omega:.12f}   err={gpp_log-2*omega:+.2e}")
    p(f"  −ρ″/ρ|_0 (G=C 读法)= {gpp_C:.12f}   err={gpp_C-2*omega:+.2e}   ← 0 点两式同值")
    p(f"  a=1: −∂²log ρ = {gpp_log_1:.12f} ; −ρ″/ρ = {gpp_C_1:.12f}   ← 偏离 0 点两式反号")
    p(f"  ψ0′/ψ0|_0 = {dpsi:+.3e} ⇒ 2(ψ′/ψ)² = {2*dpsi**2:.2e}")
    p(f"  T03 式 4(E0−V0)+2(ψ′/ψ)² = {gpp_id:.12f}（旁证）")

    a2, E2d, rho_a = m6_2d(omega, L)
    gpp_2d = dd_rich(a2, np.log(rho_a))
    var_2d = float((a2**2 * rho_a).sum() * (a2[1] - a2[0]))
    p(f"  [2D 复核 281×241] E0(2D)={E2d:.12f}  G″(边际 b)={gpp_2d:.12f} "
      f"(err {gpp_2d-2*omega:+.2e})  Var_num(a)={var_2d:.12f}")

    # [2] E[H]
    Eb_D, Eb_P = b_spec(L, bc="dirichlet"), b_spec(L, bc="periodic")
    EH_D, EH_P = omega / 2 + Eb_D[1], omega / 2 + Eb_P[1]
    p("\n[2] b 方向自由粒子数值谱")
    p(f"  Dirichlet: {np.array2string(Eb_D[:4], precision=10)}   解析 n²π²/(2L²)="
      f"{np.array2string(np.array([n**2*np.pi**2/(2*L**2) for n in range(4)]), precision=10)}")
    p(f"  periodic : {np.array2string(Eb_P[:4], precision=10)}   解析 2π²n²/L²   ="
      f"{np.array2string(np.array([2*np.pi**2*n**2/L**2 for n in range(4)]), precision=10)}")
    p(f"  E[H]=ω/2+E_b1:  periodic {EH_P:.12f}   Dirichlet {EH_D:.12f}")

    # [3] Var(G), G=a
    var_a = 1.0 / (2 * omega)
    p("\n[3] Var(G), G=a：Var(a)=1/(2ω)=" f"{var_a:.12f}  (2D 数值 {var_2d:.12f}；E0={E0a:.12f})")

    # [4] Δ_dyn
    dP, dD = EH_P - var_a, EH_D - var_a
    p("\n[4] Δ_dyn = E[H] − Var(G)")
    p(f"  periodic : {EH_P:.12f} − {var_a:.12f} = {dP:+.12f}   ⇒ {'>0' if dP>0 else '<0'}")
    p(f"  Dirichlet: {EH_D:.12f} − {var_a:.12f} = {dD:+.12f}   ⇒ {'>0' if dD>0 else '<0'}")
    p(f"  Δ_dyn − E_b1 = {dP-Eb_P[1]:+.2e}  ⇒ ω=1 时 Δ_dyn 恰等于谱隙 m（Var(a)=ω/2=E0）")

    # [5] 扫描
    p("\n[5a] 扫 ω @ L=10（G″、Var 均重算）")
    p(f"  {'ω':>8}{'G″(num)':>15}{'2ω':>9}{'E[H]':>13}{'Var_num':>13}{'Δ_dyn':>14}")
    for om in [0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0, 1.5, 2.0, 5.0, 10.0]:
        xw, E0w, pw = ho1d(om, 9.0 / np.sqrt(min(om, 1.0)))
        gw = dd_rich(xw, np.log(pw**2))
        vw = float((xw**2 * pw**2).sum() * (xw[1] - xw[0]))
        d = (om / 2 + Eb_P[1]) - vw
        p(f"  {om:8.3f}{gw:15.9f}{2*om:9.4f}{om/2+Eb_P[1]:13.9f}{vw:13.9f}{d:+14.9f}")
    fP = lambda om: (om / 2 + Eb_P[1]) - 1.0 / (2 * om)
    fD = lambda om: (om / 2 + Eb_D[1]) - 1.0 / (2 * om)
    wP, wD = brentq(fP, 1e-8, 1.0, xtol=1e-15), brentq(fD, 1e-8, 1.0, xtol=1e-15)
    p(f"  过零 ω*(periodic)={wP:.12f}（二次式核对 {np.sqrt(Eb_P[1]**2+1)-Eb_P[1]:.12f}）")
    p(f"  过零 ω*(Dirichlet)={wD:.12f}（二次式核对 {np.sqrt(Eb_D[1]**2+1)-Eb_D[1]:.12f}）")
    p(f"  样本 ω=0.5: Δ={fP(0.5):+.9f}；ω=0.9: Δ={fP(0.9):+.9f}；ω=2: Δ={fP(2.0):+.9f}")

    p("\n[5b] 扫 L @ ω=1")
    p(f"  {'L':>8}{'E_b1(=m)':>15}{'E[H]':>15}{'Var(a)':>10}{'Δ_dyn':>16}")
    for LL in [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 500.0]:
        eb = 2 * np.pi**2 / LL**2
        p(f"  {LL:8.1f}{eb:15.9f}{0.5+eb:15.9f}{0.5:10.4f}{0.5+eb-0.5:+16.9f}")

    p("\n[6] 标度：Var(a)∝1/ω、E[H]∝ω ⇒ Δ_dyn=(ω−1/ω)/2+m 非齐次，符号随 ω 翻（ω* 如上）；"
      "ω=1 的正号是单位取法的产物。")

    # [7] G 的选择敏感性（同一模型、同一 E[H]，只换 G）
    p("\n[7] G 的选择敏感性 @ ω=1,L=10（E[H]=0.697392088 固定）")
    for name, var in [("G = a        (任务口径)", 1 / (2 * omega)),
                      ("G = ω a      (归一化尺度)", omega**2 / (2 * omega)),
                      ("G = 2ω a = −∂log ρ (score)", (2 * omega) ** 2 / (2 * omega)),
                      ("G = a² (q11 式第二项, Var=E[a⁴]−E[a²]²)", 3 / (4 * omega**2) - (1 / (2 * omega)) ** 2)]:
        p(f"  {name:42s} Var={var:10.6f}  Δ_dyn={EH_P-var:+12.6f}")
    p(f"  旁注：Dirichlet 的隙 m = E_b[1]−E_b[0] = {Eb_D[1]-Eb_D[0]:.12f}（不是 E_b[1]；b 基态在 π²/(2L²)）")

    txt = "\n".join(o)
    print(txt)
    open("p04_m6_numeric.log", "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
