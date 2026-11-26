"""O05 — Var(G_v+Phi'|a=0) 的 beta-标度：单 plaquette SU(2) 精确一维求积。

模型（约定同 H1A-EXACT.md §2.1 / w7/na_plaquette.py）：
  xi in R^3, |xi|<2pi;  S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2)*u),
  u = xi_3/rho;  平坦图 dU = J dX, J = sin^2(rho/2)/rho^2, Phi = -log J
  v = e_3, 条件化 <v,xi> = xi_3 = 0  =>  切片上 rho=r, u=0, 且 Phi' == 0 (逐点精确)
  切片上精确（Hermite/BLAS 无关，纯解析）：
     G   = d_3 S          = beta*sin(th/2)*sin(r/2)/r
     H   = d_3^2 S        = beta*cos(th/2)*sin(r/2)/(2r)
     P'' = d_3^2 Phi      = 2/r^2 - cot(r/2)/r        (3D 径向: Phi_r/r)
  切片测度 w(r)dr = 2*pi*r*e^{-S}*J dr = 2*pi*sin^2(r/2)/r * e^{-S} dr
  故 Var(G_v+Phi'|0) = Var(G|0), E[H+Phi''|0] = E[H|0] + E[P''|0]。
只写 o05_*.py|.log（禁 git commit）。
"""
import numpy as np

NPTS = 400001
BETAS = [0.5, 1.0, 2.5, 6.0, 10.0, 20.0, 50.0, 200.0]
THETAS = [0.0, np.pi / 2, 3 * np.pi / 4, np.pi]
TAB_BETAS = [1.0, 2.5, 6.0, 10.0]


def slice_stats(beta, theta, n=NPTS):
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    r = np.linspace(1e-9, 2 * np.pi - 1e-9, n)
    Sr, Cr = np.sin(r / 2), np.cos(r / 2)
    S = beta * (1.0 - Cr * c)
    w = 2 * np.pi * Sr**2 / r * np.exp(-S)              # 切片权重 (J 已含)
    Z = np.trapezoid(w, r)
    G = beta * s * Sr / r                                # 力
    H = beta * c * Sr / (2 * r)                          # 条件曲率 (S_W 部分)
    P2 = 2.0 / r**2 - (Cr / Sr) / r                      # Phi''_vv
    Ga = np.trapezoid(w * G, r) / Z
    var = np.trapezoid(w * (G - Ga) ** 2, r) / Z
    EH = np.trapezoid(w * H, r) / Z
    EP2 = np.trapezoid(w * P2, r) / Z
    # |grad F|^2 = (dG/dr)^2 （切片上一维径向梯度）
    dG = np.gradient(G, r)
    Egrad2 = np.trapezoid(w * dG**2, r) / Z
    return dict(var=var, EH=EH, EP2=EP2, grad2=Egrad2, meanG=Ga, Z=Z, r=r, w=w)


def poincare_1d(r, w, M=600):
    """一维 (径向 m=0) 条件 Poincaré 常数：inf int f'^2 dw / int (f-fbar)^2 dw."""
    F = np.concatenate(([0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(r))))
    F /= F[-1]
    keep = np.concatenate(([True], np.diff(F) > 0))
    u = np.linspace(0.0, 1.0, M + 1)
    ru = np.interp(u, F[keep], r[keep])
    wu = np.interp(ru, r, w)
    Z = np.trapezoid(w, r)
    a = (0.5 * (wu[1:] + wu[:-1]) / Z) ** 2               # 中点系数 a(u) = p(u)^2, 长度 M
    dg = np.empty(M + 1)                                  # (-d/du a d/du), Neumann 自由端
    dg[0] = a[0]
    dg[M] = a[M - 1]
    for i in range(1, M):
        dg[i] = a[i - 1] + a[i]
    T = np.diag(dg) + np.diag(-a, 1) + np.diag(-a, -1)
    T *= M * M
    ev = np.linalg.eigvalsh(T)
    return float(ev[1])                                  # ev[0]≈0 (常数模)


def main():
    print("[O05] Var(G_v+Phi'|0) vs E[H+Phi''|0]  beta-标度; 单 plaquette SU(2), 精确求积")
    print(f"[O05] 网格 n={NPTS}  (r in (0,2pi) 均匀, 端点权为 0)")
    # --- 存档自检: theta=pi 应给出 Var=beta^2*V0, E[H]==0, E[Phi'']=0.2051153
    print("\n[CHK] theta=pi 对账 (H1A-EXACT §2.1: Var=beta^2*0.0110593, E[H]=0, W=0.2051153)")
    for b in [1.0, 2.5, 6.0, 10.0]:
        d = slice_stats(b, np.pi)
        print(f"  beta={b:6.2f}  Var={d['var']:.6f}  Var/beta^2={d['var']/b**2:.7f}"
              f"  E[H]={d['EH']:.3e}  E[Phi'']={d['EP2']:.7f}")
    # --- 主表
    print("\n[P1] 主表: Var, E[H], E[H+Phi''], 比值")
    print("  theta      beta      Var          E[H]        E[P'']      E[H+P'']"
          "    Var/E[H]   Var/(E[H]+E[P''])")
    for th in THETAS:
        for b in TAB_BETAS:
            d = slice_stats(b, th)
            tot = d["EH"] + d["EP2"]
            rat = d["var"] / d["EH"] if d["EH"] > 1e-12 else float("inf")
            print(f"  {th:6.4f} {b:7.2f}  {d['var']:11.6f}  {d['EH']:10.6f}"
                  f"  {d['EP2']:10.6f}  {tot:10.6f}  {rat:10.4f}  {d['var']/tot:10.4f}")
    # --- 幂指数
    print("\n[P2] 局部幂 dlog(Var)/dlog(beta) 与 dlog(E[H])/dlog(beta)")
    for th in THETAS:
        st = [slice_stats(b, th) for b in BETAS]
        ks = [np.log(st[i + 1]["var"] / st[i]["var"]) / np.log(BETAS[i + 1] / BETAS[i])
              if st[i]["var"] > 0 and st[i + 1]["var"] > 0 else float("nan")
              for i in range(len(BETAS) - 1)]
        ms = [np.log(st[i + 1]["EH"] / st[i]["EH"]) / np.log(BETAS[i + 1] / BETAS[i])
              if st[i]["EH"] > 0 and st[i + 1]["EH"] > 0 else float("nan")
              for i in range(len(BETAS) - 1)]
        print(f"  theta={th:.4f}")
        for i, b in enumerate(BETAS[:-1]):
            print(f"     beta {b:6.2f}->{BETAS[i+1]:6.2f}   k={ks[i]:8.4f}   m={ms[i]:8.4f}")
        print(f"     Var(beta->200)={st[-1]['var']:.6f}  E[H](200)={st[-1]['EH']:.6f}")
    # --- 条件 Poincare 常数 lambda 的 beta-标度
    print("\n[P3] 条件 Poincare 常数 lambda (径向 m=0, CDF 变换 + 三对角本征);"
          " 及 BL 界 (1/lambda)E|grad F|^2")
    print("  theta      beta     lambda       (1/lam)E|gradF|^2    Var        bound/Var")
    for th in [np.pi / 2, np.pi]:
        for b in [1.0, 2.5, 6.0, 10.0, 50.0]:
            d = slice_stats(b, th, n=200001)
            lam = poincare_1d(d["r"], d["w"])
            bound = d["grad2"] / lam
            print(f"  {th:6.4f} {b:7.2f}  {lam:10.5f}   {bound:16.6f}  {d['var']:10.6f}"
                  f"   {bound/d['var']:10.3f}")
    # --- 结论所需的一句: Var 的 beta 上确界 (theta 依赖, beta 无关?)
    print("\n[P4] sup_beta Var(beta,theta)（看它是否 beta-无关、是否 theta-一致）")
    for th in THETAS:
        vs = [slice_stats(b, th, n=100001)["var"] for b in [0.5, 1, 2, 5, 10, 100, 1000]]
        print(f"  theta={th:.4f}  Var 序列={['%.5f' % v for v in vs]}  sup={max(vs):.6f}")


if __name__ == "__main__":
    main()
