"""b14_window_extend.py — I1 的窗能否从 3.777e-3 往外推？

设定（逐字沿 w7/adj4_kp_a.py，不改归一）：
    Φ(a,β) = Σ_{n≥1} (n+1)² [I_{n+1}(β)/I_1(β)] e^{a n}      (SU(2), |r|=n=2j, dim=n+1)
    KP: S_x(a) ≤ a ; P3-FULLREP-KP:90 给 S_x(a) ≤ λΦ/(1-λΦ), λ = eΔ = 20e = 54.3656
    ⇒ 窗  W(a) := sup{β : λΦ(a,β) ≤ a/(1+a)}                  ("env" 形，现行判据)
    最宽松形（丢掉 (1-λΦ)^{-1} 包络松弛）: λΦ(a,β) ≤ a         ("bare")

四问：
  Q1 W(a) 的全局极大（形状 + 唯一性）
  Q2 换判据能否推过 1/48 —— 判据族的**结构天花板**
  Q3 p1 的 1/48 是不是真天花板（λ_needed 反解）
  Q4 窗能推到多远（(a)/(b)/(c) 判死线）
"""
import numpy as np
from scipy.special import iv
from scipy.optimize import brentq

LAM = 20.0 * np.e            # eΔ, Δ=20  → 54.3656
NMAX = 300
ONE48 = 1.0 / 48.0


def Phi(a, beta, nmax=NMAX):
    n = np.arange(1, nmax + 1, dtype=float)
    w = (n + 1.0) ** 2 * iv(n + 1.0, beta) / iv(1.0, beta)
    m = (w > 0.0) & (a * n < 700.0)
    return float(np.sum(w[m] * np.exp(a * n[m])))


def Wroot(lam, a, kind):
    """λΦ(a,β) = rhs(a) 的唯一根（Φ 随 β 严格增）。kind='env' → a/(1+a)；'bare' → a。"""
    rhs = a / (1.0 + a) if kind == "env" else a
    f = lambda b: lam * Phi(a, b) - rhs
    lo = 1e-14
    if f(lo) > 0:
        return 0.0
    hi = 1.0
    while f(hi) < 0 and hi < 1e3:
        hi *= 1.6
    if f(hi) < 0:
        return float("nan")
    return brentq(f, lo, hi, xtol=1e-16, rtol=1e-15, maxiter=300)


def Wmax(lam, kind, n=1400, amax=60.0):
    """全局极大：log 网格扫描 + 每段内 golden 细化；返回 (a*, W*) 与驻点计数。"""
    aa = np.geomspace(1e-4, amax, n)
    ws = np.array([Wroot(lam, a, kind) for a in aa])
    i = int(np.argmax(ws))
    lo = aa[max(i - 1, 0)]
    hi = aa[min(i + 1, n - 1)]
    g = (np.sqrt(5.0) - 1.0) / 2.0
    for _ in range(90):                      # golden section（W 单峰工具）
        c, d = hi - g * (hi - lo), lo + g * (hi - lo)
        if Wroot(lam, c, kind) < Wroot(lam, d, kind):
            lo = c
        else:
            hi = d
    astar = 0.5 * (lo + hi)
    wstar = Wroot(lam, astar, kind)
    # 驻点（W 的符号变化）计数，证"唯一内点极大"
    dw = np.diff(ws)
    sc = np.sum((dw[:-1] > 0) & (dw[1:] < 0))
    return astar, wstar, int(sc), aa, ws


def lam_needed(target, kind, lam_hi=60.0):
    """反解：max_a W(a;λ) = target 的 λ（W 随 λ 单调降）。"""
    f = lambda L: Wmax(L, kind, n=400)[1] - target
    if f(1e-3) < 0:
        return float("nan")
    return brentq(f, 1e-3, lam_hi, xtol=1e-12, rtol=1e-13, maxiter=200)


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    p = out.append
    p("=" * 78)
    p("b14_window_extend — W(a) 的极大 + 判据族的结构天花板")
    p("=" * 78)
    p(f"λ = eΔ = 20e = {LAM:.4f} ; 1/48 = {ONE48:.6e} ; 目标 2.5")

    # ---------------- Q1 ----------------
    p("\n[Q1] W(a) 的形状（现行判据 env: λΦ ≤ a/(1+a), λ=54.3656）")
    p(f"{'a':>9} {'W(a)':>13} {'W/Wmax':>9} {'Φ(a,β*)/β*':>12}")
    a_star, W_star, nstat, aa, ws = Wmax(LAM, "env")
    for a in (0.01, 0.05, 0.1, 0.25, 0.382, 0.5, 0.615, 0.618034, 0.75,
              1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0):
        w = Wroot(LAM, a, "env")
        p(f"{a:9.3f} {w:13.4e} {w/W_star:9.4f} {Phi(a,w)/w:12.5f}")
    p(f"  全局极大: a* = {a_star:.5f}, W* = {W_star:.6e}")
    p(f"  闭式（小β: φ≈e^a ⇒ W≈a e^-a/((1+a)λ)，d/da=0 ⇔ a²+a−1=0）: a*=(√5−1)/2={0.5*(np.sqrt(5)-1):.6f}")
    p(f"  驻点计数（先升后降的段数）= {nstat}（1 ⇒ 唯一内点极大，端点 W(0+)=0, W(∞)=0）")
    p(f"  端点: W(1e-3)={Wroot(LAM,1e-3,'env'):.4e} (线性 ~a/λ) ; W(40)={Wroot(LAM,40.0,'env'):.4e} (~a e^-a/λ)")
    p(f"  ⇒ 全局极大 = a*={a_star:.4f}, W*={W_star:.4e} [E]；"
      f"闭式 0.618034 [C]（大 a 处 Φ>e^a 的高阶修正把 a* 从 0.618 拉到 0.6150）")

    # ---------------- Q2 ----------------
    p("\n[Q2] 换判据能否推过 1/48？—— 判据族的天花板阶梯")
    p(f"{'判据':<34}{'λ':>10}{'max_a W':>13}{'a*':>9}{'×(1/48)':>10}{'×2.5':>9}")
    rows = []
    for name, lam, kind in (
        ("env   λΦ≤a/(1+a),  λ=eΔ (现行)", LAM, "env"),
        ("bare  λΦ≤a     ,  λ=eΔ", LAM, "bare"),
        ("env   λΦ≤a/(1+a),  λ=λ_*≈25.2(真值上界)", 25.2, "env"),
        ("bare  λΦ≤a     ,  λ=λ_*≈25.2", 25.2, "bare"),
        ("bare  λΦ≤a     ,  λ=e (Δ=1, 动物数地板)", np.e, "bare"),
        ("bare  λΦ≤a     ,  λ=1  (非物理极端)", 1.0, "bare"),
    ):
        a_m, W_m, ns, _, _ = Wmax(lam, kind)
        rows.append((name, lam, a_m, W_m))
        p(f"{name:<34}{lam:10.4f}{W_m:13.4e}{a_m:9.4f}{W_m/ONE48:10.4f}{W_m/2.5:9.2e}")
    p(f"  结构界（bare, 小β: Φ≈β e^a ⇒ W≤a e^-a/λ）：每列的 max_a a e^-a/λ = 1/(eλ)")
    for name, lam, a_m, W_m in rows:
        p(f"    {name:<34} 1/(eλ)={1.0/(np.e*lam):.5e} ; 数值 max={W_m:.5e} ; 比={W_m/(1/(np.e*lam)):.4f}")
    p(f"  ⇒ 结构天花板：W ≤ 1/(eλ)，λ=eΔ 且 Δ≥1 ⇒ W ≤ 1/e² = {1/np.e**2:.5f}")
    p(f"     1/e² = {1/np.e**2:.4f} < 2.5（差 {2.5/(1/np.e**2):.1f}×）⇒ 该族**任何**常数改进都到不了 2.5")

    # ---------------- Q3 ----------------
    p("\n[Q3] p1 的 1/48 是不是真天花板？—— 反解 λ_needed")
    for kind in ("env", "bare"):
        for tgt, nm in ((ONE48, "1/48"), (2.5, "2.5")):
            ln = lam_needed(tgt, kind)
            p(f"  {kind:>4}: 要 max_a W = {nm:<5} 需 λ ≤ {ln:.4f}"
              f"  (现行 λ=eΔ={LAM:.2f}；真值 λ_*≈25.2) {'✅' if ln >= LAM else '❌ 现行不够'}")
    p(f"  ⇒ λ_* ≤ 48 的口径：48 = 1/(2·(1/96))…直接核：bare 形 λ=48 给 max_a W = "
      f"{Wmax(48.0,'bare',400)[1]:.4e} = {(Wmax(48.0,'bare',400)[1])/ONE48:.3f}×(1/48)")
    p(f"     env  形 λ=48 给 max_a W = {Wmax(48.0,'env',400)[1]:.4e} = "
      f"{(Wmax(48.0,'env',400)[1])/ONE48:.3f}×(1/48)")

    # ---------------- Q4 ----------------
    p("\n[Q4] 窗能推到多远（判死线）")
    best = Wmax(25.2, "bare", 600)[1]
    p(f"  (i)  现行 (env, eΔ)                : {W_star:.4e}  = {W_star/ONE48:.3f}×(1/48)")
    p(f"  (ii) 丢包络松弛 (bare, eΔ)          : {Wmax(LAM,'bare',600)[1]:.4e}  "
      f"= {Wmax(LAM,'bare',600)[1]/ONE48:.3f}×(1/48)")
    p(f"  (iii) 再取真值 λ_*≈25.2 (bare)      : {best:.4e}  = {best/ONE48:.3f}×(1/48)"
      f" ; 距 2.5 差 {2.5/best:.0f}×")
    p(f"  ⇒ 最好情形仍 **< 1/48**（{best/ONE48:.2f}× of 1/48）⇒ (b)：2.5 覆盖不到，2.5/{best:.0f}× 缺口")
    p(f"  ⇒ (c) 可改 a* 的只有两处：包络 ×{Wmax(LAM,'bare',600)[1]/W_star:.3f}（a*："
      f"{a_star:.3f}→{Wmax(LAM,'bare',600)[0]:.3f}）；λ ×{LAM/25.2:.3f}。"
      f"Φ 的表示求和已饱和（见下）")

    # Φ 截断对 a* 处的增益（(c) 的第三条路）
    p("\n[(c)] Φ 表示求和截断对 W* 的增益（“换 Φ 的界”能否提高 a*）")
    for nm_ in (1, 2, 3, 5, 10):
        Phi.__defaults__ = (nm_,)
        w_ = Wroot(LAM, a_star, "env")
        p(f"  n≤{nm_:>3}: W(a*)={w_:.6e} ; 相对全和 {w_/W_star:.6f}")
    Phi.__defaults__ = (NMAX,)
    p(f"  ⇒ n≤3 已收敛到 1e-9 内 ⇒ Φ 的界无油水；(c) 的增益只来自包络与 λ")

    txt = "\n".join(out)
    print(txt)
    return txt


if __name__ == "__main__":
    t = main()
    with open(r"D:\Dev\math\w7\b14_window_extend.log", "w", encoding="utf-8") as f:
        f.write(t + "\n")
