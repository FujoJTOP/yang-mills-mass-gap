"""adj4_kp_a.py — KP 判据里的 a 取多少？Phi(a) 三值 + 窗 W(a)。

约定（语料 SU(2)，e^{β cos θ} 权，c_0=2I_1(β)/β）：
    c_n = 2(n+1) I_{n+1}(β)/β ,  dim r = n+1 ,  |r| = n = 2j
    Φ(a) = Σ_{r≠0} (c_r dim r / c_0) e^{a|r|} = Σ_{n≥1} (n+1)² [I_{n+1}(β)/I_1(β)] e^{a n}

KP（per-site 形，RS8 §1.3 ≡ KP 定理在 a(γ)=a|γ| 的特例）：S_x(a) ≤ a。
P3-FULLREP-KP:90 给 S_x(a) ≤ λΦ(a)/(1-λΦ(a))，λ=eΔ=20e=54.3656。
⇒ 窗  W(a) := sup{β>0 : λΦ(a,β) ≤ a/(1+a)}。
"""
import numpy as np
from scipy.special import iv

LAM = 20.0 * np.e          # eΔ, Δ=20 (P3-FULLREP-KP §4.1)


def Phi(a, beta, nmax=300):
    n = np.arange(1, nmax + 1, dtype=float)
    w = (n + 1.0) ** 2 * iv(n + 1.0, beta) / iv(1.0, beta)
    m = (w > 0.0) & (a * n < 600.0)          # 防 0*inf 与 exp 溢出（不可见项本就 <1e-300）
    return float(np.sum(w[m] * np.exp(a * n[m])))


def W(a, lo=1e-9, hi=4.0):
    """解 λΦ(a,β)=a/(1+a)（Φ 随 β 严格增，唯一根）。a→0 时根→0。"""
    tgt = a / (1.0 + a)
    f = lambda b: LAM * Phi(a, b) - tgt
    if f(lo) > 0:
        return 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    p = out.append
    p("=" * 72)
    p("adj4_kp_a — KP 的 a 该取多少；Φ(a) 三值；窗 W(a)")
    p("=" * 72)
    p(f"λ = eΔ = 20e = {LAM:.4f} (P3-FULLREP-KP §4.1, Δ=20)")

    # ---- §1 复算 cb4 的表（定标校验）----
    p("\n[1] 定标校验（对上 CB4.md §3 的表）")
    p(f"{'β':>10} {'Φ(0)':>16} {'Φ(1)':>16} {'Φ(1)/Φ(0)':>12}")
    for b in (1 / 48, 0.018394, 2.5):
        f0, f1 = Phi(0.0, b), Phi(1.0, b)
        p(f"{b:10.6f} {f0:16.9f} {f1:16.9f} {f1/f0:12.5f}")
    p("  CB4: 1/48 → 0.02099646854 / 0.05784772974 / 2.75512 ; 2.5 → 5.050788395 / 85.18214939 / 16.8651")

    # ---- §2 Φ(a) 三值 @ β=0.018394 ----
    b = 0.018394
    p(f"\n[2] Φ(a) @ β={b}（SU(2)，含全部表示；n≤300 已收敛）")
    p(f"{'a':>7} {'Φ(a)':>14} {'Φ(a)/Φ(0)':>11} {'Φ(a)/β':>10}")
    vals = {}
    for a in (0.0, 1.0, 2.0):
        v = Phi(a, b)
        vals[a] = v
        p(f"{a:7.3f} {v:14.9f} {v/vals[0.0]:11.5f} {v/b:10.5f}")
    p(f"  Φ(0)=β+3β²/8 检验: β+3β²/8 = {b + 0.375*b*b:.9f} vs Φ(0) {vals[0.0]:.9f}")
    p(f"  Φ(2)/Φ(1) = {vals[2.0]/vals[1.0]:.5f}")

    # ---- §3 窗 W(a)：a 是自由的 ⇒ 对 a 优化 ----
    p("\n[3] 窗 W(a)：λΦ(a,β) ≤ a/(1+a) 的根")
    p(f"{'a':>7} {'W(a)':>12} {'×W(1)':>10} {'Φ(a)/β':>10}")
    W1 = W(1.0)
    grid = [0.0, 0.05, 0.1, 0.25, 0.382, 0.5, 0.618034, 0.75, 1.0, 1.5, 2.0, 3.0]
    best = (0.0, 0.0)
    for a in grid:
        w = W(a) if a > 0 else 0.0
        if w > best[1]:
            best = (a, w)
        p(f"{a:7.3f} {w:12.3e} {w/W1:10.3f} {Phi(a,b)/b:10.5f}")
    # 细扫 a* 与上界（闭式 a e^{-a}/(1+a) 的极大点 = (√5−1)/2）
    aa = np.linspace(0.2, 1.2, 201)
    ws = np.array([W(x) for x in aa])
    i = int(np.argmax(ws))
    a_star, W_star = float(aa[i]), float(ws[i])
    g = lambda x: x * np.exp(-x) / (1 + x)
    p(f"  细扫: a*={a_star:.4f}, W*={W_star:.4e} ; 闭式极大点 (√5−1)/2={0.5*(np.sqrt(5)-1):.6f}")
    p(f"  g(a)=a e^-a/(1+a): g(1)={g(1.0):.5f}, g(a*)={g(a_star):.5f} ⇒ 优化增益 {g(a_star)/g(1.0):.3f}×")
    p(f"  W(1)  = {W1:.4e}   (cb4/P6-SPEC: 3.357e−3 ; 差 {abs(W1-3.357e-3)/3.357e-3*100:.2f}%)")
    p(f"  W*    = {W_star:.4e}   (a 自由、最优)")

    # ---- §4 a=0 不是 KP 实例 ----
    p("\n[4] a=0 检验：KP 右端 a(γ)=0·|γ|=0，而左端 S_x(0) ≥ |w| > 0 ⇒ 判据恒假")
    p(f"  β→0+ : S_x(0) ≈ λΦ(0) = {LAM*vals[0.0]:.6f} > 0 = a(γ)  ⇒ 无任何 β>0 满足")
    # 语料的两个 a=0 根
    from scipy.optimize import brentq
    r1 = brentq(lambda x: LAM * Phi(0.0, x) - 1.0, 1e-9, 1.0)      # λΦ(0)=1
    r2 = 1.0 / LAM                                                  # Φ(0)=β 版
    p(f"  λΦ(0,β)=1 的根（M3C-TAIL 的判据根）= {r1:.9f}   (语料 0.018268566)")
    p(f"  1/(eΔ) = {r2:.9f}   (必要条件根 0.018394)")
    p(f"  ⇒ 这两个数是 a=0 的'活度和收敛'根，不是 KP 窗")

    # ---- §5 高报倍数 ----
    p("\n[5] 高报倍数（I1 声称 0.018269）")
    for name, w in (("W(a=1)", W1), ("W(a*)", W_star), ("W(a=2)", W(2.0))):
        p(f"  0.018269 / {name} = {0.018269/w:.3f}×   ;  0.018394 / {name} = {0.018394/w:.3f}×")

    # ---- §6 结论 ----
    p("\n[6] 判")
    p("  (1) KP: a 自由（a: 聚合物→ℝ₊，无 a≥1）；a>0 是自洽要求（右端须 ≥ 活度）。")
    p("  (2) Φ(a) 的 a = KP 的 a（同一 grading/衰减率参数），不是两个量。")
    p("  (3) a=0 使右端=0 ⇒ 0.0184 不是 KP 窗，是 a=0 的活度和收敛根。")
    p(f"  (4) 真窗：a 自由优化 a*={a_star:.3f} → {W_star:.3e}；a=1（语料未优化取法）→ {W1:.3e}。")
    p(f"  (5) 判：(a) 的数（e−3 量级，0.018269 高报 {0.018269/W_star:.1f}×）+ 理由改正；(b) 不成立；"
      f"(c) 部分（同一 a，但 a=0 处不可比）。")
    txt = "\n".join(out)
    print(txt)
    return txt


if __name__ == "__main__":
    import sys
    t = main()
    with open(r"D:\Dev\math\w7\adj4_kp_a.log", "w", encoding="utf-8") as f:
        f.write(t + "\n")
