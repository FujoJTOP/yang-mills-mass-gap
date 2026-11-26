"""l1sign_tail — (L1) 的 B=(u/(1-u))^3：逐字来源、定义域、u>1 时到底发生什么。

(X01-L1C:77-82 的推导重放；LOOP-RUNBOOK:6040-6044 的 2.329267 复现)
"""
from fractions import Fraction as F
from math import e

ED = e * 20                      # eΔ, Δ=20  (plaquette 邻接度)
K = ED                           # ū = K·u,  K=eΔ·C̄, C̄∈[1,1.0535]
UBAR = F(57, 100)                # 语料数值点 ū=0.57 (X01:97)


def B(u):
    """B = u³/(1-u)³ = (u/(1-u))³  —— X01:80 逐字"""
    return (u / (1 - u)) ** 3


def tail(A, T, N=600):
    """ε_T = Σ_{α>T} (α-1-T) A(α)  —— X01:43 (L1-c.1)"""
    return sum((a - 1 - T) * A(a) for a in range(T + 1, T + 1 + N))


def main():
    print("[1] B 的三次方逐字来源  (X01-L1C:71-82)")
    print("    |A(α)| ≤ ū^{α+1}/(1-ū)          ← Ursell树界 + 动物数 (eΔ)^|Γ|, |Γ|≥α")
    print("    ε_T = Σ_{α>T}(α-1-T)A(α),  m=α-1-T, α=m+T+1")
    print("    Σ_{m≥1} m x^m = x/(1-x)²   ⇒  ū^{T+2}/(1-ū) · ū/(1-ū)²  =  ū^{T+3}/(1-ū)³")
    print("    ⇒ |ε_T| ≤ [ū³/(1-ū)³]·ū^T     ⇒ B=(ū/(1-ū))³   ✔（三次方 = 1/(1-ū) 一次 + 1/(1-ū)² 两次）")

    A = lambda a: UBAR ** (a + 1) / (1 - UBAR)      # 饱和模型（把界取等）
    Bv = float(B(UBAR))
    r = [float(tail(A, T)) / float(UBAR) ** T for T in range(13)]
    print(f"[2] 复现 LOOP-RUNBOOK:6044 / C03:50：B(0.57)={Bv:.6f}")
    print("    T :", "  ".join(f"{T}" for T in range(13)))
    print("    ε_T/ū^T :", " ".join(f"{x:.6f}" for x in r))
    assert max(abs(x - Bv) for x in r) < 1e-9, "ratio must be constant = B"
    print(f"    ⇒ 逐位恒为 B（max|ratio-B| = {max(abs(x-Bv) for x in r):.2e}）  [T]")

    print("\n[3] B(ū) 的符号：B>0 ⟺ ū<1；B→±∞ 于 ū→1∓")
    print("    ū      1-ū        B            备注")
    rows = [(0.28315, "β=1/48, u=β/4"), (0.57, "语料数值点 X01:97"),
            (0.9, ""), (0.99, ""), (1.0, "几何级数收敛半径边"),
            (1.01, ""), (16.98926, "β=2.5, u=β/8"), (33.97852, "β=2.5, u=β/4")]
    for u, note in rows:
        b = B(u) if u != 1 else float("inf")
        print(f"    {u:<9.5f}{1-u:>9.5f}  {b:>13.6f}    {note}")
    assert B(0.9) > 0 > B(1.01)

    print("\n[4] ū>1：解析延拓 ū/(1-ū)² ≠ Σ m ū^m（部分和 vs 延拓值）")
    q = 1.7
    print(f"    ū={q}: 延拓值 ū/(1-ū)² = {q/(1-q)**2:.6f} （正！此因子不翻号）"
          f"；翻号的是 1/(1-ū) = {1/(1-q):.6f}")
    for N in (5, 10, 20, 40):
        s = sum(m * q ** m for m in range(1, N + 1))
        print(f"    N={N:<3} Σ_{{m≤N}} m·ū^m = {s:.6e}   比值/延拓 = {s/(q/(1-q)**2):+.3e}")
    print("    ⇒ 部分和正、超指数增长；延拓值 +3.469 有限正 ⇒ 延拓值不表示该和  [T]")
    print("    ⇒ B 的负号只由 |A| 界里那一个 1/(1-ū) 携带（三次方=奇次；x/(1-x)² 是平方，不翻号）")

    print("\n[5] β=2.5 处把 B<0 代进 (L1)：不等式是假，不是‘平凡真’")
    ub = 2.5 * K / 4                     # u=β/4 (FOURWALL:167), K=eΔ
    print(f"    ū(β=2.5,u=β/4) = {ub:.5f}   B = {B(ub):.6f}   (1-ū) = {1-ub:.5f}")
    A2 = lambda a: ub ** (a + 1) / (1 - ub)   # 界本身已为负 ⇒ 无意义
    print(f"    该‘界’ |A(α)|≤ū^{{α+1}}/(1-ū) 的右端在 α=0 就是负：{(1-ub) and ub/(1-ub):.4f} < 0 ⇒ 前提先垮")
    for T in (0, 1, 2, 5):
        rhs = B(ub) * ub ** T
        print(f"    T={T}: RHS = B·ū^T = {rhs:.6f} < 0 ≤ |ε_T|  ⇒ 不等式假（除 ε_T≡0）")
    print("    ⇒ 只有 ε_T≡0 才能救 ⇒ 与 ε_T≠0（周期管反例/一般情形）矛盾 ⇒ 该式在 ū>1 无内容  [T]")

    print("\n[6] 正确的 ū>1 形式：换成真半径 r = limsup|A(α)|^{1/α}（不是 |1-ū|）")
    print("    一般式：|A(α)| ≤ C·r^α  ⇒  |ε_T| ≤ [C·r²/(1-r)²]·r^T ,  r<1")
    print("    （语料的三次方 = 把 C=ū/(1-ū)、r=ū 代进去：C·r²/(1-r)² = ū³/(1-ū)³ ✔）")
    r_true = 0.57
    C = 0.5
    Br = C * r_true ** 2 / (1 - r_true) ** 2
    A3 = lambda a: C * r_true ** a
    print(f"    取 r=0.57, C=0.5（真幅度模型，与 ū=33.98 的粗界无关）：B_r = {Br:.6f} > 0")
    for T in (0, 2, 5, 10):
        num = abs(tail(A3, T))
        print(f"    T={T:<3} |ε_T| = {num:.6e}   B_r·r^T = {Br*r_true**T:.6e}   比值 = {num/(Br*r_true**T):.6f}")
    print("    ⇒ 形状不变、参数换 r ⇒ 正数、有限、带率 δ=ln(1/r)/a  [T]")
    print("    ✗ 用户候选 |ε_T| ≤ |C|·ū^T（ū>1）保留发散底 ⇒ RHS 随 T 指数上涨，恒真但零信息 [T]")

    print("\n[7] ū<1 的窗（两种 u 约定），及冻结窗是否在内")
    for scale, name in ((4, "u=β/4"), (8, "u=β/8")):
        bc = 1 / (K / scale)
        print(f"    {name}: ū<1 ⟺ u<1/K={1/K:.6f} ⟺ β<{bc:.6f} ;  β=1/48={1/48:.6f} ⇒ {'在内' if 1/48 < bc else '窗外'}")
    print(f"    冻结窗 β<1/48 ⇒ ū ≤ {K/4*(1/48):.5f}, B ≤ {B(K/4/48):.6f} > 0   [T]")


if __name__ == "__main__":
    main()
