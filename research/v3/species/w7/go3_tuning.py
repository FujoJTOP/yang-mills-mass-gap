# -*- coding: utf-8 -*-
"""
O3 攻击：检验 GV2 的"临界调谐"条件  κ_c<1  ⟺  1−ρ(n) ≲ κ_c^n
纯算术（精确有理数 + n 次根），无 MC、无随机、单线程、峰值内存 < 1 MB。

输入（GV2 §2/§3，n=1..6 精确整数）：
    B_n = 1, 20, 514, 14748, 450255, 14305816
    m_n = u^n B_n ,  u = 5/18 (SU3, β=2.5, u=β/9)  /  5/8 = 0.625 (SU2, u=β/4)
    S_n = Σ_{k≤n}(−1)^{k−1} m_k ,  Z_n = Σ_{k≤n} m_k ,  R_n = |S_n|/Z_n
    读法 A :  ρ_A(n) := |S_n|/Z_n = R_n        ⇒ κ_c(n)^n = ρ_A   （V-A 定义）
    读法 B :  ρ_B(n) := m_n/m_{n−1}            （GV2 §3.4 自己的 ρ_n，壳质量比）
    读法 C :  ρ_C     := e^{γ_ent − c₀q}       （GV2 §1.4 命题 AS 的 ρ，率比，∈(0,1)）
判据（本脚本唯一判据）：  (1−ρ)·κ_c^n  ≤ 1  ⇒ O3 在该读法下"满足"；> 1 ⇒ "不满足"
"""
from fractions import Fraction as F

B = {1: 1, 2: 20, 3: 514, 4: 14748, 5: 450255, 6: 14305816}
NS = [1, 2, 3, 4, 5, 6]
LN3 = 1.0986122886681098  # ln 3


def series(u):
    m, out = {}, {}
    Z = S = F(0)
    for n in NS:
        m[n] = u ** n * B[n]
        Z += m[n]
        S += (F(1) if n % 2 == 1 else F(-1)) * m[n]
        out[n] = dict(m=m[n], Z=Z, S=S, R=abs(S) / Z)
    return out, m


def kc(R, n):
    return float(R) ** (1.0 / n)


def block(u, label, ufrac=None):
    u = F(u) if ufrac is None else ufrac
    ser, m = series(u)
    print(f"\n===== {label} :  u = {u} = {float(u):.10f} =====")
    hdr = ("  n |        m_n |        Z_n |  R_n=|S_n|/Z_n |    kappa_c(n) "
           "| 读A: 1-R_n | 读A: (1-R)R |   rho_B=m_n/m_(n-1) | 1-rho_B | 读B: (1-rho_B)R")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    prev = None
    rows = {}
    for n in NS:
        d = ser[n]
        R = d["R"]
        k = kc(R, n)
        oneA = 1 - R
        valA = oneA * R                      # 读法 A: κ_c^n ≡ R_n
        if prev is None:
            rhoB = None
            oneB = None
            valB = None
        else:
            rhoB = m[n] / prev
            oneB = 1 - rhoB
            valB = oneB * R
        rows[n] = dict(R=R, kc=k, rhoB=rhoB, valA=valA, valB=valB)
        f = lambda x, w=12, p=8: (f"{float(x):>{w}.{p}f}" if x is not None else f"{'--':>{w}}")
        print(f"{n:3d} | {f(d['m'],14,6)} | {f(d['Z'],12,6)} | {float(R):>12.8f} | "
              f"{k:>11.8f} | {float(oneA):>9.3e} | {float(valA):>11.6e} | "
              f"{f(rhoB,17,9)} | {(f'{float(oneB):>9.3e}' if oneB is not None else '       --'):>9} | "
              f"{(f'{float(valB):>11.6e}' if valB is not None else '         --'):>11}")
        prev = m[n]
    return rows


def verdict_A(rows, label):
    print(f"  [读法 A] 判定 (1-R)·κ_c^n ≤ 1 ?   "
          + "  ".join(f"n={n}:{'OK' if rows[n]['valA'] <= 1 else 'FAIL'}({float(rows[n]['valA']):.4f})"
                      for n in NS))
    mx = max(float(rows[n]["valA"]) for n in NS)
    print(f"  [读法 A] 上确界 = {mx:.10f}   vs AM-GM 硬顶 1/4 = 0.25   ⇒ "
          f"{'恒 ≤ 1/4 < 1（与 u、n、B_n 全无关）' if mx <= 0.25 + 1e-15 else '超出 1/4 —— 读法 A 不成立'}")
    print(f"  [读法 A] 等价形式：1-ρ_A ≲ κ_c^n  ⟺  R_n ≳ 1/2 ；viol 条件 R_n < 1/2 ；"
          f"实际最小 R_n = {min(float(rows[n]['R']) for n in NS):.8f}")


def verdict_B(rows, label):
    print(f"  [读法 B] 判定 (1-ρ_B)·κ_c^n ≤ 1 ?")
    for n in NS:
        if rows[n]["rhoB"] is None:
            print(f"      n={n}: n=1 无 ρ_B（需 m_0）—— 空")
            continue
        rb = float(rows[n]["rhoB"])
        vb = float(rows[n]["valB"])
        note = "ρ_B>1 ⇒ 1-ρ_B<0 ⇒ 自动满足（零信息）" if rb > 1 else \
               ("满足" if vb <= 1 else "★ 违反（1-ρ_B>0 且 > κ_c^n）")
        print(f"      n={n}: ρ_B={rb:12.6f}  1-ρ_B={1-rb:>13.6f}  κ_c^n=R_n={float(rows[n]['R']):.8f}"
              f"  (1-ρ_B)R_n={vb:>12.6e}  {note}")


def reading_C():
    """读法 C：ρ_C = e^{γ_ent − c₀q}（命题 AS 的率比）。γ_ent=1（BianchiLemma.md:117 口径）。"""
    print("\n===== 读法 C : ρ_C = exp(γ_ent − c₀q)，γ_ent = 1 =====")
    print("  AS 的界：R_n ≥ 1 − 2ρ_C ；O3 的 RHS：1 − ρ_C ≤ κ_c^n")
    print("  目标值 |  c0q    |   ρ_C    | 1−ρ_C   | (1−ρ_C)^{1/n}  即 RHS 允许的最小 κ_c"
          "  |  AS 下界 1−2ρ_C ")
    Rows = {}
    for c0q in (1 + LN3, 2.10, 2.60, 2.80, 3.00):
        rho = 2.718281828459045 ** (1.0 - c0q)
        print(f"  c0q={c0q:6.4f} | {rho:0.8f} | {1-rho:0.8f} | "
              + "  ".join(f"n={n}:{(1-rho)**(1.0/n):.6f}" for n in (2, 6))
              + f" | 1−2ρ_C = {1-2*rho:>8.4f}")
        Rows[c0q] = rho
    # (1−ρ_C)·κ_c^n 用 GV2 实测 R_n（两种 u 口径）
    for ufrac, lbl in ((F(5, 18), "SU3 u=5/18"), (F(5, 8), "SU2 u=5/8")):
        ser, _ = series(ufrac)
        print(f"  -- {lbl} 实测 R_n 代入 (1−ρ_C)·κ_c^n：")
        for c0q in (1 + LN3, 2.60):
            rho = Rows[c0q]
            print(f"     c0q={c0q:5.4f} ρ_C={rho:.8f}: "
                  + "  ".join(f"n={n}:{(1-rho)*float(ser[n]['R']):.6f}" for n in NS[1:]))
    # 临界调谐线：1−ρ_C = κ_c^n  ⟺  c0q = γ + ln(1/(1−κ_c^n))
    print("  -- 临界调谐线（O3 取等号）所需 c0q* = γ_ent + ln(1/(1−R_n))：")
    for ufrac, lbl in ((F(5, 18), "SU3"), (F(5, 8), "SU2")):
        ser, _ = series(ufrac)
        print(f"     {lbl}: " + "  ".join(
            f"n={n}:{1.0 + __import__('math').log(1.0/(1.0-float(ser[n]['R']))):.4f}" for n in NS[1:])
            + f"   （对照 γ+ln3 = {1+LN3:.4f}，CPB-S 窗 c₀=2.60）")
    # O3 的"精确缺什么"写法的阈值：ρ ≤ 1−ε 要 ε 多大才够喂 κ_c ≥ 0.81
    print("  -- 命题 AS 的界 R_n ≥ 1−2ρ 要能推出 κ_c(n) ≥ 0.81，需要 ρ ≤ ρ_max(n)：")
    for n in NS:
        thr = (1 - 0.81 ** n) / 2
        print(f"     n={n}: 0.81^n={0.81**n:.6f}  ⇒ ρ_max = {thr:.6f}  (ε = 1−ρ ≥ {1-thr:.6f})")
    print("  ⇒ 结论：O3 行写的'ε 与 n 无关'不够；必须 ε > 1/2（最紧的一格 n=6 要 ε ≥ "
          f"{1-(1-0.81**6)/2:.6f}）。")


def geometric_tail():
    """为什么 O3 的形状是 1−ρ ≲ κ_c^n：A_k = A₀ρ^k 的几何尾模型给出 R = |1−ρ|/(1+ρ)。"""
    print("\n===== O3 形状的来源（几何尾模型，纯代数）=====")
    print("  A_k = A₀ρ^k ⇒ R = |1−ρ|/(1+ρ)。于是")
    print("  R ≤ κ^n  ⟺  |1−ρ| ≤ κ^n(1+ρ)  ⟺  (0<κ<1) 1−ρ ≤ 2κ^n/(1+κ^n) ≈ 2κ^n")
    print("  ⇒ O3 的 `≲` 后面藏的就是这个 2 倍常数：O3 ⟺ ρ ≥ 1−2κ^n/(1+κ^n)。")
    for n in (2, 6, 20, 100):
        k = 0.80
        print(f"     κ_c=0.80, n={n:3d}: 需 ρ ≥ {(1-2*k**n/(1+k**n)):.10f}  (1−ρ ≤ {2*k**n/(1+k**n):.3e})")
    # 用 GV2 §3.4 的实测 ρ_B 反算 R，核对模型
    print("  -- 核对 GV2 §3.4 的模型（ρ_B = m_n/m_{n−1} 代入 |1−ρ|/(1+ρ)）：")
    for ufrac, lbl in ((F(5, 18), "SU3"), (F(5, 8), "SU2")):
        ser, m = series(ufrac)
        for n in (2, 3, 6):
            rho = m[n] / m[n - 1]
            model = abs(1 - rho) / (1 + rho)
            print(f"     {lbl} n={n}: ρ_B={float(rho):.9f}  模型R={float(model):.8f}  "
                  f"实测R={float(ser[n]['R']):.8f}  差={abs(float(model)-float(ser[n]['R'])):.2e}")


def scan():
    """§3：能否调谐出 κ_c(n) < 0.80？—— 能，但只在 n 小处；随 n 单调回 1。"""
    print("\n===== §3 调谐扫描：κ_c(n) 对 u 与 n =====")
    print("      u | " + " | ".join(f"κ_c({n})" for n in (2, 3, 4, 5, 6))
          + " | n=2..6 是否全 <0.80 | min")
    for u in (0.02, 0.03, 0.04, 0.05, 0.06, 0.10, 0.20, 5 / 18, 0.625):
        ser, _ = series(F(u).limit_denominator(10 ** 6))
        ks = [kc(ser[n]["R"], n) for n in (2, 3, 4, 5, 6)]
        print(f"  {u:7.4f} | " + " | ".join(f"{k:8.5f}" for k in ks)
              + f" | {'是' if all(k < 0.8 for k in ks) else '否':>18} | {min(ks):.5f}")
    print("  ⇒ 构造可得 κ_c(n)<0.80，但 (i) 只在 n≤6，(ii) 随 n 单调升回 1（GV2 §3.3 的伪信号）。")
    # 反例（对 O3 的 ⟺ 本身）：找 (u,n) 使 κ_c<1 而 1−ρ > κ_c^n
    print("  -- O3 的 ⟺ 的字面反例（LHS 真、RHS 假）：")
    hits = 0
    for u1000 in range(5, 401):
        u = F(u1000, 1000)
        ser, _ = series(u)
        for n in NS:
            R = float(ser[n]["R"])
            if R < 1 and (1 - R) > R:      # 读法 A：κ_c^n = R
                print(f"     [读A] u={float(u):.3f} n={n}: κ_c={kc(ser[n]['R'], n):.6f}<1 但 "
                      f"1−ρ={1-R:.6f} > κ_c^n={R:.6f}")
                hits += 1
                if hits >= 6:
                    break
        if hits >= 6:
            break
    print("  （读法 B 的反例在 §2 单点表里：ρ_B<1 且 1−ρ_B > R_n 的格）")


def acid():
    """§4 酸测：2D（逐项全正）与自由场（ū=27.18, ε_T≡0）。"""
    print("\n===== §4 酸测 =====")
    print("  2D：逐项非负（GV2:203, `c_R≥0` 且 `M_ℓ∈ℤ_{≥0}`）⇒ z_signed=z_abs ⇒ R_n=1 精确")
    print("      ⇒ ρ_A=1 ⇒ 读A: (1−ρ)·κ_c^n = 0·1 = 0 ≤ 1  满足（且与 β、n 无关）")
    print("      ⇒ κ_c=1.000000000000（GV2:197–199，n=2,3,6,12 全同）")
    print("  自由场：c_f→1 ⇒ u→1/N=0.5，ū=eΔ·u=27.182818（GV1:22 G9 / GU4:113 CE-1）")
    print("      ⇒ 同样逐项非负 ⇒ R_n=1 ⇒ κ_c=1；O3 读A: 0 ≤ 1 满足")
    print("      ⇒ 而那里 ε_T≡0、f₁²=0 ⇒ O3 在两格都'满足' ⇒ 无判别力")
    for u in (0.5, 0.625):
        ser, m = series(F(u).limit_denominator(10))
        rb = [float(m[n] / m[n - 1]) for n in (2, 3, 6)]
        print(f"      u={u}: 读B ρ_B = {rb} ⇒ {'>1 ⇒ 自动满足' if min(rb) > 1 else '有 <1 格'}")


if __name__ == "__main__":
    print("#" * 100)
    print("# GO3 · O3「临界调谐」检验    B_n = 1, 20, 514, 14748, 450255, 14305816   (GV2 §3.1)")
    print("#" * 100)
    r3 = block(5, "SU(3) β=2.5", F(5, 18))
    verdict_A(r3, "SU3")
    verdict_B(r3, "SU3")
    r2 = block(5, "SU(2) β=2.5", F(5, 8))
    verdict_A(r2, "SU2")
    verdict_B(r2, "SU2")
    geometric_tail()
    reading_C()
    scan()
    acid()
    print("\nDONE（纯算术；无 MC）。")
