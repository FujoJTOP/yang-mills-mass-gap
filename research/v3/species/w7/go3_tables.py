# -*- coding: utf-8 -*-
"""GO3 · 输出 GO3-TUNING.md §2/§3/§4 表格用的 Markdown 行（精确有理数，供人眼核对）。
只打印，不写文件；无 MC、无随机。"""
from fractions import Fraction as F
import math

B = {1: 1, 2: 20, 3: 514, 4: 14748, 5: 450255, 6: 14305816}
NS = [1, 2, 3, 4, 5, 6]


def series(u):
    m, S, Z, out = {}, F(0), F(0), {}
    for n in NS:
        m[n] = u ** n * B[n]
        Z += m[n]
        S += (F(1) if n % 2 else F(-1)) * m[n]
        out[n] = (m[n], Z, abs(S), abs(S) / Z)
    return out, m


def table(uf, label):
    ser, m = series(uf)
    print(f"\n#### {label}  (u = {uf} = {float(uf):.10f})\n")
    print("| `n` | `m_n` | `Z_n` | `|S_n|` | `R_n=|S_n|/Z_n` | `κ_c(n)=R_n^{1/n}` "
          "| 读A `1−ρ_A` | 读A `(1−ρ_A)κ_c^n` | 读A `(1−ρ_A)/κ_c^n` "
          "| 读B `ρ_B=m_n/m_{n−1}` | `1−ρ_B` | 读B `(1−ρ_B)κ_c^n` |")
    print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    prev = None
    for n in NS:
        mn, Z, Sm, R = ser[n]
        k = float(R) ** (1 / n)
        one = 1 - R
        if prev is None:
            rb = ob = vb = None
        else:
            rb = m[n] / prev
            ob = 1 - rb
            vb = ob * R
        cell = lambda x, d=6: "--" if x is None else f"{float(x):.{d}f}"
        print(f"| {n} | {float(mn):.6f} | {float(Z):.6f} | {float(Sm):.6f} | {float(R):.8f} "
              f"| **{k:.8f}** | {float(one):.6f} | {float(one*R):.6f} | {float(one/R):.6f} "
              f"| {cell(rb,9)} | {cell(ob,6)} | {cell(vb,6) if vb is None else f'{float(vb):.6f}'} |")
        prev = m[n]


def readC():
    print("\n#### 读法 C 表：`ρ_C = e^{γ_ent−c₀q}`，`γ_ent=1`（`BianchiLemma.md:39,111,116`）\n")
    print("| `c₀q` | `ρ_C` | `1−ρ_C` | `(1−ρ_C)^{1/2}` | `(1−ρ_C)^{1/6}` | AS 下界 `1−2ρ_C` "
          "| SU3 `(1−ρ_C)/R_6` | SU2 `(1−ρ_C)/R_6` |")
    print("|---:|---:|---:|---:|---:|---:|---:|---:|")
    ser3, m3 = series(F(5, 18))
    ser2, m2 = series(F(5, 8))
    R3, R2 = float(ser3[6][3]), float(ser2[6][3])
    for c0q in (2.0, 1 + math.log(3), 2.5, 2.5986, 2.8):
        rho = math.exp(1 - c0q)
        print(f"| {c0q:.4f} | {rho:.6f} | {1-rho:.6f} | {(1-rho)**0.5:.6f} | {(1-rho)**(1/6):.6f} "
              f"| {1-2*rho:.4f} | {(1-rho)/R3:.5f} | {(1-rho)/R2:.5f} |")
    print(f"\n（R_6 = {R3:.8f} SU3 / {R2:.8f} SU2 ⇒ `(1−ρ_C)/R_n ≤ 1` 即 O3 的紧形式成立）")


def constraints():
    print("\n#### AS 的常数 / 判据阈值\n")
    print("| `k=|C₀|` | `2^{−1/k}`（AS 给的 `κ_c` 下界，ρ=1/3） | 过 `κ_c≥0.81`? |")
    print("|---:|---:|:---|")
    for k in (1, 2, 3, 4, 5, 6, 12):
        v = 2 ** (-1 / k)
        print(f"| {k} | {v:.6f} | {'是' if v >= 0.81 else '**否**'} |")
    print("\n| `n` | `3^{−1/n}`（`ρ≤1/3` ⇒ `κ_c^n≥1/3`） | 过 `κ_c≥0.81`? |")
    print("|---:|---:|:---|")
    for n in NS:
        v = 3 ** (-1 / n)
        print(f"| {n} | {v:.6f} | {'是' if v >= 0.81 else '**否**'} |")
    print("\n| `n` | `ρ_max=(1−0.81^n)/2`（AS 推 `κ_c(n)≥0.81` 所需的 ρ 上限） | `ε_min=1−ρ_max` |")
    print("|---:|---:|---:|")
    for n in NS:
        t = (1 - 0.81 ** n) / 2
        print(f"| {n} | {t:.6f} | {1-t:.6f} |")


def rhoB_threshold():
    print("\n#### 读 B 的“自动满足”适用区（`ρ_B>1 ⟺ u > B_{n−1}/B_n`）\n")
    print("| `n` | `B_{n−1}/B_n` | 读B 在 u 高于此值时自动满足 |")
    print("|---:|---:|:---|")
    for n in NS[1:]:
        print(f"| {n} | {B[n-1]/B[n]:.6f} | 是 |")
    # 读 B 的违反格（ρ_B<1 且 1−ρ_B > R_n）
    print("\n读 B 的违反格（`ρ_B<1` 且 `1−ρ_B > R_n`）：")
    hits = 0
    for u1000 in range(5, 201):
        u = F(u1000, 1000)
        ser, m = series(u)
        for n in NS[1:]:
            rb = m[n] / m[n - 1]
            if rb < 1 and (1 - rb) > ser[n][3]:
                print(f"  u={float(u):.3f} n={n}: ρ_B={float(rb):.6f} 1−ρ_B={float(1-rb):.6f} "
                      f"R_n={float(ser[n][3]):.6f}  ⇒ 违反")
                hits += 1
                break
        if hits >= 5:
            break


if __name__ == "__main__":
    table(F(5, 18), "SU(3)，`β=2.5`，`u=5/18=β/9`")
    table(F(5, 8), "SU(2)，`β=2.5`，`u=5/8=0.625=β/4`")
    readC()
    constraints()
    rhoB_threshold()
    print("\nDONE")
