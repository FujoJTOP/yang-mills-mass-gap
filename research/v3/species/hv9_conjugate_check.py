"""HV9-CONJUGATE 机检：确定性算术，零采样、零 MC。
[1] 单 plaquette：符号（-Var 还是 +Var），独立 Gauss-Legendre 求积 + 中心差分
[2] 归一化因子：广延 O=Phi 与密度 O=Phi/V 的 ±Var(O) 恒等式差一个 V
[3] 预算算术：w0=1/(4N^4) vs 窗内地板 r0/M0max（精确有理数）
[4] 见证：同一积分值下靶真/靶假两族并存（类 C 内，QI-1 逐位）
[5] 传递因子：F<1 时 F*w0 < w0（∴ 任何 F<1 的结论弱于靶）
"""
import numpy as np
from fractions import Fraction as F

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s)

# ---------- 单 plaquette：SU(2)，u = (1/2)trU = cos θ，权 ∝ e^{β cosθ} sin²θ dθ ----------
XG, WG = np.polynomial.legendre.leggauss(4000)
TH = 0.5 * np.pi * (XG + 1.0)
WT = 0.5 * np.pi * WG
U = np.cos(TH)

def mom(beta):
    """返回 (E[u], Var(u))；e^{βcosθ} 用 e^{β(u-1)} 归一（常数抵消）。"""
    d = np.exp(beta * (U - 1.0)) * np.sin(TH) ** 2
    Z = np.sum(WT * d)
    m1 = np.sum(WT * d * U) / Z
    m2 = np.sum(WT * d * U * U) / Z
    return m1, m2 - m1 * m1

say("== [1] 符号：O_lam = lam - u,  lam in {0,1,2,3} ==")
say("   beta   lam   dE/db(FD)        Var(O)          -Var(O)          FD+Var         FD-Var")
h = 1e-3
rows = []
for beta in (0.5, 1.0, 2.5, 6.0):
    m1p, _ = mom(beta + h); m1m, _ = mom(beta - h)
    _, var = mom(beta)
    dfd = -(m1p - m1m) / (2 * h)              # d/dβ E[lam - u] = -d/dβ E[u]
    for lam in (0, 1, 2, 3):
        rows.append((beta, lam, dfd, var, -var, abs(dfd + var), abs(dfd - var)))
for r in rows:
    say("   %5.2f  %2d  %+.10f  %+.10f  %+.10f  %+.2e  %+.2e" % r)
ok1 = all(r[5] < 1e-5 for r in rows) and all(r[6] > 1e-3 for r in rows)
say("   => dE/db = -Var 逐点成立 (max|FD+Var| = %.2e)；dE/db = +Var 被否 (min|FD-Var| = %.2e)"
    % (max(r[5] for r in rows), min(r[6] for r in rows)))
say("   sign_check:", "PASS" if ok1 else "FAIL")

say("")
say("== [2] 归一化因子：广延 O=lam-u 与密度 O=(lam-u)/V ==")
_, var1 = mom(2.5)
ok2 = True
for V in (1, 2, 3, 5):
    lhs_ext = -var1                      # d/dβ E[O_ext] = -Var(O_ext)
    rhs_ext = -var1
    lhs_den = -var1 / V                  # d/dβ E[O_den] = -Var(O_ext)/V
    rhs_den = -(var1 / V ** 2) * V       # -V*Var(O_den)
    chk = abs(lhs_den - rhs_den) < 1e-15
    ok2 &= chk
    say("   V=%d  Var(O_ext)=%.10f  Var(O_den)=%.10f  d/dbE[O_den]=%.10f  -V*Var(O_den)=%.10f  %s"
        % (V, var1, var1 / V ** 2, lhs_den, rhs_den, "OK" if chk else "MISMATCH"))
say("   => 逐字 'dE/db = -Var(O)' 只在广延归一成立；密度归一再乘 V（因子 V 落在这里）")
say("   norm_check:", "PASS" if ok2 else "FAIL")

say("")
say("== [3] 预算算术（精确有理数）==")
ok3 = True
for N in (2, 3, 5):
    r0 = F(1, 4) * F(1, (2 * N) ** 2)            # ¼(2N)^{-2}
    wP3_2N2 = r0 / (2 * N ** 2)                  # HZ1 的 M0 <= 2N^2
    wP3_pop = r0 / (F(N ** 2, 4))                # Popoviciu M0 <= N^2/4
    w0 = F(1, 4 * N ** 4)                        # 任务书的 w0
    ok3 &= (wP3_pop == w0) and (w0 / wP3_2N2 == 8)
    say("   N=%d  r0=%s  M0max(2N^2)=%d -> 窗内地板=%s (1/(32N^4)=%s)  M0max(N^2/4)-> %s = w0=%s  比=%s"
        % (N, r0, 2 * N ** 2, wP3_2N2, F(1, 32 * N ** 4), wP3_pop, w0, wP3_2N2 and (w0 / wP3_2N2)))
say("   => Popoviciu 收紧后 窗内地板 = 靶常数 = 1/(4N^4) 逐位相同 => 预算因子 = 1 (ln1 = 0)")
say("   budget_check:", "PASS" if ok3 else "FAIL")

say("")
say("== [4] 见证：同一积分值（M0 恒 1 => ∫M0 db = b2-b1）下靶真 / 靶假 ===")
w0 = F(1, 64)

def two_atom(w, a, b):
    """nu = w d_a + (1-w) d_b, M0=1. 返回 (w1, delta1, eps, Delta, QI1, M0)"""
    w = F(w); a = F(a); b = F(b)
    M0 = F(1); s1 = a; s2 = b; eps = s2 - s1
    mean = w * a + (1 - w) * b
    d1 = mean - s1
    Delta = F(0)                       # 次一支撑点之上无质量
    QI1 = 1 - d1 / eps + Delta / eps
    return w, d1, eps, Delta, QI1, M0

ok4 = True
for tag, w in (("A(靶真)", F(9, 10)), ("B(靶假)", F(1, 100))):
    w1, d1, eps, Delta, QI1, M0 = two_atom(w, 1, 2)
    in_C = (M0 > 0) and (eps > 0) and (0 < d1) and (d1 < eps)   # ④ 在 δ∈(0,ε) 上成立
    bud = d1 <= (1 - w0) * eps + Delta / eps                     # 预算形
    tgt = w1 >= w0
    ok4 &= (QI1 == w1) and (bud == tgt) and in_C
    say("   %s: w1=%s(=QI-1 %s) d1=%s eps=%s Delta=%s 类C=%s 预算形=%s 靶=%s"
        % (tag, w1, QI1, d1, eps, Delta, in_C, bud, tgt))
say("   ∫_{b1}^{b2} M0 db = b2-b1 两族相同（且 <O>(b) = c-b 整条剖面相同）而靶真/靶假并存")
say("   => 积分值不是 w1 的充分统计量；不是'端点数据缺'，是'字母表不含原子'")
say("   witness_check:", "PASS" if ok4 else "FAIL")

say("")
say("== [5] 传递因子：以'种子 + 因子'为骨架的证明，F<1 时必然弱于靶 ==")
ok5 = True
for k in (1, 5, 10, 20):
    Fk = 1 - F(1, 2 ** k)
    ok5 &= (Fk * w0 < w0)
    say("   F=1-2^-%2d=%s  F*w0=%s < w0=%s  (缺口 = w0/2^%d = %s)"
        % (k, Fk, Fk * w0, w0, k, w0 / 2 ** k))
say("   => 需要 F>=1 才能用 P3 种子（种子只给界不给值）；F=1 即 net-non-decay")
say("   transfer_check:", "PASS" if ok5 else "FAIL")

say("")
allok = ok1 and ok2 and ok3 and ok4 and ok5
say("ALL CHECKS PASS" if allok else "SOME CHECK FAILED")
with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
