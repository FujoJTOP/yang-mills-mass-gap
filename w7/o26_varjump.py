# -*- coding: utf-8 -*-
"""o26_varjump.py -- O26: Var(G_v+Phi'|0) 单 plaquette -> L=2 的广延性判决.
 (1) 周期 4D L=2 的 plaquette 计数; (2) Var/N_p, Delta/N_p; (3) L^2 后符号.
数据源: research/v3/species/O11-PHYS.md (单 plaquette 精确+2群MC; L=2 模式A 2群).
定义源: GPP-STRUCT:24,49,62 (eq:split / ||v||^-4); ch05:208,233-235 (Lemma B)."""
import itertools, math
from collections import Counter


def plaquettes(L):
    def canon(site, d):
        return (tuple(s % L for s in site), d)
    terms = []
    for mu, nu in itertools.combinations(range(4), 2):
        for x in itertools.product(range(L), repeat=4):
            xm = tuple((x[i] + (1 if i == mu else 0)) % L for i in range(4))
            xn = tuple((x[i] + (1 if i == nu else 0)) % L for i in range(4))
            terms.append(frozenset([canon(x, mu), canon(xm, nu), canon(x, nu), canon(xn, mu)]))
    c = Counter(terms)
    return len(terms), len(c), sorted(set(c.values()))


VAR1, D1, EK1 = 0.069121, 0.1359947, 0.0
EK2, VAR2A, VAR2B, D2, S2 = 1.780, 2.1443, 1.8579, -0.2209, 0.1198
VAR2 = 0.5 * (VAR2A + VAR2B)

print("=" * 76)
print("O26  Var-jump: 单 plaquette -> L=2 广延性判决   [E] 复算, [T] 定义")
print("=" * 76)

print("\n[1] PLAQUETTE 计数 (周期 4D, mu<nu, 所有位点)")
print(f"{'L':>3} {'项数 6L^4':>10} {'不同 link 集合':>14} {'多重度':>8} {'6L^4/2':>8}")
NP = {}
for L in (1, 2, 4):
    nt, nd, mult = plaquettes(L)
    NP[L] = (nt, nd, mult)
    print(f"{L:>3} {nt:>10} {nd:>14} {str(mult):>8} {nt//2:>8}")
n96 = NP[2][0]
print(f"  判: L=2 -> 6L^4 = {n96} 个 plaquette, 每个 link 集合多重度 {NP[2][2]}")
print("      => 96 是每一项各算一次, 无重复; '48 = 6L^4/2' 不是 plaquette 数.")
print("      (48 = 96 条 plaquette 里每一对反定向? 否: x 与 x+mu+nu 的 link 集合互不相同, 见多重度)")
print(f"  旁证 L=1: 6L^4 = 6, 不同 link 集合 {NP[1][1]} => '单 plaquette 模型' != L=1 格点 (只有 1 项)")

print(f"\n[2] Var / N_p   (单 plaquette: N_p = 1, Var = {VAR1})")
print(f"  L=2  N_p = {n96}: Var = {VAR2:.4f}  ->  Var/N_p = {VAR2/n96:.5f}")
print(f"  若硬按 N_p=48 (已否证): Var/N_p = {VAR2/48:.5f}   <- 仅作对照")
print(f"  单 plaquette:            Var/N_p = {VAR1/n96*96:.5f}")
pred = VAR1 * n96
alpha = math.log(VAR2 / VAR1) / math.log(n96)
print(f"  严格广延预期 Var = {VAR1}*{n96} = {pred:.3f} vs 实测 {VAR2:.3f} (比 {VAR2/pred:.3f})")
print(f"  幂次 Var ~ N_p^alpha: alpha = {alpha:.3f} < 1  => Var 次广延 (30x 远小于 96x)")

print(f"\n[3] Delta / N_p   Delta = E[kappa|0]-Var(F|0) = G''(0)  [T] eq:split")
print(f"  单 plaquette: Delta = {D1:+.7f}   Delta/N_p = {D1/1:+.6f}   (N_p=1)")
print(f"  L=2 (N_p={n96}): Delta = {D2:+.4f} +- {S2:.4f}   Delta/N_p = {D2/n96:+.6f}")
print(f"  L=2 (N_p=48 对照):                        Delta/N_p = {D2/48:+.6f}")
print(f"  符号: +  vs  -   =>  {'同号' if D1*D2>0 else '【变号】'}  (N_p=96 或 48 都变号)")

print(f"\n[4] 乘 L^2  (||v||^2=1/L 的 L-link 平均约定, ch05:214-228, c=1/sqrt(L))")
for L, D in ((1, D1), (2, D2)):
    print(f"  L={L}: Delta = {D:+.7f}  x L^2={L**2}  =>  {D*L**2:+.5f}")
print(f"  符号: + vs -  =>  {'同号' if D1*D2>0 else '【变号】'} (乘正因子不改符号)")
print(f"  若 o11 的 L=2 已含 L^2 => 真值 {D2:+.4f}; 未含 => {D2*4:+.4f}(=4x, 更负)")

print(f"\n[5] 两项各自的标度 (Delta = E[kappa] - Var(F) 是两个和式之差)")
print(f"  E[kappa]: 单 plaquette {EK1:+.4f} (SU(2) Cartan 切片精确 0) vs L=2 {EK2:+.4f}")
print(f"  Var(F)  : 单 plaquette {VAR1:.4f} vs L=2 {VAR2:.4f} = {VAR2/VAR1:.1f}x  (N_p 比 {n96}x)")
print("  => 两项都随体积涨; Delta 是两广延量之差 => Delta 亦广延(至 O(1) 边界修正).")
print("     Delta 跨体积不可直比, Delta/N_p 才可比.  [T] (V = beta*sum_P + sum_l Phi_l 是局域和式)")

print("\n[6] W_full / a^2 Delta_0^2 归属")
print("  o11 的 V = S_W + Phi, Phi = -log J  =>  kappa = H_vv+Phi'', F = G_v+Phi'")
print("  => Delta = E[H+Phi''] - Var(G_v+Phi') = Delta_dyn + W_geo + X = G''(0)  [T]")
print("     即 W_full = W_geo + X 【完全在】 Delta 内;  a^2 Delta_0^2 【不在】(目标是右端阈值, 非测量项)")
print("\nDONE")

