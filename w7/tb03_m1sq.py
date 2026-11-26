"""tb03 — L=1（单 plaquette, SU(2), beta=2.5）上 q11 证书值的**算术**与**方向**对账。

只做三件事，不做新物理：
  [1] S - D = G'' 对账（三个 theta）
  [2] 两种读法给数：G 读作 m1^2 -> m1 = sqrt(G)；G 读作 m1 -> m1^2 = G^2
  [3] theta 单调 + beta 窗口余量（beta=2.5 vs beta*=4.3066...）
无采样、无求积、纯 mpmath 高精度算术。
"""
import mpmath as mp

mp.mp.dps = 60
M = mp.mpf

# q11_cert.log / r10_audit_cert 的证书值（S, D, S-D）
TH = [
    ("0",    "0.729924128620548449",    "0",                      "0.729924128620548449"),
    ("pi/2", "0.5581836427245872055",   "0.0152397203708974417",  "0.5429439223536897638"),
    ("pi",   "0.2051152971230731480",   "0.0691205860810849369",  "0.1359947110419882110"),
]
BETA, BETA_STAR = M("2.5"), M("4.30660785255975887488820795413835478236042")

print("=" * 78)
print("tb03 — L=1 (single plaquette, SU(2), beta=2.5): G''=S-D, 两种读法, 方向")
print("=" * 78)

print("\n[1] G'' = S - D 对账   (证书值, q11 + r10)")
print(f"{'theta':>6} {'S':>24} {'D':>24} {'S-D(证书)':>24} {'|S-D-证书|':>10}")
for th, s, d, g in TH:
    S, D, G = M(s), M(d), M(g)
    print(f"{th:>6} {mp.nstr(S,20):>24} {mp.nstr(D,20):>24} {mp.nstr(S-D,20):>24} "
          f"{mp.nstr(abs(S-D-G),3):>10}")
    assert abs(S - D - G) < mp.mpf("1e-18"), th

print("\n[2] 两种读法（同一批数，口径不同）")
print(f"{'theta':>6} {'G''=S-D':>24} {'m1=sqrt(G'')':>22} {'m1^2 (=G''^2, 若G''是m1)':>26}")
rows = []
for th, s, d, g in TH:
    G = M(g)
    r = mp.sqrt(G)
    rows.append((th, G, r, G * G))
    print(f"{th:>6} {mp.nstr(G,19):>24} {mp.nstr(r,19):>22} {mp.nstr(G*G,19):>26}")
print("  注：G'' 在语料里的身份是 1/chi（Ward），且 1/chi >= m1^2（调和均值 >= 最小）")
print("      => 把 G'' 当 m1^2 读，它天然是 m1^2 的**上**界形，不是下界形（Q6:93, R5:86）。")

print("\n[3] theta 单调（证书值，beta=2.5）")
print(f"{'pair':>12} {'dG''':>22} {'符号':>6}")
for (a, ga, _, _), (b, gb, _, _) in zip(rows[:-1], rows[1:]):
    d = gb - ga
    print(f"{a + '->' + b:>12} {mp.nstr(d,19):>22} {'-' if d < 0 else '+':>6}")
    assert d < 0

print("\n[4] beta 窗口（theta=pi 支）")
print(f"  beta = {mp.nstr(BETA,10)} , beta* = {mp.nstr(BETA_STAR,20)}")
print(f"  beta*/beta - 1 = {mp.nstr(BETA_STAR/BETA - 1,10)}   (余量 {mp.nstr(100*(BETA_STAR/BETA-1),6)} %)")
assert BETA < BETA_STAR

print("\n[5] 幅值体检（防把 0.73 误当'大数'）")
print(f"  max G'' = {mp.nstr(rows[0][1],19)} , min G'' = {mp.nstr(rows[-1][1],19)} , "
      f"比值 = {mp.nstr(rows[0][1]/rows[-1][1],10)}")
print(f"  sqrt(G''_max) = {mp.nstr(rows[0][2],19)} ; 相差 1 个数量级内 => 无标度分离可利用")
print("\nEXIT=0")
