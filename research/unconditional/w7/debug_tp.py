import numpy as np, sys

L = 12
d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
Cz = d["Cz"]

def mi(z):
    return z if z <= L - z else z - L

print("=== 调试：W_{mu=1}(n=6, kappa=0) 与早前 w(6) 的对比 ===")
print("Cz 是否满足 Gamma(-z)=Gamma(z)? 抽查：")
for z in [(1, 0, 0, 0), (2, 3, 0, 0), (5, 0, 0, 0)]:
    neg = tuple(((-x) % L) for x in z)
    print("   Cz%s=%.6e   Cz%s=%.6e   差=%.2e" % (z, Cz[z], neg, Cz[neg], abs(Cz[z] - Cz[neg])))

print()
print("=== 对 mu=1，n=6：两种求和方式的对比 ===")
# 方式 A：z1 固定为 signed nC，其余三维在环面上遍历（每次映射到最小镜像再查 Cz）
nC = mi(6)
totA = 0.0
totA_raw = 0.0
for a in range(L):
    for b in range(L):
        for c in range(L):
            za, zb, zc = mi(a), mi(b), mi(c)
            idx = (nC % L, za % L, zb % L, zc % L)
            totA += Cz[idx]
            idx2 = (nC % L, a, b, c)
            totA_raw += Cz[idx2]
print("  方式A（逐维映射到最小镜像）: %.6e" % totA)
print("  方式A'（只固定 z1，其余原样）: %.6e" % totA_raw)

print()
print("=== 检查 Cz 是否真的只依赖最小镜像坐标 ===")
bad = 0
for z1 in [2, 5]:
    for z2 in [3]:
        ref = Cz[(z1 % L, z2 % L, 0, 0)]
        alt = Cz[((L - z1) % L, (L - z2) % L, 0, 0)]
        print("   Cz(%d,%d,0,0)=%.6e ；Cz(%d,%d,0,0)=%.6e ；差=%.2e"
              % (z1, z2, ref, L - z1, L - z2, alt, abs(ref - alt)))
        if abs(ref - alt) > 1e-12:
            bad += 1
print("  => 若差=0，说明 Cz 已按最小镜像存储，两种索引一致")

print()
print("=== 结论线索：Sum_z Cz = chi ? ===")
print("   Sum Cz = %.8f ；chi = %.8f" % (Cz.sum(), float(d["chi_q"])))
