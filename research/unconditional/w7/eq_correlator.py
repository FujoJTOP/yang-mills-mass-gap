import numpy as np, math, sys

print("=" * 78)
print("BIG STEP: 等时相关函数 Gamma_eq 的性质（对每个 L）")
print("=" * 78)
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]; vq = float(d["varq"])
    print()
    print("L = %d   Var(q) = %.8f" % (L, vq))
    # 1) 等时超立方各向同性
    e = [tuple((1 if i == k else 0) for i in range(4)) for k in range(4)]
    vals = [Cz[z] for z in e]
    print("  |z|=1 四个方向： %s" % "  ".join("%.6e" % v for v in vals))
    print("     相对散布 = %.2e" % ((max(vals) - min(vals)) / np.mean(vals)))
    pairs = [(1,1,0,0),(1,0,1,0),(1,0,0,1),(0,1,1,0),(0,1,0,1),(0,0,1,1)]
    vals2 = [Cz[z] for z in pairs]
    print("  |z|=sqrt2 六个方向：%s" % "  ".join("%.5e" % v for v in vals2))
    print("     相对散布 = %.2e" % ((max(vals2) - min(vals2)) / np.mean(vals2)))
    # 2) 沿轴衰减率（用 (n,0,0,0)）
    ax = [Cz[tuple((n if i == 0 else 0) for i in range(4))] for n in range(0, L // 2 + 1)]
    print("  沿轴 Gamma_eq(n e1)/Var(q)：%s" % "  ".join("%.3e" % (v / vq) for v in ax))
    rates = []
    for n in range(0, min(3, L // 2)):
        if ax[n] > 0 and ax[n + 1] > 0:
            rates.append(math.log(ax[n] / ax[n + 1]))
    print("     逐步衰减率：%s" % "  ".join("%.3f" % r for r in rates))
    # 3) w(0) 与 3D 求和
    r2 = np.zeros((L, L, L, L))
    idx = np.arange(L); dm = np.minimum(idx, L - idx)
    for mu in range(4):
        a = dm.reshape([L if k == mu else 1 for k in range(4)])
        r2 = r2 + a ** 2
    r2 = np.broadcast_to(r2, (L, L, L, L)).copy()
    eq = np.zeros((L, L, L))
    for z2 in range(L):
        for z3 in range(L):
            for z4 in range(L):
                zz = [0,
                      z2 if z2 <= L - z2 else z2 - L,
                      z3 if z3 <= L - z3 else z3 - L,
                      z4 if z4 <= L - z4 else z4 - L]
                eq[z2, z3, z4] = Cz[tuple((v % L) for v in zz)]
    w0 = eq.sum()
    print("  w(0) = Sum_{z_perp} Gamma_eq = %.8f  ；  w(0)/Var(q) = %.4f" % (w0, w0 / vq))
    print("  Gamma_eq 的负值个数 = %d / %d" % (int((eq < 0).sum()), eq.size))
    sys.stdout.flush()
