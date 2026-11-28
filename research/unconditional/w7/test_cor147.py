import numpy as np, sys

print("=" * 72)
print("Corollary 14.7 直接检验：Gamma(0, n e_mu) >= 0 对一切 mu, n")
print("  （这是 (TP) 的推论；若显著失败，(TP) 被否证）")
print("=" * 72)
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]
    print()
    print("L =", L)
    for mu in range(4):
        vals = []
        for n in range(L):
            idx = [0, 0, 0, 0]
            idx[mu] = n
            vals.append(Cz[tuple(idx)])
        mn = min(vals)
        negs = [(n, v) for n, v in enumerate(vals) if v < -1e-12]
        print("  mu=%d: min = %+.4e ；负值 %d 个 %s"
              % (mu, mn, len(negs), ("" if not negs else str(negs[:4]))))

print()
print("=" * 72)
print("另：数据是否具有超立方对称性？（12-plaquette 只保留 12 平面的对称）")
print("=" * 72)
for L in (10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]
    print()
    print("L =", L)
    tests = [
        ((1, 1, 0, 0), "12 平面内对角"),
        ((1, 0, 0, 0), "方向 1"),
        ((0, 1, 0, 0), "方向 2"),
        ((0, 0, 1, 0), "方向 3"),
        ((0, 0, 0, 1), "方向 4"),
        ((1, 1, 1, 0), "1,2,3"),
        ((0, 0, 1, 1), "3,4"),
    ]
    for z, lbl in tests:
        print("   Cz%s=%.6e   (%s)" % (z, Cz[z], lbl))
