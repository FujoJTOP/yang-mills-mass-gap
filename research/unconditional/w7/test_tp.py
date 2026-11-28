import numpy as np, sys

def test_tp():
    print("=" * 74)
    print("(TP) 现在是定理 => 数据必须满足 W_mu(n,kappa) >= 0 对一切 n, kappa")
    print("=" * 74)
    for L in (6, 8, 10, 12):
        d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
        Cz = d["Cz"]
        ks = 2 * np.pi * np.arange(L) / L
        worst = 1e30
        warg = None
        nviol = 0
        ntot = 0
        for mu in range(4):
            others = [m for m in range(4) if m != mu]
            for n in range(L):
                nC = n if n <= L - n else n - L
                for k1 in range(L):
                    for k2 in range(L):
                        for k3 in range(L):
                            kap = {others[0]: ks[k1], others[1]: ks[k2], others[2]: ks[k3]}
                            tot = 0.0 + 0.0j
                            for a in range(L):
                                za = a if a <= L - a else a - L
                                for b in range(L):
                                    zb = b if b <= L - b else b - L
                                    for c in range(L):
                                        zc = c if c <= L - c else c - L
                                        z = [0, 0, 0, 0]
                                        z[mu] = nC
                                        z[others[0]] = za
                                        z[others[1]] = zb
                                        z[others[2]] = zc
                                        ph = (kap[others[0]] * za
                                              + kap[others[1]] * zb
                                              + kap[others[2]] * zc)
                                        idx = tuple((zz % L) for zz in z)
                                        tot += Cz[idx] * np.exp(1j * ph)
                            re = tot.real
                            ntot += 1
                            if re < worst:
                                worst = re
                                warg = (mu, nC, k1, k2, k3)
                            if re < -1e-12:
                                nviol += 1
        print("  L=%2d: min W = %+.4e at (mu,n,k)=%s ; 违反 %d / %d"
              % (L, worst, warg, nviol, ntot))
        sys.stdout.flush()

if __name__ == "__main__":
    test_tp()
