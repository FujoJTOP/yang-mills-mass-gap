import sys, time, argparse
import numpy as np
sys.path.insert(0, "research/unconditional/w7")
import block_mc as M


def run(L, beta, sweeps, therm, every, eps, seed):
    rng = np.random.default_rng(seed)
    idxs = tuple(np.meshgrid(*[np.arange(L)] * 4, indexing="ij"))
    masks = [((sum(idxs) % 2) == p) for p in (0, 1)]
    U = M.init_u(L, rng)
    N = L ** 4
    qq = []
    at = ac = 0
    t0 = time.time()
    for sw in range(therm + sweeps):
        U, x_, y_ = M.sweep(U, beta, L, eps, rng, masks)
        at += x_; ac += y_
        if sw >= therm and (sw - therm) % every == 0:
            qf = M.qfield(U)
            # 用 clover 观测量：Q^data = Sum_x qf ；S = 6 Q^data
            Q = float(qf.sum() / 1.0)
            qq.append(Q)
    qq = np.array(qq)
    m = qq.mean()
    d = qq - m
    m2 = (d ** 2).mean()
    m3 = (d ** 3).mean()
    m4 = (d ** 4).mean()
    return dict(L=L, beta=beta, n=len(qq), mean=m, var=m2, m3=m3, m4=m4,
                skew=m3 / m2 ** 1.5 if m2 > 0 else 0.0,
                acc=at / max(ac, 1), sec=time.time() - t0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=6)
    ap.add_argument("--beta", type=float, default=2.5)
    ap.add_argument("--sweeps", type=int, default=4000)
    ap.add_argument("--therm", type=int, default=1000)
    ap.add_argument("--every", type=int, default=2)
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    r = run(a.L, a.beta, a.sweeps, a.therm, a.every, a.eps, a.seed)
    print("L=%d beta=%.4f n=%d  acc=%.3f  %.1fs" % (r["L"], r["beta"], r["n"], r["acc"], r["sec"]))
    print("  <Q>   = %.8e" % r["mean"])
    print("  Var(Q)= %.8e   ;  chi = Var/L^4 = %.8e" % (r["var"], r["var"] / r["L"] ** 4))
    print("  m3    = %.8e" % r["m3"])
    print("  m4    = %.8e" % r["m4"])
    print("  skew  = %.6f" % r["skew"])
    print("  dchi/dbeta = -6*m3/L^4 = %+.8e" % (-6.0 * r["m3"] / r["L"] ** 4))
    print("  => chi %s in beta" % ("DECREASES" if r["m3"] > 0 else "INCREASES"))
