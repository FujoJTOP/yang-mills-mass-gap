import sys, time, argparse
import numpy as np
sys.path.insert(0, "research/unconditional/w7")
import block_mc as M

def q_orientations(U):
    """return array shape (6, L,L,L,L) of q_P = 0.5*tr(U_P) for the 6 orientations"""
    L = U[0][0].shape[0]
    out = []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            P = M.qmul(M.qmul(M.qmul(U[mu], M.shift(U[nu], mu, -1)),
                              M.qdag(M.shift(U[mu], nu, -1))), M.qdag(U[nu]))
            out.append(0.5 * M.qtr2(P))
    return np.array(out)


def run(L, beta, sweeps, therm, every, eps, seed):
    rng = np.random.default_rng(seed)
    idxs = tuple(np.meshgrid(*[np.arange(L)] * 4, indexing="ij"))
    masks = [((sum(idxs) % 2) == p) for p in (0, 1)]
    U = M.init_u(L, rng)
    Qs = []           # 6 orientation sums
    for sw in range(therm + sweeps):
        U, x_, y_ = M.sweep(U, beta, L, eps, rng, masks)
        if sw >= therm and (sw - therm) % every == 0:
            Qs.append(q_orientations(U).sum(axis=(1, 2, 3, 4)))
    Qs = np.array(Qs)              # (n, 6)
    n = Qs.shape[0]
    Qc = Qs - Qs.mean(axis=0)
    Mmat = (Qc.T @ Qc) / (n * L ** 4)     # M_{AB} = Cov(Q^A,Q^B)/L^4
    return Qs.mean(axis=0), Mmat, n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=4)
    ap.add_argument("--beta", type=float, default=2.5)
    ap.add_argument("--sweeps", type=int, default=1500)
    ap.add_argument("--therm", type=int, default=400)
    ap.add_argument("--every", type=int, default=2)
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    means, Mmat, n = run(a.L, a.beta, a.sweeps, a.therm, a.every, a.eps, a.seed)
    np.set_printoptions(precision=6, suppress=False, linewidth=200)
    print("L=%d beta=%.3f n=%d" % (a.L, a.beta, n))
    print("  <Q_A>/L^4 =", means / a.L ** 4)
    print("  M_{AB} = Cov(Q^A,Q^B)/L^4 :")
    print(Mmat)
    off = Mmat[~np.eye(6, dtype=bool)]
    print("  min off-diagonal M_AB = %+.6e   (need >= 0)" % off.min())
    print("  M_AA (diagonal) =", np.diag(Mmat))
    print("  row sums Sum_B M_AB =", Mmat.sum(axis=1))
    print("  M_AA <= rowsum ?", np.all(np.diag(Mmat) <= Mmat.sum(axis=1) + 1e-12))
    print("  chi_clover = sum(M)/36 = %.8e" % (Mmat.sum() / 36.0))
    print("  max M_AA = %.8e   (bound 1)" % np.diag(Mmat).max())
