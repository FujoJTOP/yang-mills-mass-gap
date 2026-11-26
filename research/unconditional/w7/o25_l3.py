"""O25-L3: Delta(L) trend, 4D SU(2) Wilson, beta=2.5, mode A (k=0 uniform, |v|=1).

Same observable / same conditioning slice as O11-PHYS (o11_lat.py):
  V = S_W + Phi,  Phi = -log J (flat-chart Haar Jacobian),  X_link = log U_link
  kappa = d_v^2 V,  F = d_v V,  Delta = E[kappa|<v,X>=0] - Var(F|<v,X>=0) = G''(0)
  v_{l,3} = 1/sqrt(4 L^4)  (k=0 uniform z-mode, unit norm)  -> same mode at every L.

Fast path: SU(2) as unit quaternions, per-proposal work restricted to the 6
plaquettes touching the updated link; V validated against o11_physical.Vtot.

usage: py -3.12 o25_l3.py <L> <seed> <nsweep> <eps> [tag]
"""
import sys, time, math
import numpy as np
import o11_physical as o

L = int(sys.argv[1]); seed = int(sys.argv[2]); ns = int(sys.argv[3])
eps = float(sys.argv[4]); tag = sys.argv[5] if len(sys.argv) > 5 else "A"
BETA = 2.5
H = 1e-3

QC = np.array([1.0, -1.0, -1.0, -1.0])


def qmul(a, b):
    """SU(2) matrix product in the q0 I + i q.sigma chart: vector part a0 b + b0 a - a x b."""
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack((a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                     a0 * b1 + a1 * b0 - a2 * b3 + a3 * b2,
                     a0 * b2 + a2 * b0 - a3 * b1 + a1 * b3,
                     a0 * b3 + a3 * b0 - a1 * b2 + a2 * b1), axis=-1)


def qexp(x):
    r = np.sqrt((x * x).sum(-1))
    rs = np.where(r > 1e-10, r, 1.0)
    s = np.where(r > 1e-10, np.sin(r / 2) / rs, 0.5)
    return np.concatenate([np.cos(r / 2)[..., None], s[..., None] * x], axis=-1)


def phi_r(r):
    rs = np.where(r > 1e-10, r / 2.0, 1.0)
    q = np.where(r > 1e-10, np.sin(r / 2) / rs, 1.0)
    return -2.0 * np.log(np.abs(q))


def phi_scalar(r):
    if r <= 1e-10:
        return 0.0
    return -2.0 * math.log(abs(math.sin(0.5 * r) / (0.5 * r)))


def trq(q, tab, rows):
    A = q[tab[rows, 0]]; Bq = q[tab[rows, 1]]
    C = q[tab[rows, 2]] * QC; D = q[tab[rows, 3]] * QC
    return qmul(qmul(qmul(A, Bq), C), D)[..., 0]


def Vparts(x, tab, beta):
    """return (S_W part, Phi part) separately, so E[H]=<d^2 S_W> and <Phi''> can be split."""
    q = qexp(x)
    tr = trq(q, tab, np.arange(len(tab)))
    return beta * np.sum(1.0 - tr), phi_r(np.sqrt((x * x).sum(-1))).sum()


def Vtotq(x, tab, beta):
    A, Bq = Vparts(x, tab, beta)
    return A + Bq


def chain(tab, near, beta, v, nsweep, eps, rng, burn, h=H):
    nlink = v.shape[0]
    npl = len(tab)
    tabL = tab.tolist()
    nearL = [a.tolist() for a in near]
    x = rng.normal(0, 0.6, (nlink, 3))
    x -= (x * v).sum() * v
    q = qexp(x)
    r = np.sqrt((x * x).sum(-1))
    tr = trq(q, tab, np.arange(npl))
    ph = phi_r(r)
    Fs = np.empty(nsweep); Ks = np.empty(nsweep)
    KA = np.empty(nsweep); KB = np.empty(nsweep); FAa = np.empty(nsweep); FBb = np.empty(nsweep)
    nacc = 0
    tot = nsweep + burn
    for sw in range(tot):
        idx = rng.integers(0, nlink, nlink)
        dl = rng.uniform(-eps, eps, (nlink, 3))
        for i in range(nlink):
            l = int(idx[i]); d = dl[i]
            c = d[2] / nlink                      # exact projection: z -> z - c on every link
            xl = x[l] + d
            x[l] = xl
            x[:, 2] -= c
            rows = nearL[l]
            qq = qexp(x[tab[rows]])            # (6,4,3) -> (6,4,4)
            A = qq[:, 0]; Bq = qq[:, 1]; C = qq[:, 2] * QC; D = qq[:, 3] * QC
            trn = qmul(qmul(qmul(A, Bq), C), D)[:, 0]
            rl = math.sqrt(xl[0] * xl[0] + xl[1] * xl[1] + xl[2] * xl[2])
            dS = beta * float(np.sum(tr[rows] - trn)) + phi_scalar(rl) - float(ph[l])
            if dS <= 0.0 or rng.random() < math.exp(-dS):
                tr[rows] = trn
                ph[l] = phi_scalar(rl)
                nacc += 1
            else:
                x[l] = xl - d
                x[:, 2] += c
        if sw < burn:
            continue
        Ap, Bp = Vparts(x + h * v, tab, beta)
        A0, B0 = Vparts(x, tab, beta)
        Am, Bm = Vparts(x - h * v, tab, beta)
        kA = (Ap - 2 * A0 + Am) / h ** 2
        kB = (Bp - 2 * B0 + Bm) / h ** 2
        FA = (Ap - Am) / (2 * h); FB = (Bp - Bm) / (2 * h)
        Fs[sw - burn] = FA + FB
        Ks[sw - burn] = kA + kB
        KA[sw - burn] = kA; KB[sw - burn] = kB
        FAa[sw - burn] = FA; FBb[sw - burn] = FB
    out = dict(F=Fs, K=Ks, KA=KA, KB=KB, FA=FAa, FB=FBb)
    return out, nacc / (tot * nlink)


if __name__ == "__main__":
    tab, near = o.lat_tables(L)
    nlink = 4 * L ** 4
    # ---- validation of the quaternion energy against o11_physical.Vtot ----
    rngv = np.random.default_rng(4242)
    xv = rngv.normal(0, 0.5, (nlink, 3))
    dv = Vtotq(xv, tab, BETA) - o.Vtot(xv, tab, BETA)
    trv = np.abs(trq(qexp(xv), tab, np.arange(len(tab))) - o.tr_plaq(o.expi(xv), tab, np.arange(len(tab)))).max()
    v = np.zeros((nlink, 3)); v[:, 2] = 1.0 / math.sqrt(nlink)
    t0 = time.time()
    out, acc = chain(tab, near, BETA, v, ns, eps, np.random.default_rng(seed), burn=max(400, ns // 10))
    wall = time.time() - t0
    Ks, Fs = out['K'], out['F']
    d0, sg = o.jk_lat(Ks, Fs, nblk=20)
    aK = o.ac(Ks, min(300, ns // 4)); aF = o.ac(Fs, min(300, ns // 4))
    tauK = 1 + 2 * aK[1:].sum(); tauF = 1 + 2 * aF[1:].sum()
    neff = len(Ks) / max(tauK, tauF, 1.0)
    V4 = L ** 4
    print(f"L={L} seed={seed} eps={eps} nsweep={ns} tag={tag} acc={acc:.3f} "
          f"E[kappa]={Ks.mean():+.4f} varF={Fs.var():.4f} E[F]={Fs.mean():+.4f} "
          f"tauK={tauK:.1f} tauF={tauF:.1f} Neff={neff:.0f} "
          f"Delta={d0:+.4f} sig_jk={sg:.4f} ratio={d0/max(sg,1e-12):+.2f} "
          f"EK_L4={Ks.mean()/V4:+.5f} VF_L4={Fs.var()/V4:+.5f} "
          f"dE={out['KA'].mean():+.4f} dPhi={out['KB'].mean():+.4f} "
          f"vFA={out['FA'].var():.4f} vFB={out['FB'].var():.4f} "
          f"dvV={dv:.2e} dtr={trv:.2e} t={wall:.1f}s "
          f"ms_sweep={1000*wall/ns:.2f}", flush=True)
    np.save(f"o25_L{L}_{tag}_{seed}.npy", np.array([Ks, Fs, out['KA'], out['KB'], out['FA'], out['FB']]))
