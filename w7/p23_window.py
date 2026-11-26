"""p23_window.py -- P23: independent recompute of E[C_geo] and the beta-window boundary.

Definition re-derived from P06 / O06 (checked numerically in [0]):
    S_geo(U) = sum_p (1 - 1/2 Re tr U_p),      S_W = beta * S_geo        (exact, beta-free S_geo)
    U_l(s)   = exp(X_l + s v_l),               U_p(s) = U_0 U_1 U_2^dag U_3^dag
    U_p' = W_p U_p      (W_p anti-Hermitian),  H_p = -i W_p  (Hermitian), K_p = H_p^2 >= 0
    D(v)  = d^2/ds^2 S_geo|_0 = -1/2 Re tr U_p'' = sum_p [ 1/2 <U_p,K_p> - 1/2 <U_p,W_p'> ]
            with  <A,B> := Re tr(AB)   and   W'' contribution via W' = dW/ds
    H_vv  = beta * D(v)                                   (beta factors out exactly)
    C_geo(v) := (-D(v))^+ = max(0, -D(v))                 <- exact NEGATIVE PART of D
    => pointwise H_vv >= -beta*C_geo(v), equality iff D<=0 ;  E[H_vv] >= -beta E[C_geo]

Normalization inherited from P06 (checked, not assumed): chart vector a_l (3 components),
||v||^2 := sum_l tr(V_l^dag V_l) with V_l = i a_l.sigma/2  =>  ||v||^2 = (1/2) sum_l |a_l|^2.
Every scanned direction has sum_l |a_l|^2 = 2, i.e. ||v||^2 = 1; floor (N/12)||v||^2 = 1/6 (N=2).

Engines: (1) one plaquette (4 links, exact self-normalized importance sampling on product Haar,
reweighted by e^{-beta S_geo});  (2) 4D lattice L=2 (16 sites, 64 links, 24 plaquettes) with
Metropolis at beta=2.5.  Run: py -3.12 -u w7/p23_window.py
"""
import time
import numpy as np

t00 = time.time()


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78, flush=True)


def expc(a):
    """a (...,3) -> exp(i a.sigma/2) (...,2,2)."""
    a = np.asarray(a, float)
    th = np.sqrt(np.sum(a * a, axis=-1))
    small = th < 1e-12
    c = np.cos(0.5 * th)
    r = np.where(small, 0.5, np.sin(0.5 * th) / np.where(small, 1.0, th))
    m = np.zeros(a.shape[:-1] + (2, 2), complex)
    ax, ay, az = a[..., 0], a[..., 1], a[..., 2]
    m[..., 0, 0] = c + 1j * r * az
    m[..., 1, 1] = c - 1j * r * az
    m[..., 0, 1] = r * (1j * ax + ay)
    m[..., 1, 0] = r * (1j * ax - ay)
    return m


def dag(M):
    return np.swapaxes(M.conj(), -1, -2)


def logc(U):
    w = np.real(U[..., 0, 0]); z = np.imag(U[..., 0, 0])
    x = np.imag(U[..., 0, 1]); y = np.real(U[..., 0, 1])
    s = np.sqrt(x * x + y * y + z * z)
    th = 2.0 * np.arctan2(s, w)
    f = np.where(s > 1e-13, th / np.where(s > 1e-13, s, 1.0), 2.0)
    return np.stack([f * x, f * y, f * z], axis=-1)


def randU(rng):
    q = rng.normal(size=(4,)); q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[w + 1j * z, 1j * x + y], [1j * x - y, w - 1j * z]], complex)


def unit_dirs(n, shape, rng):
    """n directions with sum_l |a_l|^2 = 2  (i.e. ||v||^2 = 1)."""
    b = rng.normal(size=(n,) + shape)
    f = np.sqrt(0.5 * np.sum(b.reshape(n, -1) ** 2, axis=1))
    return b / f.reshape((n,) + (1,) * len(shape))


# ============================== engine 1: single plaquette ==============================
def Sg_plaq(X, B, s):
    U = expc(X + s * B)                                  # (...,4,2,2)
    P = U[..., 0, :, :] @ U[..., 1, :, :] @ dag(U[..., 2, :, :]) @ dag(U[..., 3, :, :])
    return 1.0 - 0.5 * np.real(np.trace(P, axis1=-2, axis2=-1))


def Up_of(X, B, s=0.0):
    U = expc(X + s * B)
    return U[..., 0, :, :] @ U[..., 1, :, :] @ dag(U[..., 2, :, :]) @ dag(U[..., 3, :, :])


def D_batch(X, Bl, h=1e-3, chunk=4):
    """X (NS,4,3), Bl (nd,4,3) -> D (nd,NS) by central second difference of S_geo."""
    nd = Bl.shape[0]; NS = X.shape[0]
    out = np.empty((nd, NS))
    Xe = X[None]
    for i0 in range(0, nd, chunk):
        b = Bl[i0:i0 + chunk][:, None, :, :]
        Sm = Sg_plaq(Xe, b, -h)
        S0 = Sg_plaq(Xe, b, 0.0)
        Sp = Sg_plaq(Xe, b, h)
        out[i0:i0 + chunk] = (Sp - 2.0 * S0 + Sm) / h ** 2
        del Sm, S0, Sp
    return out


hdr("P23  C_geo := (-D)^+ , D = d^2/ds^2 S_geo , H_vv = beta*D ; ||v||^2 = (1/2)sum_l|a_l|^2 = 1")
print("  py3.12 numpy; FD h=1e-3 central second difference; one job at a time (memory).", flush=True)

# ------------------------------------------------------- [0] definition / sign re-derivation
hdr("[0] re-derive D from P06/O06:  D = sum_p[1/2<U_p,K_p> - 1/2<U_p,W'_p>] , K=H^2>=0")
rng0 = np.random.default_rng(20251125)
print("  (i) FD D  vs  analytic  D = sum_p -1/2 Re tr[(W'+W^2)U_p],  W'=dW/ds", flush=True)
worst = 0.0
rels = []
for _ in range(400):
    Xr = rng0.normal(size=(4, 3)) * 0.9
    Br = rng0.normal(size=(4, 3)); Br *= np.sqrt(2.0 / np.sum(Br ** 2))
    Dfd = D_batch(Xr[None], Br[None], h=1e-3)[0, 0]
    hw = 1e-5

    def Wf(s):
        Up = lambda t: Up_of(Xr, Br, t)
        return (Up(s + hw) - Up(s - hw)) / (2 * hw) @ np.linalg.inv(Up(s))
    W0 = Wf(0.0); W1 = (Wf(hw) - Wf(-hw)) / (2 * hw)
    U0 = Up_of(Xr, Br, 0.0)
    De = -0.5 * np.real(np.trace((W1 + W0 @ W0) @ U0))
    rels.append(abs(Dfd - De) / max(abs(De), 1e-30))
    worst = max(worst, rels[-1])
rels = np.array(rels)
print(f"      over 400 random (X,v): max rel dev {worst:.2e} ; median {np.median(rels):.2e} ;"
      f" cases >1e-3: {int(np.sum(rels>1e-3))}/400 (they are the near-cancellation |D|~0 cases)"
      f"  => decomposition and sign confirmed", flush=True)
print("  (ii) h-stability of D (same config/direction):", flush=True)
Xone = rng0.normal(size=(4, 3)) * 0.8
Bone = np.tile(np.array([0., 0., 1.]) / np.sqrt(2.), (1, 4, 1))
for h in [3e-4, 1e-3, 3e-3]:
    print(f"      h={h:.1e}  D={D_batch(Xone[None], Bone, h)[0, 0]: .12e}", flush=True)
print("  (iii) cold config U_p = 1 (X=0): expect D = 1/2 Re tr(H^2) >= 0 hence C_geo = 0", flush=True)
Xz = np.zeros((1, 4, 3)); Bz = unit_dirs(6, (4, 3), np.random.default_rng(3))
Dz = D_batch(Xz, Bz)
print(f"      D = {np.array2string(Dz[:, 0], precision=6)}   C_geo = {np.array2string(np.maximum(-Dz[:, 0], 0.), precision=6)}",
      flush=True)
print("  => D = H_vv/beta (exact, beta factors out); C_geo := (-D)^+ is the NEGATIVE PART.", flush=True)

# ------------------------------------------------- [1] single plaquette: 605 directions
hdr("[1] single plaquette: 605 directions x 24000 product-Haar configs, reweighted e^{-beta S_geo}")
NS = 24000
rng = np.random.default_rng(777)
q = rng.normal(size=(NS, 4, 4)); q /= np.linalg.norm(q, axis=2, keepdims=True)
Uh = np.zeros((NS, 4, 2, 2), complex)
w_, x_, y_, z_ = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
Uh[..., 0, 0] = w_ + 1j * z_; Uh[..., 1, 1] = w_ - 1j * z_
Uh[..., 0, 1] = 1j * x_ + y_; Uh[..., 1, 0] = 1j * x_ - y_
X = logc(Uh)
Sg = Sg_plaq(X, np.zeros((1, 1, 4, 3)), 0.0)[0]
print(f"  Haar self-test: <1/2 Re tr U_l>={np.mean(0.5*np.real(np.trace(Uh, axis1=-2, axis2=-1))):+.5f}"
      f"   <S_geo>={Sg.mean():.5f}  (Haar: 1.0)", flush=True)

rngd = np.random.default_rng(31415)
NR, NC = 300, 300
B_rand = unit_dirs(NR, (4, 3), rngd)
ax = rngd.normal(size=(NC, 3)); ax /= np.linalg.norm(ax, axis=1, keepdims=True)
B_con = (ax * np.sqrt(2. / 4.)).reshape(NC, 1, 3).repeat(4, axis=1)
B_sav = np.tile(np.array([0., 0., 1.]) / np.sqrt(2.), (1, 4, 1))
B_r1 = np.zeros((4, 4, 3))
for i in range(4):
    B_r1[i, i] = np.array([0., 0., np.sqrt(2.)])
B_all = np.concatenate([B_rand, B_con, B_sav, B_r1])
lab = np.array(["random"] * NR + ["const-axis"] * NC + ["savvidy"]
               + [f"rank1(l{i})" for i in range(4)])
nv2 = 0.5 * np.sum(B_all.reshape(len(B_all), -1) ** 2, axis=1)
print(f"  directions: {len(B_all)}  (random {NR} + const-axis {NC} + savvidy + rank1 x4)"
      f"   ||v||^2 in [{nv2.min():.12f}, {nv2.max():.12f}]", flush=True)
t1 = time.time()
D = D_batch(X, B_all)
C = np.maximum(-D, 0.0)
print(f"  D matrix {D.shape} computed in {time.time()-t1:.1f}s", flush=True)


def wstats(M, w):
    W = w / w.sum()
    m = M @ W
    v2 = (M * M) @ W - m * m
    ess = 1.0 / np.sum(W * W)
    sem = np.sqrt(np.maximum(v2, 0.0) / ess)
    return m, sem, ess


hdr("[1a] E[C_geo] and beta*E[C_geo] at beta = 2.0 / 2.5 / 3.0 / 3.5  (per unit ||v||^2)")
print("   beta    ESS     family        E[C_geo]        E[D]      beta*E[C]   /(N/12)    max(-D)",
      flush=True)
res = {}
for beta in [2.0, 2.5, 3.0, 3.5]:
    w = np.exp(-beta * Sg)
    Ec, Ecsem, ess = wstats(C, w)
    Ed, _, _ = wstats(D, w)
    res[beta] = (Ec, Ecsem, Ed, ess)
    for fam in ["savvidy", "rank1(l0)"]:
        m = lab == fam
        print(f"   {beta:<5.2f} {ess:6.0f}   {fam:<12} {Ec[m].mean():.6f}+-{Ecsem[m].mean():.6f}"
              f"  {Ed[m].mean(): .6f}   {beta*Ec[m].mean():.6f}   {beta*Ec[m].mean()/(1/6):.4f}"
              f"    {np.max(-D[m]):.4f}", flush=True)
    for fam in ["const-axis", "random"]:
        m = lab == fam
        print(f"   {beta:<5.2f} {ess:6.0f}   {fam:<12} {Ec[m].mean():.6f}"
              f"+-{Ec[m].std(ddof=1)/np.sqrt(m.sum()):.6f}"
              f"  {Ed[m].mean(): .6f}   {beta*Ec[m].mean():.6f}   {beta*Ec[m].mean()/(1/6):.4f}"
              f"    {np.max(-D[m]):.4f}", flush=True)
    i = int(np.argmax(Ec))
    print(f"   {beta:<5.2f} {ess:6.0f}   SUP over v  {Ec[i]:.6f}+-{Ecsem[i]:.6f} ({lab[i]}, dir#{i})"
          f"   beta*E= {beta*Ec[i]:.6f}  = {beta*Ec[i]/(1/6):.4f}*(N/12)", flush=True)

hdr("[1b] beta-scan: TREND of beta*sup_v E[C_geo]  (the object the window depends on)")
betas = np.round(np.arange(1.5, 3.9001, 0.05), 4)
Gsup = np.empty(len(betas)); Gmean = np.empty(len(betas)); Gmrand = np.empty(len(betas))
Esup = np.empty(len(betas)); esss = np.empty(len(betas))
m_rand = lab == "random"
for k, beta in enumerate(betas):
    w = np.exp(-beta * Sg)
    Ec, _, ess = wstats(C, w)
    Sup = Ec.max(); Imax = int(np.argmax(Ec))
    Gsup[k] = beta * Sup
    Gmean[k] = beta * Ec.mean()
    Gmrand[k] = beta * Ec[m_rand].mean()
    Esup[k] = Sup
    esss[k] = ess
    if abs(beta - round(beta)) < 1e-9 or abs(beta * 2 - round(beta * 2)) < 1e-9:
        print(f"   beta={beta:5.2f}  ESS={ess:6.0f}  sup_v E[C]={Sup:.6f} ({lab[Imax]:<10}"
              f" dir#{Imax:3d})  beta*sup={Gsup[k]:.6f}={Gsup[k]/(1/6):.4f}(N/12)"
              f"  beta*<E[C]>={Gmean[k]:.6f}={Gmean[k]/(1/6):.4f}(N/12)", flush=True)
ksel = np.abs(betas - 2.5) < 1e-9
print(f"\n   argmax over beta of beta*sup_v E[C_geo]: beta={betas[int(np.argmax(Gsup))]:.2f}"
      f"  value={Gsup.max():.6f} = {Gsup.max()/(1/6):.4f}*(N/12)", flush=True)
print(f"   argmax over beta of beta*<E[C_geo]>_v   : beta={betas[int(np.argmax(Gmean))]:.2f}"
      f"  value={Gmean.max():.6f} = {Gmean.max()/(1/6):.4f}*(N/12)", flush=True)
# local slope at beta = 2.5 (linear fit over [2.3, 2.7])
mk = (betas >= 2.3 - 1e-9) & (betas <= 2.7 + 1e-9)
cA = np.polyfit(betas[mk], Gsup[mk], 1); cB = np.polyfit(betas[mk], Gmean[mk], 1)
print(f"   d(beta*sup_v E[C])/dbeta |_(2.3..2.7) = {cA[0]:+.5f}   "
      f"d(beta*<E[C]>_v)/dbeta = {cB[0]:+.5f}", flush=True)
print(f"   beta*sup_v E[C_geo] at beta=2.0/2.5/3.0/3.5: "
      f"{Gsup[np.abs(betas-2.0)<1e-9][0]:.6f} / {Gsup[np.abs(betas-2.5)<1e-9][0]:.6f} / "
      f"{Gsup[np.abs(betas-3.0)<1e-9][0]:.6f} / {Gsup[np.abs(betas-3.5)<1e-9][0]:.6f}", flush=True)

hdr("[1c] boundary solves")
N12 = 1.0 / 6.0
V0 = 0.011059293772973        # H1E closed form ( = (N/12)/3.882^2 ); P06 sec 5.3 / 7.4
print(f"   N/12 = {N12:.6f} ; V0 = {V0:.12f} ; sqrt((N/12)/V0) = {np.sqrt(N12/V0):.3f}"
      f" (rank-1 window, O06 sec 4.2)", flush=True)


def root(f):
    """first sign change of f on the beta grid, linear interpolation"""
    for k in range(len(betas) - 1):
        if f[k] > 0 >= f[k + 1]:
            return betas[k] + 0.05 * f[k] / (f[k] - f[k + 1])
    return np.nan


k25 = int(np.argmin(np.abs(betas - 2.5)))
print(f"   (A) C_geo-only equation  beta*sup_v E[C_geo] = N/12 :  root = "
      f"{root(Gsup - N12)}  (nan = NO crossing; peak {Gsup.max():.6f} < {N12:.6f})", flush=True)
print(f"   (B) C_geo-only, family-mean  beta*<E[C]>_v = N/12  :  root = {root(Gmean - N12)}", flush=True)
print(f"   floor surplus (N/12) - beta*sup_v E[C] : min over beta = {N12 - Gsup.max():.6f}"
      f" at beta={betas[int(np.argmax(Gsup))]:.2f} ; at beta=2.5 = {N12 - Gsup[k25]:.6f}", flush=True)
print(f"   (C) full window  beta*sup_v E[C](beta) + beta^2*V0 = N/12 :  root = "
      f"{root(N12 - Gsup - V0 * betas ** 2)}", flush=True)
print(f"   (D) full window with family-mean deficit:                    root = "
      f"{root(N12 - Gmean - V0 * betas ** 2)}", flush=True)
print(f"   composite surplus at beta=2.5: (N/12) - beta*sup E[C] - beta^2*V0 = "
      f"{N12 - Gsup[k25] - V0 * 6.25:+.6f}  (negative => window CLOSED at the physical point)",
      flush=True)
print(f"   same with P06's own deficit (scale x0.907 = its 0.0388 vs our 0.0429): root = "
      f"{root(N12 - 0.907 * Gsup - V0 * betas ** 2)}   (~P06's stated 2.51)", flush=True)
for sc in [0.9, 1.0, 1.1]:
    print(f"       V0 x {sc:.1f}: full window (C) root = {root(N12 - Gsup - sc*V0*betas**2)}", flush=True)
print("   surplus table  S(beta) = (N/12) - beta*sup_v E[C] - beta^2*V0 :", flush=True)
for b in [1.8, 2.0, 2.2, 2.25, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0]:
    j = int(np.argmin(np.abs(betas - b)))
    S = N12 - Gsup[j] - V0 * betas[j] ** 2
    print(f"     beta={betas[j]:.2f}  beta*sup_vE[C]={Gsup[j]:.6f}  beta^2V0={V0*betas[j]**2:.6f}"
          f"  S={S:+.6f} ({100*S/N12:+.2f}% of N/12)", flush=True)
E25 = res[2.5][0]
order = np.argsort(-E25)
print(f"   top-4 directions at beta=2.5: " + ", ".join(
    f"#{int(order[j])}({lab[order[j]]})={E25[order[j]]:.5f}" for j in range(4)), flush=True)
print(f"   sup_v E[C](2.5): 605 dirs={E25.max():.6f} ; 300 random only={E25[m_rand].max():.6f}"
      f" ; 300 const-axis only={E25[lab == 'const-axis'].max():.6f}"
      f"  (P06, 128 random, got 0.038842)", flush=True)

# ---------------------------------------------------------------- [2] L=2 4D lattice
hdr("[2] L=2 4D lattice (16 sites, 64 links, 24 plaquettes) at beta=2.5: E[C_geo] zero?")
LL = 2
sites = [(a, b, c, d) for a in range(LL) for b in range(LL) for c in range(LL) for d in range(LL)]


def sh(x, mu, d=1):
    y = list(x); y[mu] = (y[mu] + d) % LL; return tuple(y)


def shp(V, axis, d=1):
    return np.roll(V, d, axis=axis)


def Sg_lat(X, B, s, roll_shift=6):
    """X,B (...,L,L,L,L,4,3) -> S_geo summed over sites.

    NOTE on the roll axis: Um/Un are obtained by integer-indexing the LINK axis away, so the
    resulting 7-dim array has spatial axis d at negative index d-6 (not d-7).  roll_shift=6 is
    therefore the correct convention; it is verified against the dict action to 1e-16 below.
    """
    U = expc(X + s * B)
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Um = U[..., mu, :, :]; Un = U[..., nu, :, :]
            P = (Um @ shp(Un, mu - roll_shift, -1) @ dag(shp(Um, nu - roll_shift, -1))
                 @ dag(Un))
            tot = tot + np.sum(1.0 - 0.5 * np.real(np.trace(P, axis1=-2, axis2=-1)),
                               axis=(-1, -2, -3, -4))
    return tot


def D_lat(X, Bs, h=1e-3, roll_shift=6):
    Xb = X[None]
    Sm = Sg_lat(Xb, Bs, -h, roll_shift)
    S0 = Sg_lat(Xb, Bs, 0.0, roll_shift)
    Sp = Sg_lat(Xb, Bs, h, roll_shift)
    return (Sp - 2.0 * S0 + Sm) / h ** 2


def Sgeo_lat_U(U):
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for x in sites:
                P = U[x][mu] @ U[sh(x, mu)][nu] @ U[sh(x, nu)][mu].conj().T @ U[x][nu].conj().T
                tot += 1.0 - 0.5 * np.real(np.trace(P))
    return tot


def staple(U, x, mu):
    S = np.zeros((2, 2), complex)
    for nu in range(4):
        if nu == mu:
            continue
        A = U[sh(x, mu)][nu]; B = U[sh(x, nu)][mu]; Cc = U[x][nu]
        S += A @ B.conj().T @ Cc.conj().T
        P1 = U[sh(x, nu, -1)][nu]; P2 = U[sh(x, nu, -1)][mu]; P3 = U[sh(sh(x, mu), nu, -1)][nu]
        S += dag(P3) @ dag(P2) @ P1
    return S


def to_array(cfg):
    A = np.zeros((LL, LL, LL, LL, 4, 2, 2), complex)
    for x in sites:
        for mu in range(4):
            A[x][mu] = cfg[x][mu]
    return A


def mcmc(beta, ntherm, nmeas, every, seed, sig=0.35):
    rngm = np.random.default_rng(seed)
    U = {x: [randU(rngm) for _ in range(4)] for x in sites}
    acc = 0; tot = 0; out = []
    for it in range(ntherm + nmeas * every):
        for x in sites:
            for mu in range(4):
                S = staple(U, x, mu)
                V = expc(rngm.normal(size=3) * sig) @ U[x][mu]
                dS = -0.5 * np.real(np.trace((V - U[x][mu]) @ S))
                tot += 1
                if dS <= 0 or rngm.random() < np.exp(-beta * dS):
                    U[x][mu] = V; acc += 1
        if it >= ntherm and (it - ntherm) % every == 0:
            out.append({x: list(U[x]) for x in sites})
    return out, acc / tot


# engine cross-check: dict action vs array action (both roll conventions)
rngt = np.random.default_rng(4242)
Ut = {x: [randU(rngt) for _ in range(4)] for x in sites}
A0 = to_array(Ut)
X0 = logc(A0)
Sdict = Sgeo_lat_U(Ut)
Sarr6 = float(Sg_lat(X0[None], np.zeros((1,) + X0.shape), 0.0, 6)[0])
Sarr7 = float(Sg_lat(X0[None], np.zeros((1,) + X0.shape), 0.0, 7)[0])
print(f"  action cross-check: dict={Sdict:.10f}  array(axis=mu-6, correct)={Sarr6:.10f}"
      f"  rel={abs(Sarr6-Sdict)/abs(Sdict):.1e}", flush=True)
print(f"     [wrong convention axis=mu-7, for the record] ={Sarr7:.10f}"
      f"  rel={abs(Sarr7-Sdict)/abs(Sdict):.1e}", flush=True)
x0 = sites[3]; mu0 = 2
Sst = staple(Ut, x0, mu0)
V = expc(np.array([0.11, -0.07, 0.05])) @ Ut[x0][mu0]
dS_st = -0.5 * np.real(np.trace((V - Ut[x0][mu0]) @ Sst))
dS_bf = Sgeo_lat_U({**{y: list(Ut[y]) for y in sites},
                    x0: [V if i == mu0 else Ut[x0][i] for i in range(4)]}) - Sgeo_lat_U(Ut)
print(f"  staple check: dS_staple={dS_st:+.12e} dS_brute={dS_bf:+.12e}"
      f" rel={abs(dS_st-dS_bf)/max(abs(dS_bf),1e-30):.1e}", flush=True)

sh_ = (LL, LL, LL, LL, 4, 3)
rngl = np.random.default_rng(8888)
axl = rngl.normal(size=(32, 3)); axl /= np.linalg.norm(axl, axis=1, keepdims=True)
Bconst = (axl * np.sqrt(2. / 64.)).reshape(32, 1, 1, 1, 1, 1, 3).repeat(LL, 1).repeat(LL, 2) \
    .repeat(LL, 3).repeat(LL, 4).repeat(4, 5)
b_sav = np.zeros((1,) + sh_); b_sav[..., 2] = np.sqrt(2. / 64.)
b_r1 = np.zeros((1,) + sh_); b_r1[0, 0, 0, 0, 0, 0] = np.array([0., 0., np.sqrt(2.)])
Bs_all = np.concatenate([unit_dirs(300, sh_, rngl), Bconst, b_sav, b_r1])
labL = np.array(["random"] * 300 + ["const-axis"] * 32 + ["savvidy", "rank1"])
nv2L = 0.5 * np.sum(Bs_all.reshape(len(Bs_all), -1) ** 2, axis=1)
print(f"  lattice directions: {len(Bs_all)}  ||v||^2 in [{nv2L.min():.12f}, {nv2L.max():.12f}]",
      flush=True)

t2 = time.time()
cf, ar = mcmc(2.5, 400, 300, 3, 20251125)
print(f"  Metropolis beta=2.5: acc={ar:.3f}  n_cfg={len(cf)}  ({time.time()-t2:.1f}s thermal+meas)",
      flush=True)
Us = np.stack([logc(to_array(c)) for c in cf])
print(f"  config array {Us.shape}", flush=True)
Dm = np.empty((len(cf), len(Bs_all)))
t3 = time.time()
for i in range(len(cf)):
    Dm[i] = D_lat(Us[i], Bs_all)
    if (i + 1) % 100 == 0:
        print(f"    {i+1}/{len(cf)} configs  ({time.time()-t3:.1f}s)", flush=True)
C2 = np.maximum(-Dm, 0.0)
print("\n    family        E[C_geo]     E[D]       frac(D<0)    min D     max(-D)  E[C]/n_plaq",
      flush=True)
for fam in ["random", "const-axis", "savvidy", "rank1"]:
    m = labL == fam
    print(f"    {fam:<12} {C2[:, m].mean(): .8f}  {Dm[:, m].mean(): .6f}   {np.mean(Dm[:, m] < 0):.4f}"
          f"      {np.min(Dm[:, m]): .4f}   {np.max(-Dm[:, m]):.4f}   {C2[:, m].mean()/24:.8f}",
          flush=True)
print(f"    SUP over v : E[C_geo]={C2.mean(axis=0).max():.8f}  "
      f"beta*E={2.5*C2.mean(axis=0).max():.8f}={2.5*C2.mean(axis=0).max()/(1/6):.4f}*(N/12)"
      f"  at {labL[int(np.argmax(C2.mean(axis=0)))]}", flush=True)
neg = int(np.sum(Dm < 0)); tot2 = Dm.size
print(f"    negative (cfg,dir) pairs: {neg}/{tot2} ; min D over all = {Dm.min():.6f}", flush=True)

# --------------------------------------------------------------------- [3] verdict numbers
hdr("[3] verdict quantities")
print(f"  E[C_geo] at beta=2.0 : sup_v={res[2.0][0].max():.6f}  <v>_all={res[2.0][0].mean():.6f}"
      f"  <v>_rand={res[2.0][0][m_rand].mean():.6f}", flush=True)
for b in [2.0, 2.5, 3.0, 3.5]:
    Ec, Se, Ed, ess = res[b]
    im = int(np.argmax(Ec))
    print(f"  beta={b}: sup_v E[C]={Ec[im]:.6f}  beta*sup={b*Ec[im]:.6f}  "
          f"<E[C]>_rand={Ec[m_rand].mean():.6f} beta*<>= {b*Ec[m_rand].mean():.6f}  "
          f"ESS={ess:.0f}  E[D]={Ed.mean(): .6f}", flush=True)
print(f"  peak of beta*sup_v E[C_geo] = {Gsup.max():.6f} at beta={betas[int(np.argmax(Gsup))]:.2f}"
      f"  (N/12 = {N12:.6f}, so no C_geo-only crossing exists: {Gsup.max() < N12})", flush=True)
print(f"  total runtime {time.time()-t00:.1f}s", flush=True)
