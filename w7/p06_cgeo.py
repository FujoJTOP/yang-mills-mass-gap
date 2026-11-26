"""p06_cgeo.py -- P06 (o06 gap): the geometric negative part C_geo of H_vv/beta.

Definition fixed HERE (exact negative part; O06 Sec 4.4 only had a triangle majorant):
    S_geo(U) := sum_p (1 - (1/2) Re tr U_p)              [ = S_W / beta ; beta-free ]
    D(v)     := d^2/ds^2 S_geo(U(s))|_{s=0},  U_l(s) = exp(X_l + s v_l)   [o06 chart]
    H_vv      = beta * D(v)                              [exact, beta clean]
    C_geo(v) := (-D(v))^+ = max(0, -D(v))                [exact negative part]
  =>  H_vv >= -beta * C_geo(v)   pointwise, equality iff D(v) <= 0.

Normalization: ||v||^2 := sum_l tr(v_l^dag v_l); every scanned v has ||v||^2 = 1, so the
Phi'' floor is (N/12)||v||^2 = 1/6 (N=2) and all numbers below are per unit ||v||^2.

Engines: (1) one plaquette (4 links); (2) 4D lattice L=2 (16 sites, 64 links, 24 plaq).
Measures: product Haar ; Haar x exp(-beta S_geo)  (single plaq: exact self-normalized IS,
lattice: sequential Metropolis).  D(v) is beta-free, so one Haar sample serves all beta.
Run: py -3.12 -u w7/p06_cgeo.py
"""
import sys, time
import numpy as np

t00 = time.time()
def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78, flush=True)


def expc(a):
    """a (...,3) -> exp(i a.sigma/2)  (...,2,2).  a = chart vector, |a|<2pi."""
    a = np.asarray(a, float)
    th = np.sqrt(np.sum(a * a, axis=-1))
    small = th < 1e-12
    c = np.cos(0.5 * th)
    ratio = np.where(small, 0.5, np.sin(0.5 * th) / np.where(small, 1.0, th))
    m = np.zeros(a.shape[:-1] + (2, 2), complex)
    ax, ay, az = a[..., 0], a[..., 1], a[..., 2]
    m[..., 0, 0] = c + 1j * ratio * az
    m[..., 1, 1] = c - 1j * ratio * az
    m[..., 0, 1] = ratio * (1j * ax + ay)
    m[..., 1, 0] = ratio * (1j * ax - ay)
    return m


def dag(M):
    return np.swapaxes(M.conj(), -1, -2)


def logc(U):
    """U (...,2,2) in SU(2) -> chart a (...,3), U = expc(a), |a| in [0,2pi)."""
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
    """n random directions b with sum_l |b_l|^2 = 2  (i.e. ||v||^2 = 1)."""
    b = rng.normal(size=(n,) + shape)
    f = np.sqrt(0.5 * np.sum(b.reshape(n, -1) ** 2, axis=1))
    return b / f.reshape((n,) + (1,) * len(shape))


# ---------------- engine 1: single plaquette (4 links, chart axes (...,4,3)) -------------
def Sg_plaq(X, B, s):
    U = expc(X + s * B)                       # (...,4,2,2)
    P = U[..., 0, :, :] @ U[..., 1, :, :] @ dag(U[..., 2, :, :]) @ dag(U[..., 3, :, :])
    return 1.0 - 0.5 * np.real(np.trace(P, axis1=-2, axis2=-1))


def D_plaq(X, Bs, h=1e-3):
    Xb = X[None]
    Sm = Sg_plaq(Xb, Bs, -h); S0 = Sg_plaq(Xb, Bs, 0.0); Sp = Sg_plaq(Xb, Bs, h)
    return (Sp - 2.0 * S0 + Sm) / h ** 2


# ---------------- engine 2: 4D lattice (links indexed as (x,y,z,t,mu)) -------------------
def Sg_lat(X, B, s):
    """X,B (...,L,L,L,L,4,3) -> S_geo (...,) summed over sites only."""
    U = expc(X + s * B)                       # (...,L,L,L,L,4,2,2)
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Um = U[..., mu, :, :]; Un = U[..., nu, :, :]
            P = Um @ np.roll(Un, -1, axis=mu - 6) @ dag(np.roll(Um, -1, axis=nu - 6)) @ dag(Un)
            tot = tot + np.sum(1.0 - 0.5 * np.real(np.trace(P, axis1=-2, axis2=-1)),
                               axis=(-1, -2, -3, -4))
    return tot


def D_lat(X, Bs, h=1e-3):
    Xb = X[None]
    Sm = Sg_lat(Xb, Bs, -h); S0 = Sg_lat(Xb, Bs, 0.0); Sp = Sg_lat(Xb, Bs, h)
    return (Sp - 2.0 * S0 + Sm) / h ** 2


L = 2
sites = [(a, b, c, d) for a in range(L) for b in range(L) for c in range(L) for d in range(L)]


def sh(x, mu, d=1):
    y = list(x); y[mu] = (y[mu] + d) % L; return tuple(y)


def staple(U, x, mu):
    """staple S s.t.  dS_geo = -1/2 Re tr(dU . S)  for the link U_mu(x)."""
    S = np.zeros((2, 2), complex)
    for nu in range(4):
        if nu == mu:
            continue
        # plaquette P_mu_nu(x)      : U_mu(x) at slot 1
        A = U[sh(x, mu)][nu]; B = U[sh(x, nu)][mu]; C = U[x][nu]
        S += A @ B.conj().T @ C.conj().T
        # plaquette P_mu_nu(x-nu)   : U_mu(x)^dag at slot 3
        P1 = U[sh(x, nu, -1)][nu]; P2 = U[sh(x, nu, -1)][mu]; P3 = U[sh(sh(x, mu), nu, -1)][nu]
        S += dag(P3) @ dag(P2) @ P1
    return S


def Sgeo_lat_U(U):
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for x in sites:
                P = U[x][mu] @ U[sh(x, mu)][nu] @ U[sh(x, nu)][mu].conj().T @ U[x][nu].conj().T
                tot += 1.0 - 0.5 * np.real(np.trace(P))
    return tot


def to_array(cfg):
    A = np.zeros((L, L, L, L, 4, 2, 2), complex)
    for x in sites:
        for mu in range(4):
            A[x][mu] = cfg[x][mu]
    return A


def mcmc(beta, ntherm, nmeas, every, seed, sig=0.35):
    rng = np.random.default_rng(seed)
    U = {x: [randU(rng) for _ in range(4)] for x in sites}
    acc = 0; tot = 0; out = []
    for it in range(ntherm + nmeas * every):
        for x in sites:
            for mu in range(4):
                S = staple(U, x, mu)
                V = expc(rng.normal(size=3) * sig) @ U[x][mu]
                dS = -0.5 * np.real(np.trace((V - U[x][mu]) @ S))
                tot += 1
                if dS <= 0 or rng.random() < np.exp(-beta * dS):
                    U[x][mu] = V; acc += 1
        if it >= ntherm and (it - ntherm) % every == 0:
            out.append({x: list(U[x]) for x in sites})
    return out, acc / tot


hdr("P06  C_geo := (-D)^+ ,  D = d^2/ds^2 sum_p (1 - 1/2 Re tr U_p),  H_vv = beta*D")
print("  py3.12 numpy, FD step h=1e-3 (second difference), ||v||^2 = sum_l tr(v^dag v) = 1",
      flush=True)

# ------------------------------------------------------------------ [0] sign / definition
hdr("[0] sign check of  H_vv = beta[ sum_p 1/2<U_p,K_p> - sum_p 1/2<U_p,W'_p> ],  K_p=H_p^2>=0")
rng = np.random.default_rng(1)
Xp = rng.normal(size=(4, 3)) * 0.7
beta0 = 2.5
for tag, b3 in [("b=e3", np.array([0., 0., 1.])), ("b=e1", np.array([1., 0., 0.])),
                ("b=random", rng.normal(size=3))]:
    B = np.tile(b3, (4, 1))[None]
    h = 1e-4
    Hvv = beta0 * (Sg_plaq(Xp[None], B, h)[0] - 2 * Sg_plaq(Xp[None], B, 0.0)[0]
                   + Sg_plaq(Xp[None], B, -h)[0]) / h ** 2
    bvec = B[0]

    def Up_at(s):
        U = expc(Xp + s * bvec)
        return U[0] @ U[1] @ dag(U[2]) @ dag(U[3])

    hw = 1e-5
    Wp = lambda s: (Up_at(s + hw) - Up_at(s - hw)) / (2 * hw) @ np.linalg.inv(Up_at(s))
    W0 = Wp(0.0); W1 = (Wp(hw) - Wp(-hw)) / (2 * hw)
    Hp = -1j * W0; Hp = 0.5 * (Hp + Hp.conj().T); K = Hp @ Hp
    U0 = Up_at(0.0)
    tK = np.real(np.trace(K))
    dec = beta0 * (0.5 * np.real(np.trace(U0 @ K)) - 0.5 * np.real(np.trace(U0 @ W1)))
    ev = np.linalg.eigvalsh(K)
    print(f"  {tag:<10} FD H_vv={Hvv: .10e}  beta[1/2<U,K>-1/2<U,W'>]={dec: .10e}"
          f"  rel={abs(Hvv-dec)/max(abs(Hvv),1e-30):.1e}  eig(K)={ev[0]: .3e},{ev[1]: .3e}"
          f"  <U,K>/trK={np.real(np.trace(U0@K))/tK: .4f}", flush=True)
print("  => decomposition confirmed (rel<=1e-5); K_p>=0 always; only U_p carries the sign.",
      flush=True)

# ------------------------------------------------------- [1] single plaquette, Haar sample
hdr("[1] single plaquette: D(v) per unit ||v||^2 on product-Haar 4-link samples")
NS = 20000
q = rng.normal(size=(NS, 4, 4)); q /= np.linalg.norm(q, axis=2, keepdims=True)
Uhaar = np.zeros((NS, 4, 2, 2), complex)
w_, x_, y_, z_ = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
Uhaar[..., 0, 0] = w_ + 1j * z_; Uhaar[..., 1, 1] = w_ - 1j * z_
Uhaar[..., 0, 1] = 1j * x_ + y_; Uhaar[..., 1, 0] = 1j * x_ - y_
Xhaar = logc(Uhaar)
print(f"  Haar check: <1/2 Re tr U_l>={np.mean(0.5*np.real(np.trace(Uhaar,axis1=-2,axis2=-1))):+.5f}"
      f"  <|chart a|>={np.mean(np.linalg.norm(Xhaar,axis=-1)):.4f} (Haar in this chart peaks at pi)",
      flush=True)
Sgeo = Sg_plaq(Xhaar, np.zeros((NS, 4, 3)), 0.0)
print(f"  <S_geo>_Haar = {Sgeo.mean():.5f}", flush=True)

rngd = np.random.default_rng(7)
dirs = [("savvidy(all links e3)", np.tile(np.array([0., 0., 1.]) / np.sqrt(2.), (1, 4, 1)))]
for i in range(4):
    B = np.zeros((1, 4, 3)); B[0, i] = np.array([0., 0., np.sqrt(2.)])
    dirs.append((f"rank1(link{i} e3)", B))
ax = rngd.normal(size=(128, 3)); ax /= np.linalg.norm(ax, axis=1, keepdims=True)
dirs.append(("const-axis 128 (all links same)", (ax * np.sqrt(2. / 4.)).reshape(128, 1, 3).repeat(4, axis=1)))
dirs.append(("random 128 (per-link)", unit_dirs(128, (4, 3), rngd)))

names, Dmats = [], []
for name, B in dirs:
    nd = B.shape[0]
    D = np.empty((nd, NS))
    for i in range(nd):
        D[i] = D_plaq(Xhaar, B[i:i + 1])[0]
    names.append((name, nd)); Dmats.append(D)
Dpl = np.vstack(Dmats)
labels = np.concatenate([[n] * k for n, k in names])
print(f"  directions scanned: {Dpl.shape[0]}", flush=True)

def wstats(D, w):
    W = w / w.sum()
    return (np.maximum(-D, 0.0) @ W), (D @ W)

for beta in [2.5]:
    w = np.exp(-beta * Sgeo)
    Ec, Ed = wstats(Dpl, w)
    ess = w.sum() ** 2 / np.sum(w ** 2)
    print(f"\n  --- beta={beta}: measure exp(-beta S_geo) x Haar, ESS={ess:.0f}/{NS} ---", flush=True)
    print("    family                             E[C_geo]     E[D]=E[H_vv]/beta   max(-D)   frac(D<0)",
          flush=True)
    for n, k in names:
        m = labels == n
        print(f"    {n:<33} {Ec[m].mean(): .6f}  {Ed[m].mean(): .6f}   {np.max(-Dpl[m]):.4f}"
              f"   {np.mean(Dpl[m] < 0):.3f}", flush=True)
    imax = int(np.argmax(Ec)); mrand = labels == "random 128 (per-link)"
    print(f"    sup_v E[C_geo] = {Ec.max():.6f} at {labels[imax]!r} (dir #{imax})")
    print(f"    random-128 mean={Ec[mrand].mean():.6f} +- {Ec[mrand].std(ddof=1)/np.sqrt(128):.6f}"
          f" ;  const-axis mean={Ec[labels=='const-axis 128 (all links same)'].mean():.6f}"
          f" +- {Ec[labels=='const-axis 128 (all links same)'].std(ddof=1)/np.sqrt(128):.6f}")
    print(f"    N/12=1/6={1/6:.6f}:  E[C_geo]/(N/12)={Ec.max()/(1/6):.4f}"
          f" ; deficit beta*E[C_geo]={beta*Ec.max():.6f} = {beta*Ec.max()/(1/6):.4f}*(N/12)")
    print(f"    random-128: beta*E[C_geo]={beta*Ec[mrand].mean():.6f} ="
          f" {beta*Ec[mrand].mean()/(1/6):.4f}*(N/12)", flush=True)

# ------------------------------------------------------------- [2] L=2 lattice, 4D
hdr("[2] L=2 4D lattice (16 sites, 64 links, 24 plaquettes)")
rngt = np.random.default_rng(99)
Ut = {x: [randU(rngt) for _ in range(4)] for x in sites}
x0 = sites[3]; mu0 = 2
Sst = staple(Ut, x0, mu0)
V = expc(np.array([0.11, -0.07, 0.05])) @ Ut[x0][mu0]
dS_st = -0.5 * np.real(np.trace((V - Ut[x0][mu0]) @ Sst))
dS_bf = Sgeo_lat_U({**{y: list(Ut[y]) for y in sites},
                   x0: [V if i == mu0 else Ut[x0][i] for i in range(4)]}) - Sgeo_lat_U(Ut)
A0 = to_array(Ut)
print(f"  staple check : dS_staple={dS_st:+.12e}  dS_brute={dS_bf:+.12e}"
      f"  rel={abs(dS_st-dS_bf)/max(abs(dS_bf),1e-30):.1e}", flush=True)
print(f"  array check  : Sg_lat(arr)={Sg_lat(logc(A0)[None], np.zeros((1,)+logc(A0).shape), 0.0)[0]:.10f}"
      f"  Sgeo_lat_U(dict)={Sgeo_lat_U(Ut):.10f}", flush=True)

sh_ = (L, L, L, L, 4, 3)
rngl = np.random.default_rng(17)
axl = rngl.normal(size=(128, 3)); axl /= np.linalg.norm(axl, axis=1, keepdims=True)
Bconst = (axl * np.sqrt(2. / 64.)).reshape(128, 1, 1, 1, 1, 1, 3).repeat(L, 1).repeat(L, 2) \
    .repeat(L, 3).repeat(L, 4).repeat(4, 5)
b_sav = np.zeros((1,) + sh_); b_sav[..., 2] = np.sqrt(2. / 64.)
b_r1 = np.zeros((1,) + sh_); b_r1[0, 0, 0, 0, 0, 0] = np.array([0., 0., np.sqrt(2.)])
Bs_all = np.concatenate([unit_dirs(128, sh_, rngl), Bconst, b_sav, b_r1])
lab = np.array(["random 128"] * 128 + ["const-axis 128"] * 128
               + ["savvidy(const e3)", "rank1(one link e3)"])

def scan_lat(Us, tag):
    Xs_l = np.stack([logc(U) for U in Us])
    Dm = np.empty((len(Us), Bs_all.shape[0]))
    for i, X in enumerate(Xs_l):
        Dm[i] = D_lat(X, Bs_all)
    C = np.maximum(-Dm, 0.0)
    print(f"  --- {tag}: n_cfg={len(Us)} ---", flush=True)
    print("    family                 E[C_geo]    E[D]      frac(D<0)  max(-D)   E[C_geo]/n_plaq",
          flush=True)
    for n in ["random 128", "const-axis 128", "savvidy(const e3)", "rank1(one link e3)"]:
        m = lab == n
        print(f"    {n:<22} {C[:, m].mean(): .6f}  {Dm[:, m].mean(): .6f}"
              f"   {np.mean(Dm[:, m] < 0):.3f}    {np.max(-Dm[:, m]):.3f}    {C[:, m].mean()/24:.6f}",
              flush=True)
    Ec = C.mean(axis=0)
    print(f"    sup_v E[C_geo]={Ec.max():.6f} at {lab[int(np.argmax(Ec))]!r}"
          f" ; random mean={Ec[:128].mean():.6f} +- {Ec[:128].std(ddof=1)/np.sqrt(128):.6f}"
          f" ; const-axis mean={Ec[128:256].mean():.6f}"
          f" +- {Ec[128:256].std(ddof=1)/np.sqrt(128):.6f}", flush=True)
    print(f"    E[C_geo]/(N/12)={Ec.max()/(1/6):.4f} ; beta*E[C_geo] at beta=2.5: "
          f"{2.5*Ec.max():.6f} = {2.5*Ec.max()/(1/6):.4f}*(N/12)", flush=True)
    return Dm, C, Ec

rg = np.random.default_rng(11)
Us_h = [to_array({x: [randU(rg) for _ in range(4)] for x in sites}) for _ in range(200)]
scan_lat(Us_h, "product Haar (beta=0)")

t0 = time.time()
cf, ar = mcmc(2.5, 500, 300, 4, 4242)
Us_m = [to_array(c) for c in cf]
plq = np.mean([np.mean([0.5 * np.real(np.trace(U[x][mu] @ U[sh(x, mu)][nu]
                                        @ U[sh(x, nu)][mu].conj().T @ U[x][nu].conj().T))
                        for x in sites for mu in range(4) for nu in range(mu + 1, 4)]) for U in Us_m])
print(f"  MCMC beta=2.5: acc={ar:.3f}  n_cfg={len(Us_m)}  <1/2 Re tr U_p>={plq:.4f}"
      f"  ({time.time()-t0:.1f}s)", flush=True)
Dlat, Clat, Elat = scan_lat(Us_m, "Metropolis beta=2.5 (physical)")

hdr("[2b] L=2 lattice: beta scan (reduced direction set: savvidy, const-axis-32, rank1)")
rngb2 = np.random.default_rng(23)
ax2 = rngb2.normal(size=(32, 3)); ax2 /= np.linalg.norm(ax2, axis=1, keepdims=True)
B2 = np.concatenate([(ax2 * np.sqrt(2. / 64.)).reshape(32, 1, 1, 1, 1, 1, 3).repeat(L, 1)
                     .repeat(L, 2).repeat(L, 3).repeat(L, 4).repeat(4, 5), b_sav, b_r1])
lab2 = np.array(["const-axis"] * 32 + ["savvidy", "rank1"])
print("    beta   <plaq>   E[C](sav)  E[C](const)  E[D](sav)  E[D](const)  beta*E[C](const)  /(N/12)",
      flush=True)
for beta in [1.0, 2.5, 5.0]:
    if beta == 2.5:
        U2 = Us_m
    else:
        cf2, ar2 = mcmc(beta, 400, 200, 4, 9091)
        U2 = [to_array(c) for c in cf2]
    X2 = np.stack([logc(U) for U in U2])
    D2 = np.empty((len(U2), B2.shape[0]))
    for i, X in enumerate(X2):
        D2[i] = D_lat(X, B2)
    C2 = np.maximum(-D2, 0.0)
    pl = np.mean([np.mean([0.5 * np.real(np.trace(U[x][mu] @ U[sh(x, mu)][nu]
                                          @ U[sh(x, nu)][mu].conj().T @ U[x][nu].conj().T))
                            for x in sites for mu in range(4) for nu in range(mu + 1, 4)])
                   for U in U2])
    Ecs = C2[:, lab2 == "savvidy"].mean(); Ecc = C2[:, lab2 == "const-axis"].mean()
    Eds = D2[:, lab2 == "savvidy"].mean(); Edc = D2[:, lab2 == "const-axis"].mean()
    print(f"   {beta:<5.2f}  {pl:.4f}   {Ecs: .6f}    {Ecc: .6f}     {Eds: .6f}"
          f"   {Edc: .6f}      {beta*max(Ecs,Ecc):.6f}      {beta*max(Ecs,Ecc)/(1/6):.4f}", flush=True)

# ------------------------------------------------------------- [3] beta scaling
hdr("[3] beta scaling of E[C_geo] (single plaquette: exact reweighting of one Haar sample)")
mrand = labels == "random 128 (per-link)"; mcst = labels == "const-axis 128 (all links same)"
msav = labels == "savvidy(all links e3)"
print("    beta   ESS     E[C](savvidy) E[C](const-ax) E[C](random)  beta*E[C] const-ax   /(N/12)",
      flush=True)
for beta in [0.0, 0.5, 1.0, 2.5, 5.0, 10.0]:
    w = np.ones(NS) if beta == 0 else np.exp(-beta * Sgeo)
    Ec, Ed = wstats(Dpl, w)
    ess = w.sum() ** 2 / np.sum(w ** 2)
    print(f"   {beta:<5.2f} {ess:6.0f}   {Ec[msav].mean(): .6f}     {Ec[mcst].mean(): .6f}"
          f"    {Ec[mrand].mean(): .6f}     {beta*Ec[mcst].mean(): .6f}      "
          f"{beta*Ec[mcst].mean()/(1/6): .4f}", flush=True)

# ------------------------------------------------------------- [4] bound bookkeeping
hdr("[4] bound bookkeeping: C_geo <= C_geo^tri (O06 4.4 majorant), FD stability, analytic K")
Xone = np.random.default_rng(5).normal(size=(4, 3)) * 0.8
Bone = np.tile(np.array([0., 0., 1.]) / np.sqrt(2.), (1, 4, 1))
print("  h-stability of D (same config/direction):", flush=True)
for h in [3e-4, 1e-3, 3e-3]:
    print(f"    h={h:.1e}  D={D_plaq(Xone, Bone, h)[0]: .10e}", flush=True)

def W_of(B, X, hw=1e-5):
    def Up(s):
        U = expc(X + s * B)
        return U[0] @ U[1] @ dag(U[2]) @ dag(U[3])
    W = lambda s: (Up(s + hw) - Up(s - hw)) / (2 * hw) @ np.linalg.inv(Up(s))
    W0 = W(0.0); W1 = (W(hw) - W(-hw)) / (2 * hw)
    H = -1j * W0; H = 0.5 * (H + H.conj().T)
    return Up(0.0), H @ H, W1

rngb = np.random.default_rng(21)
Dv = np.empty(3000); Tv = np.empty(3000)
for i in range(3000):
    Xr = rngb.normal(size=(4, 3)) * 0.9
    Br = rngb.normal(size=(4, 3)); Br /= np.sqrt(0.5 * np.sum(Br ** 2))
    U0, K, W1 = W_of(Br, Xr)
    Dv[i] = D_plaq(Xr, Br[None])[0]
    Tv[i] = 0.5 * np.real(np.trace(K)) + 0.5 * abs(np.real(np.trace(U0 @ W1)))
print(f"  n=3000 random (X,v): max(-D)={np.max(-Dv):.6f}  max C_geo^tri={Tv.max():.6f}"
      f" ; C_geo<=C_geo^tri for all: {bool(np.all(np.maximum(-Dv,0)<=Tv+1e-9))}", flush=True)
print(f"  ratio sup(C_geo^tri)/sup(C_geo) = {Tv.max()/np.max(-Dv):.2f} (majorant looseness)", flush=True)
print(f"\n  total runtime {time.time()-t00:.1f}s", flush=True)
