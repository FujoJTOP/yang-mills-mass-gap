# HV2-PLAQ witness: plaquette-carrier ceiling bookkeeping
# (a) free-field row of the k_min-projected dimension-4 carrier, by exact Gaussian sampling (no lattice MC)
# (b) beta=0 exact Haar support, (c) three existing beta=2.5 SU(2) fields, (d) shell decomposition,
# (e) normalization cross-checks
import numpy as np, os, time
np.set_printoptions(precision=6, suppress=True, linewidth=200)
t0 = time.time()
def log(s): print(s, flush=True)

S2 = (4*np.pi**2-6)/9.0
log(f"[0] link-carrier anchor s^2=(4pi^2-6)/9 = {S2:.9f}")

# ------------------------------------------------------------------ [1] free-field row (Gaussian sampling)
def free_row(L, N=2, ns=3000, seed=11):
    V = L**4
    n = np.arange(L); k = 2*np.pi*n/L
    K = np.stack(np.meshgrid(k, k, k, k, indexing='ij'), axis=-1)
    kh = 2*np.sin(K/2.0); k2 = (kh*kh).sum(-1); k2[0, 0, 0, 0] = np.inf
    Proj = (np.eye(4)[None, None, None, None]-np.einsum('...m,...n->...mn', kh, kh)/k2[..., None, None])
    Kern = Proj/np.sqrt(k2)[..., None, None]
    rng = np.random.default_rng(seed)
    acck = np.zeros((L,)*4); sm = 0.0; sm2 = 0.0; ntot = 0
    for it in range(ns):
        eps = rng.normal(size=(L, L, L, L, 4))
        A = np.fft.ifftn(np.einsum('...mn,...n->...m', Kern, np.fft.fftn(eps, axes=(0, 1, 2, 3))), axes=(0, 1, 2, 3)).real
        A = np.moveaxis(A, -1, 0)
        th = np.zeros((L,)*4)
        for mu in range(4):
            for nu in range(mu+1, 4):
                Fu = (np.roll(A[nu], -1, axis=mu)-A[nu])-(np.roll(A[mu], -1, axis=nu)-A[mu])
                th += 2.0*Fu*Fu                     # ordered-pair sum = 2 * sum_{mu<nu}
        phi = th/(8.0*N)                            # 1-P = Theta/(8N), one colour
        acck += np.abs(np.fft.fftn(phi-phi.mean()))**2
        sm += phi.sum(); sm2 += (phi*phi).sum(); ntot += V
    mphi = sm/ntot
    C0 = (N*N-1)*(sm2/ntot - mphi*mphi)             # <Phi(0)^2>_c, all colours
    chi = (N*N-1)*acck/(ns*V)
    return dict(C0=C0, km=chi[1, 0, 0, 0], kh2=4*np.sin(np.pi/L)**2, chi=chi)

log("\n[1] free massless SU(2), carrier Phi = 1-P (6 planes summed) = (2/N)(1/4 tr F^2) + O(F^4)")
log("     Gaussian sampling, one colour, x(N^2-1); ns = samples")
log(" L | ns   | C_P(0)  | chi_P(kmin) | chi_til | f2=chi_til*kh2 | link free row 1/kh2 | chi_til/L^2")
rows1 = {}
for L, ns in ((4, 4000), (6, 4000), (8, 3000), (10, 1200), (12, 800)):
    d = free_row(L, 2, ns); rows1[L] = d
    ct = d['km']-d['C0']
    log(f"{L:2d} | {ns:4d} | {d['C0']:.5f} | {d['km']:.5f}   | {ct:+.5f} | {ct*d['kh2']:+10.5f}   | {1/d['kh2']:10.5f}         | {ct/L**2:+.5f}")
Ls = [L for L in sorted(rows1) if L >= 6]
e1 = np.polyfit(np.log(Ls), np.log([rows1[L]['km']-rows1[L]['C0'] for L in Ls]), 1)[0]
e2 = np.polyfit(np.log(Ls), np.log([(rows1[L]['km']-rows1[L]['C0'])*rows1[L]['kh2'] for L in Ls]), 1)[0]
log(f"[1b] d ln chi_til/d ln L = {e1:+.3f} ;  d ln f2/d ln L = {e2:+.3f}   (L=6..12)")
log("     reference link carrier (log U, A-linear): free chi = 1/khat^2 -> exponent +2 ; f2 -> exponent 0")

# ------------------------------------------------------------------ [2] beta=0 exact
log("\n[2] beta=0 (exact Haar)")
for N in (2, 3):
    log(f" N={N}: per-plane C(0)=1/N^2={1/N**2:.6f} ; 6-plane anchor 6/N^2={6/N**2:.6f} ;"
        f" s^2/C0_plane={S2*N**2:.4f} ; s^2/C0_6plane={S2*N**2/6:.4f}")
log(" two distinct plaquettes share <=1 link => >=6 links appear once => Haar 1st moment 0")
log(" => C(x)=0 for every x!=0, every N => chi_til(kmin)|beta=0 = 0 exactly")

# ------------------------------------------------------------------ [3] existing fields
def dag(A): return np.conjugate(np.swapaxes(A, -1, -2))
def plaq(U):
    d = {}
    for mu in range(4):
        for nu in range(mu+1, 4):
            W = U[mu] @ np.roll(U[nu], -1, axis=mu) @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])
            d[(mu, nu)] = np.trace(W, axis1=-2, axis2=-1).real/2.0
    return d

log("\n[3] existing beta=2.5 SU(2) fields, 1 config/volume (zero new MC)")
log(" L | <P>    | Var plane | Var 6pla | chi_P(kmin) orbit+spread | chi_til  | f2      | chi_til/L^2")
fields = {}
for L in (4, 6, 8):
    f = f"D:/Dev/math/research/v3/w7/gevp_4d_L{L}_field.npy"
    if not os.path.exists(f):
        log(f"{L:2d} | MISSING"); continue
    U = np.load(f); fields[L] = U
    pl = plaq(U)
    pbar = float(np.mean([pl[c].mean() for c in pl]))
    vp = float(np.mean([pl[c].var() for c in pl]))
    O = sum(1.0-pl[c] for c in pl); vs = float(O.var())
    F = np.fft.fftn(O-O.mean()); V = L**4
    modes = []
    for ax in range(4):
        for s in (1, -1):
            idx = [0, 0, 0, 0]; idx[ax] = s % L
            modes.append(abs(F[tuple(idx)])**2/V)
    modes = np.array(modes); cm = modes.mean(); sd = modes.std(ddof=1)/np.sqrt(len(modes))
    kh2 = 4*np.sin(np.pi/L)**2; ct = cm-vs
    log(f"{L:2d} | {pbar:.4f} | {vp:.4f}    | {vs:.4f}   | {cm:.4f} +- {sd:.4f}        | {ct:+.4f} | {ct*kh2:+.4f} | {ct/L**2:+.4f}")

# ------------------------------------------------------------------ [4] shell decomposition
log("\n[4] orbit-averaged chi_til = sum_{x!=0} w(x) C(x),  w(x) = (1/4) sum_mu cos(2pi x_mu/L)")
log(" bin   | beta=2.5 L=8 (1 cfg) | free L=8")
L = 8
U = fields[L]; pl = plaq(U); O = sum(1.0-pl[c] for c in pl); Om = O-O.mean(); V = L**4
Cfull = np.fft.ifftn(np.abs(np.fft.fftn(Om))**2).real/V
xs = np.indices((L,)*4).reshape(4, -1).T; r1 = xs.sum(1)
w = 0.25*sum(np.cos(2*np.pi*xs[:, mu]/L) for mu in range(4))
w = np.where(r1 == 0, 0.0, w)
for lo, hi in ((0, 0), (1, 1), (2, 2), (3, 4), (5, 8), (9, 999)):
    m = (r1 >= lo) & (r1 <= hi)
    Cfree = np.fft.ifftn(rows1[L]['chi']).real
    log(f" {lo}-{hi:<4d}| {(w*Cfull.reshape(-1))[m].sum():+.5f}              | {(w*Cfree.reshape(-1))[m].sum():+.5f}")
Cfree = np.fft.ifftn(rows1[L]['chi']).real
log(f" total | {(w*Cfull.reshape(-1)).sum():+.5f}              | {(w*Cfree.reshape(-1)).sum():+.5f}")
for c in [(0, 1), (0, 2), (2, 3)]:
    P = pl[c]; Pm = P-P.mean()
    log(f"   C_1(plane {c}) = {(Pm*np.roll(Pm,-1,axis=0)).mean():+.5f}  (per-plane Var={P.var():.5f})")

# ------------------------------------------------------------------ [5] normalization cross-checks
log("\n[5] cross-checks")
log(f"   identity (1/V) sum_k chi(k) = C(0): machine check on the L=8 field:")
U = fields[8]; pl = plaq(U); O = sum(1.0-pl[c] for c in pl); Om = O-O.mean(); V = 8**4
Fk = np.fft.fftn(Om); left = (np.abs(Fk)**2/V).mean(); right = Om.var()
log(f"      (1/V)sum_k |Otilde|^2/V = {left:.10f}   vs   Var(O) = {right:.10f}   diff={left-right:.2e}")
log("   free-row exponent vs the two branches of HR5B-CARRIER.md:11 (A: L^+4 or L^-2 ; B: L^+2)")

# ------------------------------------------------------------------ [6] noise projection
log("\n[6] sigma(f_P^2) = chi_P(kmin)*sqrt(2/N)/sqrt(8)*kh2  (measured beta=2.5 chi_P(kmin))")
for L, cm in ((4, 0.9337), (6, 0.9445), (8, 1.5450)):
    kh2 = 4*np.sin(np.pi/L)**2
    for Nc in (4000, 16000):
        log(f"   L={L} N={Nc}: sigma(f2)={cm*np.sqrt(2.0/Nc)/np.sqrt(8)*kh2:.4f}"
            f"   (HT4 link carrier: sigma(chi_til)=0.0144 -> 2sigma cap 0.055)")
log(f"\ndone in {time.time()-t0:.1f}s")
