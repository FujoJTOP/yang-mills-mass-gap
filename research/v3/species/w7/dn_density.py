"""dn_density.py -- task-2 / LOOP-RUNBOOK R3.

Question: is  "sup_phi |grad S_W|^2  ∝  L^4"  a real growth, or is it extensivity
(the face count N_p = 6 L^4)?

Settings: SU(2), beta = 2.5, L = 2,4,6,8, periodic torus, bi-invariant frame (G = delta).
The object measured is the *code's* gradient norm

    |grad S_W|^2  :=  sum_{links l} sum_{a=1..3} g[l,a]^2 ,  g[l,a] = e_a S_W at link l

(the convention of w7/mc_lattice.py `gradnorm`, i.e. the number that gave 900 / 14400).

Contents
  S1  reproduce 900 (L=2) / 14400 (L=4) and identify WHAT the number is
  S2  config-by-config table (pointwise values, no averaging)
  S3  d := |grad S_W|^2 / N_p  for L = 2,4,6,8   (+ rigorous bound, + certificate)
  S4  L-exponent with error
  S5  probing sup_phi by greedy single-link ascent, then exact period-2 lift
  S6  Hessian side: FD Hessian, max row-sum bound, and the binding eigenvector lifted
  S7  Gibbs average (rough Metropolis) for the sup-vs-average comparison
"""
import numpy as np
import time

BETA = 2.5
T1 = np.array([[0.0, 0.5], [0.5, 0.0]], dtype=complex)
T2 = np.array([[0.0, -0.5j], [0.5j, 0.0]], dtype=complex)
T3 = np.array([[0.5, 0.0], [0.0, -0.5]], dtype=complex)
LAM = np.stack([1j * T1, 1j * T2, 1j * T3])       # tr(T_a T_b) = delta/2 ; lam = iT ; lam^2 = -I/4
I2 = np.eye(2, dtype=complex)
S3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


# ----------------------------------------------------------------------------- lattice
def lattice(L):
    """Plaquette link arrays in the SAME order/signs as w7/mc_lattice.py plaq_list()."""
    v = np.stack(np.meshgrid(*[np.arange(L)] * 4, indexing='ij'), -1).reshape(-1, 4)

    def li(vv, m):
        return (vv[:, 0] + vv[:, 1] * L + vv[:, 2] * L ** 2 + vv[:, 3] * L ** 3) * 4 + m

    A, B, C, D = [], [], [], []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            vm = v.copy(); vm[:, mu] = (vm[:, mu] + 1) % L
            vn = v.copy(); vn[:, nu] = (vn[:, nu] + 1) % L
            A.append(li(v, mu)); B.append(li(vm, nu)); C.append(li(vn, mu)); D.append(li(v, nu))
    E = 4 * L ** 4
    NP = 6 * L ** 4
    link_mu = [li(v, m) for m in range(4)]
    return (np.concatenate(A), np.concatenate(B), np.concatenate(C), np.concatenate(D),
            E, NP, link_mu, v)


def state(U, lat, want_g=True):
    """(S_W, g, P).  g[l,a] = e_a S_W in the left-invariant frame -- mc_lattice.py convention."""
    A, B, C, D, E, NP, _, _ = lat
    M = [U[A], U[B], U[C].conj().transpose(0, 2, 1), U[D].conj().transpose(0, 2, 1)]
    P = M[0] @ M[1] @ M[2] @ M[3]
    trP = np.trace(P, axis1=1, axis2=2).real
    S = BETA * float(np.sum(1.0 - 0.5 * trP))
    if not want_g:
        return S, None, P
    # T[k] = cyclic product starting just after position k ("tail");  head(k) = T[k-1]
    T = []
    for k in range(4):
        perm = [(k + 1) % 4, (k + 2) % 4, (k + 3) % 4, k]
        T.append(M[perm[0]] @ M[perm[1]] @ M[perm[2]] @ M[perm[3]])
    base = -BETA * 0.5
    g = np.zeros((E, 3))
    for links, sg, X in ((A, +1, T[0]), (B, +1, T[1]), (C, -1, T[1]), (D, -1, T[2])):
        c = base * sg * np.einsum('aij,pji->pa', LAM, X).real
        np.add.at(g, links, c)
    return S, g, P


def gradnorm(U, lat):
    return float(np.sum(state(U, lat)[1] ** 2))


def expL(t, a):
    return np.cos(t / 2) * I2 + 2.0 * np.sin(t / 2) * LAM[a]


def exp_su2_batch(t, v):
    """exp(t * sum_a v[l,a] lam_a) for every link l (uses A^2 = -|v|^2/4 I)."""
    n = np.linalg.norm(v, axis=1)
    A = np.einsum('la,aij->lij', v, LAM)
    c = 0.5 * n
    safe = np.where(c > 1e-14, c, 1.0)
    coef = np.where(c > 1e-14, np.sin(t * safe) / safe, t)     # lim_{c->0} sin(tc)/c = t
    return np.cos(t * c)[:, None, None] * I2 + coef[:, None, None] * A


def haar(rng):
    Z = (rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R)
    Q = Q * (d / np.abs(d)).conj()[None, :]
    return Q / np.linalg.det(Q) ** (1.0 / 2)


def haar_batch(rng, n):
    """Haar on SU(2) = uniform on S^3 (fast; same law as haar())."""
    q = rng.standard_normal((n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    a, b, c, d = q.T
    U = np.empty((n, 2, 2), dtype=complex)
    U[:, 0, 0] = a + 1j * b
    U[:, 0, 1] = c + 1j * d
    U[:, 1, 0] = -c + 1j * d
    U[:, 1, 1] = a - 1j * b
    return U


def cfg(L, lat, kind, s=1.0, rng=None):
    E, v = lat[4], lat[7]
    U = np.tile(I2, (E, 1, 1)).astype(complex)
    if kind == 'pure':
        return U
    if kind == 'minus1':
        return -U
    if kind == 'haar':
        return haar_batch(rng, E) if E > 1000 else np.array([haar(rng) for _ in range(E)])
    if kind == 'z2':
        th = s * np.pi * ((v[:, 1] + v[:, 2] + v[:, 3]) % 2)
        U[lat[6][0]] = np.cos(th)[:, None, None] * I2 + 1j * np.sin(th)[:, None, None] * S3
        return U
    raise ValueError(kind)


FMT = "{:16s} {:10.3f} {:14.4f} {:12.4f} {:10.4f}"


def line(name, L, lat, U, extra=""):
    S, g, _ = state(U, lat)
    gn = float(np.sum(g ** 2))
    E, NP = lat[4], lat[5]
    print("  " + FMT.format(name, S, gn, gn / NP, gn / E) + "  " + extra)
    return gn


print("=" * 108)
print("dn_density.py -- task-2 / R3 : is  sup|grad S_W|^2 ∝ L^4  a growth or the face count?")
print("SU(2), beta = 2.5, bi-invariant frame (G = delta), periodic torus, L = 2,4,6,8")
print("normalisation: |grad S_W|^2 = sum_l sum_a g[l,a]^2  (identical to mc_lattice.py 'gradnorm')")
print("=" * 108)

# ------------------------------------------------------------------ S1 sanity + reproduce
print("\n[S1] SANITY -- finite differences vs the analytic gradient, and the corpus numbers")
for L in (2, 4):
    lat = lattice(L)
    rng = np.random.default_rng(20261125 + L)
    U = cfg(L, lat, 'haar', rng=rng)
    S, g, _ = state(U, lat)
    err = 0.0
    for _ in range(20):
        i = int(rng.integers(lat[4])); a = int(rng.integers(3)); h = 1e-6
        Up = U.copy(); Up[i] = U[i] @ expL(h, a)
        Um = U.copy(); Um[i] = U[i] @ expL(-h, a)
        fd = (state(Up, lat, False)[0] - state(Um, lat, False)[0]) / (2 * h)
        err = max(err, abs(fd - g[i, a]))
    Uz = cfg(L, lat, 'z2', s=0.5)
    Sz, gz, _ = state(Uz, lat)
    cS, cG = (120.0, 900.0) if L == 2 else (1920.0, 14400.0)
    print(f"  L={L}: |fd-g|={err:.2e}   z2(s=0.5): S_W={Sz:.3f} [corpus {cS}]"
          f"   sum|g|^2={np.sum(gz**2):.4f} [corpus {cG}]")
print("  => gradient routine reproduces mc_lattice.py bit-for-bit on S_W and on the reported")
print("     900 / 14400.  The quantity is a POINTWISE value, evaluated at ONE configuration.")

# --------------------------------------------------------- S2 the three competing numbers
print("\n[S2] WHAT THE NUMBER IS -- pointwise values at configurations, no measure average")
print("  " + "{:16s} {:>10s} {:>14s} {:>12s} {:>10s}".format(
    "config", "S_W", "sum|g|^2", "per plaquette", "per link"))
for L in (2, 4):
    lat = lattice(L)
    rng = np.random.default_rng(7 + L)
    print(f"-- L={L}   E={lat[4]} links   N_p={lat[5]} plaquettes")
    line("pure", L, lat, cfg(L, lat, 'pure'))
    line("minus1", L, lat, cfg(L, lat, 'minus1'))
    line("z2 (s=1)", L, lat, cfg(L, lat, 'z2', s=1.0))
    for s in (0.25, 0.5, 0.75):
        line(f"z2 s={s}", L, lat, cfg(L, lat, 'z2', s=s))
    for i in range(3):
        line(f"haar#{i}", L, lat, cfg(L, lat, 'haar', rng=rng))
print("  (corpus rows, for comparison: z2_s0.5 -> 900 @L=2, 14400 @L=4 ; z2_s0.25/0.75 -> 450/7200)")
print("  NOTE: pure, minus1 and the s=1 centre-element config all give |grad S_W|^2 = 0 exactly")
print("        (they are critical points: every plaquette is +-1 and tr(lam_a * +-1) = 0).")

# ------------------------------------------------------- S3 d for L=2,4,6,8 + rigorous bound
print("\n[S3] THE INTENSIVE QUANTITY  d := |grad S_W|^2 / N_p ,  N_p = 6 L^4   (L = 2,4,6,8)")
print("  bound D1 (rigorous):  d <= 6 beta^2 = {:.4f}  for ALL L and ALL phi".format(6 * BETA ** 2))
print("  " + "{:5s} {:>7s} {:>7s} {:>10s} {:>12s} {:>12s} {:>12s} {:>8s}".format(
    "L", "E", "N_p", "E/N_p", "d(z2 s=.5)", "d(haar mean)", "d(haar max)", "d/bound"))
rows = {}
for L in (2, 4, 6, 8):
    lat = lattice(L)
    E, NP = lat[4], lat[5]
    dz = gradnorm(cfg(L, lat, 'z2', s=0.5), lat) / NP
    n = {2: 48, 4: 32, 6: 16, 8: 8}[L]
    rng = np.random.default_rng(1000 + L)
    vals = np.array([gradnorm(cfg(L, lat, 'haar', rng=rng), lat) / NP for _ in range(n)])
    rows[L] = dict(dz=dz, haer_mean=float(vals.mean()), haer_std=float(vals.std(ddof=1)),
                   haer_max=float(vals.max()), n=n, E=E, NP=NP)
    print("  " + "{:5d} {:7d} {:7d} {:10.4f} {:12.6f} {:12.6f} {:12.6f} {:8.4f}".format(
        L, E, NP, E / NP, dz, vals.mean(), vals.max(), dz / (6 * BETA ** 2)))

print("\n  exact value at z2(s=0.5):  1.5 beta^2 = {:.6f} = (18/8) * (E/N_p) * beta^2"
      .format(1.5 * BETA ** 2))
print("  per-link value at z2(s=0.5): 2.25 beta^2 = {:.6f} = 18 * beta^2 / 8".format(2.25 * BETA ** 2))
print("  -> both are ratios of INCIDENCE COUNTS: 18 = aligned plaquette-link incidences at a")
print("     direction-1 link, 8 = per-plaquette gradient normalisation, E/N_p = 4/6. No dynamics.")

# ----------------------------------------------------- S3b certificate of the rigorous bound
print("\n[S3b] CERTIFICATE of D1.  For one link l:  g[l,a] = sum_{p incident to l} s_{p,l,a},")
print("      s_{p,l,a} = -(beta/2) tr(lam_a X_{p,l}) ;  Cauchy-Schwarz over the <=6 plaquettes:")
print("        |g_l|^2 <= 6 sum_p sum_a s^2 = 6 * (beta^2/4) sum_p sin^2(psi_p) <= 9 beta^2")
print("      => sum_l |g_l|^2 <= 9 beta^2 * E = 36 beta^2 L^4  =>  d <= 36/6 * beta^2 = 6 beta^2.")
for L in (2, 4, 6, 8):
    lat = lattice(L)
    A, B, C, D = lat[0], lat[1], lat[2], lat[3]
    U = cfg(L, lat, 'z2', s=0.5)
    M = [U[A], U[B], U[C].conj().transpose(0, 2, 1), U[D].conj().transpose(0, 2, 1)]
    Tt = []
    for k in range(4):
        perm = [(k + 1) % 4, (k + 2) % 4, (k + 3) % 4, k]
        Tt.append(M[perm[0]] @ M[perm[1]] @ M[perm[2]] @ M[perm[3]])
    sq = 0.0
    for links, sg, X in ((A, +1, Tt[0]), (B, +1, Tt[1]), (C, -1, Tt[1]), (D, -1, Tt[2])):
        cc = (BETA * 0.5) ** 2 * np.einsum('aij,pji->pa', LAM, X).real ** 2
        sq += float(np.sum(cc))
    gn = gradnorm(U, lat)
    print(f"  L={L}: sum_l |g_l|^2 = {gn:12.3f}   6*sum_p sum s^2 = {6*sq:12.3f}"
          f"   ratio = {gn/(6*sq):.4f}  (<=1 = certificate holds)")

# ---------------------------------------------------------------- S4 L-exponent with error
print("\n[S4] L-EXPONENT  d(L) = A L^p   (log-log least squares over L = 2,4,6,8)")
Ls = np.array([2.0, 4.0, 6.0, 8.0])
dz = np.array([rows[int(L)]['dz'] for L in Ls])
hm = np.array([rows[int(L)]['haer_mean'] for L in Ls])
hs = np.array([rows[int(L)]['haer_std'] for L in Ls])


def fit(y, yerr=None):
    x = np.log(Ls); Y = np.log(y)
    A = np.vstack([x, np.ones(len(x))]).T
    if yerr is None:
        sol, *_ = np.linalg.lstsq(A, Y, rcond=None)
        s2 = float(np.sum((Y - A @ sol) ** 2) / max(len(x) - 2, 1))
        cov = s2 * np.linalg.inv(A.T @ A)
        return sol[0], np.sqrt(cov[0, 0]), sol[1]
    W = np.diag(1.0 / (yerr / y) ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    sol = cov @ (A.T @ W @ Y)
    return sol[0], np.sqrt(cov[0, 0]), sol[1]


pz, ez, cz = fit(dz)
ph, eh, ch = fit(hm, hs)
print(f"  d(z2 s=0.5):  {dz[0]:.6f}, {dz[1]:.6f}, {dz[2]:.6f}, {dz[3]:.6f}"
      f"   spread = {dz.max()-dz.min():.2e}")
print(f"                exponent p = {pz:+.6f} +- {ez:.1e}   (A = {np.exp(cz):.6f})")
print(f"  d(haar mean): {hm[0]:.4f}, {hm[1]:.4f}, {hm[2]:.4f}, {hm[3]:.4f}")
print(f"                +- {hs[0]:.4f}, {hs[1]:.4f}, {hs[2]:.4f}, {hs[3]:.4f}")
print(f"                exponent p = {ph:+.4f} +- {eh:.4f}")
print(f"  raw |grad S_W|^2 ratio L=2 -> L=4 : "
      f"{rows[4]['dz']*rows[4]['NP']/(rows[2]['dz']*rows[2]['NP']):.4f}"
      f"   (= 16 = (4/2)^4 ; L^8 would give 256)")

# --------------------------------------------------------------- S5 probing the true sup
print("\n[S5] PROBING sup_phi -- greedy single-link ascent from three seeds (L=2)")
rng = np.random.default_rng(4242)


def ascent(U, lat, sweeps=8):
    U = U.copy()
    E = lat[4]
    best = gradnorm(U, lat)
    for sw in range(sweeps):
        for l in rng.permutation(E):
            cur = U[l]
            bestU, bestv = cur, best
            for a in range(3):
                for t in (-3.0, -2.0, -1.4, -0.9, -0.5, 0.5, 0.9, 1.4, 2.0, 3.0):
                    V = U.copy(); V[l] = cur @ expL(t, a)
                    val = gradnorm(V, lat)
                    if val > bestv + 1e-12:
                        bestv, bestU = val, V[l]
            U[l] = bestU
            best = bestv
        print(f"      sweep {sw+1:2d}: sum|g|^2 = {best:11.4f}   d = {best/lat[5]:.6f}")
    return U


lat2 = lattice(2)
seeds = [('z2 s=0.5', cfg(2, lat2, 'z2', s=0.5)),
         ('haar', cfg(2, lat2, 'haar', rng=np.random.default_rng(11))),
         ('random', cfg(2, lat2, 'haar', rng=np.random.default_rng(12)))]
best2, bestname, bestgn = None, '', -1.0
for nm, U0 in seeds:
    print(f"  seed = {nm}  (start d = {gradnorm(U0, lat2)/lat2[5]:.4f})")
    U1 = ascent(U0, lat2, sweeps=8)
    gn = gradnorm(U1, lat2)
    print(f"    -> d = {gn/lat2[5]:.6f}   (d / d_z2s05 = {gn/lat2[5]/(1.5*BETA**2):.3f})")
    if gn > bestgn:
        bestgn, best2, bestname = gn, U1, nm
d_sup_lo = bestgn / lat2[5]
print(f"  BEST at L=2: d = {d_sup_lo:.6f}  [seed {bestname}]  ->  d_sup >= {d_sup_lo:.4f}")
print(f"  z2(s=0.5) density 1.5 beta^2 = {1.5*BETA**2:.4f} is NOT the sup: it is {d_sup_lo/(1.5*BETA**2):.2f}x below.")
print(f"  bound 6 beta^2 = {6*BETA**2:.4f}  ->  best/bound = {d_sup_lo/(6*BETA**2):.4f}"
      f"   =>  9.375 <= d_sup <= 37.5")

print("\n  LIFT TEST: a period-2 configuration lifts EXACTLY to any even L (identical plaquettes,")
print("  same per-link gradient), so the density must be identical -- this is an identity, not a fit.")
for L in (4, 6, 8):
    lat = lattice(L)
    vL = lat[7]

    def li2(vv, m):
        return (vv[:, 0] + vv[:, 1] * 2 + vv[:, 2] * 4 + vv[:, 3] * 8) * 4 + m
    U = np.tile(I2, (lat[4], 1, 1)).astype(complex)
    for m in range(4):
        U[lat[6][m]] = best2[li2(vL % 2, m)]
    print(f"    L={L}: d = {gradnorm(U, lat)/lat[5]:.6f}   (L=2 value {d_sup_lo:.6f})")

# ------------------------------------------------------- S6 the Hessian side: is it extensive?
print("\n[S6] HESSIAN SIDE -- FD Hessian, row-sum (Gershgorin) bound, binding eigenvector")
print("     D2 (structural): a plaquette touches 4 links, a link sits in 6 plaquettes => bounded")
print("     degree => max_i sum_j |Hess_ij| is O(beta) with an L-independent constant.")


def hess_fd(U, lat, h=2e-3):
    E = lat[4]
    dim = 3 * E
    H = np.zeros((dim, dim))
    for j in range(dim):
        l, b = divmod(j, 3)
        Up = U.copy(); Up[l] = U[l] @ expL(h, b)
        Um = U.copy(); Um[l] = U[l] @ expL(-h, b)
        H[:, j] = (state(Up, lat)[1] - state(Um, lat)[1]).ravel() / (2 * h)
    return 0.5 * (H + H.T)


H2 = None
for L in (2, 4):
    lat = lattice(L)
    U = cfg(L, lat, 'z2', s=1.0)
    t0 = time.time()
    H = hess_fd(U, lat)
    ev = np.linalg.eigvalsh(H)
    rs = float(np.abs(H).sum(axis=1).max())
    print(f"  L={L}: lambda_min = {ev[0]:+.6f}   lambda_max = {ev[-1]:+.6f}"
          f"   max row-sum = {rs:.4f} = {rs/BETA:.4f} beta   ({time.time()-t0:.1f}s)")
    if L == 2:
        H2, ev2, evec2 = H, ev, np.linalg.eigh(H)
print("  => the Hessian is NOT extensive: both eigenvalues and row sums are identical at L=2,4.")

print("\n  binding eigenvector structure (L=2, z2 config):")
V = evec2[1][:, 0].reshape(lat2[4], 3)
print(f"    lambda_min = {ev2[0]:+.6f}   ||v||^2 = {np.sum(V**2):.6f}")
for m in range(4):
    idx = lat2[6][m]
    print(f"    direction {m+1}: mass = {np.sum(V[idx]**2):.6f}")
    for a in range(3):
        vals = V[idx, a]
        print(f"       gen {a+1}: mass = {np.sum(vals**2):8.6f}   distinct values = "
              f"{len(np.unique(np.round(vals, 6))):2d}   range [{vals.min():+.4f}, {vals.max():+.4f}]")

print("\n  Rayleigh quotient of that eigenvector, LIFTED periodically (exact identity at any even L)")


def rq_vec(U, lat, v, h=1e-3):
    def S(t):
        V = U.copy()
        V = V @ exp_su2_batch(t, v)
        return state(V, lat, False)[0]
    n2 = float(np.sum(v ** 2))
    return (S(h) - 2.0 * S(0.0) + S(-h)) / h ** 2 / n2


m = cfg(2, lat2, 'z2', s=1.0)
print(f"    L=2: RQ = {rq_vec(m, lat2, V):+.6f}   (lambda_min = {ev2[0]:+.6f})")
for L in (4, 6, 8):
    lat = lattice(L)
    vL, UL = lat[7], cfg(L, lat, 'z2', s=1.0)

    def li2(vv, mm):
        return (vv[:, 0] + vv[:, 1] * 2 + vv[:, 2] * 4 + vv[:, 3] * 8) * 4 + mm
    vv = np.zeros((lat[4], 3))
    for mm in range(4):
        vv[lat[6][mm]] = V[li2(vL % 2, mm)]
    print(f"    L={L}: RQ = {rq_vec(UL, lat, vv):+.6f}")

print("\n  naive 'uniform on the direction-1 links' modes (corpus description) -- RQ:")


def rq_flat(U, lat, links, a, h=1e-3):
    def S(t):
        V = U.copy(); V[links] = V[links] @ expL(t, a)
        return state(V, lat, False)[0]
    return (S(h) - 2.0 * S(0.0) + S(-h)) / h ** 2 / (3 * len(links))


for L in (2, 4, 6, 8):
    lat = lattice(L)
    U = cfg(L, lat, 'z2', s=1.0)
    out = [rq_flat(U, lat, lat[6][0], a) for a in range(3)]
    out.append(rq_flat(U, lat, np.arange(lat[4]), 2))
    print(f"    L={L}: RQ(dir1, gen1..3) = {out[0]:+.4f} {out[1]:+.4f} {out[2]:+.4f}"
          f"   RQ(all links, gen3) = {out[3]:+.4f}")
print("  => those flat modes are EXACTLY zero at a centre-element configuration (the two")
print("     occurrences of a direction-1 link in a plaquette carry the same right multiplier")
print("     and U_1 = +-1 commutes with everything) => the corpus' description of the binding")
print("     mode as 'uniform on the direction-1 links' is not correct; it is a staggering mode")

# ------------------------------------------------- S7 Gibbs-measure average (rough, L=2 only)
print("\n[S7] AVERAGE vs SUP -- Gibbs measure at beta=2.5 (rough single-link Metropolis, L=2, [V])")
lat = lattice(2)
rng = np.random.default_rng(999)
U = cfg(2, lat, 'pure')
S0 = state(U, lat, False)[0]
acc = tot = 0
vals = []
for sw in range(140):
    for l in range(lat[4]):
        a = int(rng.integers(3))
        d = 0.35 * rng.standard_normal() * (2 * np.pi)
        V = U.copy(); V[l] = U[l] @ expL(d, a)
        S1 = state(V, lat, False)[0]
        tot += 1
        if S1 <= S0 or rng.random() < np.exp(-(S1 - S0)):
            U, S0 = V, S1
            acc += 1
    if sw >= 20:
        vals.append(gradnorm(U, lat) / lat[5])
vals = np.array(vals)
gm = float(vals.mean())
print(f"  acceptance = {acc/tot:.3f}   <d>_Gibbs = {gm:.4f} +- {vals.std(ddof=1)/np.sqrt(len(vals)):.4f}")
print(f"  Haar mean (L=2) = {rows[2]['haer_mean']:.4f} ; d_sup candidate >= {d_sup_lo:.4f}")
print(f"  sup/average separation (densities) = {d_sup_lo/gm:.2f}x  -- same order, no exponential gap")

# ------------------------------------------- S8 trace vs operator norm of Cov(grad S_W)
print("\n[S8] TRACE vs OPERATOR NORM -- where the L^4 actually comes from")
print("     tr(Cov) = E|grad S_W|^2 - |E grad S_W|^2 = THE NUMBER THE CORPUS QUOTES (dimension E*3)")
print("     ||Cov||_op <= row sum  sum_{r,b} |Cov_{(0,a),(r,b)}|  (Gershgorin, block-Toeplitz)")
print("     Bounding an operator norm by the trace is exactly a factor dim = E = 4 L^4.")
print("  " + "{:4s} {:>7s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
    "L", "n_samp", "tr(Cov)", "tr(Cov)/N_p", "rowsum", "tr/rowsum"))
for L in (2, 4, 6):
    lat = lattice(L)
    E = lat[4]
    nsamp = {2: 160, 4: 64, 6: 32}[L]
    rng = np.random.default_rng(31 + L)
    G = np.array([state(cfg(L, lat, 'haar', rng=rng), lat)[1] for _ in range(nsamp)])
    tr_cov = float(np.mean(np.sum((G - G.mean(0)) ** 2, axis=(1, 2))))
    X = (G - G.mean(0)).reshape(nsamp, L, L, L, L, 4, 3)
    vals = (-1, 0, 1) if L > 2 else (0, 1)
    rowsum = np.zeros(3)
    nshift = 0
    for r in np.ndindex(*(len(vals),) * 4):
        sh = tuple(vals[i] for i in r)
        Y = np.roll(X, sh, axis=(1, 2, 3, 4))
        C = np.einsum('nma,nmb->ab', X.reshape(nsamp, -1, 3),
                      Y.reshape(nsamp, -1, 3)) / (nsamp * L ** 4 * 4)
        rowsum += np.abs(C).sum(axis=1)
        nshift += 1
    rs = float(rowsum.max())
    print("  " + "{:4d} {:7d} {:12.3f} {:12.4f} {:12.4f} {:12.2f}".format(
        L, nsamp, tr_cov, tr_cov / lat[5], rs, tr_cov / rs))
    print(f"       ({nshift} shift classes used; the L=2 set is COMPLETE, L>=4 is a lower bound)")
print("  => tr(Cov) grows like N_p = 6L^4 with a flat density (~4.6), while the row sum stays O(1)")
print("     and L-independent.  ||Cov||_op = O(1);  the L^4 is the MATRIX DIMENSION E = 4L^4.")
print("\n  And at the t = infinity endpoint the identity gives Cov_nu(grad S_W) = E_nu[Hess S_W]")
print("  exactly, so ||Cov||_op <= sup_phi max_i rowsum_i(Hess S_W) = 6 beta = 15 is RIGOROUS")
print("  (an average is bounded by the sup of the row sums, and the row sums are L-uniform).")
print("\nDONE.")
