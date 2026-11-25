"""mc_lattice.py -- T2 on real small lattices.  SU(2), beta=2.5, L=2 and L=4.

Question: on X = SU(2)^E, can ANY left-invariant metric shape G (the free data of route N1)
make   M_G := Hess_G S_W + Ric_G/2   positive definite, i.e. lambda_dot_0 > 0 ?

Exact Hessian, left-invariant frame:
    Hess f(e_a,e_b) = e_a e_b f - Gamma^c_ab(G) e_c f
The two pieces separate:
    e_a e_b S_W   ("D2") and  e_c S_W ("g") are G-INDEPENDENT  (computed once per configuration)
    the Koszul term -Gamma^c_ab(G) g_c and Ric_G/2 are G-dependent and BLOCK DIAGONAL
so M_G = D2 + blockdiag( -Gamma^c_ab(G) g_c[link] + Ric_G/2 ).

Configurations: pure gauge, the constant Z2 2-cocycle (U_1=(-1)^(x2+x3+x4), U_mu=(-1)^x1 for
mu>=2, so every (1,mu) plaque = -1 and every (mu,nu) plaque = +1), U=-1, Haar random, and the
interpolation angle s in A -> B.
Also recorded: sum_links |grad_link S_W|^2, the EXTENSIVE quantity behind the L-uniformity
discussion (T3).
"""
import numpy as np
import json
import scipy.sparse.linalg as spla

np.set_printoptions(precision=6, suppress=True, linewidth=170)
rng = np.random.default_rng(20261125)
BETA = 2.5
N = 2
NG = 3


def su_generators(N):
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = T[k, j] = 0.5; Ts.append(T)
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = -0.5j; T[k, j] = 0.5j; Ts.append(T)
    for l in range(1, N):
        T = np.zeros((N, N), dtype=complex)
        d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d
        Ts.append(T)
    return Ts


LAM = [1j * T for T in su_generators(N)]
I2 = np.eye(N, dtype=complex)
S3 = np.array([[1, 0], [0, -1]], complex)

_f_cache = {}


def gamma_ric(G):
    key = G.tobytes()
    if key in _f_cache:
        return _f_cache[key]
    n = NG
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            C = LAM[a] @ LAM[b] - LAM[b] @ LAM[a]
            f[:, a, b] = [-2.0 * np.trace(C @ LAM[c]).real for c in range(n)]
    Gi = np.linalg.inv(G)
    T = np.einsum('cd,dab->abc', G, f)
    Gl = 0.5 * (np.transpose(T, (1, 2, 0)) - np.transpose(T, (1, 0, 2)) - np.transpose(T, (2, 0, 1)))
    Gup = np.einsum('cd,dab->cab', Gi, Gl)
    Ric = (np.einsum('ccd,dab->ab', Gup, Gup) - np.einsum('cbd,dac->ab', Gup, Gup))
    _f_cache[key] = (Gup, Ric)
    return Gup, Ric


def plaq_list(L):
    out = []
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                for x4 in range(L):
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            x = [x1, x2, x3, x4]
                            xm = x.copy(); xm[mu] = (xm[mu] + 1) % L
                            xn = x.copy(); xn[nu] = (xn[nu] + 1) % L
                            li = lambda v, m: (v[0] * L ** 0 + v[1] * L + v[2] * L ** 2 + v[3] * L ** 3) * 4 + m
                            out.append(((li(x, mu), +1), (li(xm, nu), +1),
                                        (li(xn, mu), -1), (li(x, nu), -1)))
    return out


def mats_of(p, U):
    return [U[i] if sg > 0 else U[i].conj().T for i, sg in p]


def prod(lst):
    P = I2
    for m in lst:
        P = P @ m
    return P


def S_W(U, plaqs):
    return BETA * sum(1.0 - 0.5 * np.trace(prod(mats_of(p, U))).real for p in plaqs)


def build(U, plaqs):
    """g (E,NG) and D2 (dim,dim), both G-independent.
    e_a = left-invariant field, flow U -> U e^{t lam_a}.
      positive link U_j      : lam_a inserted immediately AFTER  M_j
      negative link U_j^{-1} : -lam_a inserted immediately BEFORE M_j
    and e_a e_b uses the order lam_a lam_b (NOT symmetrised: e_a e_b - e_b e_a = f^c_ab e_c,
    which is exactly what the Koszul term cancels)."""
    E = len(U)
    dim = E * NG
    g = np.zeros((E, NG))
    D2 = np.zeros((dim, dim))
    base = -BETA * 0.5
    for p in plaqs:
        M = mats_of(p, U)
        for k in range(4):
            ik, sgk = p[k]
            rot = M[k:] + M[:k]
            tail = prod(rot[1:] + [rot[0]])       # rest * M_k  (insert AFTER M_k)
            head = prod(rot)                      # M_k * rest  (insert BEFORE M_k)
            X = tail if sgk > 0 else head
            for a in range(NG):
                g[ik, a] += base * sgk * np.trace(LAM[a] @ X).real
            for b_ in range(NG):
                for a_ in range(NG):
                    if sgk > 0:
                        v = np.trace(LAM[a_] @ LAM[b_] @ tail).real
                    else:
                        v = np.trace(LAM[b_] @ LAM[a_] @ head).real
                    D2[ik * NG + a_, ik * NG + b_] += base * v
            for kk in range(4):
                if kk == k:
                    continue
                ik2 = p[kk][0]
                for a_ in range(NG):
                    for b_ in range(NG):
                        ins = []
                        for j in range(4):
                            ij, sj = p[j]
                            if j == k and sj < 0:
                                ins.append(-LAM[a_])
                            if j == kk and sj < 0:
                                ins.append(-LAM[b_])
                            ins.append(M[j])
                            if j == k and sj > 0:
                                ins.append(LAM[a_])
                            if j == kk and sj > 0:
                                ins.append(LAM[b_])
                        D2[ik * NG + a_, ik2 * NG + b_] += base * np.trace(prod(ins)).real
    return g, D2


def M_for(g, D2, G):
    Gup, Ric = gamma_ric(G)
    E = g.shape[0]
    M = D2.copy()
    for i in range(E):
        sl = slice(i * NG, (i + 1) * NG)
        M[sl, sl] += -np.einsum('cab,c->ab', Gup, g[i]) + 0.5 * Ric
    return 0.5 * (M + M.T)


def lam_min(M):
    if M.shape[0] <= 400:
        return float(np.linalg.eigvalsh(M)[0])
    return float(spla.eigsh(M, k=1, which='SA', maxiter=50000, tol=1e-7)[0][0])


def haar():
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R)
    Q = Q * (d / np.abs(d)).conj()[None, :]
    return Q / np.linalg.det(Q) ** (1.0 / N)


def config(L, kind, s=1.0):
    E = 4 * L ** 4
    if kind == 'pure':
        return [np.eye(N, dtype=complex) for _ in range(E)]
    if kind == 'minus1':
        return [-np.eye(N, dtype=complex) for _ in range(E)]
    if kind == 'haar':
        return [haar() for _ in range(E)]
    U = [np.eye(N, dtype=complex) for _ in range(E)]
    li = lambda v, m: (v[0] + v[1] * L + v[2] * L ** 2 + v[3] * L ** 3) * 4 + m
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                for x4 in range(L):
                    v = [x1, x2, x3, x4]
                    # constant Z2 2-cocycle with z_{1 mu} = -1 (mu >= 2): only the direction-1
                    # links are flipped, alternating in x2+x3+x4, so u1(x)u1(x+e_mu)^{-1} = -1.
                    th = s * np.pi * ((x2 + x3 + x4) % 2)
                    U[li(v, 0)] = np.cos(th) * I2 + 1j * np.sin(th) * S3
    return U


# shape families
SH_FULL = [("delta", np.eye(3))]
for a in np.linspace(-2, 2, 7):
    for b in np.linspace(-2, 2, 7):
        SH_FULL.append((f"diag({a:+.1f},{b:+.1f})", np.diag(np.exp([a, b, -a - b]))))
SH_COARSE = [("delta", np.eye(3))] + [SH_FULL[i] for i in (2, 6, 8, 24, 30, 40, 42, 44)]

def expL(t, a):
    """exact exp(t lam_a) for su(2): lam_a^2 = -I/4."""
    return np.cos(t / 2) * I2 + 2.0 * np.sin(t / 2) * LAM[a]


print("=" * 104)
print("SANITY: analytic g and D2 vs finite differences (L=2, random Haar config)")
L = 2
plaqs = plaq_list(L)
U0 = config(L, 'haar')
g, D2 = build(U0, plaqs)
h = 1e-6
err_g = 0.0
for _ in range(20):
    i = int(rng.integers(len(U0))); a = int(rng.integers(NG))
    Up = list(U0); Up[i] = U0[i] @ expL(h, a)
    Um = list(U0); Um[i] = U0[i] @ expL(-h, a)
    fd = (S_W(Up, plaqs) - S_W(Um, plaqs)) / (2 * h)
    err_g = max(err_g, abs(fd - g[i, a]))
err_d = 0.0
h2 = 1e-3
for _ in range(12):
    i = int(rng.integers(len(U0))); a, b = int(rng.integers(NG)), int(rng.integers(NG))
    def val(sa, sb):
        V = list(U0); V[i] = U0[i] @ expL(sa, a) @ expL(sb, b); return S_W(V, plaqs)
    fd = (val(h2, h2) - val(h2, 0) - val(0, h2) + val(0, 0)) / h2 ** 2
    err_d = max(err_d, abs(fd - D2[i * NG + a, i * NG + b]))
err_m = 0.0
prs = [(p[0][0], p[1][0]) for p in plaqs]           # link pairs sharing a plaquette
for _ in range(12):
    i, j = prs[int(rng.integers(len(prs)))]
    a, b = int(rng.integers(NG)), int(rng.integers(NG))
    def val2(sa, sb):
        V = list(U0); V[i] = U0[i] @ expL(sa, a); V[j] = U0[j] @ expL(sb, b)
        return S_W(V, plaqs)
    fd = (val2(h2, h2) - val2(h2, 0) - val2(0, h2) + val2(0, 0)) / h2 ** 2
    err_m = max(err_m, abs(fd - D2[i * NG + a, j * NG + b]))
print(f"  |fd - g| = {err_g:.2e}   |fd - D2_diag| = {err_d:.2e}   |fd - D2_mixed| = {err_m:.2e}"
      f"   (h={h:g} for g, {h2:g} for D2)")
Gup0, Ric0 = gamma_ric(np.eye(3))
Hh = M_for(g, D2, np.eye(3))
print(f"  asymmetry of the assembled Hessian: max|H-H^T| = {np.abs(Hh - Hh.T).max():.2e}")
print(f"  and lambda_min(Hess S_W + Ric/2) at delta = {np.linalg.eigvalsh(Hh)[0]:+.6f}")

report = {}
for L, shapes in ((2, SH_FULL), (4, SH_COARSE)):
    plaqs = plaq_list(L)
    E = 4 * L ** 4
    print("=" * 104)
    print(f"L={L}  E={E} links  dim={E*NG}  #plaques={len(plaqs)}  SU(2) beta={BETA}  "
          f"({len(shapes)} shapes)")
    cfgs = {'pure': config(L, 'pure'), 'z2': config(L, 'z2'), 'minus1': config(L, 'minus1')}
    for sd in range(2):
        cfgs[f'haar{sd}'] = config(L, 'haar')
    for s in (0.25, 0.5, 0.75):
        cfgs[f'z2_s{s}'] = config(L, 'z2', s)
    data = {}
    for name, U in cfgs.items():
        g, D2 = build(U, plaqs)
        mu = {sn: lam_min(M_for(g, D2, G)) for sn, G in shapes}
        best = max(mu, key=mu.get)
        data[name] = dict(S=float(S_W(U, plaqs)), mu_delta=float(mu["delta"]),
                          sup_shape=float(mu[best]), best_shape=best,
                          gradnorm=float(np.sum(g ** 2)))
        print(f"  {name:9s} S_W={data[name]['S']:8.3f} mu(delta)={mu['delta']:+9.4f} "
              f"sup_G mu={mu[best]:+9.4f} [{best:>16s}]  sum|gradS|^2={data[name]['gradnorm']:9.2f}")
    print(f"  ==> inf_config sup_G mu  = {min(v['sup_shape'] for v in data.values()):+.6f}"
          "   (< 0  => no metric shape can give lambda_dot_0 > 0)")
    print(f"      inf_config mu(delta) = {min(v['mu_delta'] for v in data.values()):+.6f}")
    print(f"      max_config sum|grad S_W|^2 = {max(v['gradnorm'] for v in data.values()):9.2f}"
          f"   (per link: {max(v['gradnorm'] for v in data.values())/E:7.2f})")
    report[f"L{L}"] = data

with open("mc_lattice.json", "w") as fh:
    json.dump(report, fh, indent=1)
print("\nwrote mc_lattice.json")
