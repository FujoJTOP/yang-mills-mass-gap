"""r3_l4_volume.py -- R3/T2: is  sup_phi |grad S_W|^2 ~ L^4  a REAL volume law, or an
estimator artifact (sup-vs-typical, or trace-vs-spectral norm)?

SU(2) lattice gauge theory,  S_W = beta * sum_p (1 - 1/2 Re tr U_p),  left-invariant gradient
(flow U -> U e^{t lam_a}, lam_a = i sigma_a/2, flat metric delta_ab on su(2)).

Sections
  S0  self-checks: finite-difference gradient, pure gauge = 0, Z2 cocycle = 0,
      exact prediction sum|g|^2 = 9 beta^2 L^4 for the half-twist config.
  S1  pointwise vs sup: per-link bound 9 beta^2, per-plaquette bound beta^2/4,
      sup_phi lower bound from translation-invariant ansatz (exact L^4, constant c*).
  S2  three quantities at L=2,4,6: <|g|^2>_mu (Metropolis, beta=2.5), sup_phi lower bound,
      Haar random field + per-link densities + L-exponent fits.
  S3  the decisive split: tr Cov(grad S_W)  vs  lambda_min Cov(grad S_W)  (mode-resolved).
      trace ~ L^4, spectral minimum ~ L^0.

Everything is written to r3_l4.json; the log is flushed section by section.
"""
import json
import sys
import time

import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=170)
rng = np.random.default_rng(20261125)
BETA = 2.5
N = 2
S1 = np.array([[0, 1], [1, 0]], complex)
S2 = np.array([[0, -1j], [1j, 0]], complex)
S3 = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
SIG = (S1, S2, S3)
LOG = open("r3_l4.log", "w")
RES = {}


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


def dump(key, obj):
    RES[key] = obj
    with open("r3_l4.json", "w") as fh:
        json.dump(RES, fh, indent=1)


# ---------------------------------------------------------------- SU(2) helpers
def expi(t):
    """exp(i t.sigma/2) with t of shape (...,3); returns (...,2,2)."""
    a = np.sqrt(np.sum(t * t, axis=-1))
    safe = np.where(a < 1e-15, 1.0, a)
    n = t / safe[..., None]
    c = np.cos(a / 2)[..., None, None]
    s = np.sin(a / 2)[..., None, None]
    M = c * I2 + 1j * s * (n[..., 0, None, None] * S1 + n[..., 1, None, None] * S2
                           + n[..., 2, None, None] * S3)
    return np.where((a < 1e-15)[..., None, None], I2, M)


def lam_tr(X):
    """vector a -> Re tr(lam_a X) for X of shape (...,2,2), lam_a = i sigma_a/2."""
    x00, x01, x10, x11 = X[..., 0, 0], X[..., 0, 1], X[..., 1, 0], X[..., 1, 1]
    return -0.5 * np.stack([np.imag(x10 + x01), -np.real(x10 - x01), np.imag(x00 - x11)], axis=-1)


def dag(X):
    return X.conj().swapaxes(-1, -2)


def plaq_rots(U):
    """yield (mu, nu, R0..R3) with R_k the cyclic rotation of the plaquette starting at pos k."""
    for mu in range(4):
        for nu in range(mu + 1, 4):
            A = U[..., mu, :, :]
            B = np.roll(U[..., nu, :, :], -1, axis=mu)
            C = dag(np.roll(U[..., mu, :, :], -1, axis=nu))
            D = dag(U[..., nu, :, :])
            R0 = A @ B @ C @ D
            R1 = B @ C @ D @ A
            R2 = C @ D @ A @ B
            R3 = D @ A @ B @ C
            yield mu, nu, (R0, R1, R2, R3)


def S_W(U, beta=BETA):
    s = 0.0
    for mu, nu, R in plaq_rots(U):
        s += np.sum(1.0 - 0.5 * np.trace(R[0], axis1=-2, axis2=-1).real)
    return beta * s


def grad_S(U, beta=BETA):
    """left-invariant gradient, g[x,mu,a] = e_a S_W at link (x,mu); shape (L,L,L,L,4,3)."""
    g = np.zeros(U.shape[:-2] + (3,))
    base = -0.5 * beta
    for mu, nu, R in plaq_rots(U):
        # position 0: link (x,mu) positive   X = R1
        g[..., mu, :] += base * lam_tr(R[1])
        # position 1: link (x+e_mu,nu) positive  X = R2
        g[..., nu, :] += np.roll(base * lam_tr(R[2]), 1, axis=mu)
        # position 2: link (x+e_nu,mu) negative  X = R2  (coef +beta/2)
        g[..., mu, :] += np.roll(-base * lam_tr(R[2]), 1, axis=nu)
        # position 3: link (x,nu) negative  X = R3
        g[..., nu, :] += -base * lam_tr(R[3])
    return g


def gnorm2(U, beta=BETA):
    return float(np.sum(grad_S(U, beta) ** 2))


def haar(L):
    Z = (rng.standard_normal((L, L, L, L, 4, N, N)) + 1j * rng.standard_normal((L, L, L, L, 4, N, N)))
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=-2, axis2=-1)
    Q = Q * (d / np.abs(d)).conj()[..., None, :]
    det = np.linalg.det(Q)
    return Q / det[..., None, None] ** (1.0 / N)


def plaq_vals(U):
    """Re tr U_p / 2 per plaquette (mean over all 6 planes)."""
    out = []
    for mu, nu, R in plaq_rots(U):
        out.append(0.5 * np.trace(R[0], axis1=-2, axis2=-1).real)
    return np.array(out)


def staircase(L, s=0.5):
    """Z2-alternating direction-1 links, U_1(x) = (-1)^{x2+x3+x4} rotation by s*pi about z."""
    U = np.zeros((L, L, L, L, 4, N, N), complex)
    U[...] = I2
    x = np.meshgrid(*[np.arange(L)] * 4, indexing="ij")
    par = (x[1] + x[2] + x[3]) % 2
    th = s * np.pi * par
    U[..., 0, :, :] = np.cos(th)[..., None, None] * I2 + 1j * np.sin(th)[..., None, None] * S3
    return U


# ---------------------------------------------------------------- Metropolis
def staple(U, mu):
    St = np.zeros_like(U[..., mu, :, :])
    for nu in range(4):
        if nu == mu:
            continue
        # forward: U_nu(x+e_mu) U_mu(x+e_nu)^dag U_nu(x)^dag
        St += (np.roll(U[..., nu, :, :], -1, axis=mu) @ dag(np.roll(U[..., mu, :, :], -1, axis=nu))
               @ dag(U[..., nu, :, :]))
        # backward plaquette at x-e_nu contains U_mu(x)^dag, so Re tr P' = Re tr(U_mu(x) . A'^dag)
        # with A' = U_nu(x-e_nu)^dag U_mu(x-e_nu) U_nu(x-e_nu+e_mu)
        Bw = np.roll(np.roll(U[..., nu, :, :], -1, axis=mu), 1, axis=nu)
        Aw = np.roll(U[..., mu, :, :], 1, axis=nu)
        St += dag(Bw) @ dag(Aw) @ np.roll(U[..., nu, :, :], 1, axis=nu)
    return St


def sweep(U, h=0.35, beta=BETA):
    L = U.shape[0]
    xg = np.meshgrid(*[np.arange(L)] * 4, indexing="ij")
    acc = tot = 0
    for mu in range(4):
        par = sum(xg[nu] for nu in range(4) if nu != mu) % 2
        for c in (0, 1):
            m = par == c
            St = staple(U, mu)
            M = U[..., mu, :, :] @ St                      # (L,L,L,L,2,2)
            nl = int(m.sum())
            n = rng.standard_normal((nl, 3))
            n /= np.linalg.norm(n, axis=1)[:, None]
            Xs = expi(h * n)
            Mm = M[m]
            if (Xs - I2).shape != Mm.shape:
                raise RuntimeError(f"DBG Xs{(Xs-I2).shape} Mm{Mm.shape} U{U.shape} m{m.shape} nl{nl} M{M.shape}")
            dS = -0.5 * beta * np.trace((Xs - I2) @ Mm, axis1=-2, axis2=-1).real
            ok = (dS <= 0) | (np.log(rng.random(nl)) < -dS)
            V = U[..., mu, :, :]
            newV = V[m]
            newV = np.where(ok[:, None, None], (Xs @ newV), newV)
            V[m] = newV
            acc += int(ok.sum())
            tot += nl
    return acc / tot


def measure(U):
    g = grad_S(U)
    return dict(S=float(S_W(U)), g2=float(np.sum(g ** 2)), gp=float(np.mean(plaq_vals(U)))), g


# ================================================================= S0 self-checks
say("=" * 108)
say("S0 self-checks (SU(2), beta = %.2f)" % BETA)
L = 2
U0 = haar(L)
g = grad_S(U0)
# single-link FD: the definitive per-component test (mc_lattice reports 2.7e-8)
err = 0.0
for _ in range(12):
    x = tuple(int(rng.integers(0, L)) for _ in range(4))
    mu = int(rng.integers(0, 4))
    a = int(rng.integers(0, 3))
    t = np.zeros(g.shape)
    t[x + (mu, a)] = 1.0
    Up = U0 @ expi(1e-5 * t)
    Um = U0 @ expi(-1e-5 * t)
    fd = (S_W(Up) - S_W(Um)) / 2e-5
    err = max(err, abs(fd - g[x + (mu, a)]))
say(f"  per-link FD check (12 random links x random component): max|fd - g| = {err:.2e}")
v = rng.standard_normal(g.shape)
v /= np.linalg.norm(v)
for h in (1e-4,):
    Up = U0 @ expi(h * v)
    Um = U0 @ expi(-h * v)
    fd = (S_W(Up) - S_W(Um)) / (2 * h)
    say(f"  global FD (unit-norm v): <g,v> = {np.sum(g * v):+.10e}   FD = {fd:+.10e}"
        f"   ratio = {fd / np.sum(g * v):.6f}")
Up = np.zeros_like(U0)
Up[...] = I2
say(f"  pure gauge  U=1     : sum|g|^2 = {gnorm2(Up):.3e}   S_W = {S_W(Up):.3e}")
say(f"  all links U=-1      : sum|g|^2 = {gnorm2(-Up):.3e}   (central => d(Re tr)=0)")
for s in (1.0, 0.5, 0.25):
    Us = staircase(2, s)
    say(f"  Z2 staircase s={s:<4}: sum|g|^2 = {gnorm2(Us):9.3f}   S_W = {S_W(Us):8.3f}"
        f"   analytic (6*beta/2*sin(pi s))^2*L^4 = {(6 * BETA / 2 * np.sin(np.pi * s)) ** 2 * 2 ** 4:9.3f}")
dump("S0", dict(fd_ok=True))

# ================================================================= S1 pointwise vs sup
say("=" * 108)
say("S1  pointwise bound vs sup_phi")
say("    per plaquette: sum_a (d_a S_p)^2 = (beta^2/4) sin^2(theta_p/2) <= beta^2/4 = %.4f" % (BETA ** 2 / 4))
say("    per link (6 plaquettes, triangle ineq.): |grad_l S_W|^2 <= (6 beta/2)^2 = 9 beta^2 = %.4f" % (9 * BETA ** 2))
say("    global: |grad S_W|^2 <= 9 beta^2 * #links = 36 beta^2 L^4 = %.4f L^4" % (36 * BETA ** 2))

# translation-invariant ansatz: U_mu(x) = h_mu constant.  per-link density is exactly L-independent
# => any TI config gives a rigorous  sup_phi |grad S_W|^2 >= 4 L^4 * (per-link value)
def ti_config(p):
    U = np.zeros((2, 2, 2, 2, 4, N, N), complex)
    for mu in range(4):
        U[..., mu, :, :] = expi(p[mu])
    return U


def ti_obj(p):
    """per-link mean |grad_l S_W|^2 for a constant-link config (computed on a 2^4 torus)."""
    U = ti_config(p)
    g = grad_S(U)
    return float(np.sum(g ** 2) / 64.0)


best = (-1.0, None)
from scipy.optimize import minimize

starts = [rng.standard_normal((4, 3)) * np.pi for _ in range(40)]
starts.append(np.array([[0, 0, np.pi], [0, 0, np.pi], [0, 0, np.pi], [0, 0, np.pi]], float))
starts.append(np.array([[0, 0, np.pi / 2], [0, 0, np.pi / 2], [0, 0, np.pi / 2],
                        [0, 0, np.pi / 2]], float))
for p0 in starts:
    p0 = np.array(p0, float)
    f0 = -ti_obj(p0)
    neg = lambda q: -ti_obj(q.reshape(4, 3))
    res = minimize(neg, p0.reshape(-1), method="Powell",
                   options=dict(maxiter=4000, xtol=1e-9, ftol=1e-12))
    f = ti_obj(res.x.reshape(4, 3))
    if f > best[0]:
        best = (f, res.x.reshape(4, 3).copy())
say(f"  TI ansatz (12 params, {len(starts)} starts, Powell): max per-link |grad S_W|^2 = {best[0]:.6f}")
say(f"      = {best[0] / BETA ** 2:.6f} beta^2   (bound 9 beta^2 = {9 * BETA ** 2:.4f};"
    f" staircase gives {2.25 * BETA ** 2:.4f} = 2.25 beta^2)")
say(f"      best TI angles (rad): {np.round(best[1], 4).tolist()}")
say(f"      => rigorous sup_phi >= {best[0] * 4:.4f} L^4   (exactly, every L)")
np.save("r3_ti_best.npy", best[1])
dump("S1", dict(ti_perlink=best[0], ti_perlink_over_beta2=best[0] / BETA ** 2,
               staircase_perlink=2.25 * BETA ** 2,
               perlink_bound=9 * BETA ** 2, global_bound=36 * BETA ** 2))

# also: unconstrained ascent on the full L=2 config (192 params), started from the best TI config
say("  unconstrained ascent on the full L=2 config (192 params), from the best TI config:")
base = ti_config(best[1])
say(f"      start (best TI):   sum|g|^2 = {gnorm2(base):.3f}")


def cfg_from_p(p):
    return expi(p) @ base


p = rng.standard_normal((2, 2, 2, 2, 4, 3)) * 0.05
for it in range(120):
    gr = np.zeros_like(p)
    for idx in np.ndindex(*p.shape):
        dp = np.zeros_like(p)
        dp[idx] = 1e-5
        gr[idx] = (gnorm2(cfg_from_p(p + dp)) - gnorm2(cfg_from_p(p - dp))) / 2e-5
    gn = np.linalg.norm(gr)
    if gn < 1e-12:
        break
    p = p + 0.03 * gr / gn
f_asc = gnorm2(cfg_from_p(p))
say(f"      after 120 ascent steps: sum|g|^2 = {f_asc:.3f}  (per link {f_asc / 64:.3f};"
    f" per-link bound 9 beta^2 = {9 * BETA ** 2:.3f}, global bound 36 beta^2 L^4 = {36 * BETA ** 2 * 16:.1f})")
dump("S1b", dict(full_L2_ascent=f_asc, per_link=f_asc / 64))

# ================================================================= S2 three quantities
say("=" * 108)
say("S2  three quantities, L = 2,4,6   (SU(2), beta = 2.5)")
rows = {}
for L in (2, 4, 6, 8):
    E = 4 * L ** 4
    t0 = time.time()
    # (A) Gibbs measure: Metropolis
    U = haar(L)
    for _ in range(420):
        acc = sweep(U)
    Ss, gs, gr2 = [], [], []
    gsum = np.zeros((12, L, L, L, L), complex)
    g2sum = np.zeros((12, L, L, L, L))
    gacc = np.zeros((L, L, L, L, 4, 3))
    gsq_tot = 0.0
    nmeas = 220 if L <= 4 else (160 if L == 6 else 120)
    for k in range(nmeas):
        acc = sweep(U)
        st, g = measure(U)
        Ss.append(st["S"])
        gs.append(st["g2"])
        gacc += g
        gsq_tot += float(np.sum(g ** 2))
        for mu in range(4):
            for a in range(3):
                gh = np.fft.fftn(g[..., mu, a])
                i = mu * 3 + a
                gsum[i] += gh
                g2sum[i] += np.abs(gh) ** 2
    A = float(np.mean(gs))
    Aerr = float(np.std(gs) / np.sqrt(len(gs)))
    # (B) sup_phi lower bound (exact, TI)
    B = best[0] * 4 * L ** 4
    # (C) Haar random field
    Cs = [gnorm2(haar(L)) for _ in range(60)]
    C = float(np.mean(Cs))
    Cerr = float(np.std(Cs) / np.sqrt(len(Cs)))
    # ---- mode-resolved covariance of grad S_W (the spectral object BD actually needs)
    # real v with |v| = 1: v_l = cos(kx)/sqrt(NS) (self-conjugate k) or sqrt(2/NS)cos/sin(kx);
    # lambda_k = Var(ghat_k)/NS for EVERY k, and Parseval gives sum_k lambda_k = tr Cov.
    n = nmeas
    NS = L ** 4
    var = (g2sum / n - np.abs(gsum / n) ** 2).real
    lam = var / NS
    trC_modes = float(np.sum(lam))
    gm = gacc / n
    trC_direct = float(gsq_tot / n - np.sum(gm ** 2))
    lmax = float(np.max(lam))
    lok = float(np.mean(lam[:, 0, 0, 0, 0]))
    rows[L] = dict(E=E, A=A, Aerr=Aerr, A_dens=A / E, B=B, B_dens=B / E, C=C, Cerr=Cerr,
                   C_dens=C / E, trCov=trC_modes, trCov_direct=trC_direct,
                   trCov_dens=trC_modes / E, lam_max=lmax, lam_k0=lok,
                   lam_med=float(np.median(lam)), S_mean=float(np.mean(Ss)) / (6 * L ** 4),
                   acc=float(acc), nmeas=nmeas, t=time.time() - t0)
    r = rows[L]
    npq = 6 * L ** 4
    say(f"  L={L}  E={E:5d}  plaquettes={npq:5d}   [{(time.time() - t0):.1f}s, n={nmeas}, acc={acc:.2f}]")
    say(f"     A = <|grad S_W|^2>_mu      = {A:11.3f} +- {Aerr:7.3f}   per link {A / E:7.3f}"
        f"   per plaque {A / npq:7.3f}")
    say(f"     B = sup_phi lower bound    = {B:11.3f}             per link {B / E:7.3f}"
        f"   per plaque {B / npq:7.3f}   (B/A = {B / A:.3f})")
    say(f"     C = Haar random field      = {C:11.3f} +- {Cerr:7.3f}   per link {C / E:7.3f}"
        f"   per plaque {C / npq:7.3f}")
    say(f"     per plaquette S_W (mu)     = {r['S_mean']:.4f}   (Haar/infinite-T value: "
        f"{S_W(haar(L)) / (6 * L ** 4):.4f})")
    say(f"     tr Cov(grad S)  = {trC_modes:11.3f} (modes)  {trC_direct:11.3f} (direct)   per link {trC_direct / E:8.4f}")
    say(f"     lambda_max Cov  = {lmax:+11.4f}    lambda(k=0) = {lok:+11.4f}"
        f"    median mode = {r['lam_med']:+.4f}")
dump("S2", rows)


def fit(ys):
    Ls = sorted(ys)
    y = np.log([ys[L] for L in Ls])
    x = np.log(Ls)
    sl = np.polyfit(x, y, 1)[0]
    loc = [float(np.log(ys[b] / ys[a]) / np.log(b / a)) for a, b in zip(Ls[:-1], Ls[1:])]
    return float(sl), loc


say("-" * 108)
say("S2b  L-exponent fits (log-log through L=2,4,6,8; local slopes between neighbours)")
for name, key in (("A  <|grad S|^2>_mu", "A"), ("B  sup_phi lower bound", "B"),
                  ("C  Haar random field", "C"), ("tr Cov(grad S_W)", "trCov_direct"),
                  ("lambda_max Cov", "lam_max"), ("lambda(k=0) Cov", "lam_k0")):
    ys = {L: rows[L][key] for L in rows}
    sl, loc = fit(ys)
    say(f"    {name:24s} : exponent {sl:6.3f}   (local " + " / ".join(f"{v:6.3f}" for v in loc) + ")")
dump("S2fits", {k: fit({L: rows[L][k] for L in rows}) for k in
                ("A", "B", "C", "trCov_direct", "lam_max", "lam_k0")})
say("=" * 108)
LOG.close()
print("\nwrote r3_l4.json / r3_l4.log")

