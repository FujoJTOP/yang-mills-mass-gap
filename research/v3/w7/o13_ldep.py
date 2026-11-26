#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
o13_ldep.py -- Delta(L) for 4D SU(2) Wilson at beta=2.5, constrained (xi=0) ensemble.

Object (H1C-WFULL.md / GPP-STRUCT.md notation, flat chart X_ell = log U_ell in R^3,
dU = J dX, Phi = -log J = -2 log(sin r / r)):
    V = S_W + Phi,        G = d_v S_W,   Phi' = d_v Phi,   H_sv = v^T Hess_{x} S_W v
    Delta_dyn = E[H_sv | xi=0] - Var(G | xi=0)
    W_geo     = E[Phi''_vv | xi=0] - Var(Phi' | xi=0)
    X_cov     = -2 Cov(G, Phi' | xi=0)
    W_full    = W_geo + X_cov
    G''(0)    = (E[H_sv+Phi''] - Var(G+Phi'))|xi=0 = Delta_dyn + W_full      [eq:split]
    identity (I):  W_full = E[Phi''] + Var(G) - Var(G+Phi')                    [H1C eq (I)]

Mode: v = normalised constant (k=0) mode, link direction mu=1, Lie direction sigma_1:
      v[mu=1](x) = e_1 / L^2  (all x),  v = 0 elsewhere  => ||v|| = 1, v perp. pure gauge.

Exact constrained sampling: the constraint xi(X) = <v,X> = sum_ell v_ell . x_ell is AFFINE,
so a Hamiltonian flow with force and momentum projected by P = 1 - v v^T stays on the
hyperplane {xi=0} exactly -> ordinary HMC with rank-1 projected forces, exact, no Jacobian.

Validations printed in the log: quaternion S_W vs explicit 2x2 matrices; analytic force vs
finite difference; beta=0 vs exact 3D quadrature (Haar); FD step-size convergence;
symmetry checks E[G]=E[Phi']=0; unconstrained HMC vs independent single-link Metropolis.

Run:  py -3.12 o13_ldep.py L [ntraj] [nmeas_every] [nsteps] [eps]
Labels [V].  No git commit.
"""
import sys, time, json
import numpy as np

np.seterr(all='ignore')

# ------------------------------------------------------------------ quaternions
# U = q0 * I + i q.sigma ; tr(U)/2 = q0 ; q = (cos r, sin(r)/r * x), x = log U
def qmul(p, q):
    p0 = p[..., 0:1]; pv = p[..., 1:]
    q0 = q[..., 0:1]; qv = q[..., 1:]
    return np.concatenate([p0 * q0 - (pv * qv).sum(-1, keepdims=True),
                           p0 * qv + q0 * pv - np.cross(pv, qv)], -1)

def qconj(p):
    o = p.copy(); o[..., 1:] *= -1.0; return o

def x_to_q(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    q0 = np.where(sm, 1.0 - 0.5 * r * r, np.cos(rs))
    f = np.where(sm, 1.0 - r * r / 6.0, np.sin(rs) / rs)
    return np.concatenate([q0, f * x], -1)

def sincos_series(r, sm, rs):
    f = np.where(sm, 1.0 - r * r / 6.0, np.sin(rs) / rs)                  # sin r / r
    fp = np.where(sm, -r / 3.0, (rs * np.cos(rs) - np.sin(rs)) / rs ** 2)  # (sin r/r)'
    sr = np.where(sm, r, np.sin(rs))
    cr = np.where(sm, 1.0 - 0.5 * r * r, np.cos(rs))
    return f, fp, sr, cr

# ------------------------------------------------------------------ geometry
def build_geom(L):
    n = L ** 4
    idx = np.arange(n)
    c = np.stack([idx % L, (idx // L) % L, (idx // L ** 2) % L, (idx // L ** 3) % L], 1)
    def shp(mu):
        cc = c.copy(); cc[:, mu] = (cc[:, mu] + 1) % L
        return cc[:, 0] + L * cc[:, 1] + L ** 2 * cc[:, 2] + L ** 3 * cc[:, 3]
    def shm(mu):
        cc = c.copy(); cc[:, mu] = (cc[:, mu] - 1) % L
        return cc[:, 0] + L * cc[:, 1] + L ** 2 * cc[:, 2] + L ** 3 * cc[:, 3]
    return dict(L=L, n=n, idx=idx, sp=[shp(m) for m in range(4)],
                sm=[shm(m) for m in range(4)],
                planes=[(mu, nu) for mu in range(4) for nu in range(mu + 1, 4)])

# ------------------------------------------------------------------ S_W and dS_W/dx
def force_corner(xx, m, sg, beta):
    """dS_p/dx_b for one plaquette corner, with P0 = A(r) m0 + B(r)(x.m)."""
    r = np.linalg.norm(xx, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, fp, sr, _ = sincos_series(r, sm, rs)
    mv = m[..., 1:]; m0 = m[..., 0:1]
    xd = (xx * mv).sum(-1, keepdims=True)
    return beta * (sg * f * mv + (xx / rs) * (m0 * sr + sg * fp * xd))

def sw_and_force(x, g, beta, want_force=True):
    idx, sp, n = g['idx'], g['sp'], g['n']
    q = x_to_q(x)
    tot = 0.0
    accs = []
    for mu, nu in g['planes']:
        A = q[mu]; B = q[nu][sp[mu]]
        C = qconj(q[mu][sp[nu]]); D = qconj(q[nu])
        ABC = qmul(qmul(A, B), C)
        P = qmul(ABC, D)
        tot += beta * np.sum(1.0 - P[:, 0])
        if not want_force:
            continue
        m1 = qmul(qmul(B, C), D)
        m2 = qmul(qmul(C, D), A)
        m3 = qmul(qmul(D, A), B)
        for lm, site, sg, mm in ((mu, idx, +1, m1), (nu, sp[mu], +1, m2),
                                 (mu, sp[nu], -1, m3), (nu, idx, -1, ABC)):
            accs.append((lm * n + site, force_corner(x[lm][site], mm, sg, beta)))
    F = np.zeros_like(x)
    if accs:
        fi = np.concatenate([a[0] for a in accs])
        vals = np.concatenate([a[1] for a in accs])
        Ff = F.reshape(4 * n, 3)
        for b in range(3):
            Ff[:, b] += np.bincount(fi, weights=vals[:, b], minlength=4 * n)
    return tot, F

# ------------------------------------------------------------------ Phi
def phi_of(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, _, _, _ = sincos_series(r, sm, rs)
    return -2.0 * np.log(f)

def grad_phi(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, fp, _, _ = sincos_series(r, sm, rs)
    ur = np.where(sm, -1.0 / 3.0, fp / (f * rs))     # (f'/f)/r  -> -1/3 at 0
    return -2.0 * ur * x

def hess_phi_vv(x, w):
    """w^T Hess phi(x) w  (phi = -2 log(sin r/r)); exact, series-safe at r=0."""
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, fp, sr, cr = sincos_series(r, sm, rs)
    fpp = np.where(sm, -1.0 / 3.0 + r * r / 30.0,
                   ((2.0 - rs ** 2) * sr - 2.0 * rs * cr) / rs ** 3)
    u = np.where(sm, 0.0, fp / f)
    ur = np.where(sm, -1.0 / 3.0, fp / (f * rs))
    vv = fpp / f - u * u                       # -> -1/3 at 0
    xh = x / rs
    wpar = (w * xh).sum(-1, keepdims=True)
    w2 = (w * w).sum(-1, keepdims=True)
    return -2.0 * (vv * wpar ** 2 + ur * (w2 - wpar ** 2))

# ------------------------------------------------------------------ mode v
def make_v(g):
    L = g['L']; n = g['n']
    v = np.zeros((4, n, 3))
    v[0][:, 0] = 1.0 / L ** 2          # <v,v> = L^4 * L^-4 = 1
    return v

def dotv(a, v):                          # <v, a> with ||v||=1
    return float((v * a).sum())

# ------------------------------------------------------------------ HMC
def hmc(L, beta, ntraj, nsteps, eps, v, g, seed, therm=0, meas_every=1,
        h=0.05, record=None, verbose=True, constrain=True):
    rng = np.random.default_rng(seed)
    n = g['n']
    x = rng.standard_normal((4, n, 3)) * 0.30
    if constrain:
        x -= v * dotv(x, v)
    V = lambda z: sw_and_force(z, g, beta, False)[0] + phi_of(z).sum()
    def force(z):
        f = sw_and_force(z, g, beta, True)[1] + grad_phi(z)
        if constrain:
            f = f - v * dotv(f, v)
        return f
    acc = 0; natt = 0
    t0 = time.time()
    for it in range(ntraj):
        p = rng.standard_normal((4, n, 3))
        if constrain:
            p -= v * dotv(p, v)
        x0 = x.copy()
        V0 = V(x0); K0 = 0.5 * (p * p).sum()
        f = force(x0)
        p = p + 0.5 * eps * (-f)
        x = x0 + eps * p
        for k in range(nsteps - 1):
            f = force(x)
            p = p - eps * f
            x = x + eps * p
        f = force(x)
        p = p - 0.5 * eps * f
        V1 = V(x); K1 = 0.5 * (p * p).sum()
        dH = (V1 + K1) - (V0 + K0)
        natt += 1
        if np.log(rng.random()) < -dH:
            acc += 1
        else:
            x = x0
        if it >= therm and record is not None and (it - therm) % meas_every == 0:
            record.append(measure(x, g, beta, v, h))
    if verbose:
        print("  HMC L=%d beta=%.2f: %d traj x %d steps, eps=%.4f, acc=%.3f, %.1fs"
              % (L, beta, natt, nsteps, eps, acc / max(natt, 1), time.time() - t0), flush=True)
    return x, acc / max(natt, 1)

# ------------------------------------------------------------------ measurement
def measure(x, g, beta, v, h):
    S0, F = sw_and_force(x, g, beta, True)
    Gv = dotv(F, v)
    gp = grad_phi(x)
    Ph = dotv(gp, v)
    Pip = float(hess_phi_vv(x, v).sum())
    xp = x + h * v; xm = x - h * v
    Sp = sw_and_force(xp, g, beta, False)[0]
    Sm = sw_and_force(xm, g, beta, False)[0]
    Hsv = (Sp - 2.0 * S0 + Sm) / h ** 2
    return (Gv, Ph, Pip, Hsv, S0)

# ------------------------------------------------------------------ stats
def binned(vals, nb=40):
    a = np.asarray(vals, float)
    n = a.size
    nb = min(nb, max(4, n // 8))
    k = n // nb
    bm = a[:k * nb].reshape(nb, k).mean(1)
    mean = a.mean()
    sig = bm.std(ddof=1) / np.sqrt(nb)
    return mean, sig, n, k

# ------------------------------------------------------------------ validations
def val_sw_matrices(L, beta=2.5, seed=1):
    """quaternion S_W vs explicit 2x2 complex matrices."""
    g = build_geom(L); rng = np.random.default_rng(seed)
    x = rng.standard_normal((4, g['n'], 3)) * 0.7
    q = x_to_q(x)
    S, _ = sw_and_force(x, g, beta, False)
    # explicit matrices
    U = np.zeros((4, g['n'], 2, 2), complex)
    r = np.linalg.norm(x, axis=-1); f = np.where(r < 1e-12, 1.0, np.sin(r) / np.where(r < 1e-12, 1, r))
    nrm = np.where(r < 1e-12, 1.0, r)
    for mu in range(4):
        nx = x[mu] / nrm[mu][..., None]
        c, s = np.cos(r[mu]), np.sin(r[mu])
        s = np.sin(r[mu])
        U[mu, :, 0, 0] = c + 1j * s * nx[:, 2]; U[mu, :, 0, 1] = 1j * s * nx[:, 0] + s * nx[:, 1]
        U[mu, :, 1, 0] = 1j * s * nx[:, 0] - s * nx[:, 1]; U[mu, :, 1, 1] = c - 1j * s * nx[:, 2]
    idx, sp = g['idx'], g['sp']
    Sref = 0.0
    for mu, nu in g['planes']:
        A = U[mu][idx]; B = U[nu][sp[mu]]
        C = np.conj(U[mu][sp[nu]].transpose(0, 2, 1)); D = np.conj(U[nu].transpose(0, 2, 1))
        P = A @ B @ C @ D
        Sref += beta * np.sum(1.0 - 0.5 * np.real(P[:, 0, 0] + P[:, 1, 1]))
    return abs(S - Sref) / abs(Sref)

def val_force_fd(L, beta=2.5, seed=2):
    """analytic dS_W/dx vs central finite difference."""
    g = build_geom(L); rng = np.random.default_rng(seed)
    x = rng.standard_normal((4, g['n'], 3)) * 0.5
    _, F = sw_and_force(x, g, beta, True)
    e = rng.standard_normal((4, g['n'], 3))
    e /= np.linalg.norm(e)
    hh = 1e-5
    Sp = sw_and_force(x + hh * e, g, beta, False)[0]
    Sm = sw_and_force(x - hh * e, g, beta, False)[0]
    return abs((Sp - Sm) / (2 * hh) - (F * e).sum()) / abs((F * e).sum())

def val_beta0_phi(L, seed=3):
    """beta=0 (pure product Haar): measured Phi moments vs exact radial quadrature.
    Sum over links of <v,Hess v> = single-link <e1^T Hess e1>; E[Phi'^2] similarly."""
    g = build_geom(L)
    rr = np.linspace(1e-7, np.pi - 1e-7, 400001)
    w = np.sin(rr) ** 2; w /= w.sum()
    f = np.sin(rr) / rr
    fp = (rr * np.cos(rr) - np.sin(rr)) / rr ** 2
    fpp = ((2 - rr ** 2) * np.sin(rr) - 2 * rr * np.cos(rr)) / rr ** 3
    u = fp / f; ur = u / rr; vv = fpp / f - u * u
    ex = float(np.sum(w * (-2.0) * (vv / 3.0 + ur * (2.0 / 3.0))))
    ex_var = float(np.sum(w * 4.0 * u * u) / 3.0)
    v = make_v(g); rng = np.random.default_rng(seed)
    x = rng.standard_normal((4, g['n'], 3)) * 0.3
    eps, nst, N = 0.15, 8, 2500
    rec = []; acc = 0
    for it in range(N):
        p = rng.standard_normal((4, g['n'], 3))
        x0 = x.copy(); V0 = phi_of(x0).sum(); K0 = .5 * (p * p).sum()
        f0 = grad_phi(x0); p = p - .5 * eps * f0; x = x0 + eps * p
        for k in range(nst - 1):
            p = p - eps * grad_phi(x); x = x + eps * p
        p = p - .5 * eps * grad_phi(x)
        V1 = phi_of(x).sum(); K1 = .5 * (p * p).sum()
        if np.log(rng.random()) < (V0 + K0) - (V1 + K1):
            acc += 1
        else:
            x = x0
        if it > 200 and it % 5 == 0:
            rec.append((float(hess_phi_vv(x, v).sum()), dotv(grad_phi(x), v)))
    Pp = np.array([r[0] for r in rec]); Qp = np.array([r[1] for r in rec])
    return dict(ex_phipp=ex, mc_phipp=float(Pp.mean()),
                mc_phipp_sig=float(Pp.std(ddof=1) / np.sqrt(Pp.size)),
                ex_var=ex_var, mc_var=float((Qp ** 2).mean()),
                mc_var_sig=float((Qp ** 2).std(ddof=1) / np.sqrt(Qp.size)),
                acc=acc / N)

def val_metropolis(L, beta, nsweep, seed=7, step=0.30):
    """independent single-link Metropolis (unconstrained) -> <1-tr U_p/2>."""
    g = build_geom(L); rng = np.random.default_rng(seed)
    x = rng.standard_normal((4, g['n'], 3)) * 0.3
    S = sw_and_force(x, g, beta, False)[0]
    acc = 0; nprop = 0; ps = []
    for sw in range(nsweep):
        for mu in range(4):
            for s in range(g['n']):
                old = x[mu][s].copy()
                prop = old + step * rng.standard_normal(3)
                if np.linalg.norm(prop) > 0.90 * np.pi:
                    continue
                x[mu][s] = prop
                S2 = sw_and_force(x, g, beta, False)[0]
                nprop += 1
                if np.log(rng.random()) < -(S2 - S):
                    S = S2; acc += 1
                else:
                    x[mu][s] = old
        if sw >= nsweep // 4:
            ps.append(1.0 - S / (beta * 6.0 * g['n']))
    return float(np.mean(ps)), acc / max(nprop, 1), nprop

# ------------------------------------------------------------------ derived
def derived(rec):
    a = np.array(rec, float)
    return dict(G=a[:, 0], P=a[:, 1], Pip=a[:, 2], Hsv=a[:, 3], S0=a[:, 4])

def summarize(rec, beta, g, nboot=400, seed=11):
    d = derived(rec)
    Gv, Ph, Pip, Hsv, S0 = d['G'], d['P'], d['Pip'], d['Hsv'], d['S0']
    cols = np.stack([Hsv, Gv ** 2, Pip, Ph ** 2, Gv * Ph, Gv, Ph], 1)   # per-sample
    n = cols.shape[0]; nb = min(40, max(4, n // 8)); k = max(1, n // nb); nb = n // k
    B = cols[:k * nb].reshape(nb, k, 7).mean(1)                          # bin means
    def stat(bb):
        A = bb[:, 0].mean(); Bv = bb[:, 1].mean(); C = bb[:, 2].mean()
        D = bb[:, 3].mean(); E = bb[:, 4].mean()
        d_dyn = A - Bv; w_geo = C - D; xc = -2 * E
        w_full = w_geo + xc
        gpp = (A + C) - (Bv + D + 2 * E)
        vsum = Bv + D + 2 * E
        return np.array([A, Bv, C, D, E, d_dyn, w_geo, xc, w_full, gpp, vsum,
                         bb[:, 5].mean(), bb[:, 6].mean()])
    full = stat(cols)
    rng = np.random.default_rng(seed)
    boots = np.array([stat(B[rng.integers(0, nb, nb)]) for _ in range(nboot)])
    sig = boots.std(0, ddof=1)
    names = ['E[Hsv|0]', 'Var(G|0)', 'E[Phi..|0]', 'Var(Phi\'|0)', 'Cov(G,Phi\'|0)',
             'Delta_dyn', 'W_geo', 'X_cov', 'W_full', 'Gpp(0)', 'Var(G+Phi\')',
             'E[G|0]', 'E[Phi\'|0]']
    out = {nm: (float(v), float(s)) for nm, v, s in zip(names, full, sig)}
    # tau estimate: binned sigma of Delta_dyn vs naive sigma of E[Hsv]
    naiveA = np.sqrt(cols[:, 0].var(ddof=1) / n)
    out['_tau_dyn'] = float((sig[5] / (naiveA + 1e-30)) ** 2)
    out['_n'] = int(n); out['_nb'] = int(nb); out['_k'] = int(k)
    out['_plaq'] = float(1.0 - S0.mean() / (beta * 6.0 * g['n']))
    out['_Pip_sd'] = float(Pip.std()); out['_Hsv_sd'] = float(Hsv.std())
    out['_G_sd'] = float(Gv.std())
    return out

def show(L, o, tag=''):
    print("  ---- %s L=%d : moments (mean +- boot sigma) ----" % (tag, L), flush=True)
    for k in ['E[Hsv|0]', 'Var(G|0)', "E[Phi..|0]", "Var(Phi'|0)", "Cov(G,Phi'|0)",
              'Delta_dyn', 'W_geo', 'X_cov', 'W_full', 'Gpp(0)', "Var(G+Phi')",
              "E[G|0]", "E[Phi'|0]"]:
        v, s = o[k]
        print("    %-14s = %+0.6f  +- %0.6f" % (k, v, s), flush=True)
    print("    identity (I) residual |W_full - (E[Phi..]+VarG-Var(G+Phi'))| = %.2e"
          % abs(o['W_full'][0] - (o["E[Phi..|0]"][0] + o['Var(G|0)'][0] - o["Var(G+Phi')"][0])), flush=True)
    print("    <1-trU_p/2> = %.4f   sd(Hsv)=%.4f sd(G)=%.4f   n=%d bins=%d k=%d  tau_dyn~%.1f"
          % (o['_plaq'], o['_Hsv_sd'], o['_G_sd'], o['_n'], o['_nb'], o['_k'], o['_tau_dyn']), flush=True)

# ------------------------------------------------------------------ main
def main():
    L = int(sys.argv[1])
    ntraj = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    meas_every = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    nsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    eps0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.12
    beta = 2.5
    print("=" * 78, flush=True)
    print("o13_ldep.py  SU(2) beta=%.2f  L=%d   (%d^4 lattice, %d links, %d plaquettes)"
          % (beta, L, L, 4 * L ** 4, 6 * L ** 4), flush=True)
    print("mode v: k=0 constant, direction mu=1, Lie sigma_1, v=e1/L^2, ||v||=1", flush=True)
    print("=" * 78, flush=True)
    g = build_geom(L); v = make_v(g)
    t00 = time.time()
    # ---------------- validations
    print("[val 1] S_W quaternion vs explicit 2x2 matrices : rel.dev = %.2e"
          % val_sw_matrices(L), flush=True)
    print("[val 2] analytic dS_W/dx vs central FD            : rel.dev = %.2e"
          % val_force_fd(L), flush=True)
    v0 = val_beta0_phi(L)
    print("[val 3] beta=0 Haar: E[Phi''_vv] ex=%.6f mc=%.6f+-%.6f | E[Phi'^2] ex=%.6f mc=%.6f+-%.6f | acc=%.3f"
          % (v0['ex_phipp'], v0['mc_phipp'], v0['mc_phipp_sig'], v0['ex_var'],
             v0['mc_var'], v0['mc_var_sig'], v0['acc']), flush=True)
    # ---------------- eps tuning
    eps = eps0; acc = 0.0
    for it in range(6):
        _, acc = hmc(L, beta, 12, nsteps, eps, v, g, seed=100 + it, therm=4,
                     record=None, verbose=False)
        print("  tune eps=%.4f -> acc=%.3f" % (eps, acc), flush=True)
        if 0.6 < acc < 0.9:
            break
        eps *= np.exp(np.clip((acc - 0.78) / 3.0, -0.6, 0.6))
    # ---------------- production
    rec = []
    x, acc = hmc(L, beta, ntraj, nsteps, eps, v, g, seed=2024, therm=max(60, ntraj // 5),
                 meas_every=meas_every, h=0.05, record=rec)
    o = summarize(rec, beta, g)
    show(L, o, 'constrained xi=0')
    # ---------------- FD step-size convergence (subset)
    hh = []
    for h in (0.05, 0.10, 0.20):
        sub = [measure(x, g, beta, v, h) for _ in range(1)]
        hh.append(sub[0][3])
    # reuse the first 400 stored configs? (not stored) -> report on one config only
    print("  FD-h check on one config: Hsv(h=0.05,0.10,0.20) = %.6f %.6f %.6f"
          % tuple(hh), flush=True)
    # ---------------- unconstrained HMC plaquette (cross-check of the S_W force)
    rec2 = []
    hmc(L, beta, min(ntraj, 3000), nsteps, eps, v, g, seed=777, therm=300,
        meas_every=4, h=0.05, record=rec2, constrain=False)
    o2 = summarize(rec2, beta, g)
    print("  unconstrained HMC  <1-trU_p/2> = %.4f  (constrained: %.4f)"
          % (o2['_plaq'], o['_plaq']), flush=True)
    if L <= 0:
        pm, am, np_ = val_metropolis(L, beta, nsweep=200)
        print("  single-link Metropolis (indep.) <1-trU_p/2> = %.4f +- ?  acc=%.2f  (%d prop)"
              % (pm, am, np_), flush=True)
    # ---------------- summary json
    keys = ['E[Hsv|0]', 'Var(G|0)', "E[Phi..|0]", "Var(Phi'|0)", "Cov(G,Phi'|0)",
            'Delta_dyn', 'W_geo', 'X_cov', 'W_full', 'Gpp(0)', "Var(G+Phi')"]
    js = dict(L=L, beta=beta, ntraj=ntraj, meas_every=meas_every, nsteps=nsteps, eps=eps,
              acc=float(acc), nmeas=o['_n'], plaq=o['_plaq'], tau_dyn=o['_tau_dyn'],
              moments={k: [o[k][0], o[k][1]] for k in keys},
              unconstrained_plaq=o2['_plaq'], wall_s=time.time() - t00)
    print("JSON " + json.dumps(js), flush=True)
    print("total wall %.1f s" % (time.time() - t00), flush=True)

if __name__ == '__main__':
    main()
