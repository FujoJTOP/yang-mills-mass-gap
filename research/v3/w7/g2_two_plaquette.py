#!/usr/bin/env python3
"""
g2_two_plaquette.py -- two-plaquette / two-free-link extension of the exact S2 one-plaquette result.

MODEL (SU(2), 1x2 strip, TWO free links, fixed boundary)
  Two plaquettes sharing the vertical link V; a horizontal link A is also free.
  Plaquette holonomies (boundary links chosen so the fixed products are the identity and B):
        P1 = A V ,        P2 = B V^{-1}      (A,V in SU(2) free; B in SU(2) fixed boundary)
  Wilson action   S = beta(1 - 1/2 tr P1) + beta(1 - 1/2 tr P2)
                      = beta( 2 - [1/2 tr(AV)] - [1/2 tr(B V^{-1})] ).
  Conditional measure d mu = Z^{-1} e^{-S} dA dV  (Haar).  This is mu(.|U') with fixed boundary U'.
  In variables (W=AV, V) the measure FACTORIZES: W ~ one-plaquette(mu, boundary I), V ~ one-plaquette(mu, boundary B),
  W,V independent  (Haar invariance dA = dW).  [checked numerically below]

FLAT CHART:  U = exp(xi), xi in su(2) ~ R^3, Haar dU = J(rho) d^3 xi,
  Phi := -log J = 2 log(rho/2) - 2 log sin(rho/2),   J(rho) = sin^2(rho/2)/rho^2  (prop.).
  For each free link xi_A, xi_V in R^3 (six flat coordinates total).

OBJECTS (campaign definitions, flat-chart Cartesian derivative along a direction field v):
  G_v = d_v S ,  H_vv = d_v^2 S ,  Phi'_v = d_v Phi ,  Phi''_vv = d_v^2 Phi ,
  W_geo(v) := E_mu[Phi''_vv] - Var_mu(Phi'_v) ,     G''(v) := E_mu[H_vv] - Var_mu(G_v) + W_geo(v) ,
  remainder  R(v) := G''(v) - 2 W_geo(v) = E_mu[H_vv] - Var_mu(G_v) - W_geo(v) .

DIRECTIONS
 (L) LINEAR / "linear P": v = (v_A, v_V) CONSTANT in the flat chart, i.e. xi(eps) = xi + eps v.
     This is <P^T grad S, v> with P linear (a constant linear combination of the free links).
 (C) COVARIANT: right-invariant (gauge-covariant) fluctuation, U(eps) = U exp(eps w),
     w a fixed algebra vector.  In the flat chart this is the FIELD-DEPENDENT direction field
     v(xi) = d/deps log( exp(xi) exp(eps w) )|_0,  with flat divergence  div v = -Phi'_v
     (Haar is right-invariant => the right-invariant generator is Haar-divergence-free).

THEOREM (proved in G2-TWO-PLAQUETTE.md, verified here):
  (L) For EVERY constant linear direction v and every fixed boundary B, beta>0:
        R = 0  EXACTLY  (generalized flat-chart Stein identity; e^{-S}J(.) vanishes on the chart edge).
      So G''=2W_geo, independent of block size / number of free links.  Residual here is pure MC noise.
  (C) For a field-dependent (covariant) direction the divergence correction breaks the clean identity;
      R is generally nonzero.

ACID TEST: beta=0 (pure Haar, S=0 => G=H=0 => G''=W_geo) and the flat/free group must give G''=0.

Run:  py -3.12 g2_two_plaquette.py
"""
import numpy as np

# ------------------------------------------------------------------ SU(2) unit quaternions
# U = q0 I + i(q1 s1 + q2 s2 + q3 s3),  |q|=1,  1/2 tr U = q0.
def quat_mul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0*b0 - a1*b1 - a2*b2 - a3*b3,
                     a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1,
                     a0*b3 + a1*b2 - a2*b1 + a3*b0], axis=-1)
def quat_conj(q):
    return q * np.array([1., -1., -1., -1.])

def flat_exp(t):                       # xi in R^3 (...,3) -> quaternion (...,4)
    r = np.linalg.norm(t, axis=-1); small = r < 1e-9; rs = np.where(small, 1.0, r)
    s = np.sin(r/2.0)/rs; q = np.empty(t.shape[:-1] + (4,))
    q[..., 0] = np.where(small, 1.0, np.cos(r/2.0)); q[..., 1:] = t * s[..., None]
    return q
def flat_log(q):                       # quaternion (...,4) -> (xi in R^3, rho in [0,2pi])
    q0 = np.clip(q[..., 0], -1, 1); r = 2.0*np.arccos(q0); s = np.linalg.norm(q[..., 1:], axis=-1)
    ss = np.where(s < 1e-12, 1.0, s); inv = np.where(s < 1e-12, 0.0, r/ss)
    return q[..., 1:]*inv[..., None], r

# ---- radial Phi and derivatives (cancellation-free, x = rho/2) ----
def Phi_rho(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 1.0/xr - 1.0/np.tan(xr)
    xs = x[sm];  y[sm] = xs/3 + xs**3/45 + 2*xs**5/945 + 2*xs**7/4725
    return y
def Phi_rr(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 0.5*(-1.0/xr**2 + 1.0/np.sin(xr)**2)
    xs = x[sm];  y[sm] = 1.0/6 + xs**2/30 + 2*xs**4/189 + 2*xs**6/2475
    return y
def phi_scalar(rho):                   # Phi(rho) = 2log(rho/2) - 2log sin(rho/2); Phi(0)=0
    x = rho/2.0; small = rho < 1e-3
    out = np.empty_like(rho)
    xr = x[~small]; out[~small] = 2.0*np.log(xr) - 2.0*np.log(np.sin(xr))
    xs = x[small];  out[small] = xs**2/3 + xs**4/90 + xs**6/2835
    return out
def Jrad(rho):
    x = rho/2.0; small = rho < 1e-9
    return np.where(small, 0.25, (np.sin(x)/np.where(small, 1.0, rho))**2)

def phi_grad_hess(t, v):               # analytic, constant direction v
    r = np.linalg.norm(t, axis=-1); x = r/2.0
    pr = Phi_rho(x); prr = Phi_rr(x); small = r < 1e-7; rs = np.where(small, 1.0, r)
    nhat = t/rs[..., None]; nv = np.sum(nhat*v, axis=-1); vv = np.sum(v*v, axis=-1)
    p = pr*nv
    por = np.where(small, 1.0/6.0, pr/rs)
    pp = prr*nv**2 + por*(vv - nv**2)
    pp = np.where(small, (1.0/6.0)*vv, pp)
    return p, pp

# ------------------------------------------------------------------ sampling
def samp_unif_ball(n, seed):           # uniform Lebesgue in ball radius 2pi (per link)
    rng = np.random.default_rng(seed)
    u = rng.random(n); rho = 2.0*np.pi*u**(1.0/3.0)
    ct = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2*np.pi, n); st = np.sqrt(np.clip(1-ct*ct, 0, 1))
    t = np.stack([rho*st*np.cos(ph), rho*st*np.sin(ph), rho*ct], axis=1)
    return t, rho

def boundary_quat(theta, n):
    return np.broadcast_to(np.array([np.cos(theta/2), 0., 0., np.sin(theta/2)]), (n, 4))

# ------------------------------------------------------------------ expectation core
def _agg(w, G, H, pp, ppp):
    Z = w.sum()
    def E(a): return float((w*a).sum()/Z)
    EH = E(H); EG = E(G); VG = E(G*G) - EG**2
    Ep = E(pp); Vp = E(pp*pp) - Ep**2; EPP = E(ppp)
    W = EPP - Vp; Gpp = EH - VG + W
    return dict(EH=EH, VG=VG, W=W, Gpp=Gpp, R=Gpp - 2*W, EPP=EPP, Vp=Vp, EG=EG, Ep=Ep)

# ------------------------------------------------------------------ ONE-LINK validation
def one_link_linear(beta, theta, v, n=1_500_000, seed=3):
    t, rho = samp_unif_ball(n, seed); B = boundary_quat(theta, n)
    def S_of(tq):
        return beta*(1.0 - quat_mul(flat_exp(tq), B)[..., 0])
    eps = 1e-3; va = np.broadcast_to(v, t.shape)
    Sp = S_of(t+eps*va); Sm = S_of(t-eps*va); S0 = S_of(t)
    G = (Sp-Sm)/(2*eps); H = (Sp-2*S0+Sm)/eps**2
    pp, ppp = phi_grad_hess(t, va)
    return _agg(np.exp(-S0)*Jrad(rho), G, H, pp, ppp)

# ------------------------------------------------------------------ TWO-LINK, LINEAR P
def two_link_linear(beta, theta, vA, vV, n=2_000_000, seed=5, ret_extra=False):
    tA, rA = samp_unif_ball(n, seed);      tV, rV = samp_unif_ball(n, seed+1)
    B = boundary_quat(theta, n)
    def S_of(a, b):
        W = quat_mul(flat_exp(a), flat_exp(b)); BV = quat_mul(B, quat_conj(flat_exp(b)))
        return beta*(2.0 - W[..., 0] - BV[..., 0])
    eps = 1e-3; vAa = np.broadcast_to(vA, tA.shape); vVa = np.broadcast_to(vV, tV.shape)
    Sp = S_of(tA+eps*vAa, tV+eps*vVa); Sm = S_of(tA-eps*vAa, tV-eps*vVa); S0 = S_of(tA, tV)
    G = (Sp-Sm)/(2*eps); H = (Sp-2*S0+Sm)/eps**2
    pa, ppa = phi_grad_hess(tA, vAa); pb, ppb = phi_grad_hess(tV, vVa)
    out = _agg(np.exp(-S0)*Jrad(rA)*Jrad(rV), G, H, pa+pb, ppa+ppb)
    if ret_extra:
        out['halfW'] = float((S0*0 + quat_mul(flat_exp(tA), flat_exp(tV))[..., 0]).mean())  # crude
        Wt = quat_mul(flat_exp(tA), flat_exp(tV))[..., 0]; BVt = quat_mul(B, quat_conj(flat_exp(tV)))[..., 0]
        w = np.exp(-S0)*Jrad(rA)*Jrad(rV); Z = w.sum()
        out['hW'] = float((w*Wt).sum()/Z); out['hBV'] = float((w*BVt).sum()/Z)
    return out

# ------------------------------------------------------------------ TWO-LINK, COVARIANT (right-invariant)
def two_link_covariant(beta, theta, wA, wV, n=1_200_000, seed=9, delta=2e-4):
    """Fluctuate each free link right-invariantly: A(eps)=exp(tA)exp(eps wA), V(eps)=exp(tV)exp(eps wV).
       G,H from S at the perturbed GROUP elements; Phi',Phi'' from radial Phi at the perturbed chart coord."""
    tA, rA = samp_unif_ball(n, seed); tV, rV = samp_unif_ball(n, seed+1)
    B = boundary_quat(theta, n)
    qA0 = flat_exp(tA); qV0 = flat_exp(tV)
    def right(q0, w, eps):
        return quat_mul(q0, flat_exp(eps*w[None, :]))
    def S_grp(qa, qv):
        W = quat_mul(qa, qv); BV = quat_mul(B, quat_conj(qv))
        return beta*(2.0 - W[..., 0] - BV[..., 0])
    def phi_at(q):                          # radial Phi of the chart coordinate of q
        _, r = flat_log(q); return phi_scalar(r)
    d = delta
    Ap = right(qA0, wA, d);  Am = right(qA0, wA, -d)
    Vp = right(qV0, wV, d);  Vm = right(qV0, wV, -d)
    Sp = S_grp(Ap, Vp); Sm = S_grp(Am, Vm); S0 = S_grp(qA0, qV0)
    G = (Sp-Sm)/(2*d); H = (Sp-2*S0+Sm)/d**2
    # Phi along the same right-fluctuation flow
    phA0 = phi_at(qA0); phV0 = phi_at(qV0)
    phAp = phi_at(Ap);  phAm = phi_at(Am)
    phVp = phi_at(Vp);  phVm = phi_at(Vm)
    pp  = ((phAp-phAm)/(2*d)) + ((phVp-phVm)/(2*d))
    ppp = ((phAp-2*phA0+phAm)/d**2) + ((phVp-2*phV0+phVm)/d**2)
    w = np.exp(-S0)*Jrad(rA)*Jrad(rV)
    # guard pathological antipode samples (weight ~ 0, but avoid inf*0): clip weight where phi huge
    bad = ~np.isfinite(ppp) | ~np.isfinite(pp) | (rA > 2*np.pi-1e-3) | (rV > 2*np.pi-1e-3)
    w = np.where(bad, 0.0, w); pp = np.where(bad, 0.0, pp); ppp = np.where(bad, 0.0, ppp)
    return _agg(w, G, H, pp, ppp)

# ------------------------------------------------------------------ divergence check (covariant)
def right_field_divergence(w, n=200000, seed=13, h=1e-5):
    """Numerically check div_flat v = -Phi'_v for the right-invariant field v(xi)=d/deps log(exp(xi)exp(eps w)).
       returns (div v, -Phi'_v) averaged under... just pointwise compare on random chart points (interior)."""
    rng = np.random.default_rng(seed)
    t = rng.uniform(-1, 1, (n, 3)); t *= (np.pi*(0.6))  # interior ball
    def vfield(tq):
        qp = quat_mul(flat_exp(tq), flat_exp((h)*w[None, :]))
        qm = quat_mul(flat_exp(tq), flat_exp((-h)*w[None, :]))
        tp, _ = flat_log(qp); tm, _ = flat_log(qm)
        return (tp-tm)/(2*h)
    div = np.zeros(n)
    for i in range(3):
        e = np.zeros(3); e[i] = 1.0
        div += (vfield(t + h*e) - vfield(t - h*e))[:, i]/(2*h)
    pp, _ = phi_grad_hess(t, vfield(t))
    return float(np.mean(div)), float(np.mean(-pp))

def log(*a): print(*a, flush=True)

# ================================================================== MAIN
if __name__ == "__main__":
    z = np.array([0., 0., 1.]); x = np.array([1., 0., 0.]); o = np.zeros(3)

    log("="*100)
    log("0. ONE-LINK validation (uniform-sampled finite-variance code vs S2 exact quadrature)")
    log("   S2 exact:  b=0.5 th=0 e3 W=0.0455213 Gpp=0.0910425 | b=2.5 th=pi e1 W=0.0658653 | b=2.5 th=0 W=0.1319566")
    for beta, th, nm, v in [(0.5, 0.0, 'e3', z), (0.5, np.pi, 'e3', z),
                            (2.5, 0.0, 'e3', z), (2.5, np.pi, 'e1', x)]:
        q = one_link_linear(beta, th, v)
        log(f"   b={beta:4.1f} th={th:5.3f} {nm}: E[H]-VarG={q['EH']-q['VG']:.6f} "
            f"W={q['W']:.6f} Gpp={q['Gpp']:.6f} R={q['R']:+.6f}")

    log("="*100)
    log("1. FACTORIZATION CHECK (two-link measure = mu_I(W) x mu_B(V), W=AV independent of V)")
    log("   Expect E[1/2tr AV] = E[1/2tr BV^-1] = one-plaquette value at the SAME beta, boundary th.")
    for beta, th in [(0.5, 0.0), (2.5, 0.0)]:
        q = two_link_linear(beta, th, z, o, n=600000, ret_extra=True)
        log(f"   b={beta:4.1f} th={th:.2f}: E[1/2tr AV]={q['hW']:.4f}   E[1/2tr BV^-1]={q['hBV']:.4f}")

    dirs = {
        "A3 ": (z, o), "V3 ": (o, z), "A1 ": (x, o), "V1 ": (o, x),
        "D3 ": (z/np.sqrt(2), z/np.sqrt(2)),   # <-- LINEAR P: constant combination of BOTH free links
        "D1 ": (x/np.sqrt(2), x/np.sqrt(2)),
    }
    log("="*100)
    log("2. TWO-LINK, LINEAR directions (constant v in flat chart).  THEOREM: R = G''-2W = 0 exactly.")
    log(f"   {'beta':>4} {'th':>6} {'dir':>4} {'E[H]-VarG':>10} {'W_geo':>10} {'Gpp':>10} {'R':>10}")
    for beta in (0.5, 1.0, 2.5):
        for th in (0.0, np.pi/2, np.pi):
            for nm, (vA, vV) in dirs.items():
                q = two_link_linear(beta, th, vA, vV, n=2_000_000)
                log(f"   {beta:4.1f} {th:6.3f} {nm:>4} {q['EH']-q['VG']:10.5f} {q['W']:10.5f} "
                    f"{q['Gpp']:10.5f} {q['R']:+10.5f}")
            log("")

    log("="*100)
    log("3. ACID TEST  beta=0 (pure Haar): S=0 => G=H=0 => G''=W_geo ; must be 0 for linear dirs.")
    for nm, (vA, vV) in [("A3", (z, o)), ("D3", (z/np.sqrt(2), z/np.sqrt(2)))]:
        q = two_link_linear(0.0, 0.0, vA, vV, n=2_000_000)
        log(f"   beta=0 {nm}: E[H]-VarG={q['EH']-q['VG']:+.6f}  W_geo={q['W']:+.6f}  G''={q['Gpp']:+.6f}")

    log("="*100)
    log("4. COVARIANT (right-invariant) direction: U(r)->U exp(r w).  Field-dependent => R generally != 0.")
    log("   [divergence check]  <div_flat v> = %.5f   < -Phi'_v > = %.5f  (should be equal)" %
        right_field_divergence(z))
    log(f"   {'beta':>4} {'th':>6} {'dir':>4} {'E[H]-VarG':>10} {'W_geo':>10} {'Gpp':>10} {'R':>10}")
    cdirs = {"C3": (z, z), "C1": (x, x), "C30": (z, o), "C0V": (o, z)}
    for beta in (0.5, 2.5):
        for th in (0.0, np.pi/2, np.pi):
            for nm, (wA, wV) in cdirs.items():
                q = two_link_covariant(beta, th, wA, wV, n=1_200_000)
                log(f"   {beta:4.1f} {th:6.3f} {nm:>4} {q['EH']-q['VG']:10.5f} {q['W']:10.5f} "
                    f"{q['Gpp']:10.5f} {q['R']:+10.5f}")
            log("")
    log("   ACID beta=0 covariant (G''=W_geo):")
    for nm, (wA, wV) in [("C3", (z, z)), ("C30", (z, o))]:
        q = two_link_covariant(0.0, 0.0, wA, wV, n=1_200_000)
        log(f"   beta=0 {nm}: W_geo={q['W']:+.6f}  G''={q['Gpp']:+.6f}   (E[H]-VarG={q['EH']-q['VG']:+.6f})")
    log("="*100)
    log("DONE.")
