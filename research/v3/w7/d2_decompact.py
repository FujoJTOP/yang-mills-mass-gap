#!/usr/bin/env python3
"""
d2_decompact.py -- D2: analytic L->infinity limit of the constant-mode curvature.

Question (SM-WAVE-SYNTHESIS, H1): does the constant-mode curvature G''(a0*) at beta=2.5
tend to 0 as L->infinity (H1), or to a positive L-independent limit (the beta=0 Haar floor)?

Setup tested here: the 1xL strip of sm1_blocks.py (free links X_1..X_L, L plaquettes
P_1=X_1 X_2, P_i=B_{i-1} X_i^{-1} X_{i+1} (2<=i<=L-1), P_L=B_{L-1} X_L^{-1}; boundary B=e^{i theta sigma3/2}),
SO(3)-corrected 3-component flat chart xi_l in B(2pi), Haar J(rho)=sin^2(rho/2)/rho^2 (exact sampling),
Wilson S=beta*sum_p (1 - 1/2 tr P_p).  Constant (k=0) coarse mode v=(e3,...,e3)/sqrt L, a=<v,xi>.

We sample the EXACT target measure  dmu ∝ prod_l J(rho_l) d^3 xi_l  e^{-S}
by random-walk Metropolis (no weight degeneracy), measure
   Var(a) = <a^2>_c  (bounded  <=>  G''->positive ; grows ∝ L <=>  G''->0)
and the marginal-curvature G''(a0*) = -d^2 log rho_v / da^2 at the histogram peak.
beta=0 is the exact-Haar baseline; the one-link value 1/s^2 = 9/(4 pi^2 - 6) = 0.26883 is the Haar floor.

py -3.12 d2_decompact.py   (numpy only; one job at a time)
"""
import numpy as np

def su2_exp(x):
    # x:(...,3) -> quaternion (w,x,y,z) rep of exp(i x.sigma/2)
    th = np.linalg.norm(x, axis=-1, keepdims=True)
    n = np.divide(x, th, out=np.zeros_like(x), where=th > 1e-12)
    w = np.cos(th[..., 0] / 2); v = n * np.sin(th / 2)
    return np.concatenate([w[..., None], v], axis=-1)          # (...,4)

def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0*b0 - a1*b1 - a2*b2 - a3*b3,
                     a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1,
                     a0*b3 + a1*b2 - a2*b1 + a3*b0], axis=-1)

def qdag(q): return q * np.array([1.0, -1.0, -1.0, -1.0])

def half_tr(q): return q[..., 0]          # 1/2 tr U = Re(w) = w

def Jdens(rho):
    # J(rho)=sin^2(rho/2)/rho^2 ; safe at 0 -> 1/4
    s = np.sin(rho / 2)
    out = np.where(rho > 1e-9, (s / np.maximum(rho, 1e-9)) ** 2, 0.25)
    return out

def logJ(rho):
    return np.where(rho > 1e-9, 2*np.log(np.maximum(np.sin(rho/2), 1e-300)) - 2*np.log(np.maximum(rho, 1e-300)),
                    np.log(0.25))

def action(xi, B, beta):
    # xi: (nfree=L,3);  B: list of L-1 matrices (quat);  returns beta*(L - sum 1/2 tr P)
    q = su2_exp(xi)                       # (L,4)
    P = []
    P.append(qmul(q[0], q[1]))                                  # P_1 = X1 X2
    for i in range(1, L_of(xi)-1):
        P.append(qmul(qmul(B[i-1], qdag(q[i])), q[i+1]))        # B_{i-1} X_i^-1 X_{i+1}
    P.append(qmul(B[-1], qdag(q[-1])))                          # P_L = B_{L-1} X_L^-1
    st = np.sum([half_tr(p) for p in P])
    return beta * (len(P) - st)

def L_of(xi): return xi.shape[0]

def S_full(xi, B, beta):
    q = su2_exp(np.asarray(xi))
    Ln = xi.shape[0]
    plaq = qmul(q[0], q[1])
    tot = half_tr(plaq)
    for i in range(1, Ln-1):
        plaq = qmul(qmul(B[i-1], qdag(q[i])), q[i+1]); tot = tot + half_tr(plaq)
    plaq = qmul(B[-1], qdag(q[-1])); tot = tot + half_tr(plaq)
    return beta * (Ln - tot)

def boundary(theta, Ln):
    b = su2_exp(np.array([0.0, 0.0, theta]))
    return [b.copy() for _ in range(Ln-1)]

def metropolis(Ln, beta, theta, nsweep, seed, thin=2, burn=0.2, step=0.45):
    rng = np.random.default_rng(seed)
    B = boundary(theta, Ln)
    xi = rng.uniform(-1, 1, size=(Ln, 3)) * 1.0
    logP = -S_full(xi, B, beta) + np.sum(logJ(np.linalg.norm(xi, axis=1)))
    def a_of(x): return np.einsum('ld,d->', x, np.array([0,0,1.0])) / np.sqrt(Ln)
    a_samples = []
    v = np.array([0.0, 0.0, 1.0])
    nacc = 0; ntry = 0
    nburn = int(burn * nsweep)
    for sw in range(nsweep):
        for l in range(Ln):
            prop = xi.copy()
            prop[l] = xi[l] + rng.uniform(-step, step, size=3)
            rp = np.linalg.norm(prop[l])
            if rp >= 2*np.pi:    # stay inside the flat chart ball B(2pi)
                ntry += 1; continue
            newlogP = -S_full(prop, B, beta) + np.sum(logJ(np.linalg.norm(prop, axis=1)))
            ntry += 1
            if np.log(rng.random() + 1e-300) < (newlogP - logP):
                xi = prop; logP = newlogP; nacc += 1
        if sw >= nburn and (sw % thin == 0):
            a_samples.append(a_of(xi))
    a = np.array(a_samples)
    return a, nacc/max(ntry,1)

def Gpp_from_hist(a, nbin=160):
    lo, hi = np.percentile(a, 0.5), np.percentile(a, 99.5)
    h, e = np.histogram(a, bins=nbin, range=(lo, hi), density=True)
    c = 0.5*(e[1:]+e[:-1])
    # smooth
    k = np.ones(5)/5.0
    hs = np.convolve(h, k, mode='same')
    hs = np.maximum(hs, 1e-12)
    im = int(np.argmax(hs))
    i0, i1, i2 = im-1, im, im+1
    if i0 < 0 or i2 >= len(c): return np.nan
    d2 = (np.log(hs[i2]) - 2*np.log(hs[i1]) + np.log(hs[i0])) / (c[1]-c[0])**2
    return -d2

if __name__ == "__main__":
    print("="*92)
    print("D2 -- decompactification of the constant-mode curvature (1xL strip, k=0 LONG mode)")
    print("Haar floor 1/s^2 = 9/(4 pi^2 - 6) =", 9/(4*np.pi**2-6))
    print("="*92)
    print(f"{'beta':>5}{'L':>4}{'Var(a)':>12}{'1/Var(a)':>11}{'Gpp(a0*)':>11}{'acc':>7}{'a0*':>8}")
    for beta in (0.0, 2.5):
        for Ln in (2, 3, 4, 6, 8, 10):
            N = 260000 if Ln <= 6 else 200000
            a, acc = metropolis(Ln, beta, np.pi, N, seed=12345 + int(10*beta) + Ln, thin=2)
            var = np.var(a); g = Gpp_from_hist(a); a0 = np.median(a)
            print(f"{beta:5.1f}{Ln:4d}{var:12.4f}{1/var:11.4f}{g:11.4f}{acc:7.3f}{a0:8.4f}")
    print("\nREAD: at beta=0, Var(a)->s^2=3.72 (bounded), 1/Var->0.269 = Haar floor.")
    print("If at beta=2.5 Var(a) stays bounded as L grows => G''->positive (Haar floor).")
    print("If Var(a) grows ~ L => G'' ~ 1/L -> 0 (H1).")
    print("Done. No QED.")
