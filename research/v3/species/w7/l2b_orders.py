#!/usr/bin/env python3
"""
l2b_orders.py -- L2B: strong-coupling orders of  S = E[H_vv] , D = Var(G_v)
                 and the order-matching test  D <= C(1-f1^2) S .

Model: the corpus's exactly solvable cell (see w7/ds_ratio.py, same quadrature):
  one plaquette + fixed boundary holonomy B=exp(i th sigma3/2), U=exp(i xi.sigma/2),
  S_W = beta (1 - (1/2)Re tr(U B)),  Phi = -log J,  J ~ sin^2(rho/2)/rho^2,
  measure exp(-S_W - Phi) d^3 xi,  |xi| < 2 pi.
Two splits of the exact marginal identity  G'' = S - D  (both are A2 shifted by psi):
  WF  (Wilson-force split):  S = E[H_W] + E[Phi''],  D = Var(d_v S_W)
  A2  (canonical):           S = E[H_W] + E[Phi''],  D = Var(d_v S_W + d_v Phi)
Ground truth for both: G'' = -d^2/da^2 log rho(a) at a0 (5-point FD).

Outputs
  [C1] S,D in both splits vs beta, + FD curvature (identity self-check)
  [C2] log-log slopes: D ~ beta^p , S ~ beta^q   (orders in u = beta/(2N^2))
  [C3] D/S vs the archived closed form  beta^2/beta*^2  (beta* = 4.30661)
  [C4] order matching:  (D/S) / (1-f1^2) ,  1-f1^2 = ubar^{T+3}/(1-ubar)^3 , ubar = eDelta*u
  [C5] acid tests
Run: py -3.12 -I l2b_orders.py
"""
import numpy as np

PI = np.pi
EDELTA = 54.3656365692          # e*Delta, Delta=20 (P3-FULLREP-KP:89, C01-CBAR1:88)
BETA_STAR = 4.30661             # archived zero of G'' in this cell (DS-RATIO §2.1)
A_ARCH, B_ARCH = 0.2051152971, 0.0110592938
E3 = np.array([0.0, 0.0, 1.0])


def onb(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    a = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(v, a); e1 /= np.linalg.norm(e1)
    return v, e1, np.cross(v, e1)


def slice_nodes(a0, v, ns, nph):
    v, e1, e2 = onb(v)
    R = np.sqrt(max(4.0 * PI ** 2 - a0 * a0, 1e-300))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5 * R * (xs + 1.0); wsv = 0.5 * R * ws
    ph = PI * (xp + 1.0); wpv = PI * wp
    S, PH = np.meshgrid(s, ph, indexing='ij')
    WS, WPH = np.meshgrid(wsv, wpv, indexing='ij')
    xi = (a0 * v[None, None, :]
          + (S * np.cos(PH))[..., None] * e1[None, None, :]
          + (S * np.sin(PH))[..., None] * e2[None, None, :])
    return xi.reshape(-1, 3), (WS * S * WPH).reshape(-1), v


def fields(xi, v, beta, theta):
    r = np.linalg.norm(xi, axis=1)
    q = (xi @ v) / r
    u = xi[:, 2] / r
    v3 = v[2]
    ct, st = np.cos(theta / 2.0), np.sin(theta / 2.0)
    h = 0.5 * r; c1, s1 = np.cos(h), np.sin(h)
    cth = c1 / s1
    csch2 = 1.0 / s1 ** 2
    dvh = 0.5 * q
    d2vh = 0.5 * (1.0 - q * q) / r
    dvu = (v3 - u * q) / r
    d2vu = -2.0 * q * (v3 - u * q) / r ** 2 - u * (1.0 - q * q) / r ** 2
    G = beta * (ct * s1 * dvh + st * (c1 * dvh * u + s1 * dvu))
    H = beta * (ct * (c1 * dvh ** 2 + s1 * d2vh)
                + st * (-s1 * dvh ** 2 * u + c1 * d2vh * u
                        + 2.0 * c1 * dvh * dvu + s1 * d2vu))
    P1 = q * (2.0 / r - cth)
    P2 = 2.0 * (1.0 - q * q) / r ** 2 - cth * (1.0 - q * q) / r + 0.5 * q * q * csch2
    return dict(Sw=beta * (1.0 - ct * c1 + st * s1 * u),
                Phi=2.0 * np.log(r) - 2.0 * np.log(s1),
                G=G, H=H, P1=P1, P2=P2)


def moments(a0, beta, theta=PI, v=E3, ns=400, nph=40):
    xi, w, vu = slice_nodes(a0, v, ns, nph)
    f = fields(xi, vu, beta, theta)
    lw = np.log(w) - f['Sw'] - f['Phi']
    ww = np.exp(lw); Z = ww.sum()
    E = lambda x: float((ww * x).sum() / Z)
    EG, EP1 = E(f['G']), E(f['P1'])
    return dict(rho=np.log(Z), EH=E(f['H']), EP2=E(f['P2']),
                VarG=E(f['G'] ** 2) - EG ** 2,
                VarP1=E(f['P1'] ** 2) - EP1 ** 2,
                Cov=E(f['G'] * f['P1']) - EG * EP1)


def fd_curv(a0, beta, theta=PI, da=2e-4, ns=600, nph=48):
    l = [moments(a0 + k * da, beta, theta, ns=ns, nph=nph)['rho'] for k in (-2, -1, 0, 1, 2)]
    return -(-l[0] + 16 * l[1] - 30 * l[2] + 16 * l[3] - l[4]) / (12 * da * da)


def splits(m):
    S = m['EH'] + m['EP2']
    D_wf = m['VarG']
    D_a2 = m['VarG'] + m['VarP1'] + 2.0 * m['Cov']
    return S, D_wf, D_a2


def u_of(beta, N=2):
    return beta / (2.0 * N * N)


def P(*a): print(*a, flush=True)


BETAS = [0.0, 1.0 / 384.0, 1.0 / 48.0, 0.05, 0.1, 0.5, 1.0, 2.5]

P("#" * 100)
P("# L2B-ORDERS  (S = E[H_vv], D = Var(G_v), G'' = S - D)  --  theta=pi, v=e3, a0=0")
P(f"# u = beta/(2N^2), N=2 ;  eDelta = {EDELTA} ;  archived A={A_ARCH} B={B_ARCH} beta*={BETA_STAR}")
P("#" * 100)

# ---------------------------------------------------------------- [C1]
P("\n[C1] moments vs beta (beta=0 row is the pure-Haar cell; FD = -d2/da2 log rho)")
P(f"{'beta':>10} {'u':>11} {'EH':>11} {'EPhi2':>11} {'S':>11} "
  f"{'VarGw':>11} {'VarPhi1':>11} {'Cov':>11} {'D_wf':>11} {'D_a2':>11} {'S-D_wf':>10} {'FD':>10}")
rows = []
for b in BETAS:
    m = moments(0.0, b)
    S, Dwf, Da2 = splits(m)
    fd = fd_curv(0.0, b) if b in (0.0, 1.0 / 384.0, 1.0 / 48.0, 2.5) else float('nan')
    rows.append((b, S, Dwf, Da2, m, fd))
    P(f"{b:10.6f} {u_of(b):11.3e} {m['EH']:11.7f} {m['EP2']:11.7f} {S:11.7f} "
      f"{m['VarG']:11.7f} {m['VarP1']:11.7f} {m['Cov']:11.7f} {Dwf:11.7f} {Da2:11.7f} "
      f"{S - Dwf:10.7f} {fd:10.7f}")

P("\n  identity check (|S-D_wf - FD| and |S-D_a2 - FD|, excluding beta=0 cusp row):")
for b, S, Dwf, Da2, m, fd in rows:
    if not np.isnan(fd) and b > 0:
        P(f"   beta={b:9.6f}   |wf-FD|={abs(S-Dwf-fd):.2e}   |a2-FD|={abs(S-Da2-fd):.2e}")

# ---------------------------------------------------------------- [C2]
P("\n[C2] log-log slopes over beta in (0, 0.1]  (u ~ beta, so slope = order in u)")
sel = [(b, S, Dwf, Da2) for (b, S, Dwf, Da2, m, fd) in rows if 0.0 < b <= 0.1]
lb = np.log([s[0] for s in sel])
for name, idx in (("S = E[H]+E[Phi'']", 1), ("D_wf = Var(G_W)", 2), ("D_a2 = Var(G_W+Phi')", 3)):
    ls = np.log([s[idx] for s in sel])
    p = np.polyfit(lb, ls, 1)
    P(f"   {name:26s} : D = C * beta^{p[0]:.4f}   (C = {np.exp(p[1]):.6f})  => O(u^{p[0]:.2f})")
D0 = sel[0][2]; b0 = sel[0][0]
P(f"   ratio D_wf/beta^2 at beta->0 : {D0/b0**2:.9f}   (archived B = {B_ARCH:.9f})")
P(f"   S at beta=0                  : {rows[0][1]:.9f}   (archived A = {A_ARCH:.9f})")
P(f"   archived beta* = {BETA_STAR} <=>  B/A = {B_ARCH/A_ARCH:.9f} , 1/beta*^2 = {1/BETA_STAR**2:.9f}")

# ---------------------------------------------------------------- [C3]
P("\n[C3] D/S vs closed form beta^2/beta*^2   (WF split; A2 column for the split warning)")
P(f"{'beta':>10} {'D_wf/S':>12} {'beta^2/bs^2':>13} {'rel.err':>10} {'D_a2/S':>12} {'Gpp_a2':>11}")
for b, S, Dwf, Da2, m, fd in rows:
    cf = b ** 2 / BETA_STAR ** 2
    P(f"{b:10.6f} {Dwf/S:12.7f} {cf:13.7f} {abs(Dwf/S-cf)/max(cf,1e-30) if cf>0 else 0:10.2e} "
      f"{Da2/S:12.7f} {S-Da2:11.7f}")

# ---------------------------------------------------------------- [C4]
P("\n[C4] order matching:  X(T) = (D/S) / (1-f1^2)  with  1-f1^2 = Bbar*ubar^T,")
P("     Bbar=(ubar/(1-ubar))^3,  ubar = eDelta*u  (X01-L1C:79-80; Y04:84; C03-BYPASS:96)")
for beta in (1.0 / 48.0, 1.0 / 480.0):
    u = u_of(beta); ub = EDELTA * u
    S, Dwf, _ = splits(moments(0.0, beta))
    dsf = Dwf / S
    bb = (ub / (1.0 - ub)) ** 3
    P(f"\n  beta={beta:.6f}  u={u:.6e}  ubar={ub:.6e}  D/S={dsf:.6e}  Bbar={bb:.6e}")
    P(f"  {'T':>3} {'1-f1^2 = Bbar*ubar^T':>22} {'(D/S)/(1-f1^2)':>16} {'gap exponent':>13}")
    for T in range(0, 9):
        f1d = bb * ub ** T
        ratio = dsf / f1d
        P(f"  {T:3d} {f1d:22.6e} {ratio:16.6e} {np.log(ratio)/np.log(1.0/u):13.3f}")
    # exact matching order: solve  ubar^{T+3}/(1-ubar)^3 = D/S
    Ts = np.log(dsf * (1.0 - ub) ** 3) / np.log(ub) - 3.0
    P(f"  matching T* (ubar^(T+3)/(1-ubar)^3 = D/S) :  T* = {Ts:+.4f}"
      f"   -> {'NO integer T>=0 exists' if Ts < 0 else 'T>=0 exists'}")

# ---------------------------------------------------------------- [C5]
P("\n[C5] acid tests")
P("  (i) beta=0 pure Haar: G''=%+.7f > 0  while f1^2 = 0 (no Wilson weight, no atom)"
  % (rows[0][1] - rows[0][2]))
P("      => curvature defect WITHOUT residue: the implication residue=>curvature is not an equivalence.")
P("  (ii) free massive Gaussian: H=m^2, G linear => D = Var(G) = 0 exactly, G'' = m^2 > 0, f1^2 = 0")
P("      => D/S = 0 = O(u^inf) while 1-f1^2 = 1 : any bound D <= C(1-f1^2)S is satisfied vacuously.")
P("  (iii) D_wf/S = 0 exactly at beta=0 and grows as beta^2 => D/S is controlled by beta alone,")
P("        not by the residue; the same D/S is reproduced for every f1^2.")
P("\n# end")
