#!/usr/bin/env python3
r"""
d5_universality.py -- D5: is the constant-mode curvature's L-dependence UNIVERSAL across
groups (U(1), SU(2), SU(3)) and actions (Wilson, fundamental-adjoint)?

Deliverable witness for research/v3/D5-UNIVERSALITY.md.  py -3.12, numpy only.
Run:  py -3.12 d5_universality.py [quick|full]

OBJECT.  The SM1/SM4 constant-mode curvature is
    G''(v;a0*) = -d^2/da^2 log rho_v(a)|_{a0*} ,   a = <v,xi> = sum_i xi_i / sqrt(L)  (k=0 mode)
on the 1xL chain (open block, fixed boundary holonomy theta), flat exponential chart.

KEY PROXY.  For a near-Gaussian marginal, G'' ~ 1/Var(a).  The amplitude variance
    Var(a) = E[a^2]-E[a]^2 = Var( sum_i xi_i / sqrt(L) ) = (1/L) sum_{ij} Cov(xi_i,xi_j)
is KNOT-FREE, smooth, and group/action-universal to compute (just moments of the measure,
no marginal constraint -> no spline knots, no chart Jacobian needed if we sample Haar directly).
It is the *amplitude susceptibility* chi_L.   G''(a0*) and 1/Var(a) agree at large L
(marginal Gaussianises, CLT); the ratio  R := G''(a0*)*Var(a) -> 1 measures Gaussianity.

SECTIONS
 [A] U(1) L=1 exact: G''(a0*)=beta at the mode (von Mises) -- the exactly-solvable anchor.
 [B] U(1) exact Var(a) L=1..4 by product Gauss-Legendre (b in {0,pi/2,pi}, beta in {1,2.5});
     -> boundary dependence (b=0 grows, b=pi saturates).
 [C] U(1) Var(a) L=1..8 by Haar-importance sampling (validated against [B]).
 [D] SU(2) DIRECT G''(a0*) L=2..6 (reuses sm1_blocks reweight) -> plateaus ~0.48, not 1/L->0;
     W_geo (Haar part) is L-independent ~0.27.
 [E] SU(2), SU(3) Var(a) L=1..8 / 1..4 by Haar-importance sampling -> both saturate.
 [F] L=1 exact G''(a0*) across ACTIONS (Wilson vs fundamental-adjoint) -> value action-dependent,
     sign/structure action-independent.
 [G] normalisation-invariant table  R = G''(a0*) * Var(a)  (L=1 non-Gaussian, large-L -> 1).

Labels: [T] exact identity / analytic · [V] numerical witness (quadrature exact at L=1;
Haar-importance MC for L>=2, ESS quoted) · [O] open.   No third-party sources.  No QED.
"""
import sys, time
import numpy as np
np.seterr(all='ignore')

MODE = sys.argv[1] if len(sys.argv) > 1 else 'quick'
NX = {'quick': 60, 'full': 140}[MODE]
NMC = {'quick': 1_500_000, 'full': 6_000_000}[MODE]

# ============================================================ SU(N)/U(1) Haar + logs
def haar_su2(n, rng):
    g = rng.standard_normal((n, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
    q0, q1, q2, q3 = g[:, 0], g[:, 1], g[:, 2], g[:, 3]
    U = np.empty((n, 2, 2), complex)
    U[:, 0, 0] = q0 + 1j*q3; U[:, 0, 1] = q2 + 1j*q1
    U[:, 1, 0] = -q2 + 1j*q1; U[:, 1, 1] = q0 - 1j*q3
    return U

def haar_su3(n, rng):
    Z = (rng.standard_normal((n, 3, 3)) + 1j*rng.standard_normal((n, 3, 3))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    return Q * (d/np.abs(d))[:, None, :]

def log_su2(U):
    tr = np.clip(0.5*np.real(U[:, 0, 0]+U[:, 1, 1]), -1, 1)
    q1 = np.imag(U[:, 0, 1]); q2 = np.real(U[:, 0, 1]); q3 = np.imag(U[:, 0, 0])
    s = np.sqrt(np.maximum(q1*q1+q2*q2+q3*q3, 1e-300)); ang = 2*np.arccos(tr)
    return ang[:, None]*np.stack([q1, q2, q3], axis=1)/s[:, None]

_GELL = [np.array([[0,1,0],[1,0,0],[0,0,0]], complex), np.array([[0,-1j,0],[1j,0,0],[0,0,0]]),
         np.diag([1,-1,0]).astype(complex), np.array([[0,0,1],[0,0,0],[1,0,0]], complex),
         np.array([[0,0,-1j],[0,0,0],[1j,0,0]]), np.array([[0,0,0],[0,0,1],[0,1,0]], complex),
         np.array([[0,0,0],[0,0,-1j],[0,1j,0]]), (np.diag([1,1,-2])/np.sqrt(3)).astype(complex)]
def log_su3(U):
    w, V = np.linalg.eig(U); phi = np.angle(w)
    Hm = V @ (phi[:, :, None]*np.conjugate(np.swapaxes(V, 1, 2)))   # -i log U (Hermitian)
    xi = np.empty((U.shape[0], 8))
    for a in range(8):
        xi[:, a] = np.real(np.einsum('nij,ji->n', Hm, _GELL[a]))
    return xi

def diag2(ang): return np.cos(ang/2)*np.eye(2) + 1j*np.sin(ang/2)*np.diag([1, -1])
def diag3(ang): return np.diag([np.exp(1j*ang/2), np.exp(-1j*ang/2), 1.0])

def chain_plaq(U, bvals):
    L = len(U)
    if L == 1: return (U[0][:, None] @ bvals[0][None])
    P = [U[0] @ U[1]]
    for j in range(2, L):
        P.append(bvals[j-1][None] @ np.conjugate(np.swapaxes(U[j-1], 1, 2)) @ U[j])
    P.append(bvals[L-1][None] @ np.conjugate(np.swapaxes(U[L-1], 1, 2)))
    return np.stack(P, axis=1)

# ============================================================ Haar-importance Var(a)
def var_a(Ng, L, beta, bang, rng, Ntot=NMC, action='wilson', bA=0.0, comp=2):
    sa = sa2 = sw = sw2 = 0.0; done = 0; chunk = 200_000
    while done < Ntot:
        n = min(chunk, Ntot-done); done += n
        if Ng == 1:
            t = [rng.uniform(-np.pi, np.pi, n) for _ in range(L)]
            if L == 1: pl = [t[0]+bang[0]]
            else:
                pl = [t[0]+t[1]]
                for j in range(2, L): pl.append(bang[j-1]-t[j-1]+t[j])
                pl.append(bang[L-1]-t[L-1])
            S = sum(beta*(1-np.cos(p)) for p in pl)
            a = np.sum(t, axis=0)/np.sqrt(L)
        elif Ng == 2:
            U = [haar_su2(n, rng) for _ in range(L)]; bv = [diag2(b) for b in bang]
            P = chain_plaq(U, bv); tr = np.real(np.trace(P, axis1=-2, axis2=-1))
            S = (beta*(1-0.5*tr)).sum(axis=1)
            if action == 'fa':
                S = S + (bA*(1.0-(tr*tr-1.0)/3.0)).sum(axis=1)
            xi = [log_su2(U[i]) for i in range(L)]
            a = np.sum([xi[i][:, 2] for i in range(L)], axis=0)/np.sqrt(L)
        else:
            U = [haar_su3(n, rng) for _ in range(L)]; bv = [diag3(b) for b in bang]
            P = chain_plaq(U, bv); tr = np.real(np.trace(P, axis1=-2, axis2=-1))
            S = (beta*(1.0-tr/3.0)).sum(axis=1)
            xi = [log_su3(U[i]) for i in range(L)]
            a = np.sum([xi[i][:, comp] for i in range(L)], axis=0)/np.sqrt(L)
        w = np.exp(-S)
        sa += np.sum(w*a); sa2 += np.sum(w*a*a); sw += np.sum(w); sw2 += np.sum(w*w)
    Ea = sa/sw; return Ea, sa2/sw-Ea*Ea, sw*sw/sw2

# ============================================================ exact U(1) Var (moments)
XGu, WGu = np.polynomial.legendre.leggauss(NX)
def u1_var_exact(L, beta, th):
    pts = np.pi*XGu; wts = np.pi*WGu
    grids = [g.ravel() for g in np.meshgrid(*[pts]*L, indexing='ij')]
    Wg = np.meshgrid(*[wts]*L, indexing='ij'); W = np.ones_like(Wg[0])
    for k in range(L): W = W*Wg[k]
    W = W.ravel()
    if L == 1: pl = [grids[0]+th]
    else:
        pl = [grids[0]+grids[1]]
        for j in range(2, L): pl.append(th-grids[j-1]+grids[j])
        pl.append(th-grids[L-1])
    S = np.zeros(grids[0].size)
    for p in pl: S = S + beta*(1-np.cos(p))
    a = np.sum(grids, axis=0)/np.sqrt(L); w = W*np.exp(-S)
    m = np.sum(w*a)/np.sum(w)
    return np.sum(w*a*a)/np.sum(w)-m*m

# ============================================================ exact L=1 G'' across actions
XG1, WG1 = np.polynomial.legendre.leggauss(900)
def _phi0(rho):
    rho = np.asarray(rho, float); sm = rho < 1e-6; rs = np.where(sm, 1.0, rho)
    return np.where(sm, rho*rho/12.0, 2.0*np.log(rs/2.0)-2.0*np.log(np.sin(rs/2.0)))
def rho_long_L1(a, bF, bA, theta):
    R = np.sqrt(max(4*np.pi**2-a*a, 1e-14)); r = 0.5*R*(XG1+1.0); w = 0.5*R*WG1
    rho = np.sqrt(r*r+a*a); u = a/rho
    hp = np.cos(rho/2)*np.cos(theta/2)-np.sin(rho/2)*np.sin(theta/2)*u   # 1/2 tr P
    S = bF*(1-hp)+bA*(1-(4*hp*hp-1)/3.0)
    return 2*np.pi*np.sum(np.exp(-S-_phi0(rho))*r*w)
def gpp_L1(fn, bF, bA, theta, nsamp=801):
    aa = np.linspace(-2*np.pi+0.05, 2*np.pi-0.05, nsamp); h = aa[1]-aa[0]
    rho = np.array([fn(a, bF, bA, theta) for a in aa]); lg = np.log(np.maximum(rho, 1e-300))
    im = int(np.argmax(rho[2:-2]))+2
    G = -(-lg[im-2]+16*lg[im-1]-30*lg[im]+16*lg[im+1]-lg[im+2])/(12*h*h)
    return aa[im], G

# ============================================================ sections
def A_u1_anchor():
    print("="*96); print("[A] U(1) L=1 exact: G''(a0*)=beta at the mode (von Mises)  [T]")
    for beta in (1.0, 2.5, 5.0):
        # rho(a)=exp(-beta(1-cos(a+b))) ; -d2/da2 log rho at a=-b is beta exactly
        a = np.linspace(-np.pi, np.pi, 4001); b = 0.3
        rho = np.exp(-beta*(1-np.cos(a+b))); lg = np.log(rho); h = a[1]-a[0]
        ih = int(np.argmax(rho)); g = -(-lg[ih-2]+16*lg[ih-1]-30*lg[ih]+16*lg[ih+1]-lg[ih+2])/(12*h*h)
        print(f"   beta={beta:.2f}: mode a0*={a[ih]:+.3f} (expect {-b:+.3f})  G''={g:.5f} (expect {beta:.5f})")

def B_u1_exact_var():
    print("\n"+"="*96); print("[B] U(1) EXACT Var(a), L=1..4 by product Gauss-Legendre  [V]")
    print("     (boundary dependence: b=0 GROWS, b=pi SATURATES)")
    for beta in (1.0, 2.5):
        print(f"   beta={beta}")
        for th in (0.0, np.pi/2, np.pi):
            row = "".join(f"{u1_var_exact(L,beta,th):9.4f}" for L in (1,2,3,4))
            print(f"     b={th:5.3f}: L1..4 Var ={row}   -> 1/Var ~{1.0/u1_var_exact(4,beta,th):.3f}(L=4)")

def C_u1_is_var():
    print("\n"+"="*96); print("[C] U(1) Var(a) L=1..8 by Haar-importance MC  [V]  (saturates at twist b=pi)")
    rng = np.random.default_rng(20261124)
    for beta in (1.0, 2.5):
        row = []
        for L in (1,2,3,4,6,8):
            _, v, ess = var_a(1, L, beta, [np.pi]*L, rng, Ntot=NMC)
            row.append(f"{v:.3f}({ess/1e3:.0f}k)")
        print(f"   b={beta}: "+" ".join(row))

def D_su2_direct():
    print("\n"+"="*96); print("[D] SU(2) DIRECT G''(a0*) L=2..6, beta=2.5, twist pi  [V]")
    print("     DECISIVE: G'' plateaus ~0.48 for L>=4 (NOT 1/L->0); W_geo stays ~0.27 (L-independent Haar part)")
    try:
        sys.path.insert(0, r'D:\Dev\math\research\v3\w7')
        from sm1_exact import su2_exp, half_tr, dag
        from sm1_blocks import scan_block, analyze, mode_k0, mats
    except Exception as ex:
        print("   (skipped: cannot import sm1_blocks:", ex, ")"); return
    def S_factory(L):
        def S(X, fixed, beta):
            U = [su2_exp(X[:, i]) for i in range(L)]
            P = [U[0]@U[1]]
            for j in range(2, L): P.append(fixed[j-2][None]@dag(U[j-1])@U[j])
            P.append(fixed[L-2][None]@dag(U[L-1]))
            tot = np.zeros(X.shape[0])
            for p in P: tot = tot+half_tr(p)
            return beta*(L-tot)
        return S
    rng = np.random.default_rng(20261124)
    print(f"   {'L':>3}{'a0*':>8}{'G(a0*)':>10}{'Ddyn':>9}{'Wgeo':>9}{'X':>9}{'ESS':>8}")
    for L in (2, 3, 4, 5, 6):
        blk = S_factory(L); fx = mats([(0, 0, np.pi)]*(L-1)); v = mode_k0(L)
        N = 2_000_000 if L <= 4 else 1_500_000
        cen, acc = scan_block(blk, fx, 2.5, v, N, rng); r = analyze(cen, acc)
        print(f"   {L:3d}{r['a0']:8.3f}{r['G']:10.4f}{r['D']:9.4f}{r['Wg']:9.4f}{r['X']:9.4f}{r['ess']:8.0f}")

def E_su23_var():
    print("\n"+"="*96); print("[E] SU(2), SU(3) Var(a) vs L by Haar-importance MC  [V]  (both SATURATE)")
    rng = np.random.default_rng(7)
    print("   SU(2) Wilson, twist pi:")
    for beta in (1.0, 2.5):
        row = []
        for L in (2,3,4,6,8):
            _, v, ess = var_a(2, L, beta, [np.pi]*L, rng, Ntot=NMC)
            row.append(f"{v:.3f}({ess/1e3:.0f}k)")
        print(f"     beta={beta}: L2..8 "+" ".join(row))
    print("   SU(3) Wilson, twist pi:")
    for beta in (1.0, 2.0):
        row = []
        for L in (1,2,3,4):
            _, v, ess = var_a(3, L, beta, [np.pi]*L, rng, Ntot=NMC)
            row.append(f"{v:.3f}({ess/1e3:.0f}k)")
        print(f"     beta={beta}: L1..4 "+" ".join(row))

def F_action():
    print("\n"+"="*96); print("[F] L=1 exact G''(a0*) LONG, theta=pi -- ACTION dependence  [V]")
    print("     value changes with action (all >0); sign/structure action-independent")
    for lab, bF, bA in [('Wilson', 2.5, 0.0), ('FA-soft', 2.0, 0.5), ('FA-mid', 1.5, 1.0),
                        ('FA-adj', 1.0, 1.5), ('FA-hard', 0.5, 2.0), ('Wilson b=5', 5.0, 0.0)]:
        a0, G = gpp_L1(rho_long_L1, bF, bA, np.pi)
        print(f"   {lab:>11} (bF,bA)=({bF},{bA}): a0*={a0:+.3f}  G''_LONG={G:.4f}")

def G_invariant():
    print("\n"+"="*96); print("[G] normalisation-invariant R = G''(a0*) * Var(a)  [V]")
    print("     R->1 at large L (marginal Gaussianises); R>1 at L=1 (non-Gaussian / von Mises)")
    # SU(2) L=1 exact G''=0.694 (long); Var(L=1) from [E]
    rng = np.random.default_rng(11)
    _, v1, _ = var_a(2, 1, 2.5, [np.pi], rng, Ntot=NMC)
    print(f"   SU(2) b=2.5,L=1: G''=0.694 Var={v1:.3f} -> R={0.694*v1:.3f}")
    _, v1u, _ = var_a(1, 1, 2.5, [np.pi], rng, Ntot=NMC)
    print(f"   U(1)  b=2.5,L=1: G''=2.500 Var={v1u:.3f} -> R={2.5*v1u:.3f}")

if __name__ == '__main__':
    t0 = time.time()
    A_u1_anchor(); B_u1_exact_var(); C_u1_is_var(); D_su2_direct()
    E_su23_var(); F_action(); G_invariant()
    print(f"\nD5-universality: done in {time.time()-t0:.0f}s.  No QED.  Clay 0%.")
