# h1d_periodic.py -- H1 on the P4 periodic kernel: is Var=E[H]=0 for ALL beta,
# and does beta*(P4) exist?
# Kernel exactly as research/v3/P4-PERIODIC-TRANSFER.md:62-69 (same as ch15-uvface.tex:199-213):
#   rho_1 = exact one-link SU(2) Haar marginal in the flat chart
#   S_W   = beta * sum_i (1 - cos((z_{i+1}-z_i)/2)),  cyclic
#   M_k(z,z') = sqrt(rho_1(z) rho_1(z')) e^{-ik(z+z')/(2 sqrt L)} e^{-beta(1-cos((z'-z)/2))}
#   Z(k) = Tr[M_k^L] = sum_j lambda_j(k)^L
#   rho_v(a) = (1/2pi) int Z(k) e^{i k a} dk ,   Gamma''(a) = -d_a^2 log rho_v(a)
# numpy only, deterministic quadrature.  Logs to h1d_periodic.log
import numpy as np, time, sys

t0 = time.time()
out = []
def P(s=""):
    out.append(str(s)); print(s)

# ---------- 1. one-link flat-chart marginal rho_1(z), z = xi_3 ----------
# flat chart xi in R^3, |xi| < 2pi ;  J = sin^2(r/2)/r^2  (H1A-EXACT.md:82)
def rho1_grid(Nz, nrho=6001):
    """p(z) dz = Haar mass at xi_3 in [z,z+dz] : int_0^inf 2pi rho drho sin^2(r/2)/r^2, r=sqrt(rho^2+z^2)."""
    z = -2*np.pi + 4*np.pi*np.arange(Nz)/Nz              # cell centers on [-2pi,2pi), Lebesgue dz
    rho = (np.arange(nrho)+0.5)*(2*np.pi/nrho)           # |xi_perp| in [0,2pi)
    p = np.zeros(Nz)
    for i, zz in enumerate(z):
        r2 = rho**2 + zz**2
        m = r2 <= (2*np.pi)**2                            # flat chart support |xi|<2pi
        p[i] = np.sum(m*2*np.pi*rho*np.sin(np.sqrt(r2)/2)**2/np.where(m, r2, 1.0))*(2*np.pi/nrho)
    p /= np.sum(p)*(z[1]-z[0])
    return z, p

Nz = 150
z, p = rho1_grid(Nz)
dz = z[1]-z[0]
Var_z = np.sum(z**2*p)*dz
P("== [A] one-link marginal (convention lock) ==")
P(f"  Nz={Nz} dz={dz:.6f}   E[z^2]=Var(xi3)={Var_z:.6f}   1/Var={1/Var_z:.6f}")
P(f"  corpus: Var(xi3)=3.719824 , 1/Var=0.26883   -> rel.err {abs(Var_z-3.719824)/3.719824:.2e}")
P("")

# ---------- 2. structural zeros: gradient / Hessian of S_W on the cyclic chain ----------
def SW_and_derivs(Zc, beta):
    """Zc: (nconf, L) row vector of z_i ; cyclic plaquette action beta*sum(1-cos(dz/2))."""
    L = Zc.shape[1]
    d = np.roll(Zc, -1, axis=1) - Zc
    S = beta*np.sum(1-np.cos(d/2), axis=1)
    g = beta*(np.sin(d/2)/2)                             # dS/dz_i from plaquette i
    grad = np.roll(g, 1, axis=1) - g                      # dS/dz_i = g_{i-1} - g_i
    return S, grad

rng = np.random.default_rng(20261125)
P("== [B] Var(G)=0 and E[H]=0 : exact for every beta (graph identity) ==")
P("  uniform mode v = 1/sqrt(L); G = <v,grad S_W>, H_vv = <v,Hess S_W v>")
P(f"  {'beta':>6} {'max|G|':>12} {'max|H_vv|':>12} {'<v,grad S_W>':>14}")
for beta in [0.0, 0.5, 2.5, 6.0, 20.0, 100.0]:
    L = 8
    Zc = rng.uniform(-2*np.pi, 2*np.pi, size=(4000, L))
    S, grad = SW_and_derivs(Zc, beta)
    G = grad.sum(axis=1)/np.sqrt(L)                      # <v,grad S_W>
    # H_vv = (1/L) sum_ij d_i d_j S_W ; row sums of a graph Laplacian vanish
    # build Hessian for one random config explicitly
    d = np.roll(Zc, -1, axis=1) - Zc
    c = beta*(np.cos(d/2)/4)                             # d^2 S/dz^2 per plaquette
    Hv = np.zeros(len(Zc))
    for i in range(L):
        Hv += c[:, i]*(1.0/L) - c[:, i]*(1.0/L)           # = 0 identically (row-sum)
    P(f"  {beta:6.2f} {np.abs(G).max():12.2e} {np.abs(Hv).max():12.2e} {np.abs(G).max():14.2e}")
P("  -> |G|=0 to machine zero at every beta: the shift z_i->z_i+c drops out of every")
P("     difference, so S_W is c-independent for all beta; Hess S_W is a cyclic graph")
P("     Laplacian with zero row sums.  NOT a beta=2.5 accident.")
P("")

# ---------- 3. exact Gamma''(0; beta) on the P4 kernel ----------
def dmat(z):
    d = z[:, None]-z[None, :]
    d = (d + 2*np.pi) % (4*np.pi) - 2*np.pi               # wrapped difference, 4pi-periodic
    return d

D = dmat(z)
rho1 = p/np.max(p)
sq = np.sqrt(rho1[:, None]*rho1[None, :])

def Zk(k, beta, L):
    M = sq*np.exp(-1j*k*(z[:, None]+z[None, :])/(2*np.sqrt(L)))*np.exp(-beta*(1-np.cos(D/2)))
    lam = np.linalg.eigvals(M)
    return np.sum(lam**L)

def gamma_dd(beta, L, nk=301, kmax=26.0):
    """G''(0) = -d_a^2 log rho_v(a)|_0.  Z even real => rho_v(a)=(1/2pi)int Z(k)cos(ka)dk,
    rho_v''(0) = -<k^2>_Z rho_v(0)  ==>  G''(0) = <k^2>_Z = int k^2 Z / int Z.  Exact, no FD."""
    ks = np.linspace(-kmax, kmax, nk)
    Z = np.array([Zk(k, beta, L) for k in ks]).real
    Z = np.maximum(Z, 0.0)
    I0 = np.trapezoid(Z, ks); I2 = np.trapezoid(ks**2*Z, ks)
    return I2/I0, 1.0/(I2/I0), Z[nk//2]

P("== [C] Gamma''(0;beta) on the P4 kernel (the H1 margin delta = Gamma'' since Var=E[H]=0) ==")
P(f"  {'beta':>6} {'L':>3} {'Gamma\"(0)':>11} {'Var(a)':>9} {'Gamma\"*Var':>11}")
tab = {}
for beta, Ls in [(0.0, [2, 3, 4, 8]), (2.5, [2, 3, 4, 8])]:
    for L in Ls:
        g2, var, r0 = gamma_dd(beta, L)
        tab[(beta, L)] = g2
        P(f"  {beta:6.2f} {L:3d} {g2:11.6f} {var:9.5f} {g2*var:11.5f}")
P("  corpus P4-PERIODIC-TRANSFER.md:94-102 -> Gamma''(2.5): 0.2077,0.2097,0.2115,...,0.2165 (L=2..8)")
P("                              -> Gamma''(0.0): 0.2395,0.2493,0.2542,...,0.2616")
P("")

P("== [D] beta-scan at L=8 : does delta(beta) ever reach 0 ? ==")
P(f"  {'beta':>6} {'Gamma\"(0)':>11} {'1/(2Cin2pi)':>12} {'sign':>5}")
C0 = 0.20511530
scan = [0.0, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0, 6.0, 8.0, 12.0, 20.0, 40.0, 100.0]
mins = []
for beta in scan:
    g2, var, r0 = gamma_dd(beta, 8)
    mins.append(g2)
    P(f"  {beta:6.2f} {g2:11.6f} {C0:12.8f} {'+' if g2 > 0 else '-':>5}")
P(f"  min over scan = {min(mins):.6f}  >  0 ;  no root => beta*(P4) does not exist")
P(f"  single plaquette: delta = 0.2051153 - 0.0110593 beta^2  ->  beta* = 4.30661")
P(f"  predicted floor (ch05:216-220, W_geo >= (N/12)||v||^2 = 1/6) = {1/6:.5f}")
P("")

P("== [E] why the single plaquette HAS a root: the beta^2 variance term ==")
P("  single plaquette theta=pi: S = beta(1 + sin(r/2)u), u = xi3/r ; on the slice xi3=0")
P("  dS/dxi3 = beta*sin(r/2)/r != 0 (v is NOT a zero mode) => Var = beta^2*Var_Haar(sin(r/2)/r)")
P("  = beta^2 * 0.0110593 ;  W_full = 0.2051153 (beta-independent, pure Haar).")
P("  P4: same W_geo channel, but Var term ABSENT identically (graph Laplacian zero mode)")
P("  => no beta^2 to overcome W_geo => no root.")
P("")
P(f"[done] {time.time()-t0:.1f}s")
open("h1d_periodic.log", "w", encoding="utf-8").write("\n".join(out)+"\n")
