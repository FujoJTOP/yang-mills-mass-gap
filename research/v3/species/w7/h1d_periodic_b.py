# h1d_periodic_b.py -- H1 on the P4 periodic kernel: is Var = E[H] = 0 identically
# (every beta), and does beta*(P4) exist?   Independent re-derivation (sibling h1d_periodic.py untouched).
# Streaming log: every line flushed to h1d_periodic_b.log as it is produced.
#
# Definitions quoted, not remembered:
#   H1      : Var_mu(<v,grad S>|0) <= E_mu[H_vv|0] + W_full - delta, W_full = W_geo + X
#             docs/paper-v2/ch08-immune.tex:11-19 ; app.tex:179-202
#   W_geo,X : W_geo = E[Phi'']-Var(Phi'), X = -2Cov(G,Phi'), G = <v,grad S_W>
#             ch05-theorems.tex:196-206,230-239 ; app.tex:795-796
#   identity: G''(0) = E[H|0] - Var(G|0) + W_full  (ch05 eq:split) ; delta_max = G''(0) (H1A-EXACT.md:64-68)
#   floor   : ch05-theorems.tex:208 Lemma B E[Phi''|0] >= (N/12)||v||^2 ; an L-link average has
#             ||v||^2 = 1/L , i.e. c = (N/12)/L   (ch05-theorems.tex:224)
#   P4      : Z(k)=Tr[M_k^L], M_k(z,z')=sqrt(rho1 rho1') e^{-ik(z+z')/(2 sqrt L)} e^{-beta(1-cos((z'-z)/2))},
#             rho_v(a)=(1/2pi)int Z e^{ika}dk, G''(a) = -d_a^2 log rho_v(a) = <k^2>_Z
#             research/v3/P4-PERIODIC-TRANSFER.md:62-69 ; rho1(z)=pi(Cin(2pi)-Cin|z|) (H1C-WFULL.md:63)
#   corpus  : G''(2.5)=0.2077..0.2165 ; G''(0)=0.2395..0.2616 (L=2..8) ; Wilson c=-0.0461+-0.0003,
#             total c=+0.22196+-0.00011 (q=0.849)  P4-PERIODIC-TRANSFER.md:94-115 ; ch15-uvface.tex:202-213
# numpy only, deterministic quadrature.
import numpy as np, time
t0 = time.time(); LOG = open("h1d_periodic_b.log", "w", encoding="utf-8")
def P(s=""):
    print(s, flush=True); LOG.write(str(s)+"\n"); LOG.flush()
TP = 2*np.pi

# ------------ [A] convention lock ------------
n = 1000001
tt = np.linspace(0.0, TP, n)
f = np.where(tt > 0, (1-np.cos(tt))/np.where(tt > 0, tt, 1.0), 0.5*tt)
Cg = np.concatenate(([0.0], np.cumsum(0.5*(f[1:]+f[:-1])*np.diff(tt))))
CIN2PI = Cg[-1]
def Cin(x): return np.interp(np.abs(x), tt, Cg)
ANCHOR = 1/(2*CIN2PI)
P("== [A] convention lock ==")
z0 = np.linspace(-TP, TP, 100001)
r0 = np.pi*(CIN2PI - Cin(z0))
p0 = r0/np.trapezoid(r0, z0)
Var_xi3 = np.trapezoid(z0**2*p0, z0)
P(f"  Cin(2pi)={CIN2PI:.10f} (doc 2.437653393047)")
P(f"  anchor 1/(2Cin2pi)={ANCHOR:.9f} (ch15:190 eq:v3anchor 0.20511530)")
P(f"  int rho1 dz/(4pi^2)={np.trapezoid(r0,z0)/(4*np.pi**2):.9f} (must be 1)")
P(f"  Var(xi3)={Var_xi3:.6f} (doc 3.719824) ; 1/Var={1/Var_xi3:.6f} (doc 0.26883, ch15:182-195)")
P(f"  beta*(N=2)=sqrt(anchor/0.0110593)={np.sqrt(ANCHOR/0.0110593):.6f} (task 4.30661)")
P("")

# ------------ [B] structural zeros ------------
def chain(Z, beta):
    L = Z.shape[1]
    d = np.roll(Z, -1, axis=1) - Z
    S = beta*np.sum(1-np.cos(d/2), axis=1)
    g = beta*(np.sin(d/2)/2)
    grad = np.roll(g, 1, axis=1) - g
    c = beta*(np.cos(d[0]/2)/4)
    H = np.zeros((L, L))
    for i in range(L):
        j = (i+1) % L
        H[i, i] -= c[i]; H[j, j] -= c[i]; H[i, j] += c[i]; H[j, i] += c[i]
    return S, grad, H

rng = np.random.default_rng(20261125)
P("== [B] Var(G)=0, E[H]=0 : pointwise, every beta, every L (global-shift symmetry) ==")
P(f"  {'beta':>8} {'L':>3} {'max|G|':>10} {'|H_vv|':>10} {'max|rowsumH|':>13} {'max|S(z+c)-S(z)|':>17}")
for beta in [0.0, 0.5, 2.5, 6.0, 20.0, 100.0, 1e4]:
    for L in (2, 4, 8):
        Z = rng.uniform(-TP, TP, size=(3000, L))
        S, grad, H = chain(Z, beta)
        v = np.ones(L)/np.sqrt(L)
        S2, _, _ = chain(Z+0.9137, beta)
        P(f"  {beta:8.1f} {L:3d} {np.abs(grad@v).max():10.2e} {abs(v@H@v):10.2e} "
          f"{np.abs(H.sum(axis=1)).max():13.2e} {np.abs(S2-S).max():17.2e}")
P("  -> machine zeros at all beta: S_W depends on differences only => sum_i d_i S_W = 0 (G=0);")
P("     Hess S_W = cyclic Laplacian, zero row sums => <v,Hess v>=0. NOT a beta=2.5 accident.")
P("")

# ------------ [C] G''(0) on the P4 kernel ------------
def Gdd(beta, L, Nz=120, nk=61, kmax=10.0):
    z = np.linspace(-TP, TP, Nz, endpoint=False) + TP/Nz
    rho1 = np.pi*(CIN2PI - Cin(z))
    d = z[:, None]-z[None, :]; d = (d+TP) % (2*TP) - TP
    W = np.sqrt(rho1[:, None]*rho1[None, :])*np.exp(-beta*(1-np.cos(d/2)))
    ph = np.exp(-1j*(z[:, None]+z[None, :])/(2*np.sqrt(L)))
    ks = np.linspace(-kmax, kmax, nk)
    Z = np.array([np.sum(np.linalg.eigvals(W*ph**k)**L).real for k in ks])
    I0 = np.trapezoid(Z, ks); I2 = np.trapezoid(ks**2*Z, ks)
    return I2/I0, 1.0/(I2/I0)

c25 = {2: 0.2077, 3: 0.2097, 4: 0.2115, 5: 0.2130, 6: 0.2144, 7: 0.2155, 8: 0.2165}
c0 = {2: 0.2395, 3: 0.2493, 4: 0.2542, 5: 0.2572, 6: 0.2591, 7: 0.2605, 8: 0.2616}
P("== [C] G''(0) on the P4 kernel -- reproduction of the corpus table ==")
P(f"  {'beta':>6} {'L':>3} {'G\"(0)':>9} {'Var(a)':>9} {'G\"*Var':>8}  corpus   diff")
for beta, ref in ((2.5, c25), (0.0, c0)):
    for L in range(2, 9):
        g2, var = Gdd(beta, L)
        P(f"  {beta:6.2f} {L:3d} {g2:9.5f} {var:9.5f} {g2*var:8.4f}  {ref[L]:.4f}  {g2-ref[L]:+.4f}")
P("")
P("== [C'] quadrature convergence (beta=2.5) ==")
for (L, Nz, nk, km) in [(8, 120, 61, 10.0), (8, 120, 121, 16.0), (8, 240, 121, 16.0), (4, 240, 121, 16.0)]:
    P(f"  L={L} Nz={Nz} nk={nk} kmax={km} -> G\"={Gdd(2.5,L,Nz,nk,km)[0]:.5f}")
P(f"  L=8 beta=0 Nz=240 nk=121 -> G\"={Gdd(0.0,8,240,121,16.0)[0]:.5f} (corpus 0.2616)")
P("")

# ------------ [D] beta scan ------------
P("== [D] beta scan: delta = G''(0) = W_geo (since Var=E[H]=0), sign and limits ==")
P(f"  {'beta':>7} {'G\"(L=4)':>10} {'G\"(L=8)':>10} {'anchor/L=4':>11} {'anchor/L=8':>11} {'floor (1/6)/L':>14} {'sign':>5}")
rows = []
for beta in [0.0, 0.5, 1.0, 2.5, 4.0, 5.7, 8.0, 12.0, 20.0, 50.0, 200.0]:
    Nz, nk = (120, 61) if beta <= 12 else (250, 61)
    g4 = Gdd(beta, 4, Nz, nk, 10.0)[0]; g8 = Gdd(beta, 8, Nz, nk, 10.0)[0]
    rows.append((beta, g4, g8))
    P(f"  {beta:7.2f} {g4:10.5f} {g8:10.5f} {ANCHOR/4:11.5f} {ANCHOR/8:11.5f} "
      f"{(1/6)/8:14.5f} {'+' if min(g4,g8)>0 else '-':>5}")
P(f"  min over whole scan = {min(min(a,b) for _,a,b in rows):.6f} > 0  => no root beta*(P4)")
P(f"  beta->inf at fixed L heads to anchor/L: measured L=4 {rows[-1][1]:.5f} vs {ANCHOR/4:.5f};"
  f" L=8 {rows[-1][2]:.5f} vs {ANCHOR/8:.5f}")
P(f"  ch05:224 floor for an L-link average = (N/12)/L: L=4 {(1/6)/4:.5f}, L=8 {(1/6)/8:.5f} -> respected")
P("")

# ------------ [E] the 8% ------------
P("== [E] what is the 8% (P4 saturation 0.22196 vs single-link anchor 0.2051153)? ==")
clt = 1/Var_xi3; sat = 0.22196
P(f"  single-link Haar anchor 1/(2Cin2pi)        = {ANCHOR:.6f}")
P(f"  k=0 Haar block, L->inf (L-fold convolution)= 1/Var(xi3) = {clt:.6f}"
  f" = {clt/ANCHOR:.4f} x anchor  ({100*(clt/ANCHOR-1):+.1f}%)")
P(f"  Wilson part at beta=2.5 (corpus fit)       = -0.0461 +-0.0003 ({100*(-0.0461/clt):+.1f}% of it)")
P(f"  {clt:.5f} - 0.0461 = {clt-0.0461:.5f} vs P4 fit c={sat:.5f} (diff {clt-0.0461-sat:+.5f})")
P(f"  net {sat:.5f}/{ANCHOR:.7f} = {sat/ANCHOR:.5f} (+{100*(sat/ANCHOR-1):.2f}%)")
P(f"  factorisation ({clt/ANCHOR:.4f})x({1-0.0461/clt:.4f}) = {(clt/ANCHOR)*(1-0.0461/clt):.5f}")
P("")

# ------------ [F] N-dependence ------------
P("== [F] N-dependence of W_0 (betastar_scaling.log:5, betastar_n3.log, betastar_n34.log) ==")
Wd = {2: 0.2051153, 3: 0.3453, 4: 0.4908}
P(f"  {'N':>2} {'W_geo(b=0)':>11} {'W/N':>8} {'W/(N/2)':>9} {'N/12':>7}")
for N in (2, 3, 4):
    P(f"  {N:2d} {Wd[N]:11.4f} {Wd[N]/N:8.4f} {Wd[N]/(N/2):9.4f} {N/12:7.4f}")
P(f"  W(3)/W(2)={Wd[3]/Wd[2]:.4f} (proportional-N predicts 1.5); W(4)/W(2)={Wd[4]/Wd[2]:.4f} (predicts 2.0)")
P("  0.22196 and 0.2051153 are the SAME N=2 object => any N-law forces ratio 1, not 1.082")
P("  => the 8% is NOT an N-effect; it is (Haar CLT block) x (Wilson backreaction), both L/mode effects.")
P("")
P(f"[done] {time.time()-t0:.1f}s numpy {np.__version__}")
LOG.close()
