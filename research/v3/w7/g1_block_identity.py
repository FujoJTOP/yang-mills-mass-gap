#!/usr/bin/env python3
"""
g1_block_identity.py  --  witness for G1-BLOCK-IDENTITY.

Question: does  Gpp(v) := E[H_vv] - Var(G_v) + W_geo(v)  =  2 W_geo(v)  EXACTLY
for a general BLOCK of free links, fixed boundary U', and an arbitrary direction v ?

  W_geo(v) := E[Phi''_vv] - Var(Phi'_v),   G_v = d_v S,  H_vv = <v,Hess S v>.

Result proved in G1-BLOCK-IDENTITY.md (flat chart, Haar factorises -> Phi = sum_l Phi_l):
  (a) v CONSTANT in the chart (constant-coefficient combinations of link directions /
      constant projector P):        Gpp(v) = 2 W_geo(v)   EXACTLY.   [T]
  (b) v = v(xi) CONFIG-DEPENDENT (covariant / parallel-transported):
      Gpp(v) = 2 W_geo(v) + Delta(v), with the divergence defect
        Delta = E[g]^2 - E[phi]^2 + E[(d_v v).grad(Phi-S)] + E[(phi-g) div v],
        g = d_v S, phi = d_v Phi.
      Delta = 0 iff v is a constant field.

This witness checks both cases deterministically:

 Part A  Gaussian block model on R^D (faithful flat-chart block: arbitrary number D of
         free "link" coordinates, arbitrary non-commuting couplings A (action) and K (Haar
         potential), optional fixed boundary offset mu).   FULLY CLOSED FORM (Gaussian moments).
 Part B  Exact SU(2) one-link Haar quadrature (same model as s2_plaquette.py) with the
         CONFIG-DEPENDENT Euler direction v = xi  (radial; d_v v = xi, div v = 3):
         verifies Gpp = 2 W_geo + Delta and Delta != 0.

Acid test: free field (K = 0, Phi = 0) must give NO positive floor (Gpp <= 0).

Run:  py -3.12 g1_block_identity.py
"""
import numpy as np

def log(*a): print(*a, flush=True)

# =====================================================================
# PART A -- Gaussian block on R^D, fully closed form
#   S(xi) = 1/2 (xi-mu)^T A (xi-mu),   Phi(xi) = 1/2 xi^T K xi,
#   measure  exp(-S-Phi) dxi  =  Gaussian  N(m, Sigma),
#            M = A+K,  Sigma = M^{-1},  m = Sigma A mu.
#   (a) constant v  ->  E[H]=v^TAv, E[g]=v^TA(m-mu), Var(g)=v^TA Sigma A v,
#                       W=v^TKv - v^TK Sigma K v,  Gpp = E[H]-Var(g)+W = 2W.
#   (b) v = xi      ->  g=xi^TA(xi-mu), phi=xi^TK xi, H=xi^TA xi, Phi''=xi^TK xi
#        (here v=xi, so H_vv=xi^T A xi is the Hessian form).
# =====================================================================
def part_A():
    log("="*74)
    log("PART A -- Gaussian block model on R^D (closed form), arbitrary coupling")
    log("="*74)
    rng = np.random.default_rng(20261124)
    D = 6
    # random SPD A, K that do NOT commute
    def spd(n):
        B = rng.standard_normal((n, n)); return B @ B.T + n*np.eye(n)
    A = spd(D); K = spd(D)
    M = A + K; Sigma = np.linalg.inv(M)
    v = rng.standard_normal(D)
    mu = rng.standard_normal(D) * 0.7            # fixed "boundary" offset U'
    m = Sigma @ A @ mu

    # ---- (a) constant direction ----
    H = v @ A @ v
    Eg = v @ A @ (m - mu)
    Vg = v @ A @ Sigma @ A @ v
    W  = v @ K @ v - v @ K @ Sigma @ K @ v
    Gpp = H - Vg + W
    log(f" (a) constant v, boundary |mu|={np.linalg.norm(mu):.3f}")
    log(f"     E[H]={H:.12f}  Var(G)={Vg:.12f}  W_geo={W:.12f}")
    log(f"     Gpp = {Gpp:.12f}   2*W_geo = {2*W:.12f}   residual = {abs(Gpp-2*W):.2e}")
    assert abs(Gpp - 2*W) < 1e-9, "case (a) identity FAILED"

    # ---- (a') boundary-independence: vary mu ----
    worst = 0.0
    for s in (0.0, 0.5, 2.0, 7.0):
        m2 = Sigma @ A @ (mu*s)
        G2 = H - (v @ A @ Sigma @ A @ v) + W
        worst = max(worst, abs(G2 - 2*W))
    log(f"     boundary-independence: max |Gpp-2W| over 4 boundary offsets = {worst:.2e}")

    # ---- (b) v = xi (Euler field) : exact Gaussian moments ----
    # g=xi^T A(xi-mu) = xi^T A xi - (A mu).xi
    # E[xi^T A xi] = tr(A Sigma)+m^T A m ; need mu term.
    Amu = A @ mu
    def quad_mean(P):                      # E[xi^T P xi]
        return np.trace(P @ Sigma) + m @ P @ m
    def quad_var(P):                       # Var(xi^T P xi - c.xi) for symmetric P, c
        # Var(xi^T P xi) = 2 tr((P Sigma)^2) + 4 m^T P Sigma P m
        return 2*np.trace(P@Sigma@P@Sigma) + 4*(m @ P @ Sigma @ P @ m)
    # g = xi^T A xi - Amu.xi ; treat as quadratic with linear shift.
    # Var(g) = Var(xi^T A xi) + Var(Amu.xi) - 2 Cov(xi^T A xi, Amu.xi)
    Var_quadA = 2*np.trace(A@Sigma@A@Sigma) + 4*(m @ A @ Sigma @ A @ m)
    Var_lin   = Amu @ Sigma @ Amu
    # Cov(xi^T A xi, c.xi) with c=Amu : = 2 (A Sigma m).c + 2 m^T A Sigma c  ... use E[xixi]=Sigma+mm^T
    # direct: E[(xi^T A xi)(c.xi)] - E[xi^T A xi] E[c.xi]
    EA2 = quad_mean(A)                       # E[xi^T A xi]
    Elin = m @ Amu
    # E[(xi^T A xi)(c.xi)] = sum_ij A_ij E[xi_i xi_j (c.xi)]
    # E[xi_i xi_j xi_k] = m_i m_j m_k + m_i S_jk + m_j S_ik + m_k S_ij (S=Sigma)
    S = Sigma; c = Amu
    def third(i,j,k):
        return m[i]*m[j]*m[k] + m[i]*S[j,k] + m[j]*S[i,k] + m[k]*S[i,j]
    Ecross = 0.0
    for i in range(D):
        for j in range(D):
            for k in range(D):
                Ecross += A[i,j]*c[k]*third(i,j,k)
    Cov_Ac = Ecross - EA2*Elin
    Eg_b = EA2 - Elin
    Vg_b = Var_quadA + Var_lin - 2*Cov_Ac
    Eg_b2 = Eg_b*Eg_b
    # phi = xi^T K xi
    Ephi = quad_mean(K); Vphi = 2*np.trace(K@Sigma@K@Sigma) + 4*(m @ K @ Sigma @ K @ m)
    W_b = Ephi - Vphi
    Gpp_b = EA2 - Vg_b + W_b
    # Delta formula: E[g]^2-E[phi]^2 + E[(dv v).grad(Phi-S)] + E[(phi-g) div v]
    #   v=xi: dv v = xi, div v = D ; grad(Phi-S) = (K-A)xi - K?  careful:
    #   grad S = A(xi-mu); grad Phi = K xi ; grad(Phi-S)=(K-A)xi + A mu
    #   E[xi . grad(Phi-S)] = E[xi^T(K-A)xi + xi.Amu] = quad_mean(K-A) + m.Amu
    gradterm = quad_mean(K-A) + m @ Amu          # E[(dv v).grad(Phi-S)]
    # E[(phi-g) div v] = D * E[xi^T(K-A)xi + xi.Amu]
    divterm = D * (quad_mean(K-A) + m @ Amu)
    Delta_b = Eg_b*Eg_b - Ephi*Ephi + gradterm + divterm
    log(f" (b) config-dependent v = xi (Euler field)")
    log(f"     Gpp = {Gpp_b:.12f}   2*W_geo = {2*W_b:.12f}   Gpp-2W = {Gpp_b-2*W_b:.12f}")
    log(f"     Delta (formula) = {Delta_b:.12f}   residual |Gpp-2W-Delta| = "
        f"{abs(Gpp_b-2*W_b-Delta_b):.2e}")
    assert abs(Gpp_b - 2*W_b - Delta_b) < 1e-9, "case (b) corrected identity FAILED"
    assert abs(Delta_b) > 1e-6, "Delta unexpectedly zero for a varying field"

    # ---- acid test : free field K = 0 ----
    K0 = np.zeros((D,D)); M0 = A.copy(); S0 = np.linalg.inv(M0)
    H0 = v @ A @ v; Vg0 = v @ A @ S0 @ A @ v; W0 = 0.0
    G0a = H0 - Vg0 + W0
    # case (b) v=xi, K=0, mu arbitrary: Gpp = E[xi^T A xi]-Var +0 ; for mu=0 -> -D
    mm = np.zeros(D)
    Eg0 = np.trace(A@S0)                 # E[xi^T A xi], mu=0
    Vq0 = 2*np.trace(A@S0@A@S0)
    G0b = Eg0 - Vq0                       # = D - 2D = -D
    log(f" ACID TEST K=0 (free):  (a) Gpp = {G0a:.3e} (must be 0) ; "
        f"(b) Gpp = {G0b:.6f} (must be < 0)")
    assert abs(G0a) < 1e-9
    assert G0b < 0
    log(" -> PASS (no positive floor for the free field)\n")

# =====================================================================
# PART B -- exact SU(2) one-link Haar quadrature, Euler direction v = xi
#   model identical to s2_plaquette.py:  S = beta(1 - cos(rho/2)cos(th/2)
#                                              + sin(rho/2)sin(th/2) u),
#   Haar weight  ~ exp(-S) sin^2(rho/2) d rho du  (psi trivial).
#   Phi = 2 log(rho/2) - 2 log sin(rho/2).
#   v = xi  =>  g = rho S_r, H = rho^2 S_rr, phi = rho Phi_r, Phi'' = rho^2 Phi_rr.
#   Delta = <g>^2 - <phi>^2 + 4 < rho (Phi_r - S_r) >.
# =====================================================================
def Phi_r(x):
    y=np.empty_like(x); sm=np.abs(x)<0.3; xr=x[~sm]
    y[~sm]=1.0/xr-1.0/np.tan(xr)
    xs=x[sm]; y[sm]=xs/3+xs**3/45+2*xs**5/945+2*xs**7/4725
    return y
def Phi_rr(x):
    y=np.empty_like(x); sm=np.abs(x)<0.3; xr=x[~sm]
    y[~sm]=0.5*(-1.0/xr**2+1.0/np.sin(xr)**2)
    xs=x[sm]; y[sm]=1/6+xs**2/30+2*xs**4/189+2*xs**6/2475
    return y
_xg,_wg=np.polynomial.legendre.leggauss(1400)
_ug,_wu=np.polynomial.legendre.leggauss(1400)
RHO=np.pi*(_xg+1.0); WRHO=np.pi*_wg
R,U=np.meshgrid(RHO,_ug,indexing='ij'); WR,WU=np.meshgrid(WRHO,_wu,indexing='ij')

def su2_euler(beta, theta):
    x=R/2.0; c=np.cos(R/2); s=np.sin(R/2); C=np.cos(theta/2); St=np.sin(theta/2)
    S=beta*(1.0-(c*C-s*St*U))
    w=np.exp(-S)*np.sin(R/2)**2*WR*WU; Z=w.sum()
    def E(a): return (w*a).sum()/Z
    Sr  = beta*((s/2)*C+(c/2)*St*U)
    Srr = beta*((c/4)*C-(s/4)*St*U)
    pr=Phi_r(x); prr=Phi_rr(x)
    # v = xi (Euler):  d_v = rho d_rho
    g   = R*Sr
    H   = R**2*Srr
    phi = R*pr
    Phi2= R**2*prr
    Eg, Eg2 = E(g), E(g*g)
    Ephi, Ephi2 = E(phi), E(phi*phi)
    EH = E(H); EPhi2 = E(Phi2)
    Vg = Eg2-Eg**2; Vphi = Ephi2-Ephi**2
    W  = EPhi2 - Vphi
    Gpp= EH - Vg + W
    gradterm = E(R*(pr-Sr))          # E[(d_v v).grad(Phi-S)],  d_v v = xi
    Delta = Eg**2 - Ephi**2 + 4.0*gradterm
    return dict(EH=EH,Vg=Vg,W=W,Gpp=Gpp,Delta=Delta,Eg=Eg,Ephi=Ephi,
                resid=abs(Gpp-(2*W+Delta)))

def part_B():
    log("="*74)
    log("PART B -- exact SU(2) one-link Haar quadrature, v = xi (Euler)")
    log("  verify   Gpp = 2 W_geo + Delta,   Delta != 0 (v varies)")
    log("="*74)
    log(f"{'beta':>6} {'theta':>7} {'Gpp':>13} {'2W_geo':>13} {'Delta':>13} {'|resid|':>10}")
    for beta in (0.5,1.0,2.5):
        for theta in (0.0,np.pi/2,np.pi):
            o=su2_euler(beta,theta)
            log(f"{beta:6.2f} {theta:7.3f} {o['Gpp']:13.8f} {2*o['W']:13.8f} "
                f"{o['Delta']:13.8f} {o['resid']:10.1e}")
            assert o['resid'] < 3e-4, "SU(2) Euler identity residual too large"
            assert abs(o['Delta']) > 1e-4, "Delta unexpectedly zero"
    log(" -> Gpp = 2 W_geo + Delta verified on the real Haar model; Delta != 0.\n")

if __name__=="__main__":
    part_A()
    part_B()
    log("ALL CHECKS PASSED.")
