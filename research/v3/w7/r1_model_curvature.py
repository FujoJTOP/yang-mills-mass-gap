#!/usr/bin/env python3
"""
r1_model_curvature.py -- RECONCILIATION of the G1/G3 object vs the true coarse curvature G''
for (star), in the exactly solvable SU(2) one-plaquette / one-free-link model of S2.

Model (S2):  S(U;B)=beta(1-1/2 tr(UB)),  dmu(.|B) proportional to e^{-S} dU (Haar),
  U=exp(xi), xi in R^3, rho=|xi| in [0,2pi],  Phi=-log J = 2log(rho/2)-2log sin(rho/2).

Three objects are separated and cross-checked:
  (U)  G''_0  := E_mu[H]-Var_mu(G) + W_geo,mu   (UNCONDITIONED, no X)  = 2 W_geo,mu   [G1/G3/S2]
  (Xtr)X_mu   := -2 Cov_mu(G,Phi'),   Delta_dyn,mu=E_mu[H]-Var_mu(G),  W_geo,mu=E_mu[Phi'']-Var_mu(Phi')
  (C)  G''(a0):= -d^2/da^2 log rho_v(a)|_{a0} = E_a[d_v^2 V|a0] - Var_a(d_v V|a0)      [true (star)]
                = Delta_dyn,a + W_geo,a + X_a   (all moment conditioned on a=<v,xi>=a0)

We condition on the coarse field a0 along v=e3 and show, at theta=pi (B->-I) where the
slice measure is beta-independent Haar, the EXACT decomposition:
      G''(e3; theta=pi, a0=0) = E_Haar[Phi'(r)/r]  -  beta^2 Var_Haar[sin(r/2)/r]
                              =   W_geo,a         +   Delta_dyn,a          (X_a = 0),
with Delta_dyn,a = -beta^2 Var (the Wilson force variance, Savvidy), W_geo,a=E_Haar[Phi'/r]>0.
Run:  py -3.12 r1_model_curvature.py
"""
import numpy as np

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

# ---------------- (U) UNCONDITIONED: 3D quadrature in (rho,u) (as S2/G4) ----------------
_xg,_wg = np.polynomial.legendre.leggauss(800)
RHO=np.pi*(_xg+1.0); WRHO=np.pi*_wg
_ug,_wu=np.polynomial.legendre.leggauss(800)
R,U=np.meshgrid(RHO,_ug,indexing='ij'); WR,WU=np.meshgrid(WRHO,_wu,indexing='ij')

def uncond(beta,theta):
    x=R/2.0
    c=np.cos(R/2); s=np.sin(R/2); C=np.cos(theta/2); St=np.sin(theta/2)
    S=beta*(1.0-(c*C-s*St*U))
    w=np.exp(-S)*np.sin(R/2)**2*WR*WU
    Z=w.sum(); E=lambda a:(w*a).sum()/Z
    Sr=beta*((s/2)*C+(c/2)*St*U); Su=beta*s*St
    Srr=beta*((c/4)*C-(s/4)*St*U); Sru=beta*(c/2)*St; Suu=0.0
    pr=Phi_rho(x); prr=Phi_rr(x)
    # v=e3
    G3=U*Sr+(1-U**2)/R*Su
    H3=(Srr*U**2+2*Sru*U*(1-U**2)/R+Suu*((1-U**2)/R)**2+Sr*(1-U**2)/R+Su*(-3*U*(1-U**2)/R**2))
    P3=U*pr; PP3=U**2*prr+(1-U**2)/R*pr
    EG=E(G3); VG=E(G3*G3)-EG**2; EH=E(H3)
    EP=E(P3); VP=E(P3*P3)-EP**2; EPP=E(PP3)
    Ddyn=EH-VG; Wg=EPP-VP; X=-2*(E(G3*P3)-EG*EP)
    return Ddyn,Wg,X

# ---------------- (C) TRUE coarse curvature: marginal of a=xi_3 ----------------
_rg,_rgw=np.polynomial.legendre.leggauss(1600)
def marginal_e3(a0,beta,theta,n=1600):
    Rmax=np.sqrt(max(4*np.pi**2-a0*a0,1e-12))
    xg,wg=np.polynomial.legendre.leggauss(n)
    r=0.5*Rmax*(xg+1.0); wr=0.5*Rmax*wg
    rho=np.sqrt(r*r+a0*a0); u=a0/rho
    h=beta*(1.0-np.cos(rho/2)*np.cos(theta/2))+(2*np.log(rho/2)-2*np.log(np.sin(rho/2)))
    g=beta*np.sin(rho/2)*np.sin(theta/2)
    V=h+g*u
    return 2*np.pi*np.sum(np.exp(-V)*r*wr)

def Gpp_fd(beta,theta,a0=0.0,da=1e-3):
    vals=[marginal_e3(a0+k*da,beta,theta) for k in (-2,-1,0,1,2)]
    logs=np.array([np.log(v) for v in vals])
    fpp=(-logs[0]+16*logs[1]-30*logs[2]+16*logs[3]-logs[4])/(12*da*da)
    return -fpp

# ---------------- (C) exact slice decomposition at a0=0, v=e3 ----------------
# On the slice xi_3=0: rho=r, u=0, d3^2 = h'(r)/r, d3 = g(r)/r, Phi'_3=0.
# Slice weight w = e^{-h} * r dr  (e^{-h}=e^{-S}J already includes Haar J).
def slice_decomp(beta,theta,n=1600):
    xg,wg=np.polynomial.legendre.leggauss(n)
    r=np.pi*(xg+1.0); wr=np.pi*wg
    h=beta*(1.0-np.cos(r/2)*np.cos(theta/2))+(2*np.log(r/2)-2*np.log(np.sin(r/2)))
    # d/dr h at a0=0 (=hpr): action part + Phi part
    S_pr=beta*0.5*np.sin(r/2)*np.cos(theta/2)
    Ph_pr=Phi_rho(r/2)
    g=beta*np.sin(r/2)*np.sin(theta/2)
    w=np.exp(-h)*r*wr
    Z=w.sum(); E=lambda a:(w*a).sum()/Z
    # H_vv = S33 = S_pr/r ; Phi''=Ph_pr/r ; G_v = g/r ; Phi'_3 = 0
    H=S_pr/r; PH=Ph_pr/r; G=g/r
    EG=E(G); VG=E(G*G)-EG**2
    EPH=E(PH); VPH=0.0; EPHp=E(PH*PH)-EPH**2  # Phi'_3=0 so Var(Phi')=0
    Ddyn_a=E(H)-VG
    Wg_a=EPH-VPH
    X_a=-2.0*(E(G*0.0)-EG*0.0)   # = -2Cov(G,Phi'_3)=0 exactly
    Gpp=Ddyn_a+Wg_a+X_a
    # also G'' = E[h'/r]-Var(g/r) directly
    hp=S_pr+Ph_pr
    Gpp_direct=E(hp/r)-VG
    return dict(Gpp=Gpp,Gpp_direct=Gpp_direct,Ddyn_a=Ddyn_a,Wg_a=Wg_a,X_a=X_a,
                VG=VG,EPH=EPH,VarPhi=E(PH*PH)-EPH**2)

def log(*a): print(*a,flush=True)
N=2; NFLOOR=N/12.0

log("="*84)
log("RECONCILIATION r1 -- SU(2) one-plaquette model; object (U)=G1/G3 vs (C)=true coarse G''")
log(f"N/12 = {NFLOOR:.6f}")
log("="*84)

log("\n(A) theta=pi (B->-I), v=e3, a0=0 : EXACT slice decomposition  G''=Ddyn_a+Wg_a+X_a")
log(f"{'beta':>6} {'Gpp_exact':>12} {'Gpp_fd':>12} {'Ddyn_a':>11} {'Wg_a':>11} {'X_a':>8} {'2W_uncond':>11} {'Gpp0':>10}")
for beta in (0.0,1.0,2.5,4.0,4.09,5.0,10.0,20.0):
    d=slice_decomp(max(beta,1e-9),np.pi)
    fd=Gpp_fd(max(beta,1e-6),np.pi)
    Dd,Wg,X=uncond(max(beta,1e-9),np.pi)
    log(f"{beta:6.2f} {d['Gpp']:12.6f} {fd:12.6f} {d['Ddyn_a']:11.6f} {d['Wg_a']:11.6f} {d['X_a']:8.1e} {2*Wg:11.6f} {Dd+Wg:10.6f}")
log("   Gpp0 = Ddyn_mu+W_geo_mu (=2W_geo_mu, the G1/G3 object) ; 2W_uncond = same.")

log("\n   beta=0 CHECK vs ch05 theorem G''>=N/12 and S2's W_geo=0:")
d0=slice_decomp(1e-9,np.pi); Dd0,Wg0,X0=uncond(1e-9,np.pi)
log(f"     true (C) G''(e3;pi,0) = {d0['Gpp']:.6f}   ( >= N/12={NFLOOR:.4f} ?  {d0['Gpp']>=NFLOOR} )")
log(f"     decomposed: Ddyn_a={d0['Ddyn_a']:.2e}  Wg_a={d0['Wg_a']:.6f}  X_a={d0['X_a']:.2e}")
log(f"     UNCOND (U): W_geo_mu={Wg0:.2e}  G''_0=2W_geo_mu={2*Wg0:.2e}   <-- S2's 'W_geo=0 at beta=0'")

log("\n(B) theta=0 (B=I): unconditioned vs conditioned (should both be positive, G''!=2W_geo)")
log(f"{'beta':>6} {'G''_C(slice)':>13} {'2W_mu':>10} {'Ddyn_a':>10} {'Wg_a':>10}")
for beta in (0.5,2.5,5.0):
    d=slice_decomp(beta,0.0); Dd,Wg,X=uncond(beta,0.0)
    log(f"{beta:6.2f} {d['Gpp']:13.6f} {2*Wg:10.6f} {d['Ddyn_a']:10.6f} {d['Wg_a']:10.6f}")

log("\n(C) UNCONDITIONED identity check:  Delta_dyn,mu = W_geo,mu ; X_mu = -2 W_geo,mu ; sum = 0")
log(f"{'beta':>6} {'th':>6} {'Ddyn_mu':>11} {'Wg_mu':>11} {'X_mu':>11} {'Ddyn+Wg+X':>11}")
for beta in (0.5,2.5,5.0):
    for th in (0.0,np.pi/2,np.pi):
        Dd,Wg,X=uncond(beta,th)
        log(f"{beta:6.2f} {th:6.3f} {Dd:11.6f} {Wg:11.6f} {X:11.6f} {Dd+Wg+X:11.2e}")

log("\n(D) EXACT beta-law at theta=pi,e3,a0=0 (theta=pi => slice measure beta-independent Haar):")
b=np.array([0.0,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0])
G=np.array([slice_decomp(max(t,1e-9),np.pi)['Gpp'] for t in b])
A=np.polyfit(b*b,G,1)
log(f"     G''(beta) = {A[1]:.6f} + ({A[0]:.6f})*beta^2   [expect 0.2057 - 0.0123 beta^2]")
log(f"     zero crossing beta* = {np.sqrt(-A[1]/A[0]):.4f}   [expect ~4.09]")
log(f"     values: " + "  ".join(f"{t:g}:{g:+.4f}" for t,g in zip(b,G)))
log("\nClay 0%. No QED.")
