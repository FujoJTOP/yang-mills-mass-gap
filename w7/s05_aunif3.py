import numpy as np
TW=np.pi*2; B=2.5; TH=np.pi/2; C=np.cos(TH/2); T=np.sin(TH/2); K0=0.5*C/T
xx,wg=np.polynomial.legendre.leggauss(6000); tt=0.5*(xx+1); w0=0.5*wg
D1C=-T/8-C**2/(8*T); D3C=T/480; D1Q=T/320+C**2/(192*T); D3Q=-11*T/107520
def dc(r):
    r=np.asarray(r,float); sm=r<0.02; rr=np.where(sm,1.0,r); s=np.sin(0.5*rr); co=np.cos(0.5*rr)
    a1=3*T*(co/(2*rr**2)-s/rr**3)-(C**2/(4*T))*s/rr
    a3=T*(3*s/rr**5-1.5*co/rr**4-s/(4*rr**3))
    return np.where(sm,D1C+D1Q*r**2,a1), np.where(sm,D3C+D3Q*r**2,a3)
def wf(r):
    r=np.asarray(r,float); o=np.empty_like(r); m=r!=0; o[m]=np.sin(0.5*r[m])/r[m]; o[~m]=0.5; return o
def wpf(r):
    rr=np.where(r==0,1.0,r); return (0.5*rr*np.cos(0.5*rr)-np.sin(0.5*rr))/rr**2
print(f"{'a':>7} {'E[G]':>13} {'Var(G)':>13} {'E[Delta]':>14} {'beta|ED|':>13} {'M(a)':>13} {'gap':>12}")
for av in [0.0,0.1,0.25,0.5,0.75,1.0,1.2,1.25,1.2886819904,1.3,1.5]:
    lo=abs(av); rho=lo+(TW-lo)*tt; wt=(TW-lo)*w0
    dens=wf(rho)*np.exp(-B*(1-C*np.cos(0.5*rho)+T*av*wf(rho))); p=wt*dens; p/=p.sum()
    G=B*((C/2)*(av/rho)*np.sin(0.5*rho)+T*(wf(rho)+(av*av/rho)*wpf(rho)))
    d1,d3=dc(rho); D=av*d1+av**3*d3
    EG=np.sum(p*G); VG=np.sum(p*(G-EG)**2); ED=np.sum(p*D); M=K0*EG-VG
    print(f"{av:7.5f} {EG:+.6e} {VG:.6e} {ED:+.6e} {B*abs(ED):.6e} {M:+.6e} {M-B*abs(ED):+.6e}")
