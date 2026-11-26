import numpy as np, math
from numpy.polynomial.legendre import leggauss
B=2.5
def run(th,N=3000,M=600):
    ct=math.cos(th/2); st=math.sin(th/2)
    x,wx=leggauss(N); y,wy=leggauss(M)
    psi=0.5*math.pi*(x+1); wpsi=0.5*math.pi*wx
    u=y; wu=wy
    PS,UU=np.meshgrid(psi,u,indexing="ij"); W=np.outer(wpsi,wu)
    rho=2*PS; s2=np.sin(rho/2); c2=np.cos(rho/2)
    S=B*(1-ct*c2+st*s2*UU)
    w=W*np.sin(PS)**2*np.exp(-S)
    S_r=B*(ct*s2/2+st*c2*UU/2); S_rr=B*(ct*c2/4-st*s2*UU/4)
    S_u=B*st*s2; S_ur=B*st*c2/2
    G=S_r*UU+S_u*(1-UU**2)/rho
    H=(S_rr*UU**2+S_ur*UU*(1-UU**2)/rho+S_r*(1-UU**2)/rho
       +S_ur*UU*(1-UU**2)/rho-3*S_u*UU*(1-UU**2)/rho**2)
    Ph_r=-c2/s2+2/rho; Ph_rr=0.5/s2**2-2/rho**2
    P1=Ph_r*UU; P2=Ph_rr*UU**2+Ph_r*(1-UU**2)/rho
    F=G+P1; Z=w.sum(); E=lambda q:(w*q).sum()/Z
    return dict(EF=E(F),Ek=E(H+P2),EF2=E(F*F),VarF=E(F*F)-E(F)**2,EG2=E(G**2),EGp1=E(G*P1),Ep1s=E(P1**2),EH=E(H),Ep2=E(P2))
for th in (0.0,math.pi/2,0.7*math.pi,math.pi):
    d=run(th)
    print("th=%.3fpi EF=%+.2e Ek=%+.9f EF2=%+.9f  res=%+.2e | EH=%+.6f Ep2=%+.6f VarF=%+.6f EG2=%+.6f 2EGp1=%+.6f Ep1s=%+.6f"%(
        th/math.pi,d["EF"],d["Ek"],d["EF2"],d["Ek"]-d["EF2"],d["EH"],d["Ep2"],d["VarF"],d["EG2"],2*d["EGp1"],d["Ep1s"]))
