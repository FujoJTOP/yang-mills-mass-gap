import numpy as np
def aij(R1,R2,r1,r2):
    a11=r1/(2*np.pi*R1); a22=r2/(2*np.pi*R2)
    a12=r1/(2*np.pi*R2)+r2/(2*np.pi*R1)
    return a11,a12,a22
print("det check  det = -(r1R1-r2R2)^2/(16 pi^2 R1^2 R2^2)")
for (R1,R2,r1,r2) in [(1.,2.,1.,1.),(1.,2.,1.,0.),(1.,1.,1.,1.),(3.,2.,2.,3.),(1.,1.,0.,1.)]:
    a11,a12,a22=aij(R1,R2,r1,r2)
    det=a11*a22-a12**2/4
    pred=-(r1*R1-r2*R2)**2/(16*np.pi**2*R1**2*R2**2)
    print(f" R1={R1},R2={R2},rho={r1},{r2}: det={det: .10e} pred={pred: .10e}  Delta={a12**2-4*a11*a22: .10e}")
# diagonal Riccati: u' = 2u^2 + u/R + K , u=-f'/f>=0
def blow(R0,u0,K,Rmax=200.0,n=400000):
    R=R0;u=u0;h=(Rmax-R0)/n
    for _ in range(n):
        u=u+0.5*h*((2*u*u+u/R+K)+0)  # explicit euler, detect overflow
        R+=h
        if u>1e6: return R
    return None
for K in [0.0,0.1,-0.1]:
    print("K=",K,"  blow-up R =",blow(1.0,0.05,K))
