import numpy as np
from scipy.special import iv
def Rcheck(beta,J=200,lams=None):
    a=np.arange(J+1); t=2*iv(a+1,beta)/beta
    r0=iv(2,beta)/iv(1,beta)
    worst=0.0; wl=None
    for lam in lams:
        M=iv(np.abs(a[:,None]-a[None,:]),lam)-iv(a[:,None]+a[None,:]+2,lam)
        T=M@np.diag(t)@M; T=(T+T.T)/2
        ev=np.sort(np.linalg.eigvalsh(T))[::-1]
        r=ev[1]/ev[0]
        if r>worst: worst=r; wl=lam
    return r0,worst,wl
lams=[0.001,0.01,0.05,0.1,0.2,0.5,0.8,1.0,1.5,2.0,3.0,4.0,6.0,8.0,12.0,16.0,24.0,32.0]
print("theorem (R) check across beta :  is max_lambda t2/t1 <= I2/I1 ?")
print("%8s %20s %20s %10s %10s"%("beta","I2/I1","max t2/t1","at lambda","margin"))
for beta in (0.2,1.0,2.5,5.0,10.0):
    r0,w,wl=Rcheck(beta,lams=lams)
    print("%8.1f %20.10f %20.10f %10.3f %10s"%(beta,r0,w,wl,"OK" if w<=r0 else "VIOLATION"), flush=True)
