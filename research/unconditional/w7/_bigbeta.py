import numpy as np
from scipy.special import ive
def Rcheck(beta,J=200):
    a=np.arange(J+1)
    td=ive(a+1,beta)/ive(1,beta)      # t_j/t_0 ; the common e^beta factor cancels in the ratio
    r0=td[1]
    worst=0.0; wl=None
    for lam in (0.001,0.01,0.1,0.5,1.0,2.0,4.0,8.0,16.0,32.0):
        M=ive(np.abs(a[:,None]-a[None,:]),lam)*np.exp(lam)-ive(a[:,None]+a[None,:]+2,lam)*np.exp(lam) if False else None
        # M_f entries: I_{|a-b|}(lam) - I_{a+b+2}(lam); use ive*exp(lam) -> prone to overflow; use direct iv for moderate lam
        from scipy.special import iv
        M=iv(np.abs(a[:,None]-a[None,:]),lam)-iv(a[:,None]+a[None,:]+2,lam)
        T=M@np.diag(td)@M; T=(T+T.T)/2
        ev=np.sort(np.linalg.eigvalsh(T))[::-1]
        r=ev[1]/ev[0]
        if r>worst: worst=r; wl=lam
    return r0,worst,wl
print("obstacle 3 check: does (R) survive at large beta where K is astronomically large?")
print("%10s %22s %22s %10s %10s"%("beta","I2/I1 = t1/t0","max t2/t1","at lam","verdict"))
for beta in (10.0,20.0,50.0,100.0,200.0,400.0):
    r0,w,wl=Rcheck(beta)
    print("%10.1f %22.16f %22.16f %10.3f %10s"%(beta,r0,w,wl,"OK" if w<=r0 else "VIOLATION"), flush=True)
