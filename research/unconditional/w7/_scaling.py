import numpy as np, time
from scipy.special import iv
from scipy.linalg import eigh
beta=2.5
def q2(lam,J=400):
    a=np.arange(J+1); t=2*iv(a+1,beta)/beta
    M=iv(np.abs(a[:,None]-a[None,:]),lam)-iv(a[:,None]+a[None,:]+2,lam)
    T=M@np.diag(t)@M; T=(T+T.T)/2
    ev,V=eigh(T,subset_by_index=[J-1,J])
    n=4000
    th=(np.arange(n)+0.5)*np.pi/n
    S=np.sin(np.outer(th,np.arange(1,J+2)))
    out=[]
    for k in (0,1):
        p=S@V[:,k]
        out.append(np.sum(th**2*p**2)/np.sum(p**2))
    return out
print("scaling of <theta^2>_i in lam   (does it decay like 1/lam or slower?)")
print("%8s %14s %14s %14s %14s"%("lam","<th2>_1","lam*<th2>_1","<th2>_2","ratio"))
t0=time.time()
for lam in (8.0,16.0,32.0,64.0,128.0):
    q1,q2v=q2(lam)
    print("%8.1f %14.8e %14.8f %14.8e %14.6f"%(lam,q1,lam*q1,q2v,q2v/q1), flush=True)
print("total %.1fs"%(time.time()-t0))
