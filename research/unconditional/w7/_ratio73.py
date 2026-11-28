import numpy as np, time
from scipy.special import iv
from scipy.linalg import eigh
beta=2.5
def ratio(lam,J=1200):
    a=np.arange(J+1); t=2*iv(a+1,beta)/beta
    M=iv(np.abs(a[:,None]-a[None,:]),lam)-iv(a[:,None]+a[None,:]+2,lam)
    T=M@np.diag(t)@M; T=(T+T.T)/2
    ev,V=eigh(T,subset_by_index=[J-1,J])
    n=6000
    th=(np.arange(n)+0.5)*np.pi/n
    S=np.sin(np.outer(th,np.arange(1,J+2)))
    q=[]
    for k in (0,1):
        p=S@V[:,k]
        q.append(np.sum(th**2*p**2)/np.sum(p**2))
    return q[0]/q[1]     # = <th2>_2 / <th2>_1
print("large-lam limit of <theta^2>_2 / <theta^2>_1   (rank-1 x e^{-x^2/2}, rank-2 x^3 e^{-x^2/2} predicts 7/3 = 2.333333)")
print("%8s %18s %14s"%("lam","ratio","dev from 7/3"))
t0=time.time()
for lam in (32.0,64.0,128.0,256.0,512.0):
    r=ratio(lam)
    print("%8.1f %18.12f %+14.8f"%(lam,r,r-7.0/3.0), flush=True)
print("total %.1fs"%(time.time()-t0))
