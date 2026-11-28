import numpy as np, functools
from scipy.special import ive
from scipy.linalg import expm
bt=2.5
print("=== 1-link model with jm = max(20, 4*beta_s) to avoid truncation ===")
print()
rows=[]
for bs in (0.5,1,2,4,8,16,32,64,128):
    jm=int(max(20, 4*bs))
    nJ=jm+1
    Ap=np.zeros((nJ,nJ))
    for k in range(nJ-1):
        Ap[k,k+1]=0.5; Ap[k+1,k]=0.5
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    Mf=expm(bs*Ap/2)
    w0=Mf[:,0]
    sc=float(w0 @ (td*w0))
    T=(Mf @ np.diag(td) @ Mf)/sc
    ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
    r=ev[1]/ev[0]
    rows.append((bs,ev[0],ev[1],r,-np.log(r)))
    print("beta_s=%7.2f jm=%4d : t1=%.8f t2=%.8f r=%.10f m a=%.8f"%(bs,jm,ev[0],ev[1],r,-np.log(r)), flush=True)
b=np.array([x[0] for x in rows]); ma=np.array([x[4] for x in rows])
mask=b>=2
p,c=np.polyfit(np.log(b[mask]),np.log(ma[mask]),1)
print()
print("fit m a ~ C beta_s^p (beta_s>=2) : p=%.6f  C=%.6f"%(p,np.exp(c)), flush=True)
print("m a / beta_s^0.5:", np.round(ma/np.sqrt(b),6), flush=True)
