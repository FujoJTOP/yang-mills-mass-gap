import numpy as np
from scipy.special import ive
from scipy.linalg import expm
bt=2.5
print("=== 1-link model: check m a ~ 1/jm^2 at fixed beta_s ===")
print()
for bs in (0.5,1.0,2.0,4.0):
    print("beta_s=%.1f:"%bs, flush=True)
    rows=[]
    for jm in (10,20,40,80,160):
        nJ=jm+1
        Ap=np.zeros((nJ,nJ))
        idx=np.arange(nJ-1)
        Ap[idx,idx+1]=0.5; Ap[idx+1,idx]=0.5
        td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
        Mf=expm(bs*Ap/2)
        w0=Mf[:,0]
        sc=float(w0 @ (td*w0))
        T=(Mf @ np.diag(td) @ Mf)/sc
        ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
        r=ev[1]/ev[0]
        ma=-np.log(r)
        rows.append((jm,ev[0],ev[1],r,ma))
        print("   jm=%4d : t1=%.8f t2=%.8f r=%.8f m a=%.8f  m a*jm^2=%.6f"%(jm,ev[0],ev[1],r,ma,ma*jm*jm), flush=True)
    print()
