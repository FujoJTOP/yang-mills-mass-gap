import numpy as np
from scipy.special import ive
from scipy.linalg import expm
print("=== is m a related to the t_j profile? scan beta_t at fixed beta_s ===")
print()
jm=40
nJ=jm+1
Ap=np.zeros((nJ,nJ))
idx=np.arange(nJ-1)
Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0
bs=1.0
Mf=expm(bs*Ap/4)
for bt in (0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0):
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    w0=Mf[:,0]
    sc=float(w0 @ (td*w0))
    T=(Mf @ np.diag(td) @ Mf)/sc
    ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
    r=ev[1]/ev[0]
    ma=-np.log(r)
    # candidate: -ln(t_1) etc
    print("beta_t=%.2f : m a=%.8f   -ln(t1)=%.8f  -ln(t2)/2=%.8f  t1=%.6f"%(bt,ma,-np.log(td[1]) if td[1]>0 else np.inf,-np.log(td[2])/2 if td[2]>0 else np.inf,td[1]), flush=True)
