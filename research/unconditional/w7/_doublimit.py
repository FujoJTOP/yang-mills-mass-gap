import numpy as np
from scipy.special import ive
from scipy.linalg import expm
print("=== 1-link model in the double limit beta_s = beta_t = beta -> inf ===")
print("    K diag = t_j = I_{j+1}(b)/I_1(b) ~ 1 - j(j+2)/(2b)")
print()
jm=80; nJ=jm+1
Ap=np.zeros((nJ,nJ)); idx=np.arange(nJ-1); Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0
for b in (1,2,4,8,16,32,64,128,256,512):
    bb=float(b)
    td=np.array([ive(j+1,bb)/ive(1,bb) for j in range(nJ)])
    Mf=expm(bb*Ap/4)
    w0=Mf[:,0]; sc=float(w0 @ (td*w0))
    T=(Mf @ np.diag(td) @ Mf)/sc
    ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
    r=ev[1]/ev[0]; ma=-np.log(r)
    print("beta=%-6g : t1=%.10f  t2=%.10f  r=%.10f  m a=%.10f"%(b,ev[0],ev[1],r,ma), flush=True)
print()
print("=> if m a -> 0 as beta -> inf, the continuum limit exists for this model")
