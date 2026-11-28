import numpy as np
from scipy.special import ive
from scipy.linalg import expm
print("=== refine m a_inf in the double limit (larger jm, larger beta) ===")
print()
jm=400; nJ=jm+1
Ap=np.zeros((nJ,nJ)); idx=np.arange(nJ-1); Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0
res=[]
for b in (64,128,256,512,1024,2048,4096):
    bb=float(b)
    td=np.array([ive(j+1,bb)/ive(1,bb) for j in range(nJ)])
    Mf=expm(bb*Ap/4)
    w0=Mf[:,0]; sc=float(w0 @ (td*w0))
    T=(Mf @ np.diag(td) @ Mf)/sc
    ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
    r=ev[1]/ev[0]; ma=-np.log(r)
    res.append((b,r,ma))
    print("beta=%-6g : r=%.12f  m a=%.12f"%(b,r,ma), flush=True)
print()
# extrapolate r_inf from r(beta) = r_inf + c/beta
import numpy as np
bs=np.array([x[0] for x in res],dtype=float)
rs=np.array([x[1] for x in res])
mask=bs>=256
p=np.polyfit(1/bs[mask], rs[mask], 1)
print("fit r = r_inf + c/beta  (beta>=256):  r_inf = %.12f   c = %.6f"%(p[1],p[0]))
r_inf=p[1]
print("m a_inf = -ln(r_inf) = %.12f"%(-np.log(r_inf)))
print()
print("candidate closed forms:")
for name,val in [("1/(2pi)",1/(2*np.pi)),("1/(2sqrt(3))",1/(2*np.sqrt(3))),
                 ("1/(3sqrt(2))",1/(3*np.sqrt(2))),("e^{-2}",np.exp(-2)),
                 ("1/(4sqrt(3))",1/(4*np.sqrt(3))),("sqrt(2)-1",np.sqrt(2)-1)]:
    print("   %-12s = %.12f   (ratio to r_inf: %.6f)"%(name,val,val/r_inf))
print()
print("m a_inf vs ln(2)=%.12f, ln(pi)=%.12f, pi/sqrt(2)-... "%(np.log(2),np.log(np.pi)))
