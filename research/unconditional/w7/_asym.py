import numpy as np
from scipy.special import ive
from scipy.linalg import expm
bt=2.5
print("=== asymptotic analysis of m a(beta_s) with jm=8*beta_s ===")
print()
rows=[]
for bs in (4,8,16,32,64,128,256):
    jm=int(8*bs)
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
    rows.append((bs,ev[0],ev[1],r,ma))
    print("beta_s=%7.1f jm=%5d : t1=%.10f t2=%.10f m a=%.10f"%(bs,jm,ev[0],ev[1],ma), flush=True)
b=np.array([x[0] for x in rows],dtype=float)
ma=np.array([x[4] for x in rows])
print()
# fits
for lo,hi in [(0,len(b)),(2,len(b)),(3,len(b))]:
    sl=slice(lo,hi)
    p,c=np.polyfit(np.log(b[sl]),np.log(ma[sl]),1)
    print("fit over beta_s in [%g,%g] : m a ~ beta_s^%.5f  (C=%.5f)"%(b[lo],b[-1],p,np.exp(c)), flush=True)
print()
# log fit: m a = A ln(beta_s) + B
A,B=np.polyfit(np.log(b),ma,1)
print("log fit : m a = %.6f ln(beta_s) + %.6f"%(A,B), flush=True)
print("residuals of log fit:", np.round(ma-(A*np.log(b)+B),6), flush=True)
# sqrt fit: m a = A sqrt(beta_s) + B
A2,B2=np.polyfit(np.sqrt(b),ma,1)
print("sqrt fit: m a = %.6f sqrt(beta_s) + %.6f"%(A2,B2), flush=True)
print("residuals of sqrt fit:", np.round(ma-(A2*np.sqrt(b)+B2),6), flush=True)
