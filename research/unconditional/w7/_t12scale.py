import numpy as np, functools
from scipy.special import ive
from scipy.linalg import expm
bt=2.5
jm=20
nJ=jm+1
Ap=np.zeros((nJ,nJ))
for k in range(nJ-1):
    Ap[k,k+1]=0.5; Ap[k+1,k]=0.5
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
Kd=np.diag(td)
print("=== 1-link model: scaling of t1, t2 with beta_s ===")
print()
bs_list=[2,4,8,16,32,64,128,256,512]
t1s=[]; t2s=[]
for bs in bs_list:
    Mf=expm(bs*Ap/2)
    w0=Mf[:,0]
    sc=float(w0 @ (td*w0))
    T=(Mf @ Kd @ Mf)/sc
    ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
    t1s.append(ev[0]); t2s.append(ev[1])
    print("beta_s=%6.1f : t1=%.8f t2=%.8f  t1/beta_s=%.6f  t2=%.6f"%(bs,ev[0],ev[1],ev[0]/bs,ev[1]), flush=True)
t1s=np.array(t1s); t2s=np.array(t2s); b=np.array(bs_list,dtype=float)
p1,c1=np.polyfit(np.log(b),np.log(t1s),1)
p2,c2=np.polyfit(np.log(b),np.log(t2s),1)
print()
print("fit t1 ~ C beta_s^p1 : p1=%.6f  C1=%.6f"%(p1,np.exp(c1)), flush=True)
print("fit t2 ~ C beta_s^p2 : p2=%.6f  C2=%.6f"%(p2,np.exp(c2)), flush=True)
print()
print("t1/beta_s^1.5:", np.round(t1s/b**1.5,6), flush=True)
print("t1/beta_s^2  :", np.round(t1s/b**2,6), flush=True)
print("t2 * sqrt(beta_s):", np.round(t2s*np.sqrt(b),6), flush=True)
print("t2 -> ? (limit of t2)", flush=True)
