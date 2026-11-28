import numpy as np
from scipy.special import ive
from scipy.linalg import expm
jm=60; nJ=jm+1
Ap=np.zeros((nJ,nJ))
idx=np.arange(nJ-1)
Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0
bt=2.5
print("=== effective 2x2 analysis of T = M_f K M_f / scale ===")
print()
for bs in (0.5,1.0,2.0,4.0):
    Mf=expm(bs*Ap/4)
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    w0=Mf[:,0]
    sc=float(w0 @ (td*w0))
    Tfull=(Mf @ np.diag(td) @ Mf)/sc
    # restrict to {|0>,|1>}
    T2=Tfull[:2,:2]
    ev2=np.linalg.eigvalsh(0.5*(T2+T2.T))[::-1]
    evf=np.linalg.eigvalsh(0.5*(Tfull+Tfull.T))[::-1]
    print("beta_s=%.1f : full t1=%.8f t2=%.8f r=%.8f m a=%.8f"%(bs,evf[0],evf[1],evf[1]/evf[0],-np.log(evf[1]/evf[0])), flush=True)
    print("          2x2  t1=%.8f t2=%.8f r=%.8f m a=%.8f"%(ev2[0],ev2[1],ev2[1]/ev2[0],-np.log(ev2[1]/ev2[0])), flush=True)
    print("          T2 =", np.round(T2,6), flush=True)
    # also the 3x3
    T3=Tfull[:3,:3]
    ev3=np.linalg.eigvalsh(0.5*(T3+T3.T))[::-1]
    print("          3x3  t1=%.8f t2=%.8f r=%.8f m a=%.8f"%(ev3[0],ev3[1],ev3[1]/ev3[0],-np.log(ev3[1]/ev3[0])), flush=True)
    print()
