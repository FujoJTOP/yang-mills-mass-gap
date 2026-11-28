import numpy as np, functools
from scipy.special import ive
from scipy.linalg import expm
bt=2.5
print("=== 1-link model: scaling of m a with beta_s ===")
print()
for jm in (20,):
    nJ=jm+1
    Ap=np.zeros((nJ,nJ))
    for k in range(nJ-1):
        Ap[k,k+1]=0.5; Ap[k+1,k]=0.5
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    Kd=np.diag(td)
    rows=[]
    for bs in (0.25,0.5,1,2,4,8,16,32,64,128):
        Mf=expm(bs*Ap/2)
        w0=Mf[:,0]
        sc=float(w0 @ (td*w0))
        T=(Mf @ Kd @ Mf)/sc
        ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
        r=ev[1]/ev[0]
        rows.append((bs,ev[0],ev[1],r,-np.log(r)))
        print("beta_s=%7.2f : t1=%.8f t2=%.8f r=%.10f m a=%.8f"%(bs,ev[0],ev[1],r,-np.log(r)), flush=True)
    print()
    bss=np.array([x[0] for x in rows]); mas=np.array([x[4] for x in rows])
    # fit m a = C * beta_s^p on the large-beta range
    mask=bss>=2
    p,C=np.polyfit(np.log(bss[mask]),np.log(mas[mask]),1)
    print("fit m a = C * beta_s^p for beta_s>=2 :")
    print("   p = %.6f   C = %.6f"%(p,np.exp(C)), flush=True)
    # also check m a / beta_s^(1/2)
    print("   m a / sqrt(beta_s):", np.round(mas[mask]/np.sqrt(bss[mask]),6), flush=True)
    print("   m a / beta_s:      ", np.round(mas[mask]/bss[mask],6), flush=True)
