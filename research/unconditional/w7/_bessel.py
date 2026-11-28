import numpy as np
from scipy.special import ive, i0e, i1e
from scipy.linalg import expm
bt=2.5
jm=60
nJ=jm+1
Ap=np.zeros((nJ,nJ))
idx=np.arange(nJ-1)
Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0   # adjacency, entries 1
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
print("=== test: M_f = exp(beta_s * A/4) has Bessel kernel I_|j-k|(beta_s/2) ===")
print()
bs=1.0
Mf=expm(bs*Ap/4)
# compare with I_|j-k|(beta_s/2) for the interior
print("M_f entries vs I_|j-k|(beta_s/2):")
for j in range(0,5):
    for k in range(0,5):
        pred=ive(abs(j-k),bs/2)*np.exp(bs/2)   # ive = exponentially scaled, so I = ive*exp
        print("   M_f[%d,%d]=%.10f   I_|j-k|(bs/2)*e^{bs/2}=%.10f   (bs/2=%.2f)"%(j,k,Mf[j,k],pred,bs/2), flush=True)
    print()
print("   -> if diagonal matches e^{bs/2} and off-diagonal matches I_d(bs/2), the identification holds")
print()
print("check: M_f[0,0]=%.10f  e^{0.25}*I_0(0.5)? =%.10f"%(Mf[0,0], np.exp(0.25)*i0e(0.5)*np.exp(0.5)))
