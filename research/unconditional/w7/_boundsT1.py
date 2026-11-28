import numpy as np
from scipy.special import ive
from scipy.linalg import expm
print("=== verify: the round-240 bounds ARE T1 (lower) and trivial (upper) ===")
print()
print("m a = -ln(r) with r = t2/t1 (single-chain (R))")
print()
print("LOWER: m a >= -ln t1  <=>  -ln(t2/t1) >= -ln t1  <=>  t2/t1 <= t1  <=>  t2 <= t1^2")
print("       t2 = I3/I1, t1 = I2/I1  =>  I3/I1 <= I2^2/I1^2  =>  I1 I3 <= I2^2   <-- T1 !")
print()
print("UPPER: m a <= -ln t2  <=>  -ln(t2/t1) <= -ln t2  <=>  t2/t1 >= t2  <=>  t1 <= 1")
print("       t1 = I2/I1 <= 1  (Bessel ratio monotonic)   <-- trivial")
print()
bt=2.5
jm=60; nJ=jm+1
Ap=np.zeros((nJ,nJ))
idx=np.arange(nJ-1)
Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
print("numerical check of t2 <= t1^2 and t1 <= 1 (beta_t=2.5):")
print("  t1 = %.15f   t1^2 = %.15f"%(td[1],td[1]*td[1]), flush=True)
print("  t2 = %.15f   t2 <= t1^2 ?  %s"%(td[2],td[2]<=td[1]*td[1]), flush=True)
print("  t1 <= 1 ?  %s"%(td[1]<=1.0), flush=True)
print()
bs=1.0
Mf=expm(bs*Ap/4)
w0=Mf[:,0]
sc=float(w0 @ (td*w0))
T=(Mf @ np.diag(td) @ Mf)/sc
ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
r=ev[1]/ev[0]
ma=-np.log(r)
print("measured m a = %.15f"%ma, flush=True)
print("  -ln t1 = %.15f   (lower bound)"%(-np.log(td[1])), flush=True)
print("  -ln t2 = %.15f   (upper bound)"%(-np.log(td[2])), flush=True)
print("  lower <= m a <= upper ?  %s"%((-np.log(td[1])<=ma) and (ma<=-np.log(td[2]))), flush=True)
print()
print("=> the two-sided bound is FULLY PROVEN (lower = T1, upper = trivial)")
