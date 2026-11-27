import numpy as np
from collections import defaultdict
N=4
ints=[(b,t) for t in range(N+1) for b in range(t+1)]
Pm=4*N+2
shape=(N+1,)
dp={}
for (b,t) in ints:
    if t-b+1>N: continue
    A=np.zeros(shape); A[t-b+1]+=1
    dp[(b,t,1 if b==0 else 0)]=A
print("init keys:", {k:np.nonzero(v)[0].tolist() for k,v in dp.items()})
out=np.zeros(shape)
for k,A in dp.items():
    if k[2]==1: out+=A
for step in range(1,N):
    ndp=defaultdict(lambda: np.zeros(shape))
    for (b,t,h),A in dp.items():
        nz=np.nonzero(A)[0]
        if len(nz)==0: continue
        nmax=nz.max()
        for (b2,t2) in ints:
            if max(b,b2)>min(t,t2): continue
            dn=t2-b2+1
            if dn>N or nmax+dn>N: continue
            h2=1 if (h==1 or b2==0) else 0
            ndp[(b2,t2,h2)][dn:]+=A[:N+1-dn]
    dp=ndp
    print("step",step,"states",{k:np.nonzero(v)[0].tolist() for k,v in list(dp.items())[:8]})
    for k,A in dp.items():
        if k[2]==1: out+=A
print("out",out)
