import numpy as np, math
from collections import defaultdict

def counts(N,K):
    ints=[(b,t) for t in range(N+1) for b in range(t+1) if b<=K]
    shape=(N+1,)
    dp={}
    for (b,t) in ints:
        if t-b+1>N: continue
        A=np.zeros(shape); A[t-b+1]+=1
        dp[(b,t,1 if b==0 else 0)]=A
    out=np.zeros(shape)
    for k,A in dp.items():
        if k[2]==1: out+=A
    for step in range(1,N):
        ndp=defaultdict(lambda: np.zeros(shape))
        for (b,t,h),A in dp.items():
            for (b2,t2) in ints:
                if max(b,b2)>min(t,t2): continue
                dn=t2-b2+1
                if dn>N: continue
                h2=1 if (h==1 or b2==0) else 0
                ndp[(b2,t2,h2)][dn:]+=A[:N+1-dn]
        dp=ndp
        for k,A in dp.items():
            if k[2]==1: out+=A
    return out

N=22
for K in [0,1,2,3,4,6,22]:
    a=counts(N,K)
    r=[int(round(a[n])) for n in range(1,N+1)]
    ratios=[r[n-1]/r[n-2] for n in range(3,N+1)]
    print("K=%2d  a_n(n=%d)=%d   ratio(last)=%.6f   a_%d^(1/%d)=%.6f"%(K,N,r[-1],ratios[-1],N,N,r[-1]**(1.0/N)))
    if K<=2:
        print("    seq:",r[:10])
print()
print("e =",math.e)
