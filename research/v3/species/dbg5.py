import numpy as np, math
from collections import defaultdict
def pair_curv(I,J):
    b,t=I; b2,t2=J; c=0
    for y in range(min(b,b2)-1,max(t,t2)+2):
        k=((b<=y-1<=t)+(b2<=y-1<=t2)+(b<=y<=t)+(b2<=y<=t2))
        if k==3: c+=1
    return c
N=30
ints=[(b,t) for t in range(N+1) for b in range(t+1)]
def newstate():
    return [np.zeros(N+1) for _ in range(3)]
dp={}
for (b,t) in ints:
    if t-b+1>N: continue
    A=newstate(); A[0][t-b+1]=1; A[1][t-b+1]=(2*(t-b+1)+2); A[2][t-b+1]=0
    dp[(b,t,1 if b==0 else 0)]=A
out=[np.zeros(N+1) for _ in range(3)]
def acc(out,dp):
    for (b,t,h),A in dp.items():
        if h==1:
            for i in range(3): out[i]+=A[i]
acc(out,dp)
for step in range(1,N):
    ndp=defaultdict(newstate)
    for (b,t,h),A in dp.items():
        for (b2,t2) in ints:
            if max(b,b2)>min(t,t2): continue
            dn=t2-b2+1
            if dn>N: continue
            ov=min(t,t2)-max(b,b2)+1
            dpp=2*dn+2-2*ov
            dc=pair_curv((b,t),(b2,t2))
            h2=1 if (h==1 or b2==0) else 0
            B=ndp[(b2,t2,h2)]
            B[0][dn:]+=A[0][:N+1-dn]
            B[1][dn:]+=A[1][:N+1-dn]+dpp*A[0][:N+1-dn]
            B[2][dn:]+=A[2][:N+1-dn]+dc*A[0][:N+1-dn]
    dp=ndp
    acc(out,dp)
print(" n   C_n      <p>    p/n    <curv>  curv/n")
for n in range(1,N+1):
    c=out[0][n]
    if c<1: continue
    print("%2d %10d  %6.2f %.4f  %6.2f  %.4f"%(n,round(c),out[1][n]/c,out[1][n]/c/n,out[2][n]/c,out[2][n]/c/n))
