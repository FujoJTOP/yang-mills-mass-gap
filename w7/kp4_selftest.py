import math
def phi_int(adj,n):
    size=1<<n; H=bytearray(size); H[0]=1
    for S in range(1,size):
        low=S&-S; i=low.bit_length()-1; rest=S^low
        H[S]=1 if (adj[i]&rest)==0 and H[rest] else 0
    F=[0]*size
    for S in range(1,size):
        low=S&-S; rest=S^low; tot=H[S]; t=rest
        while True:
            blk=low|t
            if blk!=S: tot-=F[blk]*H[S^blk]
            if t==0: break
            t=(t-1)&rest
        F[S]=tot
    return F[size-1]
def cliquesum(a,b):
    n=a+b-1; rows=[0]*n
    E=[(i,j) for i in range(a) for j in range(i+1,a)]
    E+=[(i,j) for i in range(a-1,a+b-1) for j in range(i+1,a+b-1)]
    for (i,j) in E: rows[i]|=1<<j; rows[j]|=1<<i
    return rows,n
for (a,b) in [(3,3),(4,4),(5,5),(6,6),(3,4),(5,6)]:
    rows,n=cliquesum(a,b)
    print(f"K{a} 0-sum K{b}: n={n} phi={phi_int(rows,n)}  product={math.factorial(a-1)*math.factorial(b-1)}")
