"""kp4_chain.py -- extend the KP4 search: beam n=9,10; structure of the n=11 two-bundle graph."""
import itertools, time, sys

D = 4
T0 = time.time()
def log(*a): print(*a, flush=True)
def unit(i):
    v=[0]*D; v[i]=1; return tuple(v)
def add(u,v): return tuple(u[k]+v[k] for k in range(D))
DIRS=[(a,b) for a in range(D) for b in range(a+1,D)]

def phi_int(adj,n):
    size=1<<n; H=bytearray(size); H[0]=1
    for S in range(1,size):
        low=S&-S; i=low.bit_length()-1; rest=S^low
        H[S]= 1 if (adj[i]&rest)==0 and H[rest] else 0
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

def build_pool(R):
    sites=sorted(x for x in itertools.product(range(-R,R+1),repeat=D) if sum(map(abs,x))<=R)
    plaq=[(x,a,b) for x in sites for (a,b) in DIRS]
    idx={p:i for i,p in enumerate(plaq)}; P=len(plaq); linkmap={}
    for i,(x,a,b) in enumerate(plaq):
        for L in ((x,a),(add(x,unit(a)),b),(x,b),(add(x,unit(b)),a)):
            linkmap.setdefault(L,[]).append(i)
    nb=[set() for _ in range(P)]
    for L,ps in linkmap.items():
        for i in range(len(ps)):
            for j in range(i+1,len(ps)):
                nb[ps[i]].add(ps[j]); nb[ps[j]].add(ps[i])
    return plaq,idx,linkmap,nb,P

R=2*9
plaq,idx,linkmap,nb,P=build_pool(R)
log(f"pool R={R}: {P} plaquettes")
p0=idx[((0,0,0,0),0,1)]

def rows_of(vs,nb):
    n=len(vs); rows=[0]*n
    for i,v in enumerate(vs):
        r=0
        for j,w in enumerate(vs):
            if j!=i and w in nb[v]: r|=1<<j
        rows[i]=r
    return rows

def beam(nmax,width,startp=None):
    cur={frozenset([startp if startp is not None else p0])}
    out={}
    for n in range(2,nmax+1):
        t1=time.time(); nxt=set()
        for s in cur:
            cand=set()
            for v in s: cand|=nb[v]
            for c in (cand-s): nxt.add(s|{c})
        scored=[]
        for s in nxt:
            vs=sorted(s); ph=abs(phi_int(rows_of(vs,nb),n)); scored.append((ph,s))
        scored.sort(key=lambda t:-t[0])
        cur={s for _,s in scored[:width]}
        out[n]=(scored[0][0], scored[0][1], len(nxt))
        log(f"  n={n}: {len(nxt)} candidates -> best |phi|={scored[0][0]} ({scored[0][0]**(1/n):.6f})  [{time.time()-t1:.1f}s]")
    return out

log("[H] beam n=2..10 (width 120), start p0")
out=beam(10,120)
log("\n[H] best-at-n table")
for n,(m,s,c) in out.items():
    log(f"  n={n:>2}  max|phi|>={m:>8}  m^(1/n)={m**(1/n):.6f}")

log("\n[I] the n=11 two-bundle graph: structure + exact phi")
L0=((0,0,0,0),0); L1=((0,0,0,0),1)
B0=sorted(linkmap[L0]); B1=sorted(linkmap[L1])
sh=set(B0)&set(B1)
vs=sorted(set(B0)|set(B1))
pos={v:i for i,v in enumerate(vs)}
half0=[pos[v] for v in B0 if v in pos]
cut=[v for v in vs if v in sh]
log(f"  |B0|={len(B0)} |B1|={len(B1)} shared={len(sh)} n={len(vs)}")
# edges crossing between the two halves (excluding the shared vertex)
h0=[pos[v] for v in B0]; h1=[pos[v] for v in B1]
cross=0; crosspairs=[]
for i in h0:
    for j in h1:
        if j!=i and j in nb[vs[i]]:
            if vs[i] not in sh and vs[j] not in sh:
                cross+=1; crosspairs.append((plaq[vs[i]],plaq[vs[j]]))
log(f"  cross edges between the two 5-sets (cut vertex removed): {cross}")
for a,b in crosspairs[:8]: log(f"    {a} -- {b}")
rows=rows_of(vs,nb)
log(f"  |E|={sum(bin(r).count('1') for r in rows)//2}  phi^U={phi_int(rows,len(vs))}  "
    f"|phi|^(1/11)={abs(phi_int(rows,len(vs)))**(1/11):.6f}")
log(f"  block-product (cut vertex) would be 120^2=14400;  actual/14400 = {abs(phi_int(rows,len(vs)))/14400:.3f}")

log("\n[J] cut-vertex test: is the shared plaquette a cut vertex? (5-set adjacency)")
res=set(h0)&set(h1)
log(f"  shared index {res}, deg inside set = {bin(rows[list(res)[0]]).count('1')}")
log(f"\nelapsed {time.time()-T0:.1f}s")
