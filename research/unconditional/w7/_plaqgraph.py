import itertools, numpy as np
# 4D hypercubic: link U(x,mu) at x=0, mu=0.
# plaquettes containing it: P(x,mu,nu) and P(x-nu,mu,nu) for nu != mu
mu=0
links_in_plaq={}
def L(kind,site,d,m): return (kind,site,d,m)
# represent a link by (site, direction) up to orientation
def lk(site,d): return (site,d)
for nu in (1,2,3):
    # P(0, 0, nu) : links U(0,0), U(e0,nu), U(enu,0)^-1, U(0,nu)^-1
    P1=[lk((0,0,0,0),0), lk((1,0,0,0),nu), lk((0,1 if nu==1 else 0,1 if nu==2 else 0,1 if nu==3 else 0),0), lk((0,0,0,0),nu)]
    # P(0-nu, 0, nu): links U(-enu,0), U(-enu+e0,nu), U(0,0)^-1, U(-enu,nu)^-1
    neg=[0,0,0,0]; neg[nu]=-1
    s2=list(neg); s2[0]=1
    P2=[lk(tuple(neg),0), lk(tuple(s2),nu), lk((0,0,0,0),0), lk(tuple(neg),nu)]
    links_in_plaq[("+",nu)]=P1
    links_in_plaq[("-",nu)]=P2
verts=set()
for k,v in links_in_plaq.items():
    for l in v: verts.add(l)
verts.discard(lk((0,0,0,0),0))
verts=list(verts)
idx={v:i for i,v in enumerate(verts)}
n=len(verts)
print("number of other links coupled to U(0,0): %d"%n)
A=np.zeros((n,n))
for k,v in links_in_plaq.items():
    vs=[l for l in v if l!=lk((0,0,0,0),0)]
    for a,b in itertools.combinations(vs,2):
        A[idx[a],idx[b]]=1; A[idx[b],idx[a]]=1
ev=np.linalg.eigvalsh(A)
print("adjacency matrix size %d x %d ; degrees = %s"%(n,n,sorted(set(A.sum(1).astype(int)))))
print("eigenvalues (sorted) =", np.round(ev,6))
print("lambda_max = %.6f"%ev[-1])
print()
print("=> 6 disjoint K_3's would give lambda_max = 2")
