import itertools, numpy as np
def graph_for(d):
    mu=0
    def lk(site,dr): return (tuple(site),dr)
    links_in_plaq={}
    for nu in range(1,d):
        # P(0, 0, nu)
        e0=[0]*d; e0[0]=1
        en=[0]*d; en[nu]=1
        P1=[lk([0]*d,0), lk(e0,nu), lk(en,0), lk([0]*d,nu)]
        # P(0-nu, 0, nu)
        neg=[0]*d; neg[nu]=-1
        s2=list(neg); s2[0]=1
        P2=[lk(neg,0), lk(s2,nu), lk([0]*d,0), lk(neg,nu)]
        links_in_plaq[("+",nu)]=P1
        links_in_plaq[("-",nu)]=P2
    zero=lk([0]*d,0)
    verts=set()
    for v in links_in_plaq.values():
        for l in v: verts.add(l)
    verts.discard(zero); verts=list(verts); idx={v:i for i,v in enumerate(verts)}
    n=len(verts); A=np.zeros((n,n))
    for v in links_in_plaq.values():
        vs=[l for l in v if l!=zero]
        for a,b in itertools.combinations(vs,2):
            A[idx[a],idx[b]]=1; A[idx[b],idx[a]]=1
    return n,A
for d in (2,3,4,5):
    n,A=graph_for(d)
    ev=np.linalg.eigvalsh(A)
    degs=sorted(set(A.sum(1).astype(int)))
    print("d=%d :  vertices=%3d  plaquettes/link=%d  degrees=%s  lambda_max=%.6f"%(d,n,2*(d-1),degs,ev[-1]))
print()
print("=> the slope dr/dbeta_s|_0 = (u/2) lambda_max :")
for d in (2,3,4,5):
    n,A=graph_for(d); lam=np.linalg.eigvalsh(A)[-1]
    print("     d=%d : lambda_max=%.4f  slope = %.4f u"%(d,lam,lam/2))
