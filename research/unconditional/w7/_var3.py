import numpy as np, functools, scipy.sparse as sp, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=2; bs=1.0; bt=2.5
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B)
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((d,d))
    mm,nn=(-M2,-N2) if conj else (M2,N2)
    ph=((-1.0)**((M2-N2)//2)) if conj else 1.0
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2-mm,n2-nn
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,mm,j2,-m2)*d3(j3,n3,J2,nn,j2,-n2)
            if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
    return ph*A
u={}
for M2 in (-1,1):
    for N2 in (-1,1):
        u[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
        u[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
idx=[(-1,-1),(-1,1),(1,-1),(1,1)]
Ip=sp.identity(d,format='csr')
def Wilson4(a,b,c,dd):
    tot=None
    for (p,q) in idx:
        A=u[(p,q)]
        for (r,s) in idx:
            Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
            term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
            tot=term if tot is None else tot+term
    return (0.5*tot).tocsr()
W4=Wilson4(0,1,2,3)   # 4-link operator on links 0,1,2,3
# truncated basis for links 1,2,3: all j=0 or exactly one at j=1/2
basis3=[( (0,0,0),(0,0,0),(0,0,0) )]
for pos in range(3):
    for s in [b for b in B if b[0]==1]:
        t=[(0,0,0)]*3; t[pos]=s; basis3.append(tuple(t))
nb3=len(basis3)
print("nb3=%d (full %d)"%(nb3,d**3), flush=True)
# basis for (link0, links1-3): link0 in all B (d states), links 1-3 in basis3
basis04=[(c,l) for c in B for l in basis3]
n04=len(basis04)
print("n04=%d (full %d)"%(n04,d**4), flush=True)
# build the projection matrix P: (d^4) x n04
rows=[]; cols=[]; vals=[]
for j,(c,l) in enumerate(basis04):
    ii=0
    for b in (c,)+l: ii=ii*d+index[b]
    rows.append(ii); cols.append(j); vals.append(1.0)
P=sp.csr_matrix((vals,(rows,cols)),shape=(d**4,n04))
# restricted W4 on (link0, links1-3)
W4r=(P.T @ W4 @ P).tocsr()
print("W4r shape=%s nnz=%d"%(W4r.shape,W4r.nnz), flush=True)
# Two plaquettes: P1 uses (link0, links1-3); P2 uses (link0, links4-6)
# The variational basis: (c, l, r) with c in B (d), l,r in basis3 (nb3)
# W1 acts on (c,l): W4r.  W2 acts on (c,r): also W4r.
# Full restricted operator on the d*nb3*nb3 space:
# reshape to (c,l,r); W1 acts on the first two, W2 on (c) and the third (but the c index is shared!)
# W2 on (c,r) means: for fixed c, it's the same 4-link operator with the roles of (l,r) swapped
# Build: on the tensor product (c,l,r), W1 = W4r acting on (c,j) indices and identity on r
ncb=nb3
Dv=n04*nb3  # (c,l) x r
I_r=sp.identity(nb3,format='csr')
W1_full=sp.kron(W4r,I_r,format='csr')
# W2 acts on (c,r) and identity on l -> need to permute (c,l,r)->(c,r,l)
# Use a permutation matrix P_perm on the (c,l,r) space
perm=np.arange(Dv)
# index (c,l,r) -> flat = ((c*nb3)+l)*nb3 + r
# we want (c,r,l) -> flat2 = ((c*nb3)+r)*nb3 + l
perm2=np.zeros(Dv,dtype=int)
for c in range(d):
    for l in range(nb3):
        for r in range(nb3):
            perm2[((c*nb3)+l)*nb3+r] = ((c*nb3)+r)*nb3+l
Pperm=sp.csr_matrix((np.ones(Dv),(np.arange(Dv),perm2)),shape=(Dv,Dv))
W2_full=(Pperm.T @ sp.kron(W4r,I_r,format='csr') @ Pperm).tocsr()
Wtot=0.5*(W1_full+W2_full)
print("Wtot shape=%s nnz=%d"%(Wtot.shape,Wtot.nnz), flush=True)
# K operator on the variational space
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
Kdiag=np.zeros(Dv)
for c in range(d):
    kj=td[B[c][0]]
    for il,l in enumerate(basis3):
        kl=np.prod([td[b[0]] for b in l])
        for r in range(nb3):
            kr=np.prod([td[b[0]] for b in basis3[r]])
            Kdiag[((c*nb3+il)*nb3+r)]=kj*kl*kr
K=sp.diags(Kdiag)
print("K ready", flush=True)
# M_f = exp(-2 beta_s) exp(beta_s Wtot)
def Mf(v): return np.exp(-2*bs)*expm_multiply(bs*Wtot, v)
e0=np.zeros(Dv); e0[0]=1.0
w0=Mf(e0); scale=float(np.dot(w0,Kdiag*w0))
def T(v): return Mf(K*(Mf(v)))/scale
from scipy.sparse.linalg import LinearOperator, eigsh
op=LinearOperator((Dv,Dv),matvec=T,dtype=float)
ev,V=eigsh(op,k=3,which="LA",tol=1e-12,maxiter=5000)
o=np.argsort(ev)[::-1]; ev=ev[o]
print("VARIATIONAL (c in B, l,r in truncated): Dv=%d"%Dv, flush=True)
print("  t1=%.10f t2=%.10f t3=%.10f  r=t2/t1=%.10f"%(ev[0],ev[1],ev[2],ev[1]/ev[0]), flush=True)
for i in range(3):
    res=np.linalg.norm(T(V[:,i])-ev[i]*V[:,i])/abs(ev[i])
    print("   eig %d residual = %.2e"%(i,res), flush=True)
