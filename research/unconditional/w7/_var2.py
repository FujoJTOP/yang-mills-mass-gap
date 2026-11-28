import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
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
# Build the 4-link operator on the FIRST 4 links (0,1,2,3)
W4a=Wilson4(0,1,2,3)
# 3-link identities for the last 3 links
I3=sp.identity(d**3,format='csr')
# Wtot = 0.5*(W4a (x) I3 + (permuted W4b))
# The permuted W4b acts on links (0,4,5,6). Build it as a matrix:
# In the ordering (0,1,2,3,4,5,6), W4b = P^T (W4b_local (x) I_3) P where P permutes (0,4,5,6)->(0,1,2,3)
# Instead: build the full Wtot explicitly via the tensordot matvec, but only on the TRUNCATED basis.
# Truncated basis: keep each 3-link block in the low-j sectors (at most 1 link at j=1/2)
sd=[b for b in B if b[0]==0]  # j=0 sector (1 state)
e1=[b for b in B if b[0]==1]  # j=1/2 sector (4 states)
print("jm=%d: |j=0|=%d |j=1/2|=%d"%(jm,len(sd),len(e1)), flush=True)
# 3-link truncated basis: all j=0, or one link in j=1/2 (3*4 = 12 states)
basis3=[]
basis3.append((0,0,0))
for pos in range(3):
    for s in e1:
        tup=[0,0,0]; tup[pos]=s; basis3.append(tuple(tup))
print("3-link truncated basis size = %d (full %d)"%(len(basis3),d**3), flush=True)
# variational basis: |c> (j=0 only for simplicity) x |l> x |r>
basis=[]
for l in basis3:
    for r in basis3:
        basis.append((0,l,r))
print("variational basis size = %d"%len(basis), flush=True)
Nb=len(basis)
def lift(vb):
    # vb is a vector on the variational basis -> full vector
    full=np.zeros(d**7)
    for i,(c,l,r) in enumerate(basis):
        ii=index[(tuple(l)[0] if False else l[0])] if False else None
        # compute the full index
        ii=0
        for b in ( (0,0,0), l[0],l[1],l[2], r[0],r[1],r[2] ):
            ii=ii*d+index[b]
        full[ii]+=vb[i]
    return full
def project(full):
    out=np.zeros(Nb)
    for i,(c,l,r) in enumerate(basis):
        ii=0
        for b in ( (0,0,0), l[0],l[1],l[2], r[0],r[1],r[2] ):
            ii=ii*d+index[b]
        out[i]=full[ii]
    return out
print("basis ready", flush=True)
# The truncated 4-link operator on (0,1,2,3), with link 0 restricted to j=0 and links 1,2,3 to basis3
# Build the projection matrices
P0=np.zeros((d,1)); P0[0,0]=1.0
P3=np.zeros((d**3,len(basis3)))
for i,t in enumerate(basis3):
    ii=0
    for b in t: ii=ii*d+index[b]
    P3[ii,i]=1.0
# W4a restricted: P^T W4a P where P = kron(P0, P3)
P=sp.kron(sp.csr_matrix(P0),sp.csr_matrix(P3),format='csr')
W4a_tr=(P.T @ W4a @ P).tocsr()
# W4b acts on links (0,4,5,6) = (0, r-block). In the (c,l,r) basis, it acts on (c,r).
W4b_tr=(P.T @ W4a @ P).tocsr()   # same operator form, on (c,r)
print("truncated W4 size = %d, nnz=%d"%(W4a_tr.shape[0],W4a_tr.nnz), flush=True)
# full truncated operator on the (l,r) "matrix" of shape (len(basis3), len(basis3))
def matvec(v):
    V=v.reshape(len(basis3),len(basis3))
    return (W4a_tr @ V).ravel()
# For W4b acting on (c,r): the c index is j=0 (dim 1), so W4b acts only on r
# In the basis (l,r), W4b acts on r only: V -> V @ W4b_tr^T where W4b_tr acts on r
def matvec_b(v):
    V=v.reshape(len(basis3),len(basis3))
    return (W4a_tr @ V).ravel()  # same as before because c is j=0 in both
# NOTE: with c fixed to j=0, W1 and W2 act identically on the l/r blocks (symmetric)
print("NOTE: c restricted to j=0 makes W1 and W2 symmetric -> Wtot = 0.5(W4a_tr on l + W4a_tr on r)", flush=True)
