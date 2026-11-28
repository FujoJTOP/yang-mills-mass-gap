import numpy as np, functools, scipy.sparse as sp, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
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
def Wilson4_full():
    tot=None
    for (p,q) in idx:
        A=u[(p,q)]
        for (r,s) in idx:
            Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
            term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
            tot=term if tot is None else tot+term
    return (0.5*tot).tocsr()
W4=Wilson4_full()
# Truncated 4-link basis: all j=0, OR all four at j=1/2.
z=(0,0,0)
e12=[b for b in B if b[0]==1]
print("jm=%d d=%d  |e12|=%d"%(jm,d,len(e12)), flush=True)
basis4=[(z,z,z,z)]
for a in e12:
    for b in e12:
        for c in e12:
            for dd in e12:
                basis4.append((a,b,c,dd))
nb4=len(basis4)
print("truncated 4-link basis = %d (full %d)"%(nb4,d**4), flush=True)
# projection
rows=[]; cols=[]
for j,(a,b,c,dd) in enumerate(basis4):
    ii=0
    for t in (a,b,c,dd): ii=ii*d+index[t]
    rows.append(ii); cols.append(j)
P=sp.csr_matrix((np.ones(nb4),(rows,cols)),shape=(d**4,nb4))
W4r=(P.T @ W4 @ P).tocsr()
print("W4r nnz=%d (nonzero=%s)"%(W4r.nnz, W4r.nnz>0), flush=True)
# two plaquettes: W1 on (c, links1-3), W2 on (c, links4-6)
# variational space: (c, l) x r  where (c,l) is the 4-link truncated basis, r is the 3-link truncated basis
# for W2 we need (c, r) - but r is 3 links, plus c. So the 4-link truncated basis for (c,r) is the same basis4.
# The variational space must be consistent: use (c, l3) with l3 in basis3 (3-link truncated)
# But W4 needs 4 links. So we must use (c, links1-3) in basis4 AND (c, links4-6) in basis4.
# Problem: the two 4-link bases share c. The full space is (c, l3, r3) with (c,l3) in basis4 and (c,r3) in basis4.
# That's not a product. Use basis4 for the left plaquette and basis4 for the right, sharing c:
# The vectors are indexed by (i_left, i_right) with the c index consistent.
# Since basis4 always has c in B and the 3 links in {all j=0 or all j=1/2}, the c index can be any of d.
# Build the variational space: pairs (i,j) with c(i)==c(j).
import itertools
pairs=[]
for i,(a1,b1,c1,d1) in enumerate(basis4):
    for j,(a2,b2,c2,d2) in enumerate(basis4):
        if a1==a2: pairs.append((i,j))
print("pairs with matching c = %d"%len(pairs), flush=True)
# Build Pfull: (d^7) x len(pairs)
rows=[]; cols=[]
for k,(i,j) in enumerate(pairs):
    (a1,b1,c1,d1)=basis4[i]; (a2,b2,c2,d2)=basis4[j]
    ii=0
    for t in (a1,b1,c1,d1,a2,b2,c2): ii=ii*d+index[t]
    rows.append(ii); cols.append(k)
Pf=sp.csr_matrix((np.ones(len(pairs)),(rows,cols)),shape=(d**7,len(pairs)))
print("Pfull shape=%s nnz=%d"%(Pf.shape,Pf.nnz), flush=True)
# W1 acts on the first 4 links (0,1,2,3) -> in the full space it's W4 (x) I_3
# W2 acts on (0,4,5,6) -> P_perm
# Build W1_full and W2_full on the FULL space, then restrict
I3=sp.identity(d**3,format='csr')
W1_full=sp.kron(W4,I3,format='csr')
# W2_full via permutation
perm=np.arange(d**7)
# index (l0,l1,l2,l3,l4,l5,l6) -> we want (l0,l4,l5,l6,l1,l2,l3)
perm2=np.zeros(d**7,dtype=np.int64)
for l0 in range(d):
    for l1 in range(d):
        for l2 in range(d):
            for l3 in range(d):
                for l4 in range(d):
                    for l5 in range(d):
                        for l6 in range(d):
                            src=((((((l0*d+l1)*d+l2)*d+l3)*d+l4)*d+l5)*d+l6)
                            dst=((((((l0*d+l4)*d+l5)*d+l6)*d+l1)*d+l2)*d+l3)
                            perm2[dst]=src
# careful: define P such that (P v)[dst] = v[src]
Pp=sp.csr_matrix((np.ones(d**7),(np.arange(d**7),perm2)),shape=(d**7,d**7))
print("Pp built", flush=True)
Wtot_full=0.5*(W1_full + Pp.T @ W1_full @ Pp)
Wtot_r=(Pf.T @ Wtot_full @ Pf).tocsr()
print("Wtot_r shape=%s nnz=%d"%(Wtot_r.shape,Wtot_r.nnz), flush=True)
