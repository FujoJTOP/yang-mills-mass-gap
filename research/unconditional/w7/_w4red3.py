import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=2
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; dl=len(B)
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((dl,dl))
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
uops={}
for M2 in (-1,1):
    for N2 in (-1,1):
        uops[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
        uops[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
Idx=[(-1,-1),(-1,1),(1,-1),(1,1)]
tot=None
for (p,q) in Idx:
    A=uops[(p,q)]
    for (r,s) in Idx:
        Bm=uops[(q,r,'c')]; C=uops[(r,s)]; Dm=uops[(s,p,'c')]
        term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
        tot=term if tot is None else tot+term
W4=(0.5*tot).tocsr()
D4=dl**4
print("basis B for jm=%d:"%jm, flush=True)
for j in range(jm+1):
    print("  j=%d: m in %s"%(j,[b[1] for b in B if b[0]==j and b[1]==b[2]]), flush=True)
# e_J
nJ=jm+1
evecs=[]
indices=[]
for J in range(0,jm+1):
    v=np.zeros(D4); idxs=[]
    for i1,(j1,m1,n1) in enumerate(B):
        if j1!=J or m1!=n1: continue
        for i2,(j2,m2,n2) in enumerate(B):
            if j2!=J or m2!=n2: continue
            for i3,(j3,m3,n3) in enumerate(B):
                if j3!=J or m3!=n3: continue
                for i4,(j4,m4,n4) in enumerate(B):
                    if j4!=J or m4!=n4: continue
                    ii=((i1*dl+i2)*dl+i3)*dl+i4
                    v[ii]=1.0; idxs.append(ii)
    if len(idxs)==0:
        v[:]=0
    else:
        v=v/np.sqrt(len(idxs))
    evecs.append(v); indices.append(idxs)
    print("  J=%d: number of (m,m)-states = %d"%(J,len(idxs)), flush=True)
Mred=np.zeros((nJ,nJ))
for a in range(nJ):
    for b in range(nJ):
        Mred[a,b]=float(evecs[a] @ (W4 @ evecs[b]))
print()
print("Mred =", flush=True)
with np.printoptions(precision=8,suppress=True,linewidth=200):
    print(Mred, flush=True)
ev=eigsh(W4,k=1,which="LA",return_eigenvectors=False,tol=1e-12)[0]
print("exact W_max =",ev, flush=True)
# test the candidate eigenvector (0.5, 0.7071, 0.5)
c=np.array([0.5,np.sqrt(0.5),0.5])
c=c/np.linalg.norm(c)
print("Mred @ c =", Mred@c, "  c*lambda with lambda=%.6f:"%ev, ev*c, flush=True)
