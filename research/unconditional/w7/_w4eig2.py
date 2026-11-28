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
ev,V=eigsh(W4,k=1,which="LA",tol=1e-12)
v=V[:,0]
print("W_max=%.10f"%ev[0], flush=True)
# group by (j1,j2,j3,j4) and within, by (m1,n1,m2,n2,m3,n3,m4,n4)
from collections import defaultdict
groups=defaultdict(list)
for i in range(D4):
    idx=i; jj=[]; mn=[]
    for _ in range(4):
        b=B[idx%dl]; jj.append(b[0]); mn.append((b[1],b[2])); idx//=dl
    groups[tuple(jj)].append((tuple(mn),float(v[i])))
for k in sorted(groups):
    nz=[(mn,val) for (mn,val) in groups[k] if abs(val)>1e-8]
    if nz:
        print()
        print("j-tuple %s: %d nonzero"%(k,len(nz)), flush=True)
        for mn,val in nz[:20]:
            print("    (m,n) per link = %s : %.10f"%(mn,val), flush=True)
        if len(nz)>20: print("    ... (%d more)"%(len(nz)-20), flush=True)
