import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=3
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
# W_2 = (1/2) tr(u1 u2^-1)
tot=None
for (p,q) in Idx:
    A=uops[(p,q)]
    for (r,s) in Idx:
        Bm=uops[(q,r,'c')]
        term=sp.kron(A,Bm,format='csr')
        tot=term if tot is None else tot+term
W2=(0.5*tot).tocsr()
D2=dl**2
ev,V=eigsh(W2,k=1,which="LA",tol=1e-13)
v=V[:,0]
print("W_2 (2-link) jm=%d: W_max=%.10f  (cos(pi/%d)=%.10f)"%(jm,ev[0],jm+2,np.cos(np.pi/(jm+2))), flush=True)
from collections import defaultdict
groups=defaultdict(list)
for i in range(D2):
    idx=i
    b1=B[idx%dl]; b2=B[idx//dl]
    groups[(b1[0],b2[0])].append((((b1[1],b1[2]),(b2[1],b2[2]))),float(v[i]))
for k in sorted(groups):
    nz=[(mn,val) for (mn,val) in groups[k] if abs(val)>1e-8]
    if nz:
        print("j-tuple %s: %d nonzero"%(k,len(nz)), flush=True)
        for mn,val in nz[:8]:
            print("    %s : %.8f"%(mn,val), flush=True)
        # count by sign
        pos=sum(1 for _,x in nz if x>0); neg=len(nz)-pos
        print("    (%d positive, %d negative)"%(pos,neg), flush=True)
