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
# Build W_2 = (1/2) sum_{a,b} (u_1)_{ab} (u_2^{-1})_{ba}
# (u_1)_{ab}: mult_op(1, 2a, 2b)  [doubled indices; a,b in {-1/2,+1/2} -> doubled {-1,+1}]
# (u_2^{-1})_{ba} = conj(u_2_{ab}) = mult_op(1, 2a, 2b, conj=True)
# but note the conj=True version gives the conjugate matrix element in the (j,m,n) basis.
# Let's just build it directly: sum_{a,b} kron(u_ab, u_ba^{-1})
def u_ab(a,b):   # a,b in {-1,1} (doubled)
    return sp.csr_matrix(mult_op(1,a,b))
def uinv_ba(a,b): # (u^{-1})_{ba} = conj(u_{ab}) for SU(2)
    return sp.csr_matrix(mult_op(1,a,b,conj=True))
tot=None
for a in (-1,1):
    for b in (-1,1):
        A=u_ab(a,b)
        Bm=uinv_ba(a,b)
        term=sp.kron(A,Bm,format='csr')
        tot=term if tot is None else tot+term
W2=(0.5*tot).tocsr()
D2=dl**2
# hermiticity
print("jm=%d dl=%d D=%d"%(jm,dl,D2), flush=True)
Wa=W2.toarray()
print("   max|W-W^T| = %.3e"%np.abs(Wa-Wa.T).max(), flush=True)
ev=eigsh(W2,k=1,which="LA",return_eigenvectors=False,tol=1e-12)[0]
print("   W_2 W_max = %.10f   (cos(pi/(jm+2))=%.10f)"%(ev,np.cos(np.pi/(jm+2))), flush=True)
# also verify: (1/2)tr(u_1) as check = single link
U1=np.zeros((dl,dl))
for a in (-1,1): U1+=mult_op(1,a,a)
U1=0.5*U1
print("   single-link W_max = %.10f   (cos(pi/(jm+2))=%.10f)"%(np.linalg.eigvalsh(0.5*(U1+U1.T)).max(),np.cos(np.pi/(jm+2))), flush=True)
