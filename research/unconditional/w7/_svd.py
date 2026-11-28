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
tot=None
for a in (-1,1):
    for b in (-1,1):
        A=sp.csr_matrix(mult_op(1,a,b)); Bm=sp.csr_matrix(mult_op(1,a,b,conj=True))
        term=sp.kron(A,Bm,format='csr')
        tot=term if tot is None else tot+term
W2=(0.5*tot).tocsr()
D2=dl**2
ev,V=eigsh(W2,k=1,which="LA",tol=1e-13)
v=V[:,0]
M=v.reshape(dl,dl)
sv=np.linalg.svd(M,compute_uv=False)
print("W_2 max eig=%.10f"%ev[0], flush=True)
print("Schmidt singular values of the max eigenvector (reshaped dl x dl):", flush=True)
print("  ", np.round(sv,6), flush=True)
print("  number > 1e-6:", int(np.sum(sv>1e-6)), flush=True)
print("  top singular vector (left):", flush=True)
u,s,vt=np.linalg.svd(M)
ul=u[:,0]
print("   ",[(B[i],round(float(ul[i]),6)) for i in range(dl) if abs(ul[i])>1e-6][:12], flush=True)
