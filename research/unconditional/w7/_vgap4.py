import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
jm=4
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; dl=len(B)
print("jm=%d dl=%d D=%d"%(jm,dl,dl**4),flush=True)
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
    return sp.csr_matrix(ph*A)
U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
Ui={(r,s):mult_op(1,r,s,conj=True) for r in (-1,1) for s in (-1,1)}
I=sp.identity(dl,format='csr')
def pair_op(mu,nu):
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        A=(U[(p,q)] @ Ui[(r,s)]).tocsr()
        Bm=(U[(q,r)] @ Ui[(s,p)]).tocsr()
        term=None
        for l in range(4):
            M = A if l==mu else (Bm if l==nu else I)
            term=M if term is None else sp.kron(term,M,format='csr')
        acc=term if acc is None else acc+term
    return (0.5*acc).tocsr()
print("building pair ops ...",flush=True)
pairs=[(mu,nu) for mu in range(4) for nu in range(mu+1,4)]
Ps=[pair_op(mu,nu) for mu,nu in pairs]
print("pair ops built %.1fs"%(time.time()-t0),flush=True)
# Hermiticity check on P00
P0=Ps[0]
print("P hermitian? max|P-P^T| = %.3e"%(abs(P0-P0.T)).max() if hasattr(abs(P0-P0.T),'max') else 0,flush=True)
d=(P0-P0.T).tocoo()
print("P hermitian? nnz of P-P^T = %d ; max|val| = %.3e"%(d.nnz, np.abs(d.data).max() if d.nnz else 0.0),flush=True)
# accumulate V without symmetrization
V=Ps[0]
for P in Ps[1:]:
    V=(V+P).tocsr()
print("V built %.1fs nnz=%d (%0.2f GB)"%(time.time()-t0,V.nnz,V.nnz*12/1e9),flush=True)
ev=np.sort(eigsh(V,k=3,which="LA",return_eigenvectors=False,tol=1e-10))[::-1]
print("jm=4 : lambda_1=%.10f lambda_2=%.10f lambda_3=%.10f gap=%.10f gap/l1=%.8f"%(
    ev[0],ev[1],ev[2],ev[0]-ev[1],(ev[0]-ev[1])/ev[0]),flush=True)
print("total %.1fs"%(time.time()-t0),flush=True)
