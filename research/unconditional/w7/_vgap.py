import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
def make(jm):
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
    V=None
    for mu in range(4):
        for nu in range(mu+1,4):
            P=pair_op(mu,nu); V=P if V is None else (V+P).tocsr()
    return dl,(0.5*(V+V.T)).tocsr()
print("=== V spectral gap vs jm : minimal 4D lattice L=(1,1,1,1) ===",flush=True)
print("    m a ~ (lambda_1 - lambda_2) * beta_s  at large beta_s",flush=True)
print()
print("%-4s %-5s %-16s %-16s %-16s %-12s"%("jm","D","lambda_1","lambda_2","gap=l1-l2","gap/l1"),flush=True)
for jm in (1,2,3,4):
    dl,V=make(jm); D=dl**4
    k=3 if D>10 else 2
    ev=np.sort(eigsh(V,k=3,which="LA",return_eigenvectors=False,tol=1e-10))[::-1]
    print("%-4d %-5d %-16.10f %-16.10f %-16.10f %-12.8f"%(
        jm,D,ev[0],ev[1],ev[0]-ev[1],(ev[0]-ev[1])/ev[0]),flush=True)
print()
print("[%.1fs]"%(time.time()-t0),flush=True)
