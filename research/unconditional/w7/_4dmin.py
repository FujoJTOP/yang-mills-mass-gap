import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
def ops(jm):
    B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
    index={b:i for i,b in enumerate(B)}; dl=len(B)
    def mult_op(J2,M2,N2):
        A=np.zeros((dl,dl))
        for i2,(j2,m2,n2) in enumerate(B):
            m3,n3=m2-M2,n2-N2
            for j3 in range(abs(J2-j2),J2+j2+1,2):
                key=(j3,m3,n3)
                if key not in index: continue
                v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,M2,j2,-m2)*d3(j3,n3,J2,N2,j2,-n2)
                if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
        return A
    U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
    return B,dl,U
print("=== minimal genuine 4D lattice L=(1,1,1,1): 4 links, all 6 plaquette types ===",flush=True)
print("%-4s %-4s %-9s %-16s %-16s %-8s"%("jm","dl","D","W_max(V)","6*W_max(1p)","sec"),flush=True)
for jm in (1,2,3,4,5):
    B,dl,U=ops(jm)
    I=sp.identity(dl,format='csr')
    def pair_op(mu,nu):
        # (1/2) sum_{p,q,r,s}  [link mu: u_pq * uinv_rs] (x) [link nu: u_qr * uinv_sp] (x) I elsewhere
        acc=None
        for p in (-1,1):
         for q in (-1,1):
          for r in (-1,1):
           for s in (-1,1):
            A=sp.csr_matrix(U[(p,q)]) @ sp.csr_matrix(U[(s,r)].T)   # u_pq * u^{-1}_rs
            Bm=sp.csr_matrix(U[(q,r)]) @ sp.csr_matrix(U[(p,s)].T)  # u_qr * u^{-1}_sp
            term=None
            for l in range(4):
                M = A if l==mu else (Bm if l==nu else I)
                term=M if term is None else sp.kron(term,M,format='csr')
            acc=term if acc is None else acc+term
        return (0.5*acc).tocsr()
    D=dl**4
    V=None
    for mu in range(4):
        for nu in range(mu+1,4):
            P=pair_op(mu,nu); V=P if V is None else (V+P).tocsr()
    V=0.5*(V+V.T).tocsr()
    w=float(eigsh(V,k=1,which="LA",return_eigenvectors=False,tol=1e-9)[0])
    P01=pair_op(0,1)
    w1=float(eigsh(P01,k=1,which="LA",return_eigenvectors=False,tol=1e-9)[0])
    print("%-4d %-4d %-9d %-16.10f %-16.10f %-8.1f"%(jm,dl,D,w,6*w1,time.time()-t0),flush=True)
