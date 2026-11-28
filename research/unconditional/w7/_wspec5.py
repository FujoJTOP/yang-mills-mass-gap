import numpy as np, functools, scipy.sparse as sp, time, mpmath as mp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (5,):
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
    u={}
    for M2 in (-1,1):
        for N2 in (-1,1):
            u[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
            u[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
    idx=[(-1,-1),(-1,1),(1,-1),(1,1)]
    tot=None
    for (p,q) in idx:
        A=u[(p,q)]
        for (r,s) in idx:
            Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
            term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
            tot=term if tot is None else tot+term
    W4=(0.5*tot).tocsr()
    D4=dl**4
    print("jm=%d d=%d D=%d"%(jm,dl,D4), flush=True)
    t0=time.time()
    ev=eigsh(W4,k=3,which="LA",return_eigenvectors=False,tol=1e-10)
    print("   W top 3: %s  (%.0fs)"%(np.round(np.sort(ev)[::-1],8),time.time()-t0), flush=True)
    del W4
mp.mp.dps=30
print()
print("=== the closed-form conjecture: W_max(jm) = cos(pi/(jm+2)) ===")
for jm in (1,2,3,4,5,10):
    v=mp.cos(mp.pi/(jm+2))
    print("   jm=%2d : cos(pi/%d) = %s"%(jm,jm+2,mp.nstr(v,15)))

