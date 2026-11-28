import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (1,2,3,4):
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
    # single-link operator (1/2)tr(u) = (1/2) sum_a u_{aa}
    U1=np.zeros((dl,dl))
    for a in (-1,1):
        U1 += mult_op(1,a,a)
    U1=0.5*U1
    ev=np.linalg.eigvalsh(0.5*(U1+U1.T))
    print("jm=%d : single-link (1/2)tr(u) spectrum: min=%.8f max=%.8f  (cos(pi/(jm+2))=%.8f)"%(
        jm,ev.min(),ev.max(),np.cos(np.pi/(jm+2))), flush=True)
    # also the 2-link W_2 = (1/2)tr(u1 u2^-1)
    Ip=sp.identity(dl,format='csr')
    tot=None
    for M2 in (-1,1):
        for N2 in (-1,1):
            A=sp.csr_matrix(mult_op(1,M2,N2)); Ac=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
            term=sp.kron(A,Ac,format='csr')
            tot=term if tot is None else tot+term
    W2=(0.5*tot).tocsr()
    ev2=np.linalg.eigvalsh(0.5*(W2.toarray()+W2.toarray().T))
    print("     2-link W_2=(1/2)tr(u1u2^-1) spectrum: min=%.8f max=%.8f"%(ev2.min(),ev2.max()), flush=True)
