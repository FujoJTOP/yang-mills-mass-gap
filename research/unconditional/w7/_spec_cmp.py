import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (1,2,3):
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
    # single link
    U1=np.zeros((dl,dl))
    for a in (-1,1): U1+=mult_op(1,a,a)
    U1=0.5*U1
    ev1=np.linalg.eigvalsh(0.5*(U1+U1.T))
    # W_2 correct
    tot=None
    for a in (-1,1):
        for b in (-1,1):
            A=sp.csr_matrix(mult_op(1,a,b)); Bm=sp.csr_matrix(mult_op(1,a,b,conj=True))
            term=sp.kron(A,Bm,format='csr')
            tot=term if tot is None else tot+term
    W2=(0.5*tot).tocsr()
    ev2=np.linalg.eigvalsh(0.5*(W2.toarray()+W2.toarray().T))
    print("jm=%d : single-link spec (distinct): %s"%(jm,np.round(np.unique(np.round(ev1,8)),8)), flush=True)
    print("       W_2 spec (distinct):         %s"%np.round(np.unique(np.round(ev2,8)),8), flush=True)
    # compare multisets
    ok = np.allclose(np.sort(ev1),np.sort(ev2),atol=1e-8) if len(ev1)==len(ev2) else False
    if len(ev1)!=len(ev2):
        # compare after normalizing multiplicities
        u1,c1=np.unique(np.round(ev1,8),return_counts=True)
        u2,c2=np.unique(np.round(ev2,8),return_counts=True)
        print("       distinct counts: single=%s  W2=%s"%(list(c1),list(c2)), flush=True)
        ok = np.allclose(u1,u2,atol=1e-8)
    print("       distinct spectra equal? %s"%ok, flush=True)
    print()
