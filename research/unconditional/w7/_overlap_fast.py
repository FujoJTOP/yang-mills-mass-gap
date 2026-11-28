import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
print("=== FAST: restrict W2 to the j-only subspace FIRST (dl x dl), then kron (nJ x nJ) ===", flush=True)
print("%-4s %-4s %-16s %-16s %-10s %-8s"%("jm","d","W_max(V) full","W_max(V) j-only","ratio","sec"), flush=True)
for jm in (1,2,3,4,5,6):
    B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
    index={b:i for i,b in enumerate(B)}; dl=len(B); nJ=jm+1
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
            term=sp.kron(A,Bm,format='csr'); tot=term if tot is None else tot+term
    W2=(0.5*tot).tocsr()
    # full-spectrum reference (only feasible for jm<=5)
    if jm<=5:
        Wp1=sp.kron(W2, sp.identity(dl,format='csr'), format='csr')
        Wp2=sp.kron(sp.identity(dl,format='csr'), W2, format='csr')
        V=(Wp1+Wp2).tocsr()
        ev_full=float(eigsh(V,k=1,which="LA",return_eigenvectors=False,tol=1e-10)[0])
        del Wp1,Wp2,V
    else:
        ev_full=float('nan')
    # j-only: restrict W2 to the (m,m)-uniform basis (dl-dim) -> nJ x nJ
    evs=[]
    for J in range(nJ):
        v=np.zeros(dl); cnt=0
        for i,(j,m,n) in enumerate(B):
            if j==J and m==n: v[i]=1.0; cnt+=1
        evs.append(v/np.sqrt(cnt))
    En=np.array(evs)
    En2=np.kron(En,En)               # (nJ^2) x (dl^2)
    W2j=En2 @ (W2 @ En2.T)           # nJ x nJ? no: nJ^2 x nJ^2
    W2j=0.5*(W2j+W2j.T)   # nJ^2 x nJ^2
    Ij=np.eye(nJ)
    Vj=np.kron(W2j,Ij)+np.kron(Ij,W2j)   # (nJ^3) x (nJ^3), tiny
    evr=float(np.linalg.eigvalsh(0.5*(Vj+Vj.T))[-1])
    print("%-4d %-4d %-16.10f %-16.10f %-10.6f %-8.1f"%(jm,dl,ev_full,evr,evr/ev_full,time.time()-t0), flush=True)

