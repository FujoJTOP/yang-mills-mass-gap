import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def build(jm):
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
    return B,index,dl,mult_op
print("=== is W_max(V) attained inside the 'j-only' subspace? ===")
print("    j-only basis per link: e_J = uniform over (m,m) states of j=J")
print()
print("%-4s %-16s %-16s %-10s"%("jm","W_max(V) full","W_max(V) on j-only","ratio"))
for jm in (1,2,3,4,5):
    B,index,dl,mult_op=build(jm)
    tot=None
    for a in (-1,1):
        for b in (-1,1):
            A=sp.csr_matrix(mult_op(1,a,b)); Bm=sp.csr_matrix(mult_op(1,a,b,conj=True))
            term=sp.kron(A,Bm,format='csr'); tot=term if tot is None else tot+term
    W2=(0.5*tot).tocsr()
    Wp1=sp.kron(W2, sp.identity(dl,format='csr'), format='csr')
    Wp2=sp.kron(sp.identity(dl,format='csr'), W2, format='csr')
    V=(Wp1+Wp2).tocsr()
    ev_full=eigsh(V,k=1,which="LA",return_eigenvectors=False,tol=1e-11)[0]
    # j-only subspace: for each link, one vector per J (uniform over (m,m) states)
    nJ=jm+1
    evs=[]
    for J in range(nJ):
        v=np.zeros(dl); cnt=0
        for i,(j,m,n) in enumerate(B):
            if j==J and m==n: v[i]=1.0; cnt+=1
        evs.append(v/np.sqrt(cnt))
    basis=[]
    for J1 in range(nJ):
        for J2 in range(nJ):
            for J3 in range(nJ):
                basis.append(np.kron(np.kron(evs[J1],evs[J2]),evs[J3]))
    M=np.array([[float(u @ (V @ w)) for w in basis] for u in basis])
    evr=np.linalg.eigvalsh(0.5*(M+M.T))[-1]
    print("%-4d %-16.10f %-16.10f %-10.6f"%(jm,ev_full,evr,evr/ev_full), flush=True)
