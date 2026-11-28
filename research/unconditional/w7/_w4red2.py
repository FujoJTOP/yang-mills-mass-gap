import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (2,3,4):
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
    uops={}
    for M2 in (-1,1):
        for N2 in (-1,1):
            uops[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
            uops[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
    Idx=[(-1,-1),(-1,1),(1,-1),(1,1)]
    tot=None
    for (p,q) in Idx:
        A=uops[(p,q)]
        for (r,s) in Idx:
            Bm=uops[(q,r,'c')]; C=uops[(r,s)]; Dm=uops[(s,p,'c')]
            term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
            tot=term if tot is None else tot+term
    W4=(0.5*tot).tocsr()
    D4=dl**4
    # subspace: all four links at same J, all (m,m)-diagonal
    nJ=jm+1
    evecs=[]
    for J in range(0,jm+1):
        v=np.zeros(D4); cnt=0
        for i1,(j1,m1,n1) in enumerate(B):
            if j1!=J or m1!=n1: continue
            for i2,(j2,m2,n2) in enumerate(B):
                if j2!=J or m2!=n2: continue
                for i3,(j3,m3,n3) in enumerate(B):
                    if j3!=J or m3!=n3: continue
                    for i4,(j4,m4,n4) in enumerate(B):
                        if j4!=J or m4!=n4: continue
                        ii=((i1*dl+i2)*dl+i3)*dl+i4
                        v[ii]=1.0; cnt+=1
        evecs.append(v/np.sqrt(cnt))
    Mred=np.zeros((nJ,nJ))
    for a in range(nJ):
        for b in range(nJ):
            Mred[a,b]=float(evecs[a] @ (W4 @ evecs[b]))
    evW=eigsh(W4,k=1,which="LA",return_eigenvectors=False,tol=1e-12)[0]
    print("jm=%d : exact W_max=%.10f  (cos(pi/%d)=%.10f)"%(jm,evW,jm+2,np.cos(np.pi/(jm+2))), flush=True)
    with np.printoptions(precision=6,suppress=True,linewidth=200):
        print("   reduced M (all-same-J, (m,m)-diag):", flush=True)
        print(Mred, flush=True)
    ev=np.linalg.eigvalsh(0.5*(Mred+Mred.T))
    print("   reduced eigenvalues (desc):", np.round(np.sort(ev)[::-1],8), flush=True)
    print("   cos(k pi/(jm+2)):", np.round(np.cos(np.arange(1,jm+2)*np.pi/(jm+2)),8), flush=True)
    print()
