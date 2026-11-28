import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (2,3):
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
    # subspace: all four links at the SAME j, with uniform (m,n) choice over all (2j+1)^2 combos
    # basis vector for each J: e_J = (1/(2J+1)^2) sum_{m1,n1,...} |J,m1,n1>|J,m2,n2>|J,m3,n3>|J,m4,n4>
    nJ=jm+1
    evecs={}
    for J in range(0,jm+1):
        v=np.zeros(D4); cnt=0
        for i1,(a1,b1,c1) in enumerate(B):
            if a1!=J: continue
            for i2,(a2,b2,c2) in enumerate(B):
                if a2!=J: continue
                for i3,(a3,b3,c3) in enumerate(B):
                    if a3!=J: continue
                    for i4,(a4,b4,c4) in enumerate(B):
                        if a4!=J: continue
                        ii=((i1*dl+i2)*dl+i3)*dl+i4
                        v[ii]=1.0; cnt+=1
        evecs[J]=v/np.sqrt(cnt)
    Mred=np.zeros((nJ,nJ))
    for a in range(nJ):
        for b in range(nJ):
            Mred[a,b]=float(evecs[a] @ (W4 @ evecs[b]))
    print("jm=%d : reduced matrix in the 'all-same-j' subspace:"%(jm), flush=True)
    with np.printoptions(precision=8,suppress=True,linewidth=200):
        print(Mred, flush=True)
    ev=np.linalg.eigvalsh(0.5*(Mred+Mred.T))
    print("   reduced eigenvalues:", np.round(np.sort(ev)[::-1],8), flush=True)
    print("   cos(k pi/(jm+2)) k=1..jm+1:", np.round(np.sort(np.cos(np.arange(1,jm+2)*np.pi/(jm+2)))[::-1],8), flush=True)
    # exact W_max
    evW=eigsh(W4,k=1,which="LA",return_eigenvectors=False,tol=1e-12)[0]
    print("   exact W_max = %.10f"%(evW), flush=True)
