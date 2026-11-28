import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
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
    # 2-link W_2 (correct)
    tot=None
    for a in (-1,1):
        for b in (-1,1):
            A=sp.csr_matrix(mult_op(1,a,b)); Bm=sp.csr_matrix(mult_op(1,a,b,conj=True))
            term=sp.kron(A,Bm,format='csr')
            tot=term if tot is None else tot+term
    W2=(0.5*tot).tocsr()
    D2=dl**2
    bt=2.5; bs=1.0
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    js=np.array([td[b[0]] for b in B])
    Kdiag=np.ones(D2)
    for ii in range(D2):
        x=ii; p=1.0
        for _ in range(2):
            p*=js[x%dl]; x//=dl
        Kdiag[ii]=p
    Kd=sp.diags(Kdiag)
    def Mf(v): return np.exp(-2*bs)*expm_multiply(bs*W2, v)
    e0=np.zeros(D2); e0[0]=1.0
    w0=Mf(e0); scale=float(np.dot(w0,Kdiag*w0))
    def T(v): return Mf(Kd*(Mf(v)))/scale
    op=LinearOperator((D2,D2),matvec=T,dtype=float)
    ev=eigsh(op,k=2,which="LA",return_eigenvectors=False,tol=1e-12)
    ev=np.sort(ev)[::-1]
    # j-only subspace: for each (J1,J2), the uniform (m=m,n=n) over all combos
    basis=[]
    for J1 in range(0,jm+1):
        for J2 in range(0,jm+1):
            rows=[]
            for i1,(a1,b1,c1) in enumerate(B):
                if a1!=J1 or b1!=c1: continue
                for i2,(a2,b2,c2) in enumerate(B):
                    if a2!=J2 or b2!=c2: continue
                    rows.append(i1*dl+i2)
            basis.append((J1,J2,rows))
    Nred=len(basis)
    # restricted W2
    Wr=np.zeros((Nred,Nred))
    for a in range(Nred):
        va=np.zeros(D2); va[basis[a][2]]=1.0; va/=np.linalg.norm(va)
        for b in range(Nred):
            vb=np.zeros(D2); vb[basis[b][2]]=1.0; vb/=np.linalg.norm(vb)
            Wr[a,b]=float(va @ (W2 @ vb))
    # restricted K diag
    Kr=np.array([td[J1]*td[J2] for (J1,J2,_) in basis])
    # reduced T
    from scipy.linalg import expm
    Mr=np.exp(-2*bs)*expm(bs*Wr)
    e0r=np.zeros(Nred); e0r[0]=1.0
    w0r=Mr@e0r; scale_r=float(w0r @ (Kr*w0r))
    Tr=(Mr @ np.diag(Kr) @ Mr)/scale_r
    evr=np.linalg.eigvalsh(0.5*(Tr+Tr.T))[::-1]
    print("jm=%d : W_2 exact t1=%.10f t2=%.10f  r=%.10f"%(jm,ev[0],ev[1],ev[1]/ev[0]), flush=True)
    print("     j-subspace t1=%.10f t2=%.10f  r=%.10f  (Nred=%d)"%(evr[0],evr[1],evr[1]/evr[0],Nred), flush=True)
    print("     ratio t1_sub/t1_exact = %.8f"%(evr[0]/ev[0]), flush=True)
    print()
