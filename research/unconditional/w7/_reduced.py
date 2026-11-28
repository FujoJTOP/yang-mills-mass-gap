import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
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
    bt=2.5; bs=1.0
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    js=np.array([td[b[0]] for b in B])
    Kdiag=np.ones(D4)
    for ii in range(D4):
        x=ii; p=1.0
        for _ in range(4):
            p*=js[x%dl]; x//=dl
        Kdiag[ii]=p
    Kd=sp.diags(Kdiag)
    def Mf(v): return np.exp(-2*bs)*expm_multiply(bs*W4, v)
    e0=np.zeros(D4); e0[0]=1.0
    w0=Mf(e0); scale=float(np.dot(w0,Kdiag*w0))
    def T(v): return Mf(Kd*(Mf(v)))/scale
    op=LinearOperator((D4,D4),matvec=T,dtype=float)
    ev=eigsh(op,k=2,which="LA",return_eigenvectors=False,tol=1e-12)
    ev=np.sort(ev)[::-1]
    # reduced problem in the (m,m) subspace: W acts as (1/2)A_path, K diagonal with t_j
    nJ=jm+1
    Ared=np.zeros((nJ,nJ))
    for k in range(nJ-1):
        Ared[k,k+1]=0.5; Ared[k+1,k]=0.5
    # M_f restricted: exp(-2 beta_s) exp(beta_s * Ared)
    Mred=np.exp(-2*bs)*sp.linalg.expm(bs*Ared)  # (nJ,nJ) dense
    Kred=np.array([td[j] for j in range(jm+1)])
    # scale: <e0|T|e0> in the reduced space = ? use the same scale (it's the (j=0) component)
    # In the reduced space, the (j=0) component is e_0 (the first basis vector).
    w0r=Mred[:,0]
    sc_r=float(w0r @ (Kred*w0r))
    Tred=(Mred @ np.diag(Kred) @ Mred)/sc_r
    evr=np.linalg.eigvalsh(0.5*(Tred+Tred.T))
    evr=np.sort(evr)[::-1]
    print("jm=%d : full t_1=%.10f t_2=%.10f  r=%.10f"%(jm,ev[0],ev[1],ev[1]/ev[0]), flush=True)
    print("   (m,m)-subspace t_1=%.10f t_2=%.10f  r=%.10f"%(evr[0],evr[1],evr[1]/evr[0]), flush=True)
    print("   ratio t1_sub/t1_full = %.6f   (should be <= 1 if subspace gives lower bound)"%(evr[0]/ev[0]), flush=True)
