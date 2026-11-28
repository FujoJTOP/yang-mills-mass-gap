import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh, LinearOperator
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
def make(jm):
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
        return sp.csr_matrix(ph*A)
    U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
    Ui={(r,s):mult_op(1,r,s,conj=True) for r in (-1,1) for s in (-1,1)}
    # 2-link plaquette operator (same for every pair, by gauge invariance)
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        A=(U[(p,q)] @ Ui[(r,s)]).tocsr()
        Bm=(U[(q,r)] @ Ui[(s,p)]).tocsr()
        term=sp.kron(A,Bm,format='csr')
        acc=term if acc is None else acc+term
    P2=(0.5*acc).toarray()   # (dl^2)x(dl^2), dense, small
    return dl,P2
pairs=[(mu,nu) for mu in range(4) for nu in range(mu+1,4)]
for jm in (4,5):
    dl,P2=make(jm); D=dl**4
    print("jm=%d dl=%d D=%d  P2=%dx%d (%.0f MB) [%.1fs]"%(jm,dl,D,P2.shape[0],P2.shape[0],
          P2.nbytes/1e6,time.time()-t0),flush=True)
    def applyP2(x,mu,nu):
        xt=np.moveaxis(x,[mu,nu],[2,3])
        sh=xt.shape
        xt2=xt.reshape(-1,dl*dl)
        yt=(xt2 @ P2.T).reshape(sh)
        return np.moveaxis(yt,[2,3],[mu,nu])
    def Vmat(v):
        x=v.reshape([dl]*4)
        out=np.zeros_like(x)
        for mu,nu in pairs:
            out+=applyP2(x,mu,nu)
        return out.ravel()
    op=LinearOperator((D,D),matvec=Vmat,dtype=float)
    # sanity: V|0> norm
    e0=np.zeros(D); e0[0]=1.0
    v0=Vmat(e0)
    print("   V|0> : nnz=%d  norm=%.6f"%(np.count_nonzero(np.abs(v0)>1e-14),np.linalg.norm(v0)),flush=True)
    ev=np.sort(eigsh(op,k=3,which="LA",return_eigenvectors=False,tol=1e-8,maxiter=3000))[::-1]
    print("   jm=%d : l1=%.10f l2=%.10f l3=%.10f  gap=%.10f  rel=%.8f"%(
        jm,ev[0],ev[1],ev[2],ev[0]-ev[1],(ev[0]-ev[1])/ev[0]),flush=True)
    print("   [%.1fs]"%(time.time()-t0),flush=True)
