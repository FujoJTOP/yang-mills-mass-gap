import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
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
        return sp.csr_matrix(ph*A)
    U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
    Ui={(r,s):mult_op(1,r,s,conj=True) for r in (-1,1) for s in (-1,1)}
    I=sp.identity(dl,format='csr')
    def pair_op(mu,nu):
        acc=None
        for p in (-1,1):
         for q in (-1,1):
          for r in (-1,1):
           for s in (-1,1):
            A=(U[(p,q)] @ Ui[(r,s)]).tocsr()
            Bm=(U[(q,r)] @ Ui[(s,p)]).tocsr()
            term=None
            for l in range(4):
                M = A if l==mu else (Bm if l==nu else I)
                term=M if term is None else sp.kron(term,M,format='csr')
            acc=term if acc is None else acc+term
        return (0.5*acc).tocsr()
    V=None
    for mu in range(4):
        for nu in range(mu+1,4):
            P=pair_op(mu,nu); V=P if V is None else (V+P).tocsr()
    return dl,(0.5*(V+V.T)).tocsr()
print("=== trivial-rep (S_4-invariant) spectrum of V, sparse orbit method ===",flush=True)
for jm in (1,2,3):
    dl,V=build(jm); D=dl**4
    # enumerate multisets and orbit sizes
    reps={}
    for i in range(D):
        x=i; t=[]
        for _ in range(4): t.append(x%dl); x//=dl
        reps.setdefault(tuple(sorted(t)),[]).append(i)
    keys=sorted(reps.keys()); nT=len(keys)
    rowi=[];coli=[];vals=[]
    for oi,k in enumerate(keys):
        m=reps[k]; nrm=1.0/np.sqrt(len(m))
        for i in m:
            rowi.append(oi); coli.append(i); vals.append(nrm)
    Bmat=sp.csr_matrix((vals,(rowi,coli)),shape=(nT,D))
    Vt=(Bmat @ (V @ Bmat.T)).tocsr()
    Vt=(0.5*(Vt+Vt.T)).tocsr()
    kk=4 if nT>10 else 2
    ev=np.sort(eigsh(Vt,k=kk,which="LA",return_eigenvectors=False,tol=1e-11))[::-1]
    ev_full=np.sort(eigsh(V,k=3,which="LA",return_eigenvectors=False,tol=1e-10))[::-1]
    print("jm=%d dl=%d D=%d  nT=%d  nnz(Vt)=%d  [%.1fs]"%(jm,dl,D,nT,Vt.nnz,time.time()-t0),flush=True)
    print("   trivial top4 : %s"%np.round(ev,10),flush=True)
    print("   full   top3  : %s"%np.round(ev_full,10),flush=True)
    print("   lambda_1 match? %s ; trivial gap=%.10f"%(
        abs(ev[0]-ev_full[0])<1e-8, ev[0]-ev[1]),flush=True)
