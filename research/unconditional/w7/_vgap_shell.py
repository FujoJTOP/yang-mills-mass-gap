import numpy as np, functools, scipy.sparse as sp, time, itertools
from scipy.sparse.linalg import eigsh
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
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        A=(U[(p,q)] @ Ui[(r,s)]).tocsr()
        Bm=(U[(q,r)] @ Ui[(s,p)]).tocsr()
        term=sp.kron(A,Bm,format='csr')
        acc=term if acc is None else acc+term
    P2=(0.5*acc).tocsr()
    return B,index,dl,P2
pairs=[(mu,nu) for mu in range(4) for nu in range(mu+1,4)]
print("=== shell-truncated V (validate N_exc=4 = FULL against exact) ===",flush=True)
print("%-4s %-6s %-7s %-16s %-16s %-16s"%("jm","Nexc","n","lambda_1","lambda_2","gap"),flush=True)
for jm in (1,2,3,4,5):
    B,index,dl,P2=make(jm)
    z=[i for i,b in enumerate(B) if b[0]==0][0]
    nz=[i for i,b in enumerate(B) if b[0]>0]
    print("  jm=%d dl=%d nonzero-j=%d  P2 nnz=%d  [%.1fs]"%(jm,dl,len(nz),P2.nnz,time.time()-t0),flush=True)
    for Nexc in (2,3,4):
        if Nexc>4: continue
        states=[(z,)*4]
        for k in range(1,Nexc+1):
            for pos in itertools.combinations(range(4),k):
                for vals in itertools.product(nz,repeat=k):
                    s=[z]*4
                    for pp,vv in zip(pos,vals): s[pp]=vv
                    states.append(tuple(s))
        n=len(states)
        if n>2_500_000:
            print("%-4d %-6d %-7d SKIP (too big)"%(jm,Nexc,n),flush=True); continue
        idx={s:i for i,s in enumerate(states)}
        rows=[];cols=[];vals=[]
        for bi,b in enumerate(states):
            for (mu,nu) in pairs:
                sb=b[mu]*dl+b[nu]
                start,end=P2.indptr[sb],P2.indptr[sb+1]
                for ptr in range(start,end):
                    sa=P2.indices[ptr]; val=P2.data[ptr]
                    a=list(b); a[mu]=sa//dl; a[nu]=sa%dl
                    ai=idx.get(tuple(a))
                    if ai is not None:
                        rows.append(ai); cols.append(bi); vals.append(val)
        Vr=sp.csr_matrix((vals,(rows,cols)),shape=(n,n))
        Vr=(0.5*(Vr+Vr.T)).tocsr()
        kk=3 if n>10 else 2
        ev=np.sort(eigsh(Vr,k=kk,which="LA",return_eigenvectors=False,tol=1e-10))[::-1]
        print("%-4d %-6d %-7d %-16.10f %-16.10f %-16.10f"%(jm,Nexc,n,ev[0],ev[1],ev[0]-ev[1]),flush=True)
    print("  [%.1fs]"%(time.time()-t0),flush=True)
