import numpy as np, functools, scipy.sparse as sp, time, itertools
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
    return dl,V
pairs=[(mu,nu) for mu in range(4) for nu in range(mu+1,4)]
print("=== S_4 (direction-permutation) reduction of V on L=(1,1,1,1) ===",flush=True)
print("    trivial rep dim = number of orbits of S_4 on basis states",flush=True)
print()
for jm in (1,2,3):
    dl,V=build(jm); D=dl**4
    # orbit reps: sorted tuple of the 4 link indices (S_4 acts by permutation)
    reps={}
    for i in range(D):
        x=i; t=[]
        for _ in range(4): t.append(x%dl); x//=dl
        key=tuple(sorted(t))
        reps.setdefault(key,[]).append(i)
    keys=sorted(reps.keys())
    nT=len(keys)
    # trivial-rep basis: normalized orbit sums
    Bs=[]
    for k in keys:
        v=np.zeros(D)
        for i in reps[k]: v[i]=1.0
        v/=np.linalg.norm(v); Bs.append(v)
    Bm=np.array(Bs)                    # nT x D
    Vt=Bm @ (V @ Bm.T)                 # nT x nT
    Vt=0.5*(Vt+Vt.T)
    evt=np.linalg.eigvalsh(Vt)[::-1]
    ev_exact=np.sort(eigsh(V,k=3,which="LA",return_eigenvectors=False,tol=1e-10))[::-1]
    print("jm=%d dl=%d D=%d  n_trivial_orbits=%d"%(jm,dl,D,nT),flush=True)
    print("   exact top3      : %.10f %.10f %.10f"%(ev_exact[0],ev_exact[1],ev_exact[2]),flush=True)
    print("   trivial-rep top : %.10f  (matches lambda_1? %s)"%(
        evt[0], "YES" if abs(evt[0]-ev_exact[0])<1e-8 else "no"),flush=True)
    print("   trivial-rep top3: %.10f %.10f %.10f"%(evt[0],evt[1],evt[2] if nT>2 else np.nan),flush=True)
    print("   [%.1fs]"%(time.time()-t0),flush=True)
