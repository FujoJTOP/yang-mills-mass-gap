import numpy as np, functools, scipy.sparse as sp, time, itertools
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
def P2of(jm):
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
    return dl,(0.5*acc).tocsr()
def orbit_size(ms):
    from math import factorial
    from collections import Counter
    n=1
    c=Counter(ms)
    for v in c.values(): n*=factorial(v)
    return 24//n
print("=== trivial-rep gap via V|_triv = 6 P_(01)  (bypasses full-space memory) ===",flush=True)
for jm in (5,):
    dl,P2=P2of(jm)
    # multiset basis
    keys=[ms for ms in itertools.combinations_with_replacement(range(dl),4)]
    kidx={ms:i for i,ms in enumerate(keys)}
    nT=len(keys)
    Ns=np.array([orbit_size(ms) for ms in keys],dtype=float)
    # build A = P_(01)|_triv :  A[o,o'] ; iterate o'
    rows=[];cols=[];vals=[]
    P2c=P2.tocsr()
    for oi,ms in enumerate(keys):
        # distinct permutations of ms and their multiplicities in the orbit sum (multiset perm counts)
        seen={}
        for pm in set(itertools.permutations(ms)):
            # multiplicity of this ordered tuple among the 24 permutations of ms
            from collections import Counter
            c=Counter(ms); norm=1
            for v in c.values():
                from math import factorial
                norm*=factorial(v)
            seen[pm]=1.0/np.sqrt(Ns[oi])   # each distinct tuple contributes 1 (orbit sum) / sqrt(N)
        for (a_,b_,c_,d_),w in seen.items():
            sub=a_*dl+b_
            st,en=P2c.indptr[sub],P2c.indptr[sub+1]
            for ptr in range(st,en):
                sa=P2c.indices[ptr]; val=P2c.data[ptr]
                a2,b2=sa//dl,sa%dl
                # NOTE: A = P_(01) acting on the ket -> matrix element <o|P|o'> ; we build <o| P |o'> by
                # applying P to |o'> and projecting.  Here we applied P to the ket's (a,b).
                om=tuple(sorted((a2,b2,c_,d_)))
                j2i=kidx.get(om)
                if j2i is None: continue
                rows.append(j2i); cols.append(oi); vals.append(w*val/np.sqrt(Ns[j2i]))
    A=sp.csr_matrix((vals,(rows,cols)),shape=(nT,nT))
    A=A/np.sqrt(1.0)  # normalization already handled per-column via w and row 1/sqrt(N_o)
    A=(0.5*(A+A.T)).tocsr() if False else A
    # correct symmetrization: A as built is <o|P|o'> with P hermitian -> symmetrize carefully
    Ah=0.5*(A+A.T)
    Vt=(6.0*Ah).tocsr()
    kk=4 if nT>10 else 2
    ev=np.sort(eigsh(Vt,k=kk,which="LA",return_eigenvectors=False,tol=1e-10,maxiter=5000))[::-1]
    print("jm=%d dl=%d nT=%d nnz(Vt)=%d [%.1fs]"%(jm,dl,nT,Vt.nnz,time.time()-t0),flush=True)
    print("   trivial top4: %s"%np.round(ev,10),flush=True)
    print("   gap=%.10f  rel=%.8f"%(ev[0]-ev[1],(ev[0]-ev[1])/ev[0]),flush=True)

