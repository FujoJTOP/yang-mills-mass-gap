import numpy as np, functools, scipy.sparse as sp, time, itertools
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
jm=1
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; dl=len(B)
def mult_op(J2,M2,N2):
    A=np.zeros((dl,dl))
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2-M2,n2-N2
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,M2,j2,-m2)*d3(j3,n3,J2,N2,j2,-n2)
            if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
    return sp.csr_matrix(A)
Umat={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
def Uinv(r,s): return Umat[(s,r)].T.tocsr()
I=sp.identity(dl,format='csr')
Ns=2; Nd=4; NL=Ns*Nd; D=dl**NL
def shift(site,d):
    L=2 if d==0 else 1
    return (site+1)%2 if L==2 else 0
def lk(site,d): return site*Nd+d
plaq=[(lk(x,mu),lk(shift(x,mu),nu),lk(shift(x,nu),mu),lk(x,nu))
      for x in range(Ns) for mu in range(Nd) for nu in range(mu+1,Nd)]
def plaq_op(slots):
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        fac=[Umat[(p,q)],Umat[(q,r)],Uinv(r,s),Uinv(s,p)]
        onl={}
        for k in range(4): onl.setdefault(slots[k],[]).append(fac[k])
        term=None
        for l in range(NL):
            if l not in onl: M=I.copy()
            else:
                M=onl[l][0]
                for X in onl[l][1:]: M=(M@X).tocsr()
            term=M if term is None else sp.kron(term,M,format='csr')
        acc=term if acc is None else acc+term
    return (0.5*acc).tocsr()
print("=== overlap correction vs NUMBER of plaquettes (jm=1, 4D torus links, D=%d) ==="%D,flush=True)
print("%-6s %-16s %-16s %-10s"%("N_p","W_max(V)","N_p*W_max(1p)","ratio"),flush=True)
single=None
for Np in (1,2,3,4,6,8,10,12):
    V=None
    for k in range(Np):
        P=plaq_op(plaq[k]); V=P if V is None else (V+P).tocsr()
    V=0.5*(V+V.T).tocsr()
    w=float(eigsh(V,k=1,which="LA",return_eigenvectors=False,tol=1e-10)[0])
    if Np==1: single=w
    print("%-6d %-16.10f %-16.10f %-10.6f"%(Np,w,Np*single,w/(Np*single)),flush=True)
print("[%.1fs]"%(time.time()-t0),flush=True)
