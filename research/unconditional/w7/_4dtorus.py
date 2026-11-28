import numpy as np, functools, scipy.sparse as sp, time
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
def Uinv(r,s): return Umat[(s,r)].T.tocsr()   # real matrices: u^{-1}=u^T
I=sp.identity(dl,format='csr')
Ns=2; Nd=4; NL=Ns*Nd; D=dl**NL
def shift(site,d):
    L=2 if d==0 else 1
    return (site+1)%2 if L==2 else 0
def lk(site,d): return site*Nd+d
plaq=[]
for x in range(Ns):
    for mu in range(Nd):
        for nu in range(mu+1,Nd):
            plaq.append((lk(x,mu), lk(shift(x,mu),nu), lk(shift(x,nu),mu), lk(x,nu)))
print("4D torus (2,1,1,1): %d links, %d plaquettes, dl=%d, D=%d"%(NL,len(plaq),dl,D), flush=True)
def plaquette_op(slots):
    # slot 0: u ; slot 1: u ; slot 2: u^{-1} ; slot 3: u^{-1}
    acc=None
    for p in (-1,1):
        for q in (-1,1):
            for r in (-1,1):
                for s in (-1,1):
                    fac=[Umat[(p,q)], Umat[(q,r)], Uinv(r,s), Uinv(s,p)]
                    onlink={}
                    for k in range(4):
                        onlink.setdefault(slots[k],[]).append(fac[k])
                    term=None
                    for l in range(NL):
                        if l not in onlink:
                            M=I.copy()
                        else:
                            M=onlink[l][0]
                            for X in onlink[l][1:]: M=(M@X).tocsr()
                        term=M if term is None else sp.kron(term,M,format='csr')
                    acc=term if acc is None else acc+term
    return (0.5*acc).tocsr()
V=None
for sl in plaq:
    P=plaquette_op(sl)
    V=P if V is None else (V+P).tocsr()
print("built %.1fs nnz=%d"%(time.time()-t0,V.nnz), flush=True)
Vh=(0.5*(V+V.T)).tocsr()
ev=eigsh(Vh,k=3,which="LA",return_eigenvectors=False,tol=1e-10)
ev=np.sort(ev)[::-1]
print("W_max(V) = %.10f"%ev[0], flush=True)
print("top3:",np.round(ev,8), flush=True)
print("per-plaquette W_max=%.10f ; naive 12x=%.10f ; ratio=%.6f"%(
  np.cos(np.pi/3),12*np.cos(np.pi/3),ev[0]/(12*np.cos(np.pi/3))), flush=True)
print("total %.1fs"%(time.time()-t0), flush=True)

