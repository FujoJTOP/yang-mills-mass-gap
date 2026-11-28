import numpy as np, functools, scipy.sparse as sp, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
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
Ns=2; Nd=4; NL=Ns*Nd; D=dl**NL
def shift(site,d):
    L=2 if d==0 else 1
    return (site+1)%2 if L==2 else 0
def lk(site,d): return site*Nd+d
# plaquette slots: (U_mu(x), U_nu(x+mu), U_mu(x+nu)^{-1}, U_nu(x)^{-1})
# slot types: 0=U, 1=U, 2=U^{-1}, 3=U^{-1}
plaq=[(lk(x,mu),lk(shift(x,mu),nu),lk(shift(x,nu),mu),lk(x,nu))
      for x in range(Ns) for mu in range(Nd) for nu in range(mu+1,Nd)]
print("4D torus (2,1,1,1): %d links, %d plaquettes, jm=%d, D=%d"%(NL,len(plaq),jm,D),flush=True)
def plaq_op(slots):
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        # slot0=U_pq  slot1=U_qr  slot2=(U^-1)_rs  slot3=(U^-1)_sp
        fac=[U[(p,q)],U[(q,r)],Ui[(r,s)],Ui[(s,p)]]
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
V=None
for sl in plaq:
    P=plaq_op(sl); V=P if V is None else (V+P).tocsr()
V=(0.5*(V+V.T)).tocsr()
print("V built %.1fs nnz=%d"%(time.time()-t0,V.nnz),flush=True)
e0=np.zeros(D); e0[0]=1.0
Vv=V@e0
print("V|0> : nnz=%d norm=%.6f (must be >1 for a genuine operator)"%(
    np.count_nonzero(np.abs(Vv)>1e-14),np.linalg.norm(Vv)),flush=True)
w=float(eigsh(V,k=1,which="LA",return_eigenvectors=False,tol=1e-10)[0])
print("W_max(V)=%.10f   (12*0.5=%.6f)"%(w,12*0.5),flush=True)
bt=2.5
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
Kdiag=np.ones(D)
for ii in range(D):
    x=ii; pr=1.0
    for _ in range(NL): pr*=js[x%dl]; x//=dl
    Kdiag[ii]=pr
Kd=sp.diags(Kdiag)
print()
print("=== mass gap (CORRECTED, 8-link torus) ===",flush=True)
for bs in (0.25,0.5,1.0,1.5,2.0):
    def Mf(v): return expm_multiply(0.5*bs*V, v)
    w0=Mf(e0); sc=float(np.dot(w0,Kdiag*w0))
    def T(v): return Mf(Kd*(Mf(v)))/sc
    ev=np.sort(eigsh(LinearOperator((D,D),matvec=T,dtype=float),k=2,which="LA",
                     return_eigenvectors=False,tol=1e-9))[::-1]
    r=ev[1]/ev[0]
    print("beta_s=%.2f : t1=%.10f t2=%.10f r=%.10f m a=%.10f"%(bs,ev[0],ev[1],r,-np.log(r)),flush=True)
print("[%.1fs]"%(time.time()-t0),flush=True)
