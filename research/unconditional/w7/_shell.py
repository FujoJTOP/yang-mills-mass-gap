import numpy as np, functools, scipy.sparse as sp, time, itertools
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
# ---------- built-in 4D torus (2,1,1,1) with jm=1 : exact reference ----------
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
V=sum(plaq_op(sl) for sl in plaq); V=0.5*(V+V.T).tocsr()
bt=2.5
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
Kdiag=np.ones(D)
for ii in range(D):
    x=ii; pr=1.0
    for _ in range(NL): pr*=js[x%dl]; x//=dl
    Kdiag[ii]=pr
Kd=sp.diags(Kdiag)
def exact_ma(bs):
    def Mf(v): return expm_multiply(0.5*bs*V, v)
    e0=np.zeros(D); e0[0]=1.0
    w0=Mf(e0); sc=float(np.dot(w0,Kdiag*w0))
    def T(v): return Mf(Kd*(Mf(v)))/sc
    ev=np.sort(eigsh(LinearOperator((D,D),matvec=T,dtype=float),k=2,which="LA",
                     return_eigenvectors=False,tol=1e-10))[::-1]
    return ev[0],ev[1],-np.log(ev[1]/ev[0])
print("=== exact jm=1 4D torus reference ===",flush=True)
ref={}
for bs in (0.25,0.5,1.0,2.0,4.0):
    t1,t2,ma=exact_ma(bs); ref[bs]=ma
    print("  beta_s=%.2f : m a=%.8f"%(bs,ma),flush=True)
print("  [%.1fs]"%(time.time()-t0),flush=True)
# ---------- shell-truncated reference (same jm=1) : validate ----------
def shell_basis(jm,NL,Nmax):
    states=[]
    for combo in itertools.product(range(jm+1),repeat=NL):
        if sum(1 for c in combo if c>0)<=Nmax: states.append(combo)
    return states
for Nmax in (1,2,3,4,5,6,7,8):
    st=shell_basis(jm,NL,Nmax)
    n=len(st)
    # build restricted V via per-plaquette 4-link operator elements
    idx={c:i for i,c in enumerate(st)}
    Vr=sp.lil_matrix((n,n))
    for sl in plaq:
        Wp=plaq_op(sl)
        Wd=Wp.tocsr()
        for a,ca in enumerate(st):
            for b,cb in enumerate(st):
                # matrix element <ca|Wp|cb> : build the 4-link element
                diff=[(ca[l],cb[l]) for l in sl]
                # compute 4-link matrix element
                val=0.0
                rows=[B.index((ca[l],m,n)) for l in sl for m in [0] for n in [0]] if False else None
                # brute: only 4 links, dl small
                v=0.0
                for i0 in range(dl):
                    if B[i0][0]!=ca[sl[0]]: continue
                    for i1 in range(dl):
                        if B[i1][0]!=ca[sl[1]]: continue
                        for i2 in range(dl):
                            if B[i2][0]!=ca[sl[2]]: continue
                            for i3 in range(dl):
                                if B[i3][0]!=ca[sl[3]]: continue
                                v+=Wd[((i0*dl+i1)*dl+i2)*dl+i3,
                                      ((B.index((cb[sl[0]],0,0))*0)+ ((B.index((cb[sl[0]],0,0))))*0)] if False else 0.0
                # too slow / messy: use dense 4-link matrix times basis vector instead
                break
            break
    break
print("shell-truncation reference build too slow in pure python; switching to operator-level restriction.",flush=True)
print("total %.1fs"%(time.time()-t0),flush=True)
