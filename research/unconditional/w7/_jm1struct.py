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
U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
I=sp.identity(dl,format='csr')
def pair_op(mu,nu):
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        A=(U[(p,q)] @ U[(s,r)].T).tocsr()
        Bm=(U[(q,r)] @ U[(p,s)].T).tocsr()
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
V=(0.5*(V+V.T)).tocsr(); D=dl**4
print("jm=1 : dl=%d D=%d"%(dl,D),flush=True)
bt=2.5
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
Kdiag=np.ones(D)
for ii in range(D):
    x=ii; pr=1.0
    for _ in range(4): pr*=js[x%dl]; x//=dl
    Kdiag[ii]=pr
Kd=sp.diags(Kdiag)
print("u=t1=%.12f  -ln u=%.12f"%(td[1],-np.log(td[1])),flush=True)
print()
print("=== are the top-2 eigenvectors of T beta_s-independent? ===")
ref=None
for bs in (0.5,1.0,2.0,4.0):
    def Mf(v): return expm_multiply(0.5*bs*V, v)
    e0=np.zeros(D); e0[0]=1.0
    w0=Mf(e0); sc=float(np.dot(w0,Kdiag*w0))
    def T(v): return Mf(Kd*(Mf(v)))/sc
    ev,Vv=eigsh(LinearOperator((D,D),matvec=T,dtype=float),k=2,which="LA",tol=1e-13)
    o=np.argsort(ev)[::-1]; ev=ev[o]; Vv=Vv[:,o]
    v1=Vv[:,0]/np.linalg.norm(Vv[:,0]); v2=Vv[:,1]/np.linalg.norm(Vv[:,1])
    if ref is None: ref=(v1,v2)
    ov1=abs(ref[0]@v1); ov2=abs(ref[1]@v2)
    print("beta_s=%.2f : t1/t2 = %.8f / %.8f  r=%.12f  |<v1_ref|v1>|=%.10f  |<v2_ref|v2>|=%.10f"%(
        bs,ev[0],ev[1],ev[1]/ev[0],ov1,ov2),flush=True)
print()
print("=== structure of v1 at beta_s=1: which link-basis states have weight? ===")
bs=1.0
def Mf(v): return expm_multiply(0.5*bs*V, v)
e0=np.zeros(D); e0[0]=1.0
w0=Mf(e0); sc=float(np.dot(w0,Kdiag*w0))
def T(v): return Mf(Kd*(Mf(v)))/sc
ev,Vv=eigsh(LinearOperator((D,D),matvec=T,dtype=float),k=2,which="LA",tol=1e-13)
o=np.argsort(ev)[::-1]; ev=ev[o]; Vv=Vv[:,o]
v1=Vv[:,0]
nz=[(i,float(v1[i])) for i in range(D) if abs(v1[i])>1e-6]
print("nonzero components of v1: %d (of %d)"%(len(nz),D),flush=True)
for i,val in nz[:12]:
    x=i; labs=[]
    for _ in range(4):
        labs.append(B[x%dl]); x//=dl
    print("   %s : %.8f"%(labs,val),flush=True)
print("[%.1fs]"%(time.time()-t0),flush=True)
