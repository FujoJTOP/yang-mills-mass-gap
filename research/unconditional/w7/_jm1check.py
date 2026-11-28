import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
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
        A=(U[(p,q)] @ U[(s,r)].T).tocsr(); Bm=(U[(q,r)] @ U[(p,s)].T).tocsr()
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
bt=2.5
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
Kdiag=np.ones(D)
for ii in range(D):
    x=ii; pr=1.0
    for _ in range(4): pr*=js[x%dl]; x//=dl
    Kdiag[ii]=pr
Kd=sp.diags(Kdiag)
e0=np.zeros(D); e0[0]=1.0
bs=1.0
def Mf(v): return expm_multiply(0.5*bs*V, v)
print("--- direct check: is |0> an eigenvector of T? ---",flush=True)
w0=Mf(e0); sc=float(np.dot(w0,Kdiag*w0))
print("V|0> nnz = %d  (nonzero components)"%np.count_nonzero(np.abs(Mf(e0))>1e-14),flush=True)
Tv=Mf(Kd*(Mf(e0)))/sc
print("T|0> = %s"%("|0>" if np.allclose(Tv,e0,atol=1e-10) else "NOT |0>"),flush=True)
print("||T|0> - |0>|| = %.3e ; ||T|0>|| = %.6f"%(np.linalg.norm(Tv-e0),np.linalg.norm(Tv)),flush=True)
nz=np.flatnonzero(np.abs(Tv)>1e-10)[:8]
print("first nonzero comps of T|0> :",[(int(i),round(float(Tv[i]),8)) for i in nz],flush=True)
print()
print("--- is the mass-gap sector exactly 2-dim? ---",flush=True)
# apply T twice to |0> and look
v=e0.copy()
for k in range(1,6):
    v=Tv if k==1 else Mf(Kd*(Mf(v)))/sc
    print("  ||T^%d|0> - |0>|| = %.3e   <0|T^%d|0> = %.10f"%(k,np.linalg.norm(v-e0),k,float(v[0])),flush=True)
