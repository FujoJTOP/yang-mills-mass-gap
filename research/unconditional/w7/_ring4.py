import numpy as np, functools, scipy.sparse as sp, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=3; d=sum(2*j+1 for j in range(jm+1))**1
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B)
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((d,d))
    mm,nn=(-M2,-N2) if conj else (M2,N2)
    ph=((-1.0)**((M2-N2)//2)) if conj else 1.0
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2-mm,n2-nn
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,mm,j2,-m2)*d3(j3,n3,J2,nn,j2,-n2)
            if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
    return ph*A
Ip=sp.identity(d,format='csr')
def pair_op_4(L,i,j):
    tot=None
    for M2 in range(-1,2,2):
        for N2 in range(-1,2,2):
            A=sp.csr_matrix(mult_op(1,M2,N2)); Ac=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
            mats=[A if k==i else (Ac if k==j else Ip) for k in range(L)]
            term=mats[0]
            for k in range(1,L): term=sp.kron(term,mats[k],format='csr')
            tot=term if tot is None else tot+term
    return ((tot+tot.T)*0.5).tocsr()
L=4
t0=time.time()
Sq=sp.csr_matrix((d**L,d**L))
for (i,j) in [(0,1),(1,2),(2,3),(3,0)]:
    Sq=Sq+pair_op_4(L,i,j)
Sq=Sq.tocsr()
print("L=4 ring : d=%d  D=%d  nnz(Sq)=%d  (build %.0fs)"%(d,d**L,Sq.nnz,time.time()-t0), flush=True)
td=ive(np.arange(0,jm+1)+1,2.5)/ive(1,2.5)
# K = product of t_{j_x}
import itertools
Kd=np.ones(d**L)
for idx in range(d**L):
    x=idx; prod=1.0
    for _ in range(L):
        prod*=td[B[x%d][0]]; x//=d
    Kd[idx]=prod
bs=1.0
def Mf_mult(v): return np.exp(-L*bs/2)*expm_multiply((bs/2)*Sq, v)
e0=np.zeros(d**L); e0[0]=1.0
w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
print("scale = %.6e"%scale, flush=True)
def T_mult(v): return Mf_mult(Kd*Mf_mult(v))/scale
op=LinearOperator((d**L,d**L),matvec=T_mult,dtype=float)
t1=time.time()
ev,V=eigsh(op,k=3,which="LA",tol=1e-10,maxiter=3000)
o=np.argsort(ev)[::-1]; ev=ev[o]; V=V[:,o]
print("L=4 ring, beta_s=%.2f : t_1=%.8f t_2=%.8f r=%.10f  (%.0fs)"%(bs,ev[0],ev[1],ev[1]/ev[0],time.time()-t1), flush=True)
for i in range(2):
    res=np.linalg.norm(T_mult(V[:,i])-ev[i]*V[:,i])/abs(ev[i])
    print("   eig %d residual = %.2e"%(i,res))
