import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=int(sys.argv[1]) if len(sys.argv)>1 else 3
bt=float(sys.argv[2]) if len(sys.argv)>2 else 2.5
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B); D=d**3
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
def pair_op(a,b):
    tot=None
    for M2 in range(-1,2,2):
        for N2 in range(-1,2,2):
            A=sp.csr_matrix(mult_op(1,M2,N2)); Ac=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
            mats=[A if k==a else (Ac if k==b else Ip) for k in range(3)]
            term=sp.kron(sp.kron(mats[0],mats[1],format='csr'),mats[2],format='csr')
            tot=term if tot is None else tot+term
    return ((tot+tot.T)*0.5).tocsr()
Sq=(pair_op(0,1)+pair_op(1,2)+pair_op(0,2)).tocsr()
td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
Kd=np.array([td[B[i][0]]*td[B[j][0]]*td[B[k][0]] for i in range(d) for j in range(d) for k in range(d)])
def Mf_mult(v,bs):
    return np.exp(-1.5*bs)*expm_multiply((bs/2)*Sq, v)
def r_of(bs):
    def T_mult(v):
        return Mf_mult(Kd*Mf_mult(v,bs),bs)
    op=LinearOperator((D,D),matvec=T_mult,dtype=float)
    ev=eigsh(op,k=2,which='LA',return_eigenvectors=False,tol=1e-8,maxiter=500)
    ev=np.sort(ev)[::-1]
    return ev[1]/ev[0]
print("jm=%d  d=%d  D=%d  nnz(Sq)=%d  beta_t=%.1f"%(jm,d,D,Sq.nnz,bt))
print("%8s %18s %10s"%("beta_s","r","t(s)"))
for bs in (16.0,32.0):
    t0=time.time(); r=r_of(bs); dt=time.time()-t0
    print("%8.1f %18.12f %10.1f"%(bs,r,dt), flush=True)

