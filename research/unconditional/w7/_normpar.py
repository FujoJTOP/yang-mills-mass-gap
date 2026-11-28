import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=int(sys.argv[1]); bs=16.0; bt=2.5
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
hvec=np.array([B[i][0]+B[j][0]+B[k][0] for i in range(d) for j in range(d) for k in range(d)])
keep=np.where(hvec%2==0)[0]
Sq=Sq[np.ix_(keep,keep)].tocsr()
td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
Kd=np.array([td[B[i][0]]*td[B[j][0]]*td[B[k][0]] for i in range(d) for j in range(d) for k in range(d)])[keep]
scale=None
def Mf_mult(v): return np.exp(-1.5*bs)*expm_multiply((bs/2)*Sq, v)
D=len(keep); e0=np.zeros(D); e0[0]=1.0
w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
print("jm=%d D=%d scale=%.6e"%(jm,D,scale), flush=True)
def T_mult(v): return Mf_mult(Kd*Mf_mult(v))/scale
op=LinearOperator((D,D),matvec=T_mult,dtype=float)
t0=time.time()
ev,V=eigsh(op,k=4,which='LA',tol=1e-13,maxiter=20000)
dt=time.time()-t0
o=np.argsort(ev)[::-1]; ev=ev[o]; V=V[:,o]
print("eig %.0fs"%dt, flush=True)
print("   t_1=%.10f  t_2=%.10f  r=%.12f  eig3=%.6e"%(ev[0],ev[1],ev[1]/ev[0],ev[2]))
for i in range(2):
    v=V[:,i]; rv=T_mult(v); res=np.linalg.norm(rv-ev[i]*v)/max(abs(ev[i]),1e-300)
    print("   eig %d residual = %.3e"%(i,res))

