import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=int(__import__("sys").argv[1])
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
    # q on links a,b ; identity on the third
    tot=None
    for M2 in range(-1,2,2):
        for N2 in range(-1,2,2):
            A=sp.csr_matrix(mult_op(1,M2,N2)); Ac=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
            mats=[]
            for k in range(3):
                if k==a: mats.append(A)
                elif k==b: mats.append(Ac)
                else: mats.append(Ip)
            term=sp.kron(sp.kron(mats[0],mats[1],format='csr'),mats[2],format='csr')
            tot=term if tot is None else tot+term
    return ((tot+tot.T)*0.5).tocsr()
Q12=pair_op(0,1); Q23=pair_op(1,2); Q31=pair_op(0,2)
Sq=(Q12+Q23+Q31).tocsr()
print("built, nnz(Sq) = %d , D = %d"%(Sq.nnz,D))
jit=np.array([B[i][0]+B[j][0]+B[k][0] for i in range(d) for j in range(d) for k in range(d)])
odd=[abs(Sq[a,b]) for a,b in zip(*Sq.nonzero()) if (jit[a]-jit[b])%2!=0]
print("j_tot parity check: odd-|dj| entries = %d  (0 confirms the superselection)"%len(odd))
bt=2.5
td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
Kd=np.array([td[B[i][0]]*td[B[j][0]]*td[B[k][0]] for i in range(d) for j in range(d) for k in range(d)])
def Mf_mult(v,bs):
    w=expm_multiply((bs/2)*Sq, v)
    return np.exp(-1.5*bs)*w
def T_mult(v,bs):
    return Mf_mult(Kd*Mf_mult(v,bs),bs)
def r_of(bs):
    op=LinearOperator((D,D),matvec=lambda v: T_mult(v,bs),dtype=float)
    ev=eigsh(op,k=2,which='LA',return_eigenvectors=False,tol=1e-10)
    ev=np.sort(ev)[::-1]
    return ev[1]/ev[0], ev
print()
print("triangle L_s=3 (j_max=1 doubled i.e. j<=1/2) : r(beta_s) at beta_t=2.5")
print("%10s %18s"%("beta_s","r"))
for bs in (2.0,8.0,16.0):
    try:
        r,ev=r_of(bs)
        print("%10.2f %18.12f"%(bs,r), flush=True)
    except Exception as ex:
        print("%10.2f  FAILED: %s"%(bs,ex), flush=True)




