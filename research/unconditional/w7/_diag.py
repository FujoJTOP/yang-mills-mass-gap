import numpy as np, functools
from scipy.special import ive
from scipy.linalg import eigh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
def mult_op(B,index,J2,M2,N2,conj=False):
    d=len(B); A=np.zeros((d,d))
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
jm=3
B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
def r_diag(beta):
    # t_j proportional to ive(j+1,beta); normalise by ive(1,beta) so nothing overflows
    td=ive(np.arange(0,jm+1)+1,beta)/ive(1,beta)
    Kd=np.diag(np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B]))
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*ive(J2+1,beta)/beta
        for op in ops[J2]: M+=cJ*op
    T=M@Kd@M; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1]/ev[0]
print("diagonal beta_s = beta_t = beta  (the lattice-like slice)")
print("%8s %18s %18s %18s"%("beta","r","gap 1-r","I2/I1 (single-link)"))
for beta in (1.0,2.0,4.0,8.0,16.0,32.0,64.0,128.0):
    r=r_diag(beta)
    u=ive(2,beta)/ive(1,beta)
    print("%8.1f %18.10f %18.10f %18.10f"%(beta,r,1-r,u), flush=True)
