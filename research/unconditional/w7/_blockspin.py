import numpy as np, functools
from scipy.special import iv
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
def twolink(bt,bs):
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    Kd=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    M*=np.exp(-bs)
    T=M@np.diag(Kd)@M
    return (T+T.T)/2
def singlelink(bt):
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    return t1d[1]/t1d[0]
print("block-spin test:  T' = partial trace over the 2nd link;  compare r' = t2'/t1' with r (two-link) and I2/I1")
print("%8s %8s %14s %14s %14s"%("beta_t","beta_s","r (2-link)","r (blocked)","I2/I1"))
for bt in (1.0,2.5,8.0):
    for bs in (0.0,1.0,4.0):
        T=twolink(bt,bs)
        ev=np.sort(np.linalg.eigvalsh(T))[::-1]
        r2=ev[1]/ev[0]
        # partial trace over the second link: reshape to d x d x d x d and trace the 2nd pair
        Ta=T.reshape(d,d,d,d)          # (a1,b1,a2,b2) with row=(a1,a2), col=(b1,b2)
        Tb=np.einsum("ijik->jk",Ta.reshape(d,d,d,d).transpose(0,2,1,3).reshape(d*d,d*d).reshape(d,d,d,d)) if False else None
        T4=T.reshape(d,d,d,d)          # row index = i1*d+i2 -> (i1,i2); col = j1*d+j2 -> (j1,j2)
        Tp=np.einsum("ijik->jk",T4)
        Tp=(Tp+Tp.T)/2
        ev2=np.sort(np.linalg.eigvalsh(Tp))[::-1]
        rb=ev2[1]/ev2[0]
        print("%8.1f %8.1f %14.8f %14.8f %14.8f"%(bt,bs,r2,rb,singlelink(bt)), flush=True)
