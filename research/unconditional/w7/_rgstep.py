import numpy as np, functools
from scipy.special import iv
from scipy.optimize import brentq
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
def Mf_single(bs):
    A=np.zeros((d,d))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        O=np.zeros((d,d))
        for M2 in range(-J2,J2+1,2):
            for N2 in range(-J2,J2+1,2):
                O+=mult_op(B,index,J2,M2,N2,False)
        A+=cJ*O
    return A
def r_model(bp,bt):
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    M=Mf_single(bp)
    T=M@np.diag(t1d)@M; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1]/ev[0]
bt=2.5
print("test: is Tr_2(T) of the L=1 two-link model a single-link T' = M(b')K_1 M(b') ?")
print("%10s %16s %16s %16s"%("beta_s","r (2-link)","r (blocked)","fitted beta'"))
for bs in (0.5,1.0,2.0,4.0):
    M2=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M2+=cJ*op
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    Kd=np.diag(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B]))
    T=M2@Kd@M2; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]; r2=ev[1]/ev[0]
    T4=T.reshape(d,d,d,d); Tp=np.einsum("ijik->jk",T4); Tp=(Tp+Tp.T)/2
    ev2=np.sort(np.linalg.eigvalsh(Tp))[::-1]; rb=ev2[1]/ev2[0]
    # fit beta' by matching r'
    try:
        bp=brentq(lambda x: r_model(x,bt)-rb,1e-6,20.0)
    except Exception:
        bp=float("nan")
    print("%10.1f %16.10f %16.10f %16.8f"%(bs,r2,rb,bp), flush=True)
