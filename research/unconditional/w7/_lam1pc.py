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
bt=2.5
t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
Kd=np.diag(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B]))
def t1v(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    T=M@Kd@M; T=(T+T.T)/2
    return np.linalg.eigvalsh(T)[-1]
t0=t1v(0.0)
print("check:  ln(t1/t1(0))  vs the two predictions  1.19268*bs^2/2  and  8.256*bs^2/2")
print("%10s %20s %20s %20s"%("beta_s","ln(t1/t1(0))/bs^2","pred A (1.19/2)","pred B (8.26/2)"))
for bs in (1e-3,3e-3,1e-2,3e-2,1e-1):
    L=np.log(t1v(bs)/t0)/bs**2
    print("%10.0e %20.8f %20.8f %20.8f"%(bs,L,1.19268/2/1.0,8.256/2), flush=True)
print()
print("(if the column tends to 0.596, the fit 1.19268 is right; if to 4.13, the formula 8.256 is right)")
