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
def r_of(bt,bs):
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    Kd=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    M*=np.exp(-bs)
    T=M@np.diag(Kd)@M; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1]/ev[0]
print("leading effect of V:  dr/d beta_s at beta_s = 0")
for bt in (1.0,2.5,8.0):
    r0=r_of(bt,0.0)
    h=1e-3
    rp=r_of(bt,h); rm=r_of(bt,-h) if False else r0
    slope=(rp-r0)/h
    print("   beta_t=%.1f : I2/I1 = %.10f   slope = %.8f"%(bt,r0,slope))
print()
print("and the V=1 case is exact for ANY number of links:")
print("   T = K_1^{tensor n}  =>  eigenvalues are products t_{j1}t_{j2}...t_{jn}")
print("   top two:  t_0^n  and  t_0^{n-1} t_1   =>  r = t_1/t_0 = I2/I1  ,  independent of n")
