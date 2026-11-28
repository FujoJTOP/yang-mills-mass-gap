import numpy as np, functools
from scipy.special import iv
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
def Mf_of(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        for op in ops[J2]: M+=cJ*op
    return M*np.exp(-bs)
print("check:  <e_a, (I - M_f) e_a> / beta_s  ->  <e_a, Q e_a>   as beta_s -> 0")
print("       analytic target for a=0 (the j=0 x j=0 state):  1")
print("%12s %18s %18s"%("beta_s","(I-Mf)_00/beta_s","offdiag check"))
for bs in (1e-3,1e-4,1e-5):
    Mf=Mf_of(bs)
    v=(1.0-Mf[0,0])/bs
    print("%12.0e %18.10f %18.3e"%(bs,v,np.abs(Mf[0,1])/bs))
print()
print("also: c_1/c_0 at small beta_s, compare 2(J+1)I_{J+1}/beta  with the e^{beta cos} expansion")
for bs in (1e-2,1e-3):
    print("   beta_s=%.0e : c_1 = %.6e   (2*2*I_2/bs = %.6e)"%(bs,2*2*iv(2,bs)/bs,2*2*iv(2,bs)/bs))
