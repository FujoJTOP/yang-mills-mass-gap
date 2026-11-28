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
def Mf(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    return M
jtot=np.array([B[b1][0]+B[b2][0] for b1 in B for b2 in B])
for bs in (0.5,1.0,2.0):
    M=Mf(bs)
    print("beta_s = %.1f :  max |(M_f)_{i,i'}|  for  |j_tot - j_tot'| = n" % bs)
    for n in range(0,7):
        m=0.0
        for i in range(D):
            for k in range(D):
                if abs(jtot[i]-jtot[k])==n:
                    m=max(m,abs(M[i,k]))
        print("     n=%d : %.10e"%(n,m))
    print()
print("=> factorial decay  ~ beta_s^n/n!  would give ratios  beta_s/(n+1):")
bs=1.0
for n in range(0,6):
    print("     beta_s=1 :  n!=%d   beta_s^n/n! = %.6e"%(n,bs**n/np.math.factorial(n)))
