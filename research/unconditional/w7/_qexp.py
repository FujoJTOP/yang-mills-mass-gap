import numpy as np, functools
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
Q12=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Q12+=np.kron(mult_op(B,index,1,M2,N2,False),mult_op(B,index,1,M2,N2,True))
Q12=(Q12+Q12.T)/2
def state(i1,i2):
    v=np.zeros(D); v[i1*d+i2]=1.0; return v
i0=0                     # (0,0,0)
i1=B.index((1,1,1))      # a j=1 (doubled) state
e1=state(i0,i0); e2=state(i1,i0)
print("check the key number <q_{12}> for e_1 = (j=0)x(j=0) and e_2 = (j=1)x(j=0):")
print("   <e1,Q12 e1> = %+.10e"%float(e1@(Q12@e1)))
print("   <e2,Q12 e2> = %+.10e"%float(e2@(Q12@e2)))
print()
print("   <e1,Q12^2 e1> = %+.10e"%float(e1@(Q12@(Q12@e1))))
print("   <e2,Q12^2 e2> = %+.10e"%float(e2@(Q12@(Q12@e2))))
print()
print("=> d ln t_i / d beta_s |_0 = 2 <e_i, sum_P q_P  e_i>   (with the e^{+beta q} normalisation)")
print("   For L_s=1: sum_P = Q12 only.")
print("   dln t_1 = %.6f (should be 0)"%(2*float(e1@(Q12@e1))))
print("   dln t_2 = %.6f (round 164 measured 0.5)"%(2*float(e2@(Q12@e2))))
