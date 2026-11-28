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
O1=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        O1+=np.kron(mult_op(B,index,1,M2,N2,False),mult_op(B,index,1,M2,N2,True))
O1=(O1+O1.T)/2
O2=np.zeros((D,D))
for M2 in range(-2,3,2):
    for N2 in range(-2,3,2):
        O2+=np.kron(mult_op(B,index,2,M2,N2,False),mult_op(B,index,2,M2,N2,True))
O2=(O2+O2.T)/2
import sys; bt=float(sys.argv[1])
td=2*iv(np.arange(0,jm+1)+1,bt)/bt
tvec=np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B])
Kd=np.diag(tvec)
M1=0.5*O1; M2=0.25*O2
T1=M1@Kd+Kd@M1; T1=(T1+T1.T)/2
T2=M2@Kd+2*(M1@Kd@M1)+Kd@M2; T2=(T2+T2.T)/2
ev=np.sort(tvec)[::-1]; t1=ev[0]; t2=ev[1]
deg=[i for i in range(D) if abs(tvec[i]-t2)<1e-10*t1]
W=np.eye(D)[:,deg]
A=W.T@T1@W
wa,V=eigh((A+A.T)/2)
k=int(np.argmax(wa))
x=W@V[:,k]
print("t_1=%.10f  t_2=%.10f   lambda_2 = %.10f"%(t1,t2,wa[k]/t2))
print()
term_a=2*float(x@(M2@(Kd@x)))/t2
term_b=2*float(x@(M1@(Kd@(M1@x))))/t2
print("  <x,T_2 x>/t_2  = 2<x,M_2 K x>/t_2 + 2<x,M_1 K M_1 x>/t_2")
print("     2<x,M_2 K x>/t_2   = %+.10f"%term_a)
print("     2<x,M_1 K M_1 x>/t_2 = %+.10f"%term_b)
S=0.0
for j in range(D):
    v=np.dot(np.eye(D)[j],T1@x)
    if abs(v)<1e-13: continue
    if abs(tvec[j]-t2) < 1e-9*t1: continue
    S+=2*v*v/(t2*(t2-tvec[j]))
print("     sum_{j outside D} 2<e_j,T_1 x>^2/(t_2(t_2-t_j)) = %+.10f"%S)
print()
tot=term_a+term_b+S
print("  => (t_2''/t_2) = %.10f"%tot)
print("  => lambda_2' = that - lambda_2^2 = %.10f   (measured ~0.597)"%(tot-(wa[k]/t2)**2))


