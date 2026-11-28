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
Q=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Q+=np.kron(mult_op(B,index,1,M2,N2,False),mult_op(B,index,1,M2,N2,True))
Q=(Q+Q.T)/2
bt=2.5
t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
tvec=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
Kd=np.diag(tvec)
M1=0.5*Q
e1=np.zeros(D); e1[0]=1.0
t1=tvec[0]; tj=tvec[np.argmax(np.abs(M1@e1))]
print("t_1 = %.10f   t_j (the reached state) = %.10f   u^2 = %.10f"%(t1,tj,tj/t1))
print()
c1=2*float(e1@(M1@(Kd@(M1@e1))))/t1
print("CORRECT:  2<e1,M1 K M1 e1>/t1 = %.10f   (note: <M1 K M1> = t_j ||M1 e1||^2, NOT t_1 ||.||^2)"%c1)
print("  check: t_j/t_1 * ||M1 e1||^2 * 2 = %.10f"%((tj/t1)*(2*float(e1@(M1@(M1@e1))))) )
print()
# the sum term
v=M1@e1
S=0.0
for j in range(D):
    if abs(v[j])<1e-14: continue
    T1ej=(tvec[j]+t1)*v[j]
    S+=2*T1ej**2/(t1*(t1-tvec[j]))
print("sum term = %.10f"%S)
print("total    = %.10f     (measured lambda_1' = 1.19268)"%(c1+S))
