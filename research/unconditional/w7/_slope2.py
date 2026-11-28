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
h=1e-6
M=np.zeros((D,D)); M1=np.zeros((D,D))
for J2 in range(0,jm+1):
    cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,h)/h)
    cJp=((1.0 if J2==0 else 2*(J2+1)*iv(J2+1,h)/h)-(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,-h)/(-h)))/(2*h)
    for op in ops[J2]:
        M+=cJ*op; M1+=cJp*op
T=M@Kd@M; T=(T+T.T)/2
T1=M1@Kd@M+M@Kd@M1; T1=(T1+T1.T)/2
ev,V=eigh(T); idx=np.argsort(ev)[::-1]
t1v=ev[idx[0]]; t2v=ev[idx[1]]
print("beta_t=2.5, beta_s=0 :  t1 = %.10f  t2 = %.10f  (t2 is 3-fold degenerate)"%(t1v,t2v))
deg=[i for i in idx if abs(ev[i]-t2v)<1e-8*abs(t2v)]
print("   degenerate indices: %d"%len(deg))
W=V[:,deg]
A=W.T@T1@W
wa=np.sort(np.linalg.eigvalsh((A+A.T)/2))[::-1]
print("   T1 restricted to the degenerate subspace: eigenvalues =", " ".join("%.8f"%v for v in wa))
x1=V[:,idx[0]]
s1=np.dot(x1,T1@x1)/t1v
print("   <x1,T1 x1>/t1 = %.10f"%s1)
print()
print("   first-order slope for EACH degenerate direction:  (lambda_i/t2 - s1)")
for v in wa:
    print("      %.10f"%((v/t2v)-s1))
