import numpy as np, functools
from scipy.special import ive
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=3
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((d,d))
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
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(J2,M2,N2),mult_op(J2,M2,N2,conj=True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
bt=2.5
td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
tvec=np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B])
t1v=tvec[0]; u=ive(2,bt)/ive(1,bt)
target=u*t1v
print("t_1 = %.10f   u*t_1 = %.10f   (target for the trial bound)"%(t1v,target))
print()
print("%8s %14s %22s %14s"%("beta_s","ORTHOG trial","(val)/(u t_1)","bound holds?"))
for bs in (16.0,32.0,64.0,128.0,256.0):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    Kd=np.diag(tvec)
    T=M@Kd@M; T=(T+T.T)/2
    deg=[i for i in range(D) if abs(tvec[i]-tvec[1])<1e-12*t1v]
    W=np.eye(D)[:,deg]
    A=W.T@T@W; A=(A+A.T)/2
    wa=np.linalg.eigvalsh(A)
    mx=wa[-1]
    print("%8.2f %14.10f %22.10f %14s"%(bs,mx,mx/target,"YES" if mx>=target-1e-12 else "no"), flush=True)


