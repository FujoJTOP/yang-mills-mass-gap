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
hvec=np.array([b1[0]+b2[0] for b1 in B for b2 in B])
u=ive(2,bt)/ive(1,bt); ut1=u*tvec[0]
print("t_1=%.10f  u t_1=%.10f"%(tvec[0],ut1))
print()
for bs in (1.0,4.0):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    Kd=np.diag(tvec)
    T=M@Kd@M; T=(T+T.T)/2
    deg1=[i for i in range(D) if abs(hvec[i]-1)<1e-9]
    W=np.eye(D)[:,deg1]
    A=W.T@T@W; A=(A+A.T)/2
    wa,V=np.linalg.eigh(A)
    psi=W@V[:,-1]
    w=M@psi
    norm2=float(np.dot(w,w))
    val=float(np.dot(w,tvec*w))
    print("beta_s=%.2f : trial val=%.10f   ||M psi||^2=%.10f   val/(u t_1)=%.8f"%(bs,val,norm2,val/ut1))
    print("   h-decomposition of ||M psi||^2:")
    for h in sorted(set(hvec.tolist())):
        m=hvec==h
        wn=float(np.sum(w[m]**2)); vv=float(np.sum(tvec[m]*w[m]**2))
        if wn>1e-14:
            print("      h=%2d : weight=%.8e   t*weight=%.8e   t range=[%.4f,%.4f]"%(h,wn,vv,tvec[m].min(),tvec[m].max()))
    print()
