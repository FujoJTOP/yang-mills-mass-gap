import numpy as np, functools
from scipy.special import iv, ive
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
def pieces(bs,bt):
    M=np.zeros((D,D)); M1=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        # derivative of c_J at beta_s: c_J = 2(J+1) I_{J+1}(bs)/bs ; use numerical derivative for simplicity
        h=1e-6
        cJp=( (1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs+h)/(bs+h)) - (1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs-h)/(bs-h)) )/(2*h)
        for op in ops[J2]:
            M+=cJ*op; M1+=cJp*op
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    Kd=np.diag(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B]))
    T=M@Kd@M; T=(T+T.T)/2
    T1=(M1@Kd@M+M@Kd@M1); T1=(T1+T1.T)/2
    return T,Kd,T1
print("first-order formula:  dr/dbeta_s|_0 = (1/t2)<e2,T1 e2> - (1/t1)<e1,T1 e1>")
print("at beta_s=0, e_i are basis vectors and T1 = M1 K + K M1, so <e,T1 e> = 2 t_e <e,M1 e>")
for bt in (1.0,2.5,8.0):
    T,Kd,T1=pieces(1e-8,bt)
    ev,V=eigh(T)
    idx=np.argsort(ev)[::-1]
    # top eigenspace (watch degeneracy)
    vals=ev[idx]
    print("   beta_t=%.1f : top 4 eigenvalues of T = %s"%(bt," ".join("%.6e"%v for v in vals[:4])))
    # first-order shift using the exact eigenvectors of T at beta_s=0 (basis vectors)
    x1=V[:,idx[0]]; x2=V[:,idx[1]]
    s=np.dot(x2,T1@x2)/ev[idx[1]]-np.dot(x1,T1@x1)/ev[idx[0]]
    print("        first-order slope = %.10f"%s)
