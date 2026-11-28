import numpy as np, functools, time
from scipy.special import ive
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
def makef(jm,Jmax):
    B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
    ops={}
    for J2 in range(0,min(Jmax,jm)+1):
        ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
                 for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
    def f(bt,bs):
        td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
        Kd=np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B])
        M=np.zeros((D,D))
        for J2 in range(0,min(Jmax,jm)+1):
            if J2==0: cJ=1.0
            elif bs==0.0: continue
            else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
            for op in ops[J2]: M+=cJ*op
        T=M@np.diag(Kd)@M; T=(T+T.T)/2
        ev=np.sort(np.linalg.eigvalsh(T))[::-1]
        return -np.log(ev[1]/ev[0])
    return f,D
print("truncation check: does Phi change with jmax2 ?")
print("%8s %8s %10s %18s %18s"%("jmax2","Jmax","D","Phi(beta_s=16)","Phi(beta_s=32)"))
for jm,Jmax in ((3,3),(4,3),(4,4)):
    t0=time.time()
    f,D=makef(jm,Jmax)
    a=f(1024.0,16.0); b=f(2048.0,16.0); p16=2*b-a
    a2=f(1024.0,32.0); b2=f(2048.0,32.0); p32=2*b2-a2
    print("%8d %8d %10d %18.10f %18.10f   [%.0fs]"%(jm,Jmax,D,p16,p32,time.time()-t0), flush=True)
