import numpy as np, functools
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
def build(jm):
    B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
    ops={}
    for J2 in range(0,jm+1):
        ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
                 for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
    return B,d,D,ops
def makef(jm):
    B,d,D,ops=build(jm)
    def Mf_of(bs):
        M=np.zeros((D,D))
        # normalise by c_0 = e^{-bs} so nothing overflows
        for J2 in range(0,jm+1):
            if J2==0: cJ=1.0
            elif bs==0.0: continue
            else: cJ=2*(J2+1)*ive(J2+1,bs)/bs      # = (2(J+1)I_{J+1}/bs) * e^{-bs}
            for op in ops[J2]: M+=cJ*op
        return M
    def f(bt,bs):
        td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
        Kd=np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B])
        Mf=Mf_of(bs)
        T=Mf@np.diag(Kd)@Mf; T=(T+T.T)/2
        ev=np.sort(np.linalg.eigvalsh(T))[::-1]
        return -np.log(ev[1]/ev[0])
    return f
print("Phi(beta_s) at large beta_s  (jmax2 = 3, i.e. j<=3/2)")
f3=makef(3)
print("%10s %18s %18s %18s"%("beta_s","f(1024)","f(2048)","Phi=2f(2048)-f(1024)"))
prev=None
for bs in (4.0,8.0,16.0,32.0,64.0,128.0):
    a=f3(1024.0,bs); b=f3(2048.0,bs); phi=2*b-a
    r=("" if prev is None else "   ratio %.6f"%(phi/prev))
    print("%10.1f %18.10f %18.10f %18.10f%s"%(bs,a,b,phi,r), flush=True)
    prev=phi
