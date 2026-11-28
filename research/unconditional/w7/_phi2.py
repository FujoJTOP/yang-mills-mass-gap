import numpy as np, functools, time
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
def Mf_of(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    return M*np.exp(-bs)
def neglogr(bt,bs):
    t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
    Kd=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
    Mf=Mf_of(bs)
    T=Mf@np.diag(Kd)@Mf; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return -np.log(ev[1]/ev[0])
print("two-link:  f(beta_t) = -ln r  at fixed beta_s,  and its 1/beta_t expansion")
for bs in (0.0,2.0):
    print()
    print("  beta_s = %.1f"%bs)
    bts=np.array([8.0,16.0,32.0,64.0,128.0,256.0])
    fs=np.array([neglogr(b,bs) for b in bts])
    print("     beta_t : "+" ".join("%14.4f"%b for b in bts))
    print("     f      : "+" ".join("%14.8f"%v for v in fs))
    # fit f = sum_{k=1}^{4} ck / beta^k
    A=np.vstack([1.0/bts**k for k in range(1,5)]).T
    coef,_,_,_=np.linalg.lstsq(A,fs,rcond=None)
    print("     c1..c4 : "+" ".join("%14.8f"%c for c in coef))
    print("     beta_t * c1 = %.8f   (free-field would be 3/2 = 1.5, lambda-ratio value)"%(coef[0]))
