import numpy as np, functools, time
from scipy.linalg import expm
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
from scipy.special import iv
jm=3
B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
beta_t=2.5
t1d=2*iv(np.arange(0,jm+1)+1,beta_t)/beta_t
Kd=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
def Dmat(jj,q):
    w=q[0]
    ms=np.arange(jj,-jj-1,-1); n=len(ms)
    Jz=np.diag(ms.astype(float))
    Jp=np.zeros((n,n))
    for i in range(n-1):
        m=ms[i]; Jp[i,i+1]=np.sqrt(jj*(jj+1)-m*(m-1))
    Jm=Jp.T.conj(); Jx=(Jp+Jm)/2; Jy=(Jp-Jm)/(2j)
    th=np.arccos(np.clip(w,-1,1)); s=np.sqrt(max(1-w*w,0.0))
    if s<1e-12: return np.eye(n,dtype=complex)
    nx,ny,nz=q[1]/s,q[2]/s,q[3]/s
    return expm(1j*th*(nx*Jx+ny*Jy+nz*Jz))
rng=np.random.default_rng(1)
N=60000
qs=[(rng.normal(size=4),rng.normal(size=4)) for _ in range(N)]
qs=[(a/np.linalg.norm(a),b/np.linalg.norm(b)) for a,b in qs]
qinv=[np.clip(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3],-1,1) for a,b in qs]
vals_cache=None
def build(bs):
    Mf=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        for op in ops[J2]: Mf+=cJ*op
    Mf*=np.exp(-bs)
    T=Mf@np.diag(Kd)@Mf; T=(T+T.T)/2
    ev,V=np.linalg.eigh(T); idx=np.argsort(ev)[::-1]
    return ev,V,idx
# precompute the D-matrices once (they do not depend on beta_s)
t0=time.time()
D1=[[Dmat(j/2,q[0]) for j in range(0,jm+1)] for q in qs]
D2=[[Dmat(j/2,q[1]) for j in range(0,jm+1)] for q in qs]
print("precomputed D matrices in %.1fs"%(time.time()-t0), flush=True)
def E_of(vec):
    num=0.0; den=0.0
    for s in range(N):
        v1=D1[s]; v2=D2[s]
        p=0j
        for i1,(j1,m1,n1) in enumerate(B):
            a=v1[j1][(j1-m1)//2,(j1-n1)//2]
            if a==0: continue
            base=i1*d
            for i2,(j2,m2,n2) in enumerate(B):
                p += vec[base+i2]*a*v2[j2][(j2-m2)//2,(j2-n2)//2]
        num += abs(p)**2*(1-qinv[s]); den += abs(p)**2
    return num/den
print()
print("%8s %14s %14s %14s %14s %14s"%("beta_s","E1","E2","-2r(E2-E1)","r\x27(FD)","r"))
for bs in (0.8,1.2,1.6,2.0,3.0,4.0):
    ev,V,idx=build(bs); h=0.2
    x1=V[:,idx[0]]; x2=V[:,idx[1]]; r=ev[idx[1]]/ev[idx[0]]
    evp,_,idp=build(bs+h); evm,_,idm=build(bs-h)
    rp=((evp[idp[1]]/evp[idp[0]])-(evm[idm[1]]/evm[idm[0]]))/(2*h)
    E1=E_of(x1); E2=E_of(x2)
    print("%8.2f %14.8f %14.8f %14.8f %14.8f %14.8f"%(bs,E1,E2,-2*r*(E2-E1),rp,r), flush=True)
print("total %.1fs"%(time.time()-t0))
