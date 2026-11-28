import numpy as np, functools, time
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
def build(bs):
    Mf=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        for op in ops[J2]: Mf+=cJ*op
    Mf*=np.exp(-bs)
    T=Mf@np.diag(Kd)@Mf; T=(T+T.T)/2
    ev,V=np.linalg.eigh(T); idx=np.argsort(ev)[::-1]
    return ev,V,idx
# --- Monte Carlo evaluation of  E = <Psi, (1-q) Psi> / <Psi,Psi>  on SU(2)^2
# SU(2) as unit quaternion (w,x,y,z); Haar = uniform on S^3.
def Dmat(jj, q):
    # Wigner D for SU(2): D^j_{mn}(u) = <jm|u|jn>,  u = w + i(x sx + y sy + z sz)
    w,x,y,z = q
    from scipy.linalg import expm
    sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
    H = x*sx + y*sy + z*sz     # traceless hermitian
    # u = w I + i H ;  D = exp(i * theta * n.sigma)  with u = cos th + i sin th n.sigma
    # so D^j(u) = exp(i * (th) * n . J) ; use J = sigma/2 for j=1/2
    # general j: build Jx,Jy,Jz in the |jm> basis
    ms = np.arange(jj, -jj-1, -1)
    n = len(ms)
    Jz = np.diag(ms.astype(float))
    Jp = np.zeros((n,n))
    for i in range(n-1):
        m = ms[i]
        Jp[i,i+1] = np.sqrt(jj*(jj+1) - m*(m-1))
    Jm = Jp.T.conj()
    Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/(2j)
    th = np.arccos(np.clip(w,-1,1))
    s = np.sqrt(max(1-w*w,0.0))
    if s < 1e-12:
        return np.eye(n, dtype=complex)
    nx,ny,nz = x/s, y/s, z/s
    return expm(1j*th*(nx*Jx + ny*Jy + nz*Jz))
def rand_q(rng):
    v=rng.normal(size=4); v/=np.linalg.norm(v); return v
rng=np.random.default_rng(0)
N=3000
qs=[(rand_q(rng),rand_q(rng)) for _ in range(N)]
# precompute nothing; evaluate on the fly
def E_mc(x, bs):
    ev,V,idx=build(bs)
    vec=V[:,idx[0]]
    num=0.0; den=0.0
    for (q1,q2) in qs:
        psi=0j
        for i,(j1,m1,n1) in enumerate(B):
            pass
        # faster: group by (j,m,n)
        vals1={}; vals2={}
        for j in range(0,jm+1):
            Dj=Dmat(j/2,q1)
            vals1[j]=Dj
        for j in range(0,jm+1):
            vals2[j]=Dmat(j/2,q2)
        for i,(j1,m1,n1) in enumerate(B):
            psi += vec[i]*vals1[j1][(j1-m1)//2,(j1-n1)//2]
        for i,(j2,m2,n2) in enumerate(B):
            psi += 0
        # two-link: Psi = sum over basis of the product space
        p=0j
        for i1,(j1,m1,n1) in enumerate(B):
            for i2,(j2,m2,n2) in enumerate(B):
                p += vec[i1*d+i2]*vals1[j1][(j1-m1)//2,(j1-n1)//2]*vals2[j2][(j2-m2)//2,(j2-n2)//2]
        # q = tr(q1 q2^{-1})/2
        w1,x1,y1,z1=q1; w2,x2,y2,z2=q2
        # q2^{-1} = conj
        ww= w1*w2 + x1*x2 + y1*y2 + z1*z2
        qv= -x1*w2 + w1*x2 + y1*z2 - z1*y2
        qw= -y1*w2 - x1*z2 + w1*y2 + z1*x2
        qr= -z1*w2 + x1*y2 - y1*x2 + w1*z2
        qq=ww
        qq= np.clip(ww,-1,1)
        num += abs(p)**2*(1-qq); den += abs(p)**2
    return num/den
t0=time.time()
bs=2.0
ev,V,idx=build(bs)
vecs=[V[:,idx[0]],V[:,idx[1]]]
print("MC estimate of E_1, E_2 at beta_s=2.0   (N=%d samples)"%N)
for k in range(2):
    vec=vecs[k]
    num=0.0; den=0.0
    for (q1,q2) in qs:
        vals1=[Dmat(j/2,q1) for j in range(0,jm+1)]
        vals2=[Dmat(j/2,q2) for j in range(0,jm+1)]
        p=0j
        for i1,(j1,m1,n1) in enumerate(B):
            a=vals1[j1][(j1-m1)//2,(j1-n1)//2]
            for i2,(j2,m2,n2) in enumerate(B):
                p += vec[i1*d+i2]*a*vals2[j2][(j2-m2)//2,(j2-n2)//2]
        ww=np.clip(q1[0]*q2[0]+q1[1]*q2[1]+q1[2]*q2[2]+q1[3]*q2[3],-1,1)
        num += abs(p)**2*(1-ww); den += abs(p)**2
    print("   E_%d = %.8f    [%.1fs]"%(k+1,num/den,time.time()-t0), flush=True)
