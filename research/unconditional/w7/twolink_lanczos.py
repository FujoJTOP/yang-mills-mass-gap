import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import functools, time, sys

# O2: two-link r(J_out) via matrix-free Lanczos.
#   T = M_f K M_f ,  K = K1 (x) K1  diagonal,  M_f = sum_{J,M,N} c_J  A_{JMN} (x) Ac_{JMN}
#   (A (x) B) v  =  (A @ V @ B.T).ravel()   with V = v.reshape(d,d)
# We only need the top two eigenvalues.

@functools.lru_cache(maxsize=None)
def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1,2), Rational(j2,2), Rational(j3,2),
                             Rational(m1,2), Rational(m2,2), Rational(m3,2)), 25))

def make_basis(jmax2):
    return [(j,m,n) for j in range(0,jmax2+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]

def mult_op(B, index, J2, M2, N2, conj=False):
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

def build_terms(jmax2, betap, Jchar=None):
    B=make_basis(jmax2); index={b:i for i,b in enumerate(B)}; d=len(B)
    if Jchar is None: Jchar=jmax2
    terms=[]   # list of (c, A, Ac)
    for J2 in range(0,Jchar+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,betap)/betap)
        for M2 in range(-J2,J2+1,2):
            for N2 in range(-J2,J2+1,2):
                terms.append((cJ, mult_op(B,index,J2,M2,N2,False), mult_op(B,index,J2,M2,N2,True)))
    return B, terms, d

def lanczos_top2(matvec, n, k=60, seed=0):
    rng=np.random.default_rng(seed)
    q=rng.normal(size=n); q/=np.linalg.norm(q)
    Q=np.zeros((k+1,n)); Q[0]=q
    al=np.zeros(k); be=np.zeros(k)
    for j in range(k):
        w=matvec(Q[j])
        al[j]=Q[j]@w
        w=w-al[j]*Q[j]
        if j>0: w=w-be[j-1]*Q[j-1]
        # reorthogonalize once
        w=w-Q[:j+1].T@(Q[:j+1]@w)
        be[j]=np.linalg.norm(w)
        if be[j]<1e-14: k=j+1; break
        Q[j+1]=w/be[j]
    Tk=np.diag(al[:k])+np.diag(be[:k-1],1)+np.diag(be[:k-1],-1)
    ev=np.sort(np.linalg.eigvalsh(Tk))[::-1]
    return ev

def run(jmax2, beta, betap, Kchar=None, k=50):
    t0=time.time()
    B, terms, d = build_terms(jmax2, betap, Kchar)
    n=d*d
    tk=np.array([2*iv(j+1,beta)/beta for j in range(jmax2+1)])
    Kd=np.array([tk[b1[0]]*tk[b2[0]] for b1 in B for b2 in B])
    coef=np.exp(-betap)
    nterms=len(terms)
    def Mf(v):
        V=v.reshape(d,d); out=np.zeros((d,d))
        for (c,A,Ac) in terms:
            out += c*(A@V@Ac.T)
        return (coef*out).ravel()
    def T(v):
        return Mf(Kd*Mf(v))
    ev=lanczos_top2(T, n, k=k)
    print("J_out=%.1f  d=%d  n=%d  terms=%d  t1=%.8e  t2=%.8e  r=%.10f  [%.1fs]"
          % (jmax2/2, d, n, nterms, ev[0], ev[1], ev[1]/ev[0], time.time()-t0), flush=True)
    return ev[1]/ev[0]

if __name__=="__main__":
    beta=2.5
    print("two-link r(J_out), beta=beta'=2.5, I2/I1=%.10f"%(iv(2,beta)/iv(1,beta)))
    for jmax2 in (2,3,4):
        run(jmax2, beta, beta, Kchar=min(6,jmax2))
