import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (2,3,4):
    B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
    index={b:i for i,b in enumerate(B)}; dl=len(B)
    def mult_op(J2,M2,N2,conj=False):
        A=np.zeros((dl,dl))
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
    U1=np.zeros((dl,dl))
    for a in (-1,1): U1+=mult_op(1,a,a)
    W=0.5*U1
    # j-only subspace: for each J, the uniform vector over the (2J+1) states (m=n? or all?)
    # from the closed form, the max eigenvector is constant over the (2J+1) m=n states
    nJ=jm+1
    evecs=[]
    for J in range(0,jm+1):
        v=np.zeros(dl); cnt=0
        for i,(j,m,n) in enumerate(B):
            if j==J and m==n: v[i]=1.0; cnt+=1
        evecs.append(v/np.sqrt(cnt))
    # reduced W on the j-only subspace
    Wr=np.zeros((nJ,nJ))
    for a in range(nJ):
        for b in range(nJ):
            Wr[a,b]=float(evecs[a] @ (W @ evecs[b]))
    ev=np.linalg.eigvalsh(0.5*(Wr+Wr.T))
    print("jm=%d : reduced W (j-only) eigenvalues:"%jm, flush=True)
    print("   ",np.round(np.sort(ev)[::-1],8), flush=True)
    print("   cos(k pi/(jm+2)):", np.round(np.cos(np.arange(1,jm+2)*np.pi/(jm+2)),8), flush=True)
    print("   exact W_max = %.10f   (cos(pi/(jm+2))=%.10f)"%(
        np.linalg.eigvalsh(0.5*(W+W.T)).max(), np.cos(np.pi/(jm+2))), flush=True)
    # full spectra for comparison
    evf=np.linalg.eigvalsh(0.5*(W+W.T))
    print("   FULL W spectrum (distinct):", np.round(np.unique(np.round(evf,8)),8), flush=True)
    print()
