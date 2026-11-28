import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (4,5):
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
    for a in (-1,1):
        U1 += mult_op(1,a,a)
    U1=0.5*U1
    # build the (m,m)-diagonal subspace: for each j, the vector with m=n and all m in that sector
    # Actually the eigenvector was: component at state (j,m,m) with m=-j,-j+2,...,j  all EQUAL.
    # So define the normalized basis vector e_j = (1/sqrt(2j+1)) sum_m |j,m,m>
    evecs={}
    for jj in range(0,jm+1):
        v=np.zeros(dl)
        cnt=0
        for i,(j,m,n) in enumerate(B):
            if j==jj and m==n:
                v[i]=1.0; cnt+=1
        v=v/np.sqrt(cnt)
        evecs[jj]=v
    nJ=jm+1
    M=np.zeros((nJ,nJ))
    for jj in range(nJ):
        for kk in range(nJ):
            M[jj,kk]=float(evecs[jj] @ U1 @ evecs[kk])
    print("jm=%d : M_{jj'} in the (m,m)-diagonal basis:"%jm, flush=True)
    with np.printoptions(precision=8,suppress=True,linewidth=200):
        print(M, flush=True)
    ev,V=np.linalg.eigh(M)
    print("   M eigenvalues:", np.round(ev,8), flush=True)
    print("   cos(k pi/(jm+2)):", np.round(np.cos(np.arange(1,jm+2)*np.pi/(jm+2)),8), flush=True)
