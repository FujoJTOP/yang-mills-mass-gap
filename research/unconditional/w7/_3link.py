import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def build(jm):
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
    u={}
    for M2 in (-1,1):
        for N2 in (-1,1):
            u[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
            u[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
    return B,index,dl,u
Idxx=[(-1,-1),(-1,1),(1,-1),(1,1)]
def Wilson(L,a,b,c,dd,u,dl):
    tot=None
    for (p,q) in Idxx:
        A=u[(p,q)]
        for (r,s) in Idxx:
            Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
            mats=[A,Bm,C,Dm]
            # place the 4-link operator on links a,b,c,dd of an L-link space
            # build via kron with identities
            I=sp.identity(dl,format='csr')
            term=None
            for k in range(L):
                if k==a: M=A
                elif k==b: M=Bm
                elif k==c: M=C
                elif k==dd: M=Dm
                else: M=I
                term=M if term is None else sp.kron(term,M,format='csr')
            tot=term if tot is None else tot+term
    return (0.5*tot).tocsr()
# 3-link "star": plaquette on (0,1,2,3) is 4 links -> not 3. Use a 3-link chain?
# Instead: 2 plaquettes sharing 2 links? Simplest 3-link: one plaquette (0,1,2,3) is 4 links.
# For 3 links we need a closed chain of 3 links: plaquette (0,1,2,0) - degenerate.
# Let's do the 4-link with periodic identification: links 0,1,2,3 with 0=3 -> 3 links.
# Build a 3-link W with the trace over a 3-step loop: (1/2)tr(u1 u2^-1 u3) is not gauge invariant.
# Skip: go directly to the 4-link + measure m a for the 2-plaquette chain with jm=4.
print("3-link gauge-invariant Wilson loop needs even number of links (SU(2) invariant).")
print("Skipping; next: 2-plaquette chain with small jm and larger beta_s.")
