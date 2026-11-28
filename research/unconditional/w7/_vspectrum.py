import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=2
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; dl=len(B)
print("jm=%d, dl=%d"%(jm,dl))
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
Idxx=[(-1,-1),(-1,1),(1,-1),(1,1)]
I=sp.identity(dl,format='csr')
def W4_placed(links):
    # links = tuple of 4 link indices (in the 4-link star)
    tot=None
    for (p,q) in Idxx:
        mats=[u[(p,q)], None, None, None]
        for (r,s) in []: pass
    return None
# build the 4-link star W on links (0,1,2,3)
tot=None
for (p,q) in Idxx:
    A=u[(p,q)]
    for (r,s) in Idxx:
        Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
        term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
        tot=term if tot is None else tot+term
W4=(0.5*tot).tocsr()
# 2-plaquette star: 7 links (0,1,2,3) and (3,4,5,6) sharing link 3
# build W_p1 on links 0,1,2,3 ; W_p2 on links 3,4,5,6
def place(W4op, positions, L):
    # W4op acts on 4 links; embed into L links at 'positions'
    # build by iterating over the 4^4 tensor structure
    # simpler: use the fact that W4op is a sparse matrix over d^4; embed via kron with identity
    pass
print("Direct 7-link build is expensive; instead use the shared-link structure.")
print("For jm=2, the 7-link space is d^7 = 14^7 = 1.05e8 -> infeasible.")
print("For jm=1, d=5, d^7 = 78125 -> feasible.")
jm=1
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
# 4-link plaquette operator on 4 links
tot=None
for (p,q) in Idxx:
    A=u[(p,q)]
    for (r,s) in Idxx:
        Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
        term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
        tot=term if tot is None else tot+term
W4=(0.5*tot).tocsr()
# embed W4 into 7 links: links p1=(0,1,2,3), p2=(3,4,5,6)
L=7; D=dl**L
Iop=sp.identity(dl,format='csr')
def embed(W4op, positions):
    # W4op acts on 4 links. Build kron in the order 0..L-1.
    # Decompose W4op as sum over basis of 4-link operators
    # Simpler: reshape W4op to (dl^4, dl^4) and use kron chain
    terms=[]
    # We'll build by iterating over the 4-link indices via the sparse structure
    W4d=W4op.tocoo()
    from collections import defaultdict
    groups=defaultdict(list)
    for i,j,v in zip(W4d.row,W4d.col,W4d.data):
        groups[(i,j)].append(v)
    # build sparse via kron with identity
    Mc=W4op.tocsr()
    # build the L-link operator as: kron of factors, with W4op on 4 positions
    # We need a systematic construction: 
    # represent W4op = sum_k A_k^{(1)} (x) A_k^{(2)} (x) A_k^{(3)} (x) A_k^{(4)}
    # via SVD or just direct: for embedding we can use tensor reshaping.
    # Use the reshape approach: the full operator on L links = 
    #   full = kron(Id_before, W4op_as_4link, Id_after) but reordered.
    pass
    return None
print("Embedding W4 into 7 links needs a tensor reshuffle; doing it directly.")
# Direct construction: iterate over the full 7-link basis is too slow; but W4op is 5^4=625 dim, full=78125.
# The operator V = W_p1 + W_p2 on 78125 dim, sparse.
# Build by explicit kron: for links in order 0..6, W4op on (0,1,2,3) = kron(W4op, I, I, I)
# BUT kron(A,B) treats A as acting on the FIRST factor. So:
# full = kron(kron(W4op, I), I) acts on links (0,1,2,3) = first 4 factors -> yes
n1=dl**4; n2=dl**3
# W_p1 on links 0..3: kron(W4op, I_3) with I_3 = identity on links 4,5,6
Wp1=sp.kron(W4, sp.identity(dl**3,format='csr'), format='csr')
# W_p2 on links 3..6: this is kron(I_3, W4op) acting on the LAST 4 factors
Wp2=sp.kron(sp.identity(dl**3,format='csr'), W4, format='csr')
V=(Wp1+Wp2).tocsr()
print("V = W_p1 + W_p2 : dim=%d nnz=%d"%(V.shape[0],V.nnz))
ev=eigsh(V,k=3,which="LA",return_eigenvectors=False,tol=1e-10)
ev=np.sort(ev)[::-1]
print("top eigenvalues of V (jm=1, d=5, 7 links):", np.round(ev,10))
print("cos(pi/3) = %.10f ; 2cos(pi/3) = %.10f"%(np.cos(np.pi/3),2*np.cos(np.pi/3)))
print("single-plaquette W_max (jm=1) = %.10f"%eigsh(W4,k=1,which="LA",return_eigenvectors=False,tol=1e-12)[0])
