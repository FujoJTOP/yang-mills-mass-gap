import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=2; bs=1.0; bt=2.5
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
idx=[(-1,-1),(-1,1),(1,-1),(1,1)]
tot=None
for (p,q) in idx:
    A=u[(p,q)]
    for (r,s) in idx:
        Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
        term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
        tot=term if tot is None else tot+term
W4=(0.5*tot).tocsr()
D4=dl**4
print("dl=%d D4=%d"%(dl,D4), flush=True)
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
# FULL K diagonal
Kfull=np.ones(D4)
for ii in range(D4):
    x=ii; p=1.0
    for _ in range(4):
        p*=js[x%dl]; x//=dl
    Kfull[ii]=p
Kd_full=sp.diags(Kfull)
def Mf_full(v): return np.exp(-2*bs)*expm_multiply(bs*W4, v)
e0=np.zeros(D4); e0[0]=1.0
w0=Mf_full(e0); scale=float(np.dot(w0,Kfull*w0))
def T_full(v): return Mf_full(Kd_full*(Mf_full(v)))/scale
op=LinearOperator((D4,D4),matvec=T_full,dtype=float)
ev=eigsh(op,k=3,which="LA",tol=1e-13,maxiter=5000,return_eigenvectors=False)
ev=np.sort(ev)[::-1]
print("SINGLE plaquette EXACT jm=2: t1=%.10f t2=%.10f r=%.10f"%(ev[0],ev[1],ev[1]/ev[0]), flush=True)
# VARIATIONAL with the 257-state basis (all j=0 or all four j=1/2)
z=(0,0,0); e12=[b for b in B if b[0]==1]
basis4=[(z,z,z,z)]
for a in e12:
    for b in e12:
        for c in e12:
            for dd in e12:
                basis4.append((a,b,c,dd))
nb4=len(basis4)
rows=[]; cols=[]
for j,(a,b,c,dd) in enumerate(basis4):
    ii=0
    for t in (a,b,c,dd): ii=ii*dl+index[t]
    rows.append(ii); cols.append(j)
P=sp.csr_matrix((np.ones(nb4),(rows,cols)),shape=(D4,nb4))
W4r=(P.T@W4@P).tocsr()
K_r=Kfull[[0]+[0]*0]  # placeholder
idxfull=[]
for j,(a,b,c,dd) in enumerate(basis4):
    ii=0
    for t in (a,b,c,dd): ii=ii*dl+index[t]
    idxfull.append(ii)
Kdiag_r=Kfull[idxfull]
def Mf_r(v): return np.exp(-2*bs)*expm_multiply(bs*W4r, v)
e0r=np.zeros(nb4); e0r[0]=1.0
w0r=Mf_r(e0r); scale_r=float(np.dot(w0r,Kdiag_r*w0r))
def T_r(v): return Mf_r(Kdiag_r*(Mf_r(v)))/scale_r
op_r=LinearOperator((nb4,nb4),matvec=T_r,dtype=float)
evr=eigsh(op_r,k=3,which="LA",tol=1e-13,maxiter=5000,return_eigenvectors=False)
evr=np.sort(evr)[::-1]
print("SINGLE plaquette VARIATIONAL (257) : t1=%.10f t2=%.10f r=%.10f"%(evr[0],evr[1],evr[1]/evr[0]), flush=True)
