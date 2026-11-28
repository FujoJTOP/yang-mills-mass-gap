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
index={b:i for i,b in enumerate(B)}; d=len(B)
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((d,d))
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
z=(0,0,0)
# Basis: all j=0, or all four at j=1/2 (256 states)
e12=[b for b in B if b[0]==1]
basis4=[(z,z,z,z)]
for a in e12:
    for b in e12:
        for c in e12:
            for dd in e12:
                basis4.append((a,b,c,dd))
nb4=len(basis4)
print("nb4=%d"%nb4, flush=True)
rows=[]; cols=[]
for j,(a,b,c,dd) in enumerate(basis4):
    ii=0
    for t in (a,b,c,dd): ii=ii*d+index[t]
    rows.append(ii); cols.append(j)
P=sp.csr_matrix((np.ones(nb4),(rows,cols)),shape=(d**4,nb4))
W4r=(P.T@W4@P).tocsr()
print("W4r nnz=%d"%W4r.nnz, flush=True)
# exact full operator
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
Kdiag_full=np.array([js[index[a]]*js[index[b]]*js[index[c]]*js[index[e_]] for (a,b,c,e_) in basis4])
Kd_full=sp.diags(Kdiag_full)
def Mf_full(v): return np.exp(-2*bs)*expm_multiply(bs*W4, v)
e0=np.zeros(d**4); e0[0]=1.0
w0=Mf_full(e0); scale_full=float(np.dot(w0,Kdiag_full*w0))
def T_full(v): return Mf_full(Kd_full*(Mf_full(v)))/scale_full
op_full=LinearOperator((d**4,d**4),matvec=T_full,dtype=float)
ev,V=eigsh(op_full,k=3,which="LA",tol=1e-13,maxiter=5000)
o=np.argsort(ev)[::-1]; ev=ev[o]
print("SINGLE plaquette EXACT jm=2: t1=%.10f t2=%.10f r=%.10f"%(ev[0],ev[1],ev[1]/ev[0]), flush=True)
# now the VARIATIONAL estimate with the 257-state basis
Kd_r=sp.diags(Kdiag_full)
def Mf_r(v): return np.exp(-2*bs)*expm_multiply(bs*W4r, v)
e0r=np.zeros(nb4); e0r[0]=1.0
w0r=Mf_r(e0r); scale_r=float(np.dot(w0r,Kdiag_full*w0r))
def T_r(v): return Mf_r(Kd_full*(Mf_r(v)))/scale_r
op_r=LinearOperator((nb4,nb4),matvec=T_r,dtype=float)
evr,Vr=eigsh(op_r,k=3,which="LA",tol=1e-13,maxiter=5000)
or_=np.argsort(evr)[::-1]; evr=evr[or_]
print("SINGLE plaquette VARIATIONAL (257 basis): t1=%.10f t2=%.10f r=%.10f"%(evr[0],evr[1],evr[1]/evr[0]), flush=True)
for i in range(min(3,nb4)):
    res=np.linalg.norm(T_r(Vr[:,i])-evr[i]*Vr[:,i])/abs(evr[i])
    print("   variational eig %d residual = %.2e"%(i,res), flush=True)

