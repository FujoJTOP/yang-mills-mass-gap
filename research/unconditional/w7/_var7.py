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
z=(0,0,0); e12=[b for b in B if b[0]==1]
# pattern3: (all j=0) or (all three j=1/2)
pat3=[(z,z,z)]
for a in e12:
    for b in e12:
        for c in e12:
            pat3.append((a,b,c))
npt=len(pat3)
print("pattern3 size = %d (full %d)"%(npt,d**3), flush=True)
# 4-link truncated basis: c in B (d), leaves in pat3
basis4=[(c,l) for c in B for l in pat3]
nb4=len(basis4)
print("nb4 = %d"%nb4, flush=True)
rows=[]; cols=[]
for j,(c,l) in enumerate(basis4):
    ii=0
    for t in (c,)+l: ii=ii*d+index[t]
    rows.append(ii); cols.append(j)
P=sp.csr_matrix((np.ones(nb4),(rows,cols)),shape=(d**4,nb4))
W4r=(P.T@W4@P).tocsr()
print("W4r nnz=%d"%W4r.nnz, flush=True)
# index map: (c,l)->i
i4={}
for j,(c,l) in enumerate(basis4): i4[(index[c],tuple(index[x] for x in l))]=j
# variational basis: (c, l, r) with l,r in pat3, c in B
basis=[]
for ci in range(d):
    for lp,li in [(p,i) for i,p in enumerate(pat3)]:
        for rp,ri in [(p,i) for i,p in enumerate(pat3)]:
            basis.append((ci,li,ri))
Nv=len(basis)
print("Nv = %d"%Nv, flush=True)
bmap={b:k for k,b in enumerate(basis)}
# W1: acts on (c,l). W4r is indexed by (c,l).
# Build W1 on the variational space: (c,l,r)->(c',l',r)
W4r_coo=W4r.tocoo()
rowsW=[]; colsW=[]; valsW=[]
for ip,ii_,v in zip(W4r_coo.row,W4r_coo.col,W4r_coo.data):
    cp,lpi=basis4[ip]; cm,lmi=basis4[ii_]
    cpi=index[cp]; cmi=index[cm]
    # find pattern index
    lpi_i=pat3.index(lpi); lmi_i=pat3.index(lmi)
    for ri in range(npt):
        k=bmap[(cmi,lmi_i,ri)]
        kp=bmap[(cpi,lpi_i,ri)]
        rowsW.append(kp); colsW.append(k); valsW.append(v)
W1v=sp.csr_matrix((valsW,(rowsW,colsW)),shape=(Nv,Nv))
# W2: acts on (c,r)
rowsW=[]; colsW=[]; valsW=[]
for ip,ii_,v in zip(W4r_coo.row,W4r_coo.col,W4r_coo.data):
    cp,rpi=basis4[ip]; cm,rmi=basis4[ii_]
    cpi=index[cp]; cmi=index[cm]
    lpi_i=pat3.index(rpi); lmi_i=pat3.index(rmi)
    for li in range(npt):
        k=bmap[(cmi,li,lmi_i)]
        kp=bmap[(cpi,li,lpi_i)]
        rowsW.append(kp); colsW.append(k); valsW.append(v)
W2v=sp.csr_matrix((valsW,(rowsW,colsW)),shape=(Nv,Nv))
print("W1v nnz=%d W2v nnz=%d"%(W1v.nnz,W2v.nnz), flush=True)
Wtot=0.5*(W1v+W2v)
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js=np.array([td[b[0]] for b in B])
Kpat=np.array([np.prod([js[index[x]] for x in p]) for p in pat3])
Kdiag=np.array([js[ci]*Kpat[li]*Kpat[ri] for (ci,li,ri) in basis])
Kd=sp.diags(Kdiag)
def Mf(v): return np.exp(-2*bs)*expm_multiply(bs*Wtot, v)
e0=np.zeros(Nv); e0[0]=1.0
w0=Mf(e0); scale=float(np.dot(w0,Kdiag*w0))
print("scale=%.6e"%scale, flush=True)
def T(v): return Mf(Kd*(Mf(v)))/scale
op=LinearOperator((Nv,Nv),matvec=T,dtype=float)
ev,V=eigsh(op,k=3,which="LA",tol=1e-12,maxiter=5000)
o=np.argsort(ev)[::-1]; ev=ev[o]
print("VARIATIONAL (c in B, leaves in pat3), Nv=%d"%Nv, flush=True)
print("  t1=%.10f t2=%.10f t3=%.10f  r=t2/t1=%.10f"%(ev[0],ev[1],ev[2],ev[1]/ev[0]), flush=True)
for i in range(3):
    res=np.linalg.norm(T(V[:,i])-ev[i]*V[:,i])/abs(ev[i])
    print("   eig %d residual = %.2e"%(i,res), flush=True)
