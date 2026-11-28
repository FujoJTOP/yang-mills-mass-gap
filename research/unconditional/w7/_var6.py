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
    for t in (a,b,c,dd): ii=ii*d+index[t]
    rows.append(ii); cols.append(j)
P=sp.csr_matrix((np.ones(nb4),(rows,cols)),shape=(d**4,nb4))
W4r=(P.T@W4@P).tocsr()
cidx=[index[b[0]] for b in basis4]
pairs=[(i,j) for i in range(nb4) for j in range(nb4) if cidx[i]==cidx[j]]
Np=len(pairs)
pmap={(i,j):k for k,(i,j) in enumerate(pairs)}
print("nb4=%d Np=%d"%(nb4,Np), flush=True)
# W4r entries
co=W4r.tocoo()
ip_arr=co.row.astype(int); ii_arr=co.col.astype(int); vv=co.data
# W1_r: (i,j)->(i',j) allowed only if cidx[i']==cidx[j]
rowsW=[]; colsW=[]; valsW=[]
for ip_,ii_,v in zip(ip_arr,ii_arr,vv):
    for j in range(nb4):
        if cidx[ip_]==cidx[j] and (ip_,j) in pmap and (ii_,j) in pmap:
            rowsW.append(pmap[(ip_,j)]); colsW.append(pmap[(ii_,j)]); valsW.append(v)
W1r=sp.csr_matrix((valsW,(rowsW,colsW)),shape=(Np,Np))
# W2_r: (i,j)->(i,j') allowed if cidx[i]==cidx[j']
rowsW=[]; colsW=[]; valsW=[]
for jp_,jj_,v in zip(ip_arr,ii_arr,vv):
    for i in range(nb4):
        if cidx[i]==cidx[jp_] and (i,jp_) in pmap and (i,jj_) in pmap:
            rowsW.append(pmap[(i,jp_)]); colsW.append(pmap[(i,jj_)]); valsW.append(v)
W2r=sp.csr_matrix((valsW,(rowsW,colsW)),shape=(Np,Np))
print("W1r nnz=%d W2r nnz=%d"%(W1r.nnz,W2r.nnz), flush=True)
Wtot=0.5*(W1r+W2r)
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
js4=np.array([td[b[0]] for b in B])
K4=np.array([js4[index[a]]*js4[index[b]]*js4[index[c]]*js4[index[e]] for (a,b,c,e) in basis4])
Kdiag=np.array([K4[i]*K4[j] for (i,j) in pairs])
Kd=sp.diags(Kdiag)
def Mf(v): return np.exp(-2*bs)*expm_multiply(bs*Wtot, v)
e0=np.zeros(Np); e0[0]=1.0
w0=Mf(e0); scale=float(np.dot(w0,Kdiag*w0))
print("scale=%.6e"%scale, flush=True)
def T(v): return Mf(Kd*(Mf(v)))/scale
op=LinearOperator((Np,Np),matvec=T,dtype=float)
ev,V=eigsh(op,k=3,which="LA",tol=1e-12,maxiter=5000)
o=np.argsort(ev)[::-1]; ev=ev[o]
print("VARIATIONAL (4-link truncated, constrained c), Np=%d"%Np, flush=True)
print("  t1=%.10f t2=%.10f t3=%.10f  r=t2/t1=%.10f"%(ev[0],ev[1],ev[2],ev[1]/ev[0]), flush=True)
for i in range(3):
    res=np.linalg.norm(T(V[:,i])-ev[i]*V[:,i])/abs(ev[i])
    print("   eig %d residual = %.2e"%(i,res), flush=True)
