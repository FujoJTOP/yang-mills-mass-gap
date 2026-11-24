#!/usr/bin/env python3
"""d2_batch.py -- batched vectorized Metropolis over many independent 1xL strips.
Measures Var(a) and the marginal curvature G''(a0*) of a=(1/sqrt L) sum_l xi_{l,3}.
Vectorized over NCHAIN; sequential over L (local plaquette updates).  py -3.12 d2_batch.py
"""
import sys, numpy as np
NCHAIN=12000

def qmul(a,b):
    a0,a1,a2,a3=a[...,0],a[...,1],a[...,2],a[...,3]
    b0,b1,b2,b3=b[...,0],b[...,1],b[...,2],b[...,3]
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0],axis=-1)
def qexp(x):
    th=np.linalg.norm(x,axis=-1,keepdims=True)
    n=np.divide(x,th,out=np.zeros_like(x),where=th>1e-12)
    return np.concatenate([np.cos(th/2), n*np.sin(th/2)],axis=-1)
def qdag(q): return q*np.array([1.,-1,-1,-1])
def logJ(r):
    s=np.sin(r/2)
    return np.where(r>1e-9,2*np.log(np.maximum(s,1e-300))-2*np.log(np.maximum(r,1e-300)),np.log(0.25))

def ptr_of(q,B,L):
    # half-trace of each plaquette P_0..P_{L-1};  q:(n,L,4)
    n=q.shape[0]; out=np.zeros((n,L))
    out[:,0]=qmul(q[:,0],q[:,1])[:,0]
    for i in range(1,L-1):
        out[:,i]=qmul(qmul(B[i-1],qdag(q[:,i])),q[:,i+1])[:,0]
    out[:,L-1]=qmul(B[L-2],qdag(q[:,L-1]))[:,0]
    return out

def local_plaqs(q,B,L,l,ql_new):
    # return the plaquette indices touching link l and their half-traces given ql_new (n,4), others = q (old)
    idx=[]
    if l==0: idx=[0]
    elif l==L-1: idx=[L-2,L-1]
    else: idx=[l-1,l]
    n=ql_new.shape[0]; out=np.zeros((n,len(idx)))
    for k,i in enumerate(idx):
        if i==0:
            out[:,k]=qmul(ql_new,q[:,1])[:,0] if l==0 else qmul(q[:,0],ql_new)[:,0]
        elif i==L-1:
            out[:,k]=qmul(B[L-2],qdag(ql_new))[:,0]
        else:
            # P_i = B[i-1] q_i^-1 q_{i+1}
            qi = ql_new if (i==l) else q[:,i]
            qip1 = ql_new if (i+1==l) else q[:,i+1]
            out[:,k]=qmul(qmul(B[i-1],qdag(qi)),qip1)[:,0]
    return idx,out

def run(L,beta,theta,nsweep,seed,step=0.5,burn_frac=0.35):
    rng=np.random.default_rng(seed)
    b=qexp(np.array([0,0,theta]))
    B=[b.copy() for _ in range(L-1)]
    X=rng.uniform(-0.7,0.7,size=(NCHAIN,L,3))
    q=qexp(X)
    ptr=ptr_of(q,B,L)
    sums=ptr.sum(axis=1)
    logJc=logJ(np.linalg.norm(X,axis=-1)).sum(axis=1)
    lp=-beta*(L-sums)+logJc
    acc=np.zeros(NCHAIN); ntry=np.zeros(NCHAIN)
    As=[]; burn=int(burn_frac*nsweep)
    for sw in range(nsweep):
        for l in range(L):
            new=X[:,l,:]+rng.uniform(-step,step,size=(NCHAIN,3))
            rn=np.linalg.norm(new,axis=-1)
            ok=rn<2*np.pi
            ql=qexp(new)
            idx,newptr=local_plaqs(q,B,L,l,ql)
            oldptr=ptr[:,idx]
            dsum=newptr.sum(axis=1)-oldptr.sum(axis=1)
            newlogJ=logJc-logJ(np.linalg.norm(X[:,l,:],axis=-1))+logJ(rn)
            dlp=-beta*(-dsum)+(newlogJ-logJc)
            accm=ok & (np.log(rng.random(NCHAIN)+1e-300)<dlp)
            # apply
            rows=np.nonzero(accm)[0]
            q[accm,l,:]=ql[accm]
            X[accm,l,:]=new[accm]
            ptr[rows[:,None],idx]=newptr[accm]
            lp[accm]+=dlp[accm]
            logJc[accm]=newlogJ[accm]
            acc+=accm; ntry+=ok
        if sw>=burn:
            As.append(np.sum(X[:,:,2],axis=1)/np.sqrt(L))
    A=np.array(As).ravel()
    return A

def gpp(A,nbin=140):
    lo,hi=np.percentile(A,0.5),np.percentile(A,99.5)
    h,e=np.histogram(A,bins=nbin,range=(lo,hi),density=True); c=0.5*(e[1:]+e[:-1])
    k=np.ones(5)/5; hs=np.maximum(np.convolve(h,k,mode='same'),1e-12)
    im=int(np.argmax(hs))
    if im<1 or im>=len(c)-1: return np.nan,np.nan,c[im]
    d2=(np.log(hs[im+1])-2*np.log(hs[im])+np.log(hs[im-1]))/(c[1]-c[0])**2
    return -d2,np.var(A),c[im]

if __name__=="__main__":
    print("1/s^2 = 9/(4pi^2-6) =",9/(4*np.pi**2-6),flush=True)
    print(f"{'beta':>5}{'L':>4}{'Var(a)':>10}{'1/Var':>9}",flush=True)
    for beta in (2.5,0.0):
        for L in (2,3,4,6,8,10,14):
            ns=430
            A=run(L,beta,np.pi,ns,seed=101+L+int(10*beta))   # A: snapshots after burn
            var=np.var(A)
            print(f"{beta:5.1f}{L:4d}{var:10.4f}{1/var:9.4f}",flush=True)
    print("READ: beta=0 Var->s^2=3.72, 1/Var->0.269.  If beta=2.5 Var stays bounded as L grows")
    print("      => G''->positive Haar floor (NOT H1).  Done. No QED.",flush=True)
