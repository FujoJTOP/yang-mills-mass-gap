#!/usr/bin/env python3
"""d2_quick.py -- fast companion: is Var(a) bounded (=> G''->positive) or ~L (=> G''->0)?
1xL strip, k=0 LONG mode a=(1/sqrt L) sum xi_{l,3}.  Local-update Metropolis, Var(a) only.
Unbuffered prints.  py -3.12 d2_quick.py
"""
import sys, numpy as np

def qmul(a,b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return np.array([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0])
def qexp(x):
    th=np.linalg.norm(x)
    if th<1e-12: return np.array([1.0,0,0,0])
    n=x/th; s=np.sin(th/2); return np.concatenate([[np.cos(th/2)], n*s])
def qdag(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def logJ(r):
    if r<1e-9: return np.log(0.25)
    return 2*np.log(max(np.sin(r/2),1e-300))-2*np.log(max(r,1e-300))

def plaq_sum(xi,B):
    L=len(xi); q=[qexp(x) for x in xi]; tot=0.0
    tot+=qmul(q[0],q[1])[0]
    for i in range(1,L-1): tot+=qmul(qmul(B[i-1],qdag(q[i])),q[i+1])[0]
    tot+=qmul(B[-1],qdag(q[-1]))[0]
    return tot

def run(L,beta,theta,nsweep,seed,step=0.5,thin=1):
    rng=np.random.default_rng(seed)
    b=qexp(np.array([0,0,theta])); B=[b.copy() for _ in range(L-1)]
    xi=[rng.uniform(-0.6,0.6,3) for _ in range(L)]
    def logp(): return -beta*(L-plaq_sum(xi,B))+sum(logJ(np.linalg.norm(x)) for x in xi)
    lp=logp(); v=np.array([0,0,1.0]); acc=0; ntry=0; samples=[]
    burn=nsweep//5
    for sw in range(nsweep):
        for l in range(L):
            old=xi[l].copy(); new=old+rng.uniform(-step,step,3)
            if np.linalg.norm(new)>=2*np.pi: ntry+=1; continue
            xi[l]=new; nl=logp(); ntry+=1
            if np.log(rng.random()+1e-300)<(nl-lp): lp=nl; acc+=1
            else: xi[l]=old
        if sw>=burn and sw%thin==0:
            samples.append(sum(x[2] for x in xi)/np.sqrt(L))
    return np.array(samples),acc/max(ntry,1)

print("beta    L    Var(a)   1/Var   acc",flush=True)
for beta in (0.0,2.5):
    for L in (2,3,4,6,8,10,14):
        a,acc=run(L,beta,np.pi,120000 if L<=8 else 90000,seed=7+L+int(10*beta),thin=3)
        print(f"{beta:4.1f}{L:4d}{np.var(a):10.4f}{1/np.var(a):8.4f}{acc:6.3f}",flush=True)
print("HAAR FLOOR 1/s^2 = 9/(4pi^2-6) =",9/(4*np.pi**2-6),flush=True)
