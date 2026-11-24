#!/usr/bin/env python3
"""d2_curv.py -- peak curvature G''(a0*)=-d^2 log rho/da^2 at the histogram peak, plus Var, few L."""
import numpy as np
from d2_batch import run

def gpp(A, nbin=120):
    A=A[::25]                       # subsample for speed
    lo,hi=np.percentile(A,0.5),np.percentile(A,99.5)
    h,e=np.histogram(A,bins=nbin,range=(lo,hi),density=True); c=0.5*(e[1:]+e[:-1])
    hs=np.maximum(np.convolve(h,np.ones(5)/5,mode='same'),1e-12)
    im=int(np.argmax(hs))
    d2=(np.log(hs[im+1])-2*np.log(hs[im])+np.log(hs[im-1]))/(c[1]-c[0])**2
    # 95% CI on peak by bootstrap
    boots=[]
    for s in range(40):
        B=A[np.random.default_rng(s).integers(0,len(A),len(A)//2)]
        hb,_=np.histogram(B,bins=nbin,range=(lo,hi),density=True)
        hb=np.maximum(np.convolve(hb,np.ones(5)/5,mode='same'),1e-12)
        j=int(np.argmax(hb)); boots.append(float(-(np.log(hb[min(j+1,len(c)-1)])-2*np.log(hb[j])+np.log(hb[max(j-1,0)]))/(c[1]-c[0])**2))
    return -d2, np.std(boots), c[im]

print(f"{'beta':>5}{'L':>4}{'Gpp(a0*)':>10}{'sd':>8}{'Var(a)':>9}{'1/Var':>8}{'a0*':>8}",flush=True)
for beta,L in [(2.5,3),(2.5,4),(2.5,8),(2.5,14),(0.0,4),(0.0,10)]:
    A=run(L,beta,np.pi,430,seed=201+L+int(10*beta))
    g,sd,a0=gpp(A); var=np.var(A)
    print(f"{beta:5.1f}{L:4d}{g:10.4f}{sd:8.4f}{var:9.4f}{1/var:8.4f}{a0:8.4f}",flush=True)
print("Done.",flush=True)
