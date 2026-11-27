import numpy as np, os
HERE=r"research/v3/species/w7"
M1=1.546; NB=20
def blocks(S, idx, Lt):
    X=S[idx].reshape(len(idx),30,Lt).astype(np.float64)
    X=X-X.mean(axis=0,keepdims=True)
    Ct=np.zeros((Lt,30,30))
    for t in range(Lt):
        Z=np.roll(X,-t,axis=2)
        Ct[t]=np.einsum('cai,cbi->ab',X,Z,optimize=True)/(len(idx)*Lt)
    return Ct
def lam(Ct,t0=0,t=1):
    G=0.5*(Ct[t0]+Ct[t0].T); w,V=np.linalg.eigh(G)
    keep=w>1e-2*w.max(); W=V[:,keep]/np.sqrt(w[keep])[None,:]
    M=W.T@(0.5*(Ct[t]+Ct[t].T))@W
    return np.sort(np.linalg.eigvalsh(0.5*(M+M.T)))[::-1], int(keep.sum())
f=os.path.join(HERE,"gh3_probe.npz")
d=np.load(f); S=d['S']; Lt=int(d['Lt']); N=S.shape[0]
print("probe: N",N,"Lt",Lt,"Ls",d['Ls'],"beta",d['beta'])
idx=np.array_split(np.arange(N),NB)
Ctb=[blocks(S,ix,Lt) for ix in idx]
Ct=np.mean(Ctb,axis=0)
del S, d
for lab,sel in (("FULL30",list(range(30))),("B2",[0,3,4])):
    l0,rank=lam(Ct[:,sel][:,:,sel],0,1)
    print(f"\n{lab} rank={rank}")
    for t in range(1,11):
        l,_=lam(Ct[:,sel][:,:,sel],0,t)
        r=l[1]/l[0] if (len(l)>1 and l[0]>0) else np.nan
        # JK
        vs=[]
        for b in range(NB):
            Cb=np.mean([Ctb[j] for j in range(NB) if j!=b],axis=0)
            lb,_=lam(Cb[:,sel][:,:,sel],0,t)
            vs.append(lb[1]/lb[0] if (len(lb)>1 and lb[0]>0) else np.nan)
        v=np.array(vs,float); m=np.nanmean(v)
        e=np.sqrt((NB-1)/NB*np.nansum((v-m)**2))
        D=-np.log(r) if r>0 else np.nan
        y=2*M1*D+D*D if np.isfinite(D) else np.nan
        sy=(2*M1+2*D)*e/m if (m>0 and np.isfinite(D)) else np.nan
        print(f"  t={t:2d} lam1={l[0]:.3e} r={r:.4f}+-{e:.4f} Delta={D:.4f} y={y:7.3f}+-{sy:6.3f}")
