import numpy as np, json
d=np.load("run3_L8_b2.5_N12k.npz"); S=d["S"]; Lt=int(d["Lt"])
def corr(Sm,Lt):
    x=Sm-Sm.mean(axis=0); C=np.zeros(Lt)
    for t in range(Lt): C[t]=np.mean(np.sum(x*np.roll(x,-t,axis=1),axis=1))/Lt
    return C
def A_at(C,Lt,m):
    ts=np.arange(1,4); g=np.exp(-m*ts)+np.exp(-m*(Lt-ts))
    return float(np.sum(g*C[ts])/np.sum(g*g))
def jack(Sm,fn,nb=20):
    N=Sm.shape[0]; idx=np.array_split(np.arange(N),nb); v=[]
    for b in range(nb):
        keep=np.concatenate([idx[j] for j in range(nb) if j!=b]); v.append(fn(corr(Sm[keep],Lt)))
    v=np.array(v); m=np.nanmean(v); e=np.sqrt((nb-1)/nb*np.nansum((v-m)**2)); return m,e
out={}
for lev,name in [(0,"raw"),(1,"TAPE1"),(2,"TAPE2")]:
    m,e=jack(S[:,lev,:],lambda C:A_at(C,Lt,1.546)/C[0])
    print(f"{name:6s} w1@1.546 = {m:.4f} +- {e:.4f}  (N={S.shape[0]})")
    out[name]=[m,e]
m0,e0=out["raw"]
print(f"raw vs 0.3423: diff={m0-0.3423:+.4f}  = {(m0-0.3423)/e0:+.2f} sigma")
print(f"raw vs 0.3423(run1 err .0151): {(m0-0.3423)/np.hypot(e0,0.0151):+.2f} sigma combined")
json.dump(out,open("run3_w1_jk.json","w"),indent=1)
