import sys; sys.argv=["x","A","1","0.9"]
import numpy as np, o11_physical as o, o24_lat as m
BETA=m.BETA; NL=64; h=1e-3
v=np.zeros((NL,3)); v[:,2]=1/np.sqrt(NL)
def SP(x):
    w,a=m.quat(x); tr=m.ht(w,a,m.C0,m.C1,m.C2,m.C3)
    return BETA*(m.NPQ-tr.sum(-1)), m.phi_r2((x*x).sum(-1)).sum(-1)
rng=np.random.default_rng(9)
Ks,Fs,acc,gacc,epsg,dom=m.run(v,1,40,0.9,0.2,np.random.default_rng([9]),burn=300)
i=int(np.argmax(np.abs(Ks[:,0])))
print("max|k| sample idx",i,"kappa",Ks[i,0])
# replay to that config is hard; instead sample a few fresh configs and decompose
for t in range(6):
    Ks2,_,_,_,_,_=m.run(v,1,1,0.9,0.2,np.random.default_rng([100+t]),burn=400)
print("--- decompose over random configs from the chain ---")
