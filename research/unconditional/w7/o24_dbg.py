import numpy as np
src = open("o24_lat.py").read().split("mode = sys.argv[1]")[0]
K = {}
exec(compile(src, "o24_lat_core", "exec"), K)
v = np.zeros((K["NL"],3)); v[:,2] = 1/np.sqrt(K["NL"])
# replicate the chain but keep the final x
run = K["run"]; quat=K["quat"]; ht=K["ht"]; phi_r2=K["phi_r2"]
import types
code = src.replace("    return Ks, Fs, acc, gacc / (ns + burn), epsg, ndom / ((ns + burn) + (ns + burn) * NL)",
                   "    return Ks, Fs, acc, gacc / (ns + burn), epsg, ndom / ((ns + burn) + (ns + burn) * NL), x")
K2 = {}
exec(compile(code, "o24_lat_core2", "exec"), K2)
Ks,Fs,acc,gacc,epsg,dom,x = K2["run"](v,1,50,0.9,0.2,np.random.default_rng([7]),burn=400)
r = np.linalg.norm(x[0],axis=1)
print("r: min",r.min(),"max",r.max(),"mean",r.mean(),"n>6",(r>6).sum(),"2pi",2*np.pi)
print("kappa sample",Ks[:,0][:5])
l = int(np.argmax(r)); print("worst link",l,"r",r[l])
phi=K2["phi_r2"]; BETA=K2["BETA"]; NPQ=K2["NPQ"]; h=1e-3
def S_P(xx):
    w,a=quat(xx); tr=ht(w,a,K2["C0"],K2["C1"],K2["C2"],K2["C3"])
    return BETA*(NPQ-tr.sum(-1)), phi((xx*xx).sum(-1)).sum(-1)
for hh in (1e-3,):
    xp=x[0:1]+hh*v; xm=x[0:1]-hh*v
    S0,P0=S_P(x[0:1]); Sp,Pp=S_P(xp); Sm,Pm=S_P(xm)
    print("kS",((Sp-2*S0+Sm)/hh**2)[0],"kPhi",((Pp-2*P0+Pm)/hh**2)[0])
    print("P0",P0,"Pp",Pp,"Pm",Pm)

