import numpy as np, functools, re
src=open("research/unconditional/w7/_lam2pb.py").read()
src=src.replace('return np.log(ev[1]/ev[0])','return np.log(ev[1]/ev[0])')
src=src.replace('ys=np.array([lnr(bt,h) for h in hs])','r0=lnr(bt,0.0); ys=np.array([lnr(bt,h)-r0 for h in hs])')
open("research/unconditional/w7/_lam2pd.py","w").write(src)
print("patched")
