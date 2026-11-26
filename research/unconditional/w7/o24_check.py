import sys
sys.argv = ["o24_lat.py", "A", "1", "0.9"]
import numpy as np, o11_physical as o, o24_lat as m
rng=np.random.default_rng(5)
x=rng.normal(0,0.9,(4,64,3))
V,tr,ph=m.full_V(x)
h=1e-3
v=np.zeros((64,3)); v[:,2]=1/np.sqrt(64)
Vp=m.full_V(x+h*v[None])[0]; Vm=m.full_V(x-h*v[None])[0]
K=(Vp-2*V+Vm)/h**2; F=(Vp-Vm)/(2*h)
for g in range(4):
    a=o.Vtot(x[g],m.tab,m.BETA); b=o.Vtot(x[g]+h*v,m.tab,m.BETA); c=o.Vtot(x[g]-h*v,m.tab,m.BETA)
    print(f"g={g} V={V[g]:+.8f}/{a:+.8f} dV={V[g]-a:+.2e} kap={K[g]:+.6f}/{(b-2*a+c)/h**2:+.6f} F={F[g]:+.6f}/{(b-c)/(2*h):+.6f}")
