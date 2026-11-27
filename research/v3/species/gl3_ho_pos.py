import numpy as np
from fractions import Fraction
from math import factorial
exec(open("gl3_ho_moments.py").read().split("# ---------------------------------------------------------------- 1. axis kappa")[0])

# leading coefficient rho_O of h_O at small p (same units as h_mn)
eps=1e-4
def rho_num(m,n):
    z0=np.exp(1j*eps); one=np.exp(0j)*np.ones(1)
    return np.real(h_mn(m,n,z0,1.0,1.0,1.0))/eps**2

print("rho_num(1,1) =",rho_num(1,1),"  rho_num(1,2) =",rho_num(1,2),"  rho_num(2,2) =",rho_num(2,2),"  rho_num(3,3) =",rho_num(3,3))

N=25
g=(2*np.pi*np.arange(N)/N)
Z=[np.exp(1j*g).reshape([N if i==k else 1 for k in range(4)]) for i in range(4)]
def minh(terms):
    h=0.0
    for (m,n),d in terms:
        h=h+ (d/rho_num(m,n))*h_mn(m,n,Z[0],Z[1],Z[2],Z[3])
    h=np.real(h); f=h.ravel(); i=int(np.argmin(f)); j=np.unravel_index(i,h.shape)
    return f[i], [round(g[k]/np.pi,3) for k in j]
def ivd(M):
    out=[]
    for m in range(1,M+1):
        v=1.0
        for l in range(1,M+1):
            if l!=m: v*= l*l/(l*l-m*m)
        out.append(v)
    return out
print()
print("corrected positivity: coefficients c=d/rho  (physical normalization h~|p|^2)")
for M in (1,2,3,4,5):
    d=ivd(M); terms=[((m,m),d[m-1]) for m in range(1,M+1)]
    print(f" M={M} k0={2*M}: min h = {minh(terms)[0]:+.9f} at {minh(terms)[1]}*pi ; c = {[round(d[m-1]/rho_num(m,m),6) for m in range(1,M+1)]}")
print()
print("actual Luescher-Weisz action (c0=5/3 on plaq, c1=-1/12 on 1x2 rect):")
print("  min h =", minh([((1,1),5/3*rho_num(1,1)),((1,2),-1/12*rho_num(1,2))]))
print("Wilson:")
print("  min h =", minh([((1,1),rho_num(1,1))]))
