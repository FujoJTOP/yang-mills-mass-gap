"""
tool_hardy_witness.py -- TOOL-HARDY
Hardy/H^2 norm of the Stieltjes resolvent boundary value as a NON-INTEGRAL atom detector.
Model KL measure: atom f1^2 at m1 (=1)  +  flat continuum rho0 on [4,8].
G(x+i0) = sum_i w_i/(s_i - x - i0); for x<4 (below continuum) Im G = 0 except the atom delta.
Test: truncated L^2 norm N(eps)=int_{eps<|x-m1|<delta}|G|^2 dx.
  atom  -> N(eps) ~ 2 f1^4 / eps   (diverges; N*eps -> 2 f1^4)
  smooth bump (Lorentzian width gamma, same mass) -> N(eps) -> finite  (pi f1^4 / gamma)
  atomless (f1^2=0) -> N(eps) -> finite  (no divergence)  => acid-clean
"""
import numpy as np
f1sq=0.5; m1=1.0; rho0=0.25   # continuum total mass = rho0*(8-4)=1
def cont(x): return rho0*(np.log(8.0-x)-np.log(4.0-x))
def G_atom(x): return f1sq/(m1-x) + cont(x)
def G_bump(x,gamma): return f1sq/(m1 - x - 1j*gamma) + cont(x)
def G_free(x): return cont(x)                      # f1^2=0
def Nnorm(G,eps,delta,n=400001):
    a=np.linspace(m1-delta,m1-eps,n); b=np.linspace(m1+eps,m1+delta,n)
    return np.trapezoid(np.abs(G(a))**2,a)+np.trapezoid(np.abs(G(b))**2,b)
delta=0.3; epss=[1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4]
print("predicted 2*f1^4 =",2*f1sq**2)
print("\n== ATOM (should diverge, N*eps -> 2 f1^4) ==")
for e in epss:
    N=Nnorm(G_atom,e,delta); print("eps=%.0e  N=%10.4f  N*eps=%8.5f"%(e,N,N*e))
print("\n== SMOOTH BUMP gamma=0.05, same mass (should stay finite) ==")
for e in epss:
    N=Nnorm(lambda x:G_bump(x,0.05),e,delta); print("eps=%.0e  N=%10.4f  N*eps=%8.5f"%(e,N,N*e))
print("\n== ATOMLESS f1^2=0 (acid test: finite, no divergence) ==")
for e in epss:
    N=Nnorm(G_free,e,delta); print("eps=%.0e  N=%10.4f  N*eps=%8.5f"%(e,N,N*e))
# recover f1^2 from the divergence coefficient of the atom case
e=1e-4; N=Nnorm(G_atom,e,delta); f1rec=np.sqrt(max(N*e,0)/2.0)
print("\nrecovered f1^2 from atom divergence rate =",f1rec," (true %.3f)"%f1sq)
# robustness: swap continuum mass, atom divergence coefficient unchanged?
print("\n== atom, different continuum mass rho0 (divergence rate should be continuum-robust) ==")
for r0 in [0.0,0.25,1.0]:
    globals()['rho0']=r0
    e=1e-4; N=Nnorm(G_atom,e,delta); print("rho0=%.2f  N*eps=%.5f"%(r0,N*e))
