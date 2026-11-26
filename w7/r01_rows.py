import numpy as np
def w(r): return np.sin(0.5*r)/r
def wp(r): return (0.5*r*np.cos(0.5*r)-np.sin(0.5*r))/r**2
def wpp(r): return -np.sin(0.5*r)/(4*r)-np.cos(0.5*r)/r**2+2*np.sin(0.5*r)/r**3
def f(a,rho,b,th):
    c,t=np.cos(0.5*th),np.sin(0.5*th); u=a/rho
    G=b*((c/2)*u*np.sin(0.5*rho)+t*(w(rho)+(a*a/rho)*wp(rho)))
    H=b*((c/2)*(np.sin(0.5*rho)/rho+(a*a/(2*rho**2))*np.cos(0.5*rho)-(a*a/rho**3)*np.sin(0.5*rho))
         +t*((3*a/rho)*wp(rho)-(a**3/rho**3)*wp(rho)+(a**3/rho**2)*wpp(rho)))
    return G,H
b,th=2.5,np.pi/2
for a in [0.3,0.75,1.2]:
    r=np.linspace(a+1e-9,2*np.pi,2000001); G,H=f(a,r,b,th)
    idx=np.where(np.sign(G[:-1])!=np.sign(G[1:]))[0]
    print("a=%.2f  G-zero rho:"%a, np.round(r[idx],6), " G at |a|:",round(float(f(a,np.array([a+1e-9]),b,th)[0][0]),6))
    for rr in [a+0.001,a+0.4,1.5,2.5,3.5,4.5,5.5,6.28]:
        g,h=f(a,np.array([rr]),b,th)
        print("    rho=%.3f  G=%+.6f  H=%+.6f  H/G=%+.6f"%(rr,g[0],h[0],h[0]/g[0]))
