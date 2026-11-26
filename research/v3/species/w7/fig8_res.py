import numpy as np
s=0.37
def G(R1,R2): return np.exp(-s*np.pi*(R1**2+R2**2))
h=1e-6
def D(f,i,R1,R2):
    e=[0,0]; e[i]=h
    return (f(R1+e[0],R2+e[1])-f(R1-e[0],R2-e[1]))/(2*h)
def D2(f,R1,R2):
    return (f(R1+h,R2)-2*f(R1,R2)+f(R1-h,R2))/h**2
def D12(f,R1,R2):
    return (f(R1+h,R2+h)-f(R1+h,R2-h)-f(R1-h,R2+h)+f(R1-h,R2-h))/(4*h*h)
def LHS(R1,R2,r1,r2):
    g=G(R1,R2)
    G1=D(G,0,R1,R2); G2=D(G,1,R1,R2)
    G11=D2(G,R1,R2); G22=(G(R1,R2+h)-2*g+G(R1,R2-h))/h**2
    G12=D12(G,R1,R2)
    a11=r1/(2*np.pi*R1); a22=r2/(2*np.pi*R2)
    a12=r1/(2*np.pi*R2)+r2/(2*np.pi*R1)
    a1=-r1/(2*np.pi*R1**2); a2=-r2/(2*np.pi*R2**2)
    L=a11*G11+a12*G12+a22*G22+a1*G1+a2*G2
    return L, L/g
print("sigma=",s)
print("case TANGENT fig8: rho=(1,0)")
for (R1,R2) in [(0.5,1.0),(1.0,1.0),(2.0,3.0),(5.0,0.7),(10.0,10.0)]:
    L,q=LHS(R1,R2,1.0,0.0)
    pred=4*np.pi*s*s*R1
    print(f" R1={R1:5} R2={R2:5}  LHS/G={q: .8e}  pred 4pi s^2 R1={pred: .8e}  ratio={q/pred: .10f}")
print("case GENERAL crossing: rho=(1,1)")
for (R1,R2) in [(0.5,1.0),(2.0,3.0)]:
    L,q=LHS(R1,R2,1.0,1.0)
    pred=4*np.pi*s*s*(R1+R2)
    print(f" R1={R1:5} R2={R2:5}  LHS/G={q: .8e}  pred={pred: .8e}  ratio={q/pred: .10f}")
print("power-law fit of LHS/G in R (tangent, R2=const):")
for R1 in [1.0,2.0,4.0,8.0]:
    L,q=LHS(R1,3.0,1.0,0.0); print(f"  R1={R1:4}  LHS/G={q:.6e}  /R1={q/R1:.6e}")
print("circle sector sigma=2 g^2 N check: is b=R-independent possible?  b=4 pi s^2 R1/(g^2N)")
print("  with s=2*g^2*N (g^2N=1): sigma=2 -> 4*pi*4*R1 = ",4*np.pi*4*1.0," per R1=1 ; grows with R1 -> no constant b")
