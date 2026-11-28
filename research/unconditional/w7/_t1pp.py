import numpy as np
from scipy.special import iv
def t1(lam,beta=2.5,J=200):
    a=np.arange(J+1); t=2*iv(a+1,beta)/beta
    M=iv(np.abs(a[:,None]-a[None,:]),lam)-iv(a[:,None]+a[None,:]+2,lam)
    T=M@np.diag(t)@M; T=(T+T.T)/2
    return np.linalg.eigvalsh(T)[-1]
def t1pp(lam,h=1e-3,beta=2.5):
    return (t1(lam+h,beta)-2*t1(lam,beta)+t1(lam-h,beta))/h**2
print("is t_1(lambda) concave in lambda?   (t_1'' < 0  <=>  I1 I3 < I2^2  in the chain)")
print("%10s %22s"%("lambda","t_1''"))
for lam in (0.01,0.1,0.5,1.0,2.0,4.0,8.0,16.0):
    print("%10.2f %22.8e"%(lam,t1pp(lam)), flush=True)
print()
print("reference: I1 I3 - I2^2 at beta=2.5  =  %.10e  (should be < 0)"%(iv(1,2.5)*iv(3,2.5)-iv(2,2.5)**2))
