import numpy as np
from scipy.linalg import eigh_tridiagonal
N=6; t=1.3
d=np.array([k*(k+2)/4.0 for k in range(N)]); off=np.ones(N-1)
P=np.diag(off,1)+np.diag(off,-1)
E,V=eigh_tridiagonal(d,t*off)
m=V.T@(P@V[:,0]); M=V.T@P@V
def align(W,ref):
    for j in range(W.shape[1]):
        if W[:,j]@ref[:,j]<0: W[:,j]*=-1
    return W
for h in (1e-3,1e-4,1e-5,1e-6):
    Ep,Vp=eigh_tridiagonal(d,(t+h)*off); Em,Vm=eigh_tridiagonal(d,(t-h)*off)
    Vp=align(Vp,V); Vm=align(Vm,V)
    mp=Vp.T@(P@Vp[:,0]); mm=Vm.T@(P@Vm[:,0])
    fd=(mp-mm)/(2*h)
    print(h, np.round(fd[1:],6))
# candidate: only the dE-driven part? test vdot0 by FD
h=1e-6
Ep,Vp=eigh_tridiagonal(d,(t+h)*off); Em,Vm=eigh_tridiagonal(d,(t-h)*off)
Vp=align(Vp,V); Vm=align(Vm,V)
vdot0=(Vp[:,0]-Vm[:,0])/(2*h)
form0=-sum(m[k]/(E[k]-E[0])*V[:,k] for k in range(1,N))
print("vdot0 err:", np.max(np.abs(vdot0-form0)))
vdot1=(Vp[:,1]-Vm[:,1])/(2*h)
form1=sum((M[k,1]/(E[1]-E[k]))*V[:,k] for k in range(N) if k!=1)
print("vdot1 err:", np.max(np.abs(vdot1-form1)))
print("m      :", np.round(m[1:],6))
print("E      :", np.round(E,6))
print("M[0,:] :", np.round(M[0,:],6))
print("diagM  :", np.round(np.diag(M),6))
