from fractions import Fraction as F
from math import factorial
K=22
def a(nu,k):
    num=F(1)
    for j in range(1,k+1):
        num*= (4*nu*nu-(2*j-1)**2)
    return F((-1)**k)*num/F(factorial(k)*8**k)
# series for I2 and I1 in z=1/beta, coefficients a_k
A2=[a(F(2),k) for k in range(K)]
A1=[a(F(1),k) for k in range(K)]
# ratio series R = A2/A1  (R_0 = 1)
R=[F(0)]*K
for n in range(K):
    s=A2[n]
    for j in range(1,n+1):
        s-= R[j]*A1[n-j]
    R[n]=s/A1[0]
# log R = log(1 + u) with u = R-1
U=[R[n] if n>0 else F(0) for n in range(K)]
U[0]=F(0)
logR=[F(0)]*K
# log(1+u) = sum_{m>=1} (-1)^{m+1} u^m/m
pw=[F(0)]*K; pw[1]=F(1)     # u^0 handled separately
# compute u^m iteratively
cur=[F(0)]*K; cur[0]=F(1)
for m in range(1,K):
    new=[F(0)]*K
    for i in range(K):
        if cur[i]==0: continue
        for j in range(1,K-i):
            if U[j]!=0: new[i+j]+=cur[i]*U[j]
    cur=new
    for i in range(K):
        logR[i]+= F((-1)**(m+1))*cur[i]/F(m)
print("exact coefficients of  -ln(I2/I1) = sum c_k / beta^k")
c=[-logR[k] for k in range(1,K)]
for k in range(1,K):
    print("   c_%-2d = %-22s  %+.10f"%(k,str(c[k-1]),float(c[k-1])))
print()
print("growth ratios |c_k/c_{k-1}|:")
for k in range(2,K):
    print("   k=%-2d : %.6f"%(k,abs(float(c[k-1]))/abs(float(c[k-2]))))
