from fractions import Fraction as F
from math import factorial
K=24
def a(nu,k):
    num=F(1)
    for j in range(1,k+1):
        num*= (4*nu*nu-(2*j-1)**2)
    return F((-1)**k)*num/F(factorial(k)*8**k)
A2=[a(F(2),k) for k in range(K)]
A1=[a(F(1),k) for k in range(K)]
R=[F(0)]*K
for n in range(K):
    s=A2[n]
    for j in range(0,n):          # j = 0..n-1  (correct recurrence)
        s-= R[j]*A1[n-j]
    R[n]=s/A1[0]
U=[R[n] for n in range(K)]; U[0]=F(0)
logR=[F(0)]*K
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
c=[-logR[k] for k in range(1,K)]
print("CORRECT: exact coefficients of -ln(I2/I1) = sum c_k/beta^k")
for k in range(1,17):
    print("   c_%-2d = %-26s %+.10f"%(k,str(c[k-1]),float(c[k-1])))
print("   ...")
print()
print("growth |c_k/c_{k-1}|  (should tend to 1/z_sing, the Borel singularity):")
for k in range(2,24):
    print("   k=%-2d : %.8f"%(k,abs(float(c[k-1]))/abs(float(c[k-2]))))
