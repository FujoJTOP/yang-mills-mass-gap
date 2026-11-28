import numpy as np
from scipy.special import iv
def d2(f,z,h=1e-3): return (f(z+h)-2*f(z)+f(z-h))/h**2
cands={
 "I1/I0"      : lambda z: iv(1,z)/iv(0,z),
 "I2/I1"      : lambda z: iv(2,z)/iv(1,z),
 "ln(I1/I0)"  : lambda z: np.log(iv(1,z)/iv(0,z)),
 "ln(I2/I1)"  : lambda z: np.log(iv(2,z)/iv(1,z)),
 "(I1/I0)^2"  : lambda z: (iv(1,z)/iv(0,z))**2,
}
print("which quantity has second derivative < 0 exactly when I1 I3 < I2^2 ?")
print()
hdr="%10s"%"z"+"".join("%18s"%k for k in cands)+"%16s"%("I1I3-I2^2")
print(hdr)
for z in (0.05,0.2,0.5,1.0,2.0,3.0,5.0,8.0):
    row="%10.2f"%z
    for k,f in cands.items():
        row+="%18.6e"%d2(f,z)
    A=iv(1,z)*iv(3,z)-iv(2,z)**2
    if abs(A)<1e-12: A=0.0
    row+="%16.6e"%A
    print(row)
print()
print("note:  I2/I1 - I1/I0  has the same sign as ... let me also check (I2/I1)' and (I1/I0)'")
for z in (0.1,0.5,1.0,3.0,8.0):
    print("   z=%5.2f : (I2/I1)\x27 = %+.6e   (I1/I0)\x27 = %+.6e"%(z,
        (iv(2,z+h:=1e-4)/iv(1,z+h:=1e-4)-0) if False else (iv(2,z+1e-4)/iv(1,z+1e-4)-iv(2,z-1e-4)/iv(1,z-1e-4))/2e-4,
        (iv(1,z+1e-4)/iv(0,z+1e-4)-iv(1,z-1e-4)/iv(0,z-1e-4))/2e-4))
