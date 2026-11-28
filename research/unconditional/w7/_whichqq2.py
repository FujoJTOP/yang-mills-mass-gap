import numpy as np
from scipy.special import iv
def d2(f,z,h=1e-3): return (f(z+h)-2*f(z)+f(z-h))/h**2
cands=[("I2/I1",lambda z: iv(2,z)/iv(1,z)),
       ("ln(I2/I1)",lambda z: np.log(iv(2,z)/iv(1,z))),
       ("I1/I0",lambda z: iv(1,z)/iv(0,z)),
       ("ln(I1/I0)",lambda z: np.log(iv(1,z)/iv(0,z))),
       ("I1^2/(I0 I2)",lambda z: iv(1,z)**2/(iv(0,z)*iv(2,z)))]
print("second derivatives vs the sign of  I1 I3 - I2^2")
print("%8s"%("z")+"".join("%16s"%k for k,_ in cands)+"%16s"%("I1I3-I2^2"))
for z in (0.05,0.2,0.5,1.0,2.0,3.0,5.0,8.0):
    row="%8.2f"%z
    for _,f in cands: row+="%16.6e"%d2(f,z)
    row+="%16.6e"%(iv(1,z)*iv(3,z)-iv(2,z)**2)
    print(row)
