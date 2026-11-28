import numpy as np
from scipy.special import ive, iv
bs=np.linspace(0.01,30,3000)
c1=4*ive(2,bs)/bs          # = 4 e^{-beta} I2(beta) / beta  (the code's normalization)
c2=6*ive(3,bs)/bs
imax=np.argmax(c1)
print("V's character coefficient  c_1(beta) = 4 e^{-beta} I2(beta)/beta :")
print("   peak at beta = %.4f   value = %.6f"%(bs[imax],c1[imax]))
print()
print("%8s %16s %16s %16s"%("beta","c1","c2","c1/c2"))
for b in (1,2,4,8,16,32,64):
    print("%8.1f %16.8f %16.8f %16.8f"%(b,4*ive(2,b)/b,6*ive(3,b)/b,(4*ive(2,b)/b)/(6*ive(3,b)/b)))
print()
print("=> if the peak of c1 is near beta=8, the 'V is strongest at beta~8' phenomenon has an explanation")
