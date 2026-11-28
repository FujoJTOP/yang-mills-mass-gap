import numpy as np
from scipy.special import ive
from scipy.optimize import brentq
f=lambda b: (4*ive(2,b)/b)/(6*ive(3,b)/b)-1.0
r=brentq(f,4.0,10.0)
print("c1/c2 = 1  at  beta = %.6f"%r)
print("  (c1 = 4 e^{-b} I2/b ,  c2 = 6 e^{-b} I3/b)")
print()
print("and the single-link gap ratio u(beta) = I2/I1 there:")
from scipy.special import iv
print("  u(%.4f) = %.6f"%(r,iv(2,r)/iv(1,r)))
print()
print("c1/c2 values:")
for b in (4,5,6,7,8,10,16,64):
    print("   beta=%-4.0f : c1/c2 = %.6f"%(b,(4*ive(2,b)/b)/(6*ive(3,b)/b)))
print()
print("large-beta limit of c1/c2 (should be I2/I1 asymptotically ~ 1 - 3/(2 beta) ratio? check):")
for b in (64,256,1024):
    print("   beta=%-6.0f : c1/c2 = %.8f   u = %.8f"%(b,(4*ive(2,b)/b)/(6*ive(3,b)/b),iv(2,b)/iv(1,b)))
