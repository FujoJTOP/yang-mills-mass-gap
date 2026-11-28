import sympy as sp
z,r=sp.symbols("z r",positive=True)
# everything expressed through r = I1/I2 and the Riccati  r' = 1 - r^2 + 3r/z
# I2' = I1 - (2/z) I2 = I2 (r - 2/z)
# P = I1 I2 = r I2^2 ;  P' = I2^2 [ r' + 2 r (r - 2/z) ]
rp = 1 - r**2 + 3*r/z
Pp_over_I2sq = sp.simplify(rp + 2*r*(r-2/z))
I2p_over_I2 = r - 2/z
a  = 2*z/r + 2                 # a = 2z I1'/I1  with I1' = I2 + I1/z
ap = sp.simplify(sp.diff(a,z)) # note: d/dz of a(r(z),z) needs chain rule
# chain rule: a = 2z/r + 2 ;  da/dz = 2/r - 2z r'/r^2
ap = sp.simplify(2/r - 2*z*rp/r**2)
print("P' = I2^2 * P1'  with P1' =", sp.simplify(Pp_over_I2sq))
print("a  =", sp.simplify(a))
print("a' =", ap)
# D = P/z - (P'(1+a) - P a')/(1+a)^2 ;  factor out I2^2
D_over_I2sq = sp.simplify( r/z - (Pp_over_I2sq*(1+a) - r*ap)/(1+a)**2 )
print()
print("D / I2^2 =", sp.factor(sp.simplify(D_over_I2sq)))
print()
print("numerator and denominator separately:")
num=sp.together(sp.simplify(D_over_I2sq)); n,d=sp.fraction(num)
print("   num =", sp.factor(n))
print("   den =", sp.factor(d))
