import sympy as sp
b=sp.symbols("beta",positive=True)
z=1/b
# large-z expansion of I_nu(z):  I_nu(z) ~ e^z/sqrt(2 pi z) * sum_k (-1)^k a_k(nu)/z^k
# a_k(nu) = (4nu^2-1)(4nu^2-9)...(4nu^2-(2k-1)^2)/(k! 8^k)
def series_I(nu,K=8):
    s=0
    for k in range(K):
        num=sp.Integer(1)
        for j in range(1,k+1):
            num*= (4*nu**2-(2*j-1)**2)
        s+= (-1)**k * num/(sp.factorial(k)*8**k) * z**k
    return s
for K in (5,7,9):
    ratio=sp.series(series_I(2,K)/series_I(1,K),z,0,K).removeO()
    f=-sp.series(sp.log(ratio),z,0,K).removeO()
    f=sp.expand(f)
    print("K=%d :  -ln(I2/I1) = "%K, sp.nsimplify(sp.series(f,z,0,K).removeO()))
    for k in range(1,K):
        c=sp.simplify(f.coeff(z,k))
        print("     c_%d = %s = %.10f"%(k,c,float(c)))
    print()
