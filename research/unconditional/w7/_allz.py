import sympy as sp
import mpmath as mp
mp.mp.dps=40
N=sp.symbols("N",integer=True,nonnegative=True)
Dc=lambda n: -2*sp.factorial(2*n+3)/(sp.factorial(n)*sp.factorial(n+1)*sp.factorial(n+3)*sp.factorial(n+4))
print("growth of |D_N| : the ratio  |D_{N+1}/D_N| = 2(2N+5)/((N+1)(N+4)(N+5))")
for n in (5,10,20,50,100,200):
    r=sp.Rational(2)*(2*n+5)/((n+1)*(n+4)*(n+5))
    print("   N=%-4d : ratio = %.6e   |D_N| = %.6e"%(n,float(r),float(abs(Dc(n)))))
print()
print("=> ratio ~ 4/N^2  =>  |D_N| ~ 4^N/(N!)^2 * const  =>  the series  sum D_N z^{2N}  has INFINITE radius")
print()
print("numerical check that the series reproduces A(z)=I1 I3 - I2^2 and that A(z)<0 for large z:")
for z in (1,5,10,20,40):
    s=mp.mpf(0)
    for n in range(0,120):
        dn=-mp.mpf(2)*mp.factorial(2*n+3)/(mp.factorial(n)*mp.factorial(n+1)*mp.factorial(n+3)*mp.factorial(n+4))
        s+=dn*mp.mpf(z)**(2*n)
    A=mp.besseli(1,z)*mp.besseli(3,z)-mp.besseli(2,z)**2
    print("   z=%-4d : series = %+.6e   exact = %+.6e   rel.diff = %.2e"%(z,s,A,abs(s-A)/abs(A)))
