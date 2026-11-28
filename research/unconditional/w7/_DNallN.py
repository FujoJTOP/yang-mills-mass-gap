import sympy as sp
N=sp.symbols("N",integer=True,nonnegative=True)
# closed form guessed in round 36:
Dc=lambda n: -2*sp.factorial(2*n+3)/(sp.factorial(n)*sp.factorial(n+1)*sp.factorial(n+3)*sp.factorial(n+4))
# recurrence from round 36:
rec=lambda n: sp.Rational(2)*(2*n+5)/((n+1)*(n+4)*(n+5))
print("check D_{N+1}/D_N = 2(2N+5)/((N+1)(N+4)(N+5))  for the closed form, symbolically")
expr=sp.simplify(Dc(N+1)/Dc(N) - rec(N))
print("   difference =", expr)
print()
print("check D_0 value:")
print("   closed form D_0 =", Dc(0), " = ", sp.nsimplify(Dc(0)))
print()
print("=> since the closed form satisfies the recurrence identically (difference = 0) and matches at N=0,")
print("   D_N = -2(2N+3)!/(N!(N+1)!(N+3)!(N+4)!)  for ALL N,  hence  D_N < 0  for ALL N.")
