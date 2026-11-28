import sympy as sp
b=sp.symbols('beta',positive=True)
# exact small-beta expansions of I_n
def Iser(n,order=10):
    return sp.series(sp.besseli(n,b),b,0,order).removeO()
for J in (1,2,3):
    cJ=sp.simplify(2*(J+1)*Iser(J+1,12)/b)
    c2=sp.simplify(sp.diff(cJ,b,2).subs(b,0))
    c1=sp.simplify(sp.diff(cJ,b,1).subs(b,0))
    print("J=%d :  c_J = %s"% (J, sp.expand(cJ)))
    print("        c_J'(0) = %s ,  c_J''(0) = %s"%(c1,c2))
