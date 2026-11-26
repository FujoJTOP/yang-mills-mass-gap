import numpy as np, sys
sys.path.insert(0,r"research\v3\w7")
import o05_var_scaling as o
print("[P3-fix] lambda = Poincare const of the 1D slice measure: a(u) = (p(u))^2, p = w/Z")
print("  theta      beta     lambda    lambda/beta   (1/lam)E|gradF|^2  Var      bound/Var  Var/beta")
for th in [np.pi/2, 3*np.pi/4, np.pi]:
    for b in [1.0,2.5,6.0,10.0,50.0]:
        d = o.slice_stats(b, th, n=200001)
        lam = o.poincare_1d(d["r"], d["w"], M=1000)
        bd = d["grad2"]/lam
        print("  %6.4f %7.2f %11.5f %10.4f %16.6f %10.6f %9.4f %10.6f"%(th,b,lam,lam/b,bd,d["var"],bd/d["var"],d["var"]/b))
