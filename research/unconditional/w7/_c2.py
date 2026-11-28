import numpy as np
r0=0.507195100047
bs=np.array([1e-5,1e-4,1e-3,1e-2,5e-2])
rr=np.array([0.507197636014,0.507220458929,0.507448610197,0.509722254818,0.519645859259])
y=(rr-r0)/r0
print("fit  (r-r0)/r0 = c1*bs + c2*bs^2   (c1 should be 0.5 exactly)")
for K in (1,2):
    A=np.vstack([bs**k for k in range(1,K+1)]).T
    coef,_,_,_=np.linalg.lstsq(A,y,rcond=None)
    pred=A@coef
    print("   K=%d : c = %s   max|resid| = %.3e"%(K," ".join("%.10f"%c for c in coef),np.max(np.abs(pred-y))))
print()
print("using only the smallest three points (cleanest second derivative):")
A=np.vstack([bs**1,bs**2]).T
coef,_,_,_=np.linalg.lstsq(A[:3],y[:3],rcond=None)
print("   c1 = %.10f   c2 = %.10f"%(coef[0],coef[1]))
print("   =>  r = r0 (1 + c1*bs + c2*bs^2),  peak at bs = -c1/(2c2) = %.4f"%(-coef[0]/(2*coef[1])))
print()
print("measured peak (round 111): beta_s ~ 2.09 ;  V has a NEGATIVE second-order effect (c2<0)")
