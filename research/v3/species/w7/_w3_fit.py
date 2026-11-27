import numpy as np
# full-pool t=1 numbers from w3_y_single_coord.py (JK, 20 blocks)
for t,(y6,s6,y8,s8) in {1:(6.2489,0.451,5.8700,0.400),
                        2:(5.9952,1.640,5.0170,1.640),
                        3:(1.4222,1.066,0.9683,1.194)}.items():
    u6,u8=1/36,1/64
    # 2-param a + c/L^2, exactly solved by 2 points
    a=(y8/u8 - y6/u6)/(1/u8-1/u6); c=(y6-y8)/(u6-u8)
    # jacobian
    da_dy6=-1/u6/(1/u8-1/u6); da_dy8=1/u8/(1/u8-1/u6)
    dc_dy6=1/(u6-u8); dc_dy8=-1/(u6-u8)
    sa=np.hypot(da_dy6*s6, da_dy8*s8); sc=np.hypot(dc_dy6*s6, dc_dy8*s8)
    corr=(da_dy6*dc_dy6*s6**2+da_dy8*dc_dy8*s8**2)/(sa*sc)
    # 1-param collapse y=c/L^2
    A6,A8=y6/u6,y8/u8; sA6,sA8=s6/u6,s8/u8
    w6,w8=1/sA6**2,1/sA8**2
    cbar=(A6*w6+A8*w8)/(w6+w8); scbar=1/np.sqrt(w6+w8)
    chi2c=(A6-cbar)**2*w6+(A8-cbar)**2*w8
    # 1-param constant
    wa,wb=1/s6**2,1/s8**2
    abar=(y6*wa+y8*wb)/(wa+wb); sabar=1/np.sqrt(wa+wb)
    chi2a=(y6-abar)**2*wa+(y8-abar)**2*wb
    ratio=y8/y6; srat=ratio*np.hypot(s8/y8,s6/y6)
    print(f"t={t}: y6={y6:.3f}+-{s6:.3f} y8={y8:.3f}+-{s8:.3f} ratio={ratio:.3f}+-{srat:.3f}")
    print(f"   a={a:.3f}+-{sa:.3f} (2sig {a-2*sa:.2f}..{a+2*sa:.2f}) c={c:.1f}+-{sc:.1f} corr={corr:.3f}")
    print(f"   collapse: cbar={cbar:.1f}+-{scbar:.1f} chi2={chi2c:.2f}/1dof  const: abar={abar:.3f}+-{sabar:.3f} chi2={chi2a:.2f}/1dof")
    print(f"   sigma(collapse excluded) = {(ratio-0.5625)/srat:.2f} ; sigma(const) = {(ratio-1.0)/srat:.2f}")
    print(f"   a/sigma_a = {a/sa:.2f}")
print()
print("template 6pi^2/L^2:", [round(6*np.pi**2/L**2,4) for L in (6,8,10,12)])
print("collapse pred y8/y6 =", (6/8)**2)
