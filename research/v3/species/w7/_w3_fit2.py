import numpy as np
variants={
 "A_fullpool_L6Lt12":(6.2420,0.4506,5.8568,0.4001),
 "B_fullpool_L6Lt20":(6.2410,0.1170,5.8568,0.4001),
 "C_B2":(7.9990,0.6574,7.3326,0.4894),
 "D_B4":(7.9136,0.6301,6.8333,0.5404),
}
for k,(y6,s6,y8,s8) in variants.items():
    u6,u8=1/36,1/64
    a=(y8/u8-y6/u6)/(1/u8-1/u6); c=(y6-y8)/(u6-u8)
    sa=np.hypot((1/u8)/(1/u8-1/u6)*s8, (1/u6)/(1/u8-1/u6)*s6)
    sc=np.hypot(s6,s8)/(u6-u8)
    ratio=y8/y6; srat=ratio*np.hypot(s8/y8,s6/y6)
    A6,A8=y6/u6,y8/u8; sA6,sA8=s6/u6,s8/u8
    w6,w8=1/sA6**2,1/sA8**2
    cb=(A6*w6+A8*w8)/(w6+w8); scb=1/np.sqrt(w6+w8)
    chi2c=(A6-cb)**2*w6+(A8-cb)**2*w8
    wa,wb=1/s6**2,1/s8**2
    ab=(y6*wa+y8*wb)/(wa+wb); sab=1/np.sqrt(wa+wb)
    chi2a=(y6-ab)**2*wa+(y8-ab)**2*wb
    print(f"{k}: y6={y6:.3f}+-{s6:.3f} y8={y8:.3f}+-{s8:.3f}")
    print(f"  ratio={ratio:.4f}+-{srat:.4f} collapse {abs(ratio-0.5625)/srat:.2f}sig | const {abs(ratio-1)/srat:.2f}sig")
    print(f"  a={a:.3f}+-{sa:.3f} [2sig {a-2*sa:.2f},{a+2*sa:.2f}] {a/sa:.1f}sig c={c:.1f}+-{sc:.1f}")
    print(f"  1p collapse c={cb:.1f}+-{scb:.1f} chi2={chi2c:.1f} z={np.sqrt(chi2c):.2f} | 1p const a={ab:.3f}+-{sab:.3f} chi2={chi2a:.2f} z={np.sqrt(chi2a):.2f}")
    print(f"  y6/1.6449={y6/1.6449:.2f} y8/0.92528={y8/0.92528:.2f}")
