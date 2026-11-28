from mpmath import mp, besseli, mpf
mp.dps = 50
print("=== general family: D_k = I_k I_{k+1} + I_{k+1} I_{k+2} - I_{k-1} I_{k+2} - I_k I_{k+3} ===")
print()
zs = ["1e-4","1e-2","0.1","0.5","1","2","5","10","29","100","1000"]
for k in (1,2,3,4):
    print("k=%d:"%k, flush=True)
    allneg=True
    for z in zs:
        zz=mpf(z)
        Ik=besseli(k,zz); Ik1=besseli(k+1,zz); Ik2=besseli(k+2,zz)
        Ikm1=besseli(k-1,zz); Ik3=besseli(k+3,zz)
        Dk = Ik*Ik1 + Ik1*Ik2 - Ikm1*Ik2 - Ik*Ik3
        if Dk>=0: allneg=False
        print("   z=%-8s  D_k=%s  %s"%(z, mp.nstr(Dk,10), "<0" if Dk<0 else ">=0 !!"), flush=True)
    print("   all < 0 ? %s"%allneg, flush=True)
    print()
