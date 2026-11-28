from mpmath import mp, besseli, mpf
mp.dps = 50
print("=== PROOF of the family: use the Bessel recurrence ===")
print()
print("Recurrence:  I_{nu-1}(z) - I_{nu+1}(z) = (2 nu / z) I_nu(z)")
print()
print("D_k = I_k I_{k+1} + I_{k+1} I_{k+2} - I_{k-1} I_{k+2} - I_k I_{k+3}")
print("    = I_k (I_{k+1} - I_{k+3}) + I_{k+2} (I_{k+1} - I_{k-1})")
print()
print("apply recurrence:")
print("  I_{k+1} - I_{k+3} = (2(k+2)/z) I_{k+2}")
print("  I_{k+1} - I_{k-1} = -(2k/z) I_k")
print()
print("D_k = I_k * (2(k+2)/z) I_{k+2} + I_{k+2} * (-(2k/z)) I_k")
print("    = (2/z) [ (k+2) - k ] I_k I_{k+2}")
print("    = (4/z) I_k I_{k+2}")
print()
print("==> D_k(z) = (4/z) I_k(z) I_{k+2}(z)  > 0  TRIVIALLY (since I_nu>0 for z>0)")
print()
print("and  D'(z) = (1/2) D_1 = (2/z) I_1(z) I_3(z) > 0")
print()
# numerical verification
print("numerical verification of D'_exact = (2/z) I_1 I_3  vs  (1/2)D_1:")
for z in ("1e-4","1e-2","0.1","1","2","10","100"):
    zz=mpf(z)
    I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz); I4=besseli(4,zz); I0=besseli(0,zz)
    Dp_formula=mpf(2)/zz*I1*I3
    Dp_expanded=(I1*I2+I2*I3-I0*I3-I1*I4)/2
    print("  z=%-8s  (2/z)I_1I_3 = %s   (1/2)D_1 = %s%s"%(z, mp.nstr(Dp_formula,15), mp.nstr(Dp_expanded,15),
          "  MATCH" if abs(Dp_formula-Dp_expanded)<abs(Dp_expanded)*mpf("1e-30") else "  DIFF"), flush=True)
print()
# and the family
print("numerical verification of D_k = (4/z) I_k I_{k+2}:")
for k in (1,2,3,4):
    for z in ("0.1","1","10"):
        zz=mpf(z)
        Ik=besseli(k,zz); Ik1=besseli(k+1,zz); Ik2=besseli(k+2,zz)
        Ikm1=besseli(k-1,zz); Ik3=besseli(k+3,zz)
        Dk=Ikm1*0
        Dk = Ik*Ik1 + Ik1*Ik2 - Ikm1*Ik2 - Ik*Ik3
        pred=mpf(4)/zz*Ik*Ik2
        print("  k=%d z=%-5s  D_k=%-22s (4/z)I_kI_{k+2}=%-22s %s"%(k,z,mp.nstr(Dk,12),mp.nstr(pred,12),
              "MATCH" if abs(Dk-pred)<abs(pred)*mpf("1e-30") else "DIFF"), flush=True)
