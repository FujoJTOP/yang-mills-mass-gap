import mpmath as mp
mp.mp.dps=30
I=mp.besseli
def Phi(b): return b*mp.e**b/(2*I(1,b))-1
def root(lam): return mp.findroot(lambda b: lam*Phi(b)-1, mp.mpf("0.02"))
ED=20*mp.e
Lt=lambda D: mp.mpf(D)**D/(mp.mpf(D-1))**(D-1)
print("Lt(14)=",mp.nstr(Lt(14),10)," Lt(20)=",mp.nstr(Lt(20),10))
for lam,nm in [(Lt(14),"Lt(14)=36.69"),(mp.mpf("13.90"),"lam_*>=13.90"),(mp.mpf("25.2"),"lam_*=25.2"),(mp.mpf("47.6270568100"),"47.627"),(ED,"eD")]:
    r=root(lam); print(f"  lam={nm:16s} beta_max={mp.nstr(r,10)}  x(1/48)={mp.nstr(r*48,7)}  x0.625={mp.nstr(mp.mpf('0.625')/r,7)}")
a=[1,20,458,11132,280175,7218768]
print("a_n:",a)
print("ratios:",[mp.nstr(mp.mpf(a[i+1])/a[i],6) for i in range(5)])
print("nth roots:",[mp.nstr(mp.mpf(a[i])**(mp.mpf(1)/(i+1)),6) for i in range(6)])
print("1/Phi0(1/48)=",mp.nstr(1/Phi(mp.mpf(1)/48),12),"  1/Phi0(0.0208)=",mp.nstr(1/Phi(mp.mpf('0.0208')),12))
print("1/Phi0(0.625)=",mp.nstr(1/Phi(mp.mpf('0.625')),12))
print("eD/(1/Phi0(1/48))=",mp.nstr(ED/(1/Phi(mp.mpf(1)/48)),10))
print("u(0.0183)=",mp.nstr((I(1,mp.mpf('0.0183'))-I(3,mp.mpf('0.0183')))/(2*I(1,mp.mpf('0.0183'))/mp.mpf('0.0183')),12))
print("2.446918311/0.018268566=",mp.nstr(mp.mpf("2.446918311")/mp.mpf("0.018268566"),8))
print("0.0183 in 2D: m=3g2/8=0.375, sigma=0.375 -> need a=",0.375/(4*0.375))
