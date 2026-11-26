import numpy as np, math
exec(open("research/v3/w7/q12_sos.py").read().split("def Phi_d3")[0].split("print(\"=\"*100)")[0])
# FD test of Hess[2,2] at a point with xi3 != 0
beta, th = 2.5, math.pi
for xi0 in (np.array([0.9,0.4,0.75]), np.array([1.3,-0.6,0.2]), np.array([0.1,0.05,-1.1])):
    t = 1e-5
    sp = grad_hess_S(xi0+t*np.array([0,0,1.]), beta, th)[0][0]
    sm = grad_hess_S(xi0-t*np.array([0,0,1.]), beta, th)[0][0]
    s0 = grad_hess_S(xi0, beta, th)[0][0]
    fd = (sp-2*s0+sm)/t**2
    an = grad_hess_S(xi0, beta, th)[2][0,2,2]
    print(f"xi={xi0}  FD d33S={fd:+.9f}  H[2,2]={an:+.9f}  rel={abs(fd-an)/max(abs(fd),1e-12):.2e}")
