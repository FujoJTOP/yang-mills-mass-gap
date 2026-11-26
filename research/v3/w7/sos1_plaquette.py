# sos1_plaquette.py -- P3 "SOS on the single plaquette": is H-c ever a sum of squares?
# Convention = q11_cert.py (SU(2) single plaquette, theta in [0,pi], beta=2.5 anchor):
#   U = exp(xi.sigma/2), rho=|xi|, h=rho/2=u,  c=cos(theta/2), s=sin(theta/2)
#   S_W = beta(1 - cos h * c + sin h * s * xi3/rho),  J=4 sin^2 h/rho^2
#   slice xi3=0 (a0=<v,xi>=0, v=e3):  w(u)=sin^2 u/u,  u in (0,pi)
#   G=d3 S_W = beta s sin u/(2u),  H=d33 S_W = (beta c/2) sin u/(2u),  Phi''=1/(2u^2)-cot u/(2u)
#   G''(0) = S - D,  S=E[H+Phi''],  D=beta^2 Var(g+Phi'), g=sin u/(2u)
# Sections: [0] anchors  [A] scalar H  [B] full 3x3 Hessian (+ FD check)
#           [C] SOS tests  [D] 1D Brascamp-Lieb (the route that DOES work)  [E] verdict
import mpmath as mp
mp.mp.dps = 40
B = mp.mpf('2.5')

def w(u):   return mp.sin(u)**2/u
def gg(u):  return mp.sin(u)/(2*u)              # g(u) = sin u/(2u)
def dgg(u): return (u*mp.cos(u)-mp.sin(u))/(2*u*u)
def phipp(u): return 1/(2*u*u) - mp.cos(u)/(2*u*mp.sin(u))   # = 2/sin^2 - 1/u^2
def mom(f):
    num = mp.quad(lambda u: f(u)*w(u), [0, mp.pi])
    return num/mp.quad(w, [0, mp.pi])

Z = mp.quad(w, [0, mp.pi])
NS = mp.mpf(1)/4                                # q11: Phi''w = -(1/4)d/du[sin^2 u/u^2]  => NS=1/4
S  = NS/Z
Eg, Eg2 = mom(gg), mom(lambda u: gg(u)**2)
V  = Eg2 - Eg**2
D  = B**2*V
GP = S - D
bst = mp.sqrt(S/V)
print("[0] anchors (indep. mpmath, theta=pi)")
print("    Z          =", mp.nstr(Z, 17))
print("    S=1/(4Z)   =", mp.nstr(S, 19), "   (q11 0.20511529712307314798)")
print("    Var(g)     =", mp.nstr(V, 17))
print("    D(b=2.5)   =", mp.nstr(D, 19), "   (q11 0.069120586081084936946)")
print("    G''=S-D    =", mp.nstr(GP, 19), "   (q11 0.13599471104198821104)")
print("    beta*      =", mp.nstr(bst, 12), " (q11 4.3066078525597...)")
print("    beta=4.4:  G'' =", mp.nstr(S - mp.mpf('4.4')**2*V, 8), "  <0  while H>=0 unchanged")

# ---------------- [A] the scalar H that q11 uses (= d33 S_W, one direction only)
print("\n[A] q11 scalar H(u) = (beta c/2) g(u) = beta cos(theta/2) sin(u)/(2u)")
for th, t in [('0', mp.mpf(0)), ('pi/4', mp.pi/4), ('pi/2', mp.pi/2),
              ('3*pi/4', 3*mp.pi/4), ('pi', mp.pi)]:
    c = mp.cos(t/2)
    us = [mp.pi*mp.mpf(k)/400 for k in range(1, 400)]
    vals = [B*c*gg(u)/2 for u in us]
    print("    theta=%-7s c=%-22s min=%s  max=%s" % (th, mp.nstr(c, 8),
          mp.nstr(min(vals), 8), mp.nstr(max(vals), 8)))
c_pi = mp.cos(mp.pi/2)
print("    at theta=pi: c = cos(pi/2) =", mp.nstr(c_pi, 5), " -> set exactly 0 => H == 0 identically")
print("    range of H over u in (0,pi): [0, beta*c/4];  H >= 0 for ALL theta in [0,pi]  [T]")

# ---------------- [B] the honest operator stand-in: full 3x3 Hessian of S_W on the slice
def SW(x, y, z, beta, c, s):
    r = mp.sqrt(x*x + y*y + z*z)
    h = r/2
    return beta*(1 - mp.cos(h)*c + mp.sin(h)*s*z/r)

def Hess_FD(u, th, h=mp.mpf('1e-4')):
    c, s, beta = mp.cos(th/2), mp.sin(th/2), B
    x0 = mp.matrix([2*u, 0, 0])
    f = lambda v: SW(v[0], v[1], v[2], beta, c, s)
    M = mp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            if i == j:
                ep = mp.matrix(x0); ep[i] += h
                em = mp.matrix(x0); em[i] -= h
                M[i, j] = (f(ep) - 2*f(x0) + f(em))/h**2
            else:
                pp, pm, mpp, mmm = [mp.matrix(x0) for _ in range(4)]
                pp[i] += h; pp[j] += h; pm[i] += h; pm[j] -= h
                mpp[i] -= h; mpp[j] += h; mmm[i] -= h; mmm[j] -= h
                M[i, j] = (f(pp) - f(pm) - f(mpp) + f(mmm))/(4*h**2)
    return M

def Hess_exact(u, th):
    c, s, beta = mp.cos(th/2), mp.sin(th/2), B
    fp = mp.cos(u)/(4*u) - mp.sin(u)/(4*u*u)          # s-part: d/drho[sin h / rho]
    return mp.matrix([[beta*c*mp.cos(u)/4, 0, beta*s*fp],
                      [0, beta*c*mp.sin(u)/(4*u), 0],
                      [beta*s*fp, 0, beta*c*mp.sin(u)/(4*u)]])

u0, th0 = mp.mpf(1), mp.pi/2
A, F = Hess_exact(u0, th0), Hess_FD(u0, th0)
dev = max(abs(A[i, j] - F[i, j]) for i in range(3) for j in range(3))
print("\n[B] 3x3 Hessian on the slice: analytic vs finite differences (u=1, theta=pi/2, h=1e-4)")
print("    max|H_analytic - H_FD| =", mp.nstr(dev, 6), "  (validate the formula)")
print("    Hess_exact =", [[mp.nstr(A[i, j], 6) for j in range(3)] for i in range(3)])

def spec(u, th, add_phipp=False):
    # H has the exact block form [[a,0,h],[0,d,0],[h,0,d]] => closed-form spectrum
    M = Hess_exact(u, th)
    a, d, h = M[0, 0], M[1, 1], M[0, 2]
    if add_phipp:
        a, d = a + phipp(u), d + phipp(u)
    mid = (a + d)/2
    rad = mp.sqrt(((a - d)/2)**2 + h*h)
    return sorted([d, mid - rad, mid + rad])

print("\n    spectra over u-grid (401 pts), beta=2.5:")
print("    %-8s %-9s %-22s %-22s %s" % ("theta", "u*", "min eig Hess S_W", "min eig Hess S_W+Phi''", "frac<0"))
for th in [mp.mpf(0), mp.pi/2, mp.pi]:
    us = [mp.pi*mp.mpf(k)/400 for k in range(1, 400)]
    s1 = [(spec(u, th)[0], u) for u in us]
    s2 = [(spec(u, th, True)[0], u) for u in us]
    m1 = min(s1); m2 = min(s2)
    f1 = sum(1 for v, _ in s1 if v < 0)/len(us)
    f2 = sum(1 for v, _ in s2 if v < 0)/len(us)
    print("    %-8s %-9s %-22s %-22s %.2f / %.2f" % (
        "0" if th == 0 else ("pi/2" if th == mp.pi/2 else "pi"),
        mp.nstr(m1[1], 4), mp.nstr(m1[0], 8), mp.nstr(m2[0], 8), f1, f2))

# ---------------- [C] SOS tests
print("\n[C] SOS: can H - c = sum A_i^* A_i  (<=> H - c >= 0) ?")
us = [mp.pi*mp.mpf(k)/400 for k in range(1, 400)]
Hpi = [B*mp.cos(mp.pi/2)*gg(u)/2 for u in us]                 # = 0
for cv in ['1e-3', '1e-12']:
    cv = mp.mpf(cv)
    sup = max(Hpi) - cv
    print("    theta=pi: sup_u (H - c) =", mp.nstr(sup, 6), " for c =", mp.nstr(cv, 4), " => NEGATIVE (not SOS)")
print("    => at theta=pi, H == 0 pointwise, so H - c = -c < 0 identically: literal SOS is impossible.")
A_ = [mp.sqrt(B*mp.cos(mp.pi/2)/2*gg(u)) if B*mp.cos(mp.pi/2)*gg(u)/2 > 0 else mp.mpf(0) for u in us]
print("    trivial square exists for theta<pi: A(u)=sqrt(beta c/2 * sin u/(2u)), A^2 = H exactly;")
print("    but ANY nonneg function is a square => zero information (vacuous, measure-blind).")
print("    SUFFICIENCY TEST: at beta=4.4, G'' =", mp.nstr(S - mp.mpf('4.4')**2*V, 8),
      "< 0 while H >= 0 pointwise for every beta => H-positivity CANNOT imply G''>=c. [T by counterexample]")

# ---------------- [D] what DOES work on the single plaquette: 1D Brascamp-Lieb
print("\n[D] 1D Brascamp-Lieb on the slice measure w du (log-concavity check)")
phi2 = lambda u: 2/mp.sin(u)**2 - 1/u**2
usf = [mp.pi*mp.mpf(k)/2000 for k in range(1, 2000)]
rmin = min(phi2(u) for u in usf)
print("    phi''=2/sin^2-1/u^2 :  min over (0,pi) =", mp.nstr(rmin, 10), " at u=",
      mp.nstr(min(usf, key=phi2), 5), "  => >0 : log-concave  [T]")
Egp2 = mom(lambda u: dgg(u)**2)
bl = Egp2/rmin
print("    E[g'^2] =", mp.nstr(Egp2, 10), "   BL bound E[g'^2]/rho_min =", mp.nstr(bl, 10))
print("    true Var(g) =", mp.nstr(V, 10), "  =>  Var <= BL  [%s]" % ("holds" if V <= bl else "VIOLATED"))
for bb in ['2.5', '3.0']:
    bb = mp.mpf(bb)
    print("    beta=%-4s provable margin  S - beta^2*BL = %s" % (bb, mp.nstr(S - bb*bb*bl, 8)))
print("    analytic window from BL: beta <", mp.nstr(mp.sqrt(S/bl), 8), "  (vs exact beta* =", mp.nstr(bst, 8), ")")
print("    sharpened pointwise version: E[g'^2/phi''] =", mp.nstr(mom(lambda u: dgg(u)**2/phi2(u)), 10),
      " => Var <=", mp.nstr(mom(lambda u: dgg(u)**2/phi2(u)), 10), " [C]")
print("    NOTE: this BL is a *1D marginal* fact (u only). The 3D density is NOT log-concave:")
print("          Hess S_W + Phi'' has min eig", mp.nstr(min(spec(u, mp.pi, True)[0] for u in usf), 8),
      "at theta=pi => 3D/4D BL prerequisite FALSE.")

print("\n[E] verdict: SOS fails on the single plaquette in every reading; see SOS1.md")
