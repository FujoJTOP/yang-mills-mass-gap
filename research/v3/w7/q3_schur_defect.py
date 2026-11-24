#!/usr/bin/env python3
r"""
q3_schur_defect.py  --  Q3 witness: integrate out the unstable mode exactly and bound the
Schur complement (H1 = chi = Var(xi_{k->0}) <= 1/lambda  <=>  G''_{k=0} >= lambda > 0).

Deterministic, py -3.12 -I, numpy only (scipy.special.sici for the strip Haar marginal).
Exact Gauss-Legendre quadrature / exact linear algebra; no Monte Carlo.

PARTS
 (1) HAYNSWORTH INERTIA LEMMA [T]+[V].  For a symmetric block H=[[a,c^T],[c,K]]:
        det H = det K * Sigma ,  Sigma = a - c^T K^{-1} c ,
        In(H) = In(K) (+) In(Sigma)   (so  n_-(H) = n_-(K) + n_-(Sigma)).
     => Sigma>0  <=>  the negative modes of H are ALL in the fibre K.
        K>0 & Sigma<0 => the negativity is irreducible (a coarse-direction effect); no
        fibre-mode excision can fix it.  K indefinite & n_-(H)=n_-(K) => Sigma>0 (excision works).
 (2) ONE-PLAQUETTE SU(2) (S2/G4 exact model) [T]+[V].
        V(xi)=h(rho)+g(rho)u, u=xi3/rho, h=beta(1-cos(rho/2)cos(th/2))+Phi(rho),
        g=beta sin(rho/2)sin(th/2), Phi=2log(rho/2)-2log sin(rho/2); coarse v=e3.
        Per-r Hessian blocks; measure-averaged S=E[alpha], D=Var(f), G''=S-D (=A-B beta^2 at th=pi),
        Sigma(r), E[Sigma], eta=D/S, and the inertia of K and H.
        theta=0  : K indefinite (Savvidy radial), c_r=0 => Sigma=alpha>0  => excision restores positivity.
        theta=pi : K>0 everywhere, Sigma<0 at beta=6 (collective mode) => excision CANNOT help.
        NOTE: rho(a) is tilted at th=pi (g u term), so G''<0 is a log-convex shoulder, NOT a
              symmetric double well; the barrier reading is a property of the (tilted) marginal.
 (3) 1xL STRIP [T]+[V].  Exact constant-mode marginal curvature G''(0) via the periodic
        Cartan-chain transfer matrix Z(k)=Tr[M(k)^L] (P2/P4), beta=2.5, L=2..8.  S_W is EXACTLY
        flat on the constant mode (shift symmetry), so its Hessian block is Haar-only; the Gaussian
        flat-point Schur complement is Phi''(0)=1/6.  Both positive & L-uniform: no Savvidy
        instability on the strip constant mode.
 (4) RANK-ONE DEFECT LEMMA [T]+[V].  Exact identity G''=S-D, chi=1/(S(1-eta)), eta=D/S.
        Schur+Brascamp-Lieb (K>0): G'' >= E[Sigma]=:rho_rest , hence chi <= 1/rho_rest
        <= (1+eta)/rho_rest whenever rho_rest>0.  Verified on the one-plaquette numbers and a
        Gaussian+quartic toy (with BL slack shown).  The binding hypothesis is rho_rest>0
        (= E<K^{-1}c,c> < S), NOT eta<1.
 (5) ACID TEST [T].  Free massless field V=1/2 k^2 x^2: G''=rho_rest=k^2 -> 0 as k->0; chi->inf.
        No false positive floor; the positivity hypothesis rho_rest>0 correctly fails.

Run:  py -3.12 -I q3_schur_defect.py
Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)
np.seterr(invalid='ignore', divide='ignore')

# =====================================================================================
# shared one-plaquette SU(2) model (S2/G4 flat-chart, boundary B=exp(i th sig3/2))
# =====================================================================================
def Phi(r):        return 2.0*np.log(r/2.0) - 2.0*np.log(np.sin(r/2.0))
def Phip(r):
    x = r/2.0
    return 1.0/x - 1.0/np.tan(x)
def Phipp(r):
    x = r/2.0
    return 0.5*(-1.0/x**2 + 1.0/np.sin(x)**2)
def hfun(r,b,th):  return b*(1.0-np.cos(r/2.0)*np.cos(th/2.0)) + Phi(r)
def hp(r,b,th):    return b*0.5*np.sin(r/2.0)*np.cos(th/2.0) + Phip(r)
def hpp(r,b,th):   return b*0.25*np.cos(r/2.0)*np.cos(th/2.0) + Phipp(r)
def gfun(r,b,th):  return b*np.sin(r/2.0)*np.sin(th/2.0)
def gp(r,b,th):    return b*0.5*np.cos(r/2.0)*np.sin(th/2.0)

def blocks(r,b,th):
    """Hessian entries at slice xi3=0, transverse radius r, basis (radial rhat, tangential, e3):
       H = [[K_rad, 0, c_r],[0, K_tang, 0],[c_r, 0, alpha]].
       alpha = d^2_e3 V ; c_r = e3-radial coupling ; K_rad=h'', K_tang=h'/r ; c_tang=0."""
    alpha = hp(r,b,th)/r
    cr    = gp(r,b,th)/r - gfun(r,b,th)/r**2
    Krad  = hpp(r,b,th)
    Ktang = hp(r,b,th)/r
    return alpha, cr, Krad, Ktang

def H3(alpha,cr,Krad,Ktang):
    return np.array([[Krad,0.0,cr],[0.0,Ktang,0.0],[cr,0.0,alpha]])

# =====================================================================================
def part1_haynsworth():
    print("="*92)
    print("(1) HAYNSWORTH INERTIA LEMMA  [T]+[V]:  In(H)=In(K)+In(Sigma),  det H=det K*Sigma")
    print("="*92)
    rng = np.random.default_rng(20261124)
    print("  (a) K>0 & Sigma<0  => n_-(H)=1  (negativity irreducible; fibre-mode excision cannot help)")
    print("  (b) K indefinite, n_-(H)=n_-(K) => Sigma>0  (fibre-mode excision restores positivity)")
    ok_det = ok_in = True
    for trial in range(6):
        m = 1; nf = 2; n = m+nf
        A = rng.standard_normal((n,n)); H = A+A.T
        if trial % 2 == 0:                      # make K PD
            Q,_ = np.linalg.qr(rng.standard_normal((nf,nf)))
            K = Q@np.diag([0.5+rng.random(), 0.5+rng.random()])@Q.T
            K = 0.5*(K+K.T)
        else:                                   # make K indefinite (one negative)
            Q,_ = np.linalg.qr(rng.standard_normal((nf,nf)))
            K = Q@np.diag([-0.6-rng.random(), 0.5+rng.random()])@Q.T
            K = 0.5*(K+K.T)
        c = rng.standard_normal(nf); a = rng.standard_normal()
        Hb = np.block([[np.array([[a]]), c.reshape(1,nf)],[c.reshape(nf,1), K]])
        Sigma = a - c @ np.linalg.solve(K, c)
        evH = np.linalg.eigvalsh(Hb); evK = np.linalg.eigvalsh(K)
        nH = int((evH < -1e-10).sum()); nK = int((evK < -1e-10).sum())
        nS = int((Sigma < -1e-10).sum())
        ok_det &= abs(np.linalg.det(Hb) - np.linalg.det(K)*Sigma) < 1e-9
        ok_in  &= (nH == nK + nS)
    print(f"     det H = det K * Sigma  (all trials): {ok_det}")
    print(f"     n_-(H) = n_-(K) + n_-(Sigma)  (all trials): {ok_in}")
    # explicit demonstrations
    K1 = np.diag([1.0,1.0]); c1 = np.array([2.0,0.0]); a1 = 0.5   # Sigma = 0.5-4 = -3.5<0
    H1 = np.block([[np.array([[a1]]), c1.reshape(1,2)],[c1.reshape(2,1),K1]])
    print(f"     demo K={np.diag(K1)} PD, c={c1}, a={a1}: Sigma={a1-c1@np.linalg.solve(K1,c1):+.3f}, "
          f"evals(H)={np.linalg.eigvalsh(H1)}  -> one negative in H, none in K  (IRREDUCIBLE)")
    K2 = np.diag([-1.0,1.0]); c2 = np.array([0.0,0.3]); a2 = 2.0
    H2 = np.block([[np.array([[a2]]), c2.reshape(1,2)],[c2.reshape(2,1),K2]])
    print(f"     demo K={np.diag(K2)} indefinite, c={c2}, a={a2}: Sigma={a2-c2@np.linalg.solve(K2,c2):+.3f}>0, "
          f"evals(H)={np.linalg.eigvalsh(H2)}  -> negatives of H sit in K  (EXCISABLE)")
    print()

# =====================================================================================
def moments(b,th,n=2200):
    """measure-averaged S=E[alpha], D=Var(f), G''=S-D, E[Sigma]; weight e^{-h} r dr (G4)."""
    xg,wg = np.polynomial.legendre.leggauss(n)
    r = np.pi*(xg+1.0); wr = np.pi*wg
    h = hfun(r,b,th); w = np.exp(-h)*r*wr; Z = w.sum()
    f = gfun(r,b,th)/r                 # f = d_{e3}V at slice
    alpha = hp(r,b,th)/r               # alpha = d^2_{e3}V at slice
    S = (w*alpha).sum()/Z
    D = (w*f*f).sum()/Z - ((w*f).sum()/Z)**2
    cr = gp(r,b,th)/r - gfun(r,b,th)/r**2
    Krad = hpp(r,b,th)
    Sig = alpha - np.where(np.abs(Krad)>1e-9, cr**2/Krad, 0.0)
    ESig = (w*Sig).sum()/Z
    return S, D, S-D, ESig

def part2_one_plaquette():
    print("="*92)
    print("(2) ONE-PLAQUETTE SU(2):  Hessian inertia, Schur complement Sigma, exact G''  [T]+[V]")
    print("="*92)
    print("  Exact marginal curvature (th=pi, v=e3):  G'' = A - B beta^2,  A=0.20511530, B=0.01105929")
    print("  (S2/G4; negative for beta>4.30661).  Below: S=E[alpha], D=Var(f), G''=S-D, E[Sigma], eta=D/S.")
    print(f"  {'beta':>5} {'th':>6} {'S':>9} {'D':>9} {'Gpp':>9} {'E[Sigma]':>9} {'eta=D/S':>8} {'A-Bb^2':>9}")
    for b in (2.5,6.0):
        for th in (0.0,np.pi):
            S,D,G,ES = moments(b,th)
            AB = 0.20511530 - 0.01105929*b*b
            print(f"  {b:5.1f} {th:6.3f} {S:9.5f} {D:9.5f} {G:9.5f} {ES:9.5f} {D/S:8.4f} {AB:9.5f}")
    print("  => th=0 : D=0 (no boundary force), G''=S=alpha>0, E[Sigma]=S>0.")
    print("     th=pi: G''=A-B b^2 reproduces the paper exactly; E[Sigma]<0 at beta=6; eta=D/S>1.")

    print("\n  Hessian inertia along the slice (basis (radial,tangential,e3)); c_r=0 <=> tangential ex-axis:")
    print(f"  {'beta':>5} {'th':>6} {'r':>5} {'K_rad':>8} {'K_tang':>8} {'alpha':>8} {'c_r':>7}"
          f" {'n-(K)':>6} {'n-(H)':>6} {'Sigma':>9}")
    for b in (2.5,6.0):
        for th in (0.0,np.pi):
            for r in (1.0,2.0,3.0):
                alpha,cr,Krad,Ktang = blocks(r,b,th)
                Hb = H3(alpha,cr,Krad,Ktang)
                nH = int((np.linalg.eigvalsh(Hb) < -1e-10).sum())
                nK = int((np.linalg.eigvalsh(np.diag([Krad,Ktang])) < -1e-10).sum())
                Sig = alpha - (cr**2/Krad if abs(Krad)>1e-9 else 0.0)
                print(f"  {b:5.1f} {th:6.3f} {r:5.1f} {Krad:8.3f} {Ktang:8.3f} {alpha:8.3f} {cr:7.3f}"
                      f" {nK:6d} {nH:6d} {Sig:9.4f}")
    print("  => th=0 : c_r=0, Sigma=alpha>0; H indefinite only where K is (n-(H)=n-(K)): EXCISABLE.")
    print("     th=pi: K>0 (n-(K)=0) yet H indefinite & Sigma<0 at beta=6 (n-(H)=n-(Sigma)=1): a")
    print("            COLLECTIVE (radial,e3) mode; there is NO fibre negative mode to excise.")

    # scan: is K ever indefinite at th=pi?  and where Sigma<0?
    for b,th in ((6.0,np.pi),(6.0,0.0)):
        rs = np.linspace(0.05, 2*np.pi-0.05, 600)
        nK = nS = nH = 0
        for r in rs:
            alpha,cr,Krad,Ktang = blocks(r,b,th)
            Sig = alpha - (cr**2/Krad if abs(Krad)>1e-9 else 0.0)
            Ksig = np.array([Krad,Ktang]); Hb = H3(alpha,cr,Krad,Ktang)
            if (Ksig<0).any(): nK+=1
            if Sig<0: nS+=1
            if (np.linalg.eigvalsh(Hb)<0).any(): nH+=1
        print(f"     beta={b}, th={th:.3f}: over {len(rs)} r-values -> K-indef {nK}, Sigma<0 {nS}, H-indef {nH}")
    print()

# =====================================================================================
def part3_strip():
    print("="*92)
    print("(3) 1xL STRIP (periodic Cartan chain):  exact constant-mode Schur-complement curvature [T]+[V]")
    print("="*92)
    from scipy.special import sici
    EULER = 0.5772156649015328606; TWO_PI = 2.0*np.pi
    def Cin(x):
        x=np.asarray(x,float); _,ci=sici(x)
        return np.where(x>0, EULER+np.log(np.where(x>0,x,1.0))-ci, 0.0)
    def rho1(z):
        z=np.asarray(z,float); az=np.abs(z)
        return np.where(az<=TWO_PI, np.pi*(Cin(TWO_PI)-Cin(np.minimum(az,TWO_PI))), 0.0)
    class PeriodicChain:
        def __init__(self,L,beta,nz=120):
            self.L=L; xg,wg=np.polynomial.legendre.leggauss(nz)
            self.z=TWO_PI*xg; self.w=TWO_PI*wg; self.r=rho1(self.z); self.sw=np.sqrt(self.r*self.w)
            d=self.z[:,None]-self.z[None,:]
            self.M0=np.exp(-beta*(1.0-np.cos(0.5*d)))*self.sw[:,None]*self.sw[None,:]
            self.half=(self.z[:,None]+self.z[None,:])/(2.0*np.sqrt(L))
        def Z(self,k):
            k=np.atleast_1d(np.asarray(k,float)); out=np.empty(k.shape,complex)
            for i,kk in enumerate(k):
                out[i]=np.trace(np.linalg.matrix_power(self.M0*np.exp(-1j*kk*self.half), self.L))
            return out
    def gpp0(Zfun,kw=30.,K=50.,dk=0.05,delta=0.08,power=6):
        kj=np.arange(-int(round(K/dk)),int(round(K/dk))+1)*dk
        Z=np.asarray(Zfun(kj),complex)*np.exp(-(kj/kw)**power)
        def rho(a): return (dk/(2*np.pi))*np.real(np.exp(1j*np.outer(a,kj))@Z)
        xs=delta*np.array([-2,-1,0,1,2]); rr=rho(xs); lg=np.log(np.maximum(rr,1e-300))
        return -(-lg[0]+16*lg[1]-30*lg[2]+16*lg[3]-lg[4])/(12*delta*delta)
    print("  S_W is EXACTLY flat on the constant mode (uniform shift) => its Hessian block on the")
    print("  constant direction is Haar-only : alpha_const = E[Phi''], coupling c_const=0 at flat pt.")
    print("  Gaussian flat-point Schur complement  Sigma0 = Phi''(0) = 1/6 = %.6f > 0 (beta-indep)." % (1/6))
    print(f"  {'L':>3} {'G''(0) exact':>14} {'>1/6?':>6}")
    vals={}
    for L in range(2,9):
        g=gpp0(PeriodicChain(L,2.5).Z); vals[L]=g
        print(f"  {L:3d} {g:14.5f} {('yes' if g>1/6 else 'NO'):>6}")
    # fit c + A q^L
    try:
        from scipy.optimize import curve_fit
        Ls=np.array(list(vals.keys()),float); y=np.array([vals[L] for L in vals])
        f=lambda L,c,A,q: c+A*q**L
        p,_=curve_fit(f,Ls,y,p0=[0.22,2.0,0.6],maxfev=400000)
        print(f"  fit c+A q^L: c={p[0]:.5f} (saturating, L-independent)  q={p[2]:.3f}")
    except Exception:
        pass
    print("  => strip constant-mode Schur curvature is POSITIVE, saturating (~0.22), L-uniform.")
    print("     No Savvidy instability on the constant mode: S_W flat there; negatives live in the")
    print("     gapped (4pi^2/L^2) non-constant fibre.  Exact exceeds the Gaussian 1/6 (measure gain).")
    print()

# =====================================================================================
def part4_defect_lemma():
    print("="*92)
    print("(4) RANK-ONE DEFECT LEMMA  [T]+[V]:  chi <= (1+eta)/rho_rest  when rho_rest>0")
    print("="*92)
    print("  Exact variance identity:  G'' = S - D,  S=E[alpha]=bare stiffness, D=Var(f)=force variance.")
    print("  Gaussian-LOCAL: 1/G'' = 1/(S(1-eta)), eta=D/S.  (In general Var(a) != 1/G'';")
    print("  the rigorous ceiling is BL: Var(a) <= E[1/G''] <= 1/rho_* , rho_* = essinf_a0 rho_rest.)")
    print("  Schur + Brascamp-Lieb (fibre K>0 log-concave):  G''(a0) >= E[Sigma|a0] =: rho_rest(a0),")
    print("  because BL gives D=Var(f) <= E<K^{-1}c,c> => rho_rest = S - E<K^{-1}c,c> <= S - D = G''.")
    print("  Hence chi <= E[1/G''] <= 1/rho_* <= (1+eta)/rho_*  whenever rho_*>0.  [T]")
    print("  SHARP form chi <= 1/rho_*; (1+eta) is cosmetic slack.  Binding hypothesis: rho_*>0.\n")
    print("  one-plaquette check (th=pi):")
    print(f"  {'beta':>5} {'S':>8} {'D':>8} {'G''':>8} {'E[Sigma]=rho_rest':>18} {'1/rho_rest':>11} {'(1+eta)/rho':>12} {'1/G(loc)':>9}")
    for b in (2.5,6.0):
        S,D,G,ES = moments(b,np.pi)
        eta=D/S
        one = 1.0/ES if abs(ES)>1e-9 else float('nan')
        loc = 1.0/G if abs(G)>1e-9 else float('nan')
        print(f"  {b:5.1f} {S:8.4f} {D:8.4f} {G:8.4f} {ES:18.5f} {one:11.3f} {(1+eta)*one:12.3f} {loc:9.3f}")
    print("  => beta=2.5: rho_rest=+0.055>0 -> Var(a)<=E[1/G'']<=1/rho_rest=18.0 (and <=(1+eta)/rho).")
    print("     beta=6.0: rho_rest=-0.657<0 -> HYPOTHESIS FAILS; lemma does not apply; G''<0.")

    print("\n  Gaussian+quartic toy  V(x,y)=1/2 A x^2 + b x y + 1/2 K y^2 + (lam/4) y^4 :")
    print("  alpha=A (const), c=d_x d_y V=b, fibre Hessian K(y)=K+3 lam y^2.  BL says")
    print("  G''=A-Var(f) >= E[Sigma] = A - b^2 E[1/K(y)] =: rho_rest.  The gap is the BL slack.")
    print("  Verify G''>=rho_rest and that small eta=D/S does NOT make rho_rest>0 (the slack")
    print("  b^2(E[1/K(y)]-Var(y)) can dominate).")
    def toy(A,b,K,lam,x0=0.0,n=120001):
        # fibre weight at x0:  w(y)=exp(-1/2 K y^2 - 1/4 lam y^4)  (the b x0 y term is a tilt; absorb)
        s = np.sqrt(1.0/lam) if lam>0 else 1.0
        L = max(8.0, 10.0*s, 10.0/np.sqrt(K) if K>0 else 8.0)
        y = np.linspace(-L, L, n); dy = y[1]-y[0]
        def fibre_w(x):
            V = 0.5*K*y*y + 0.25*lam*y**4 + b*x*y
            e = -V; e -= e.max(); return np.exp(e)
        def lrho(x):
            V = 0.5*A*x*x + b*x*y + 0.5*K*y*y + 0.25*lam*y**4   # NO max-shift (keep x-dependence)
            return np.log(np.trapezoid(np.exp(-V), y, dx=dy))
        h=0.02
        ls=[lrho(x0+k*h) for k in range(-2,3)]
        G = -(-ls[0]+16*ls[1]-30*ls[2]+16*ls[3]-ls[4])/(12*h*h)
        # conditional fibre moments at x0
        w=fibre_w(x0); Z=np.trapezoid(w,y,dx=dy)
        m1=np.trapezoid(w*y,y,dx=dy)/Z; y2=np.trapezoid(w*y*y,y,dx=dy)/Z
        Vary = y2 - m1*m1
        Ky = K + 3.0*lam*y*y
        invKbar = np.trapezoid(w/Ky, y, dx=dy)/Z          # E[1/K(y)]
        S = A; D = b*b*Vary                                # S=E[alpha], D=Var(f)
        rho_rest = A - b*b*invKbar                         # E[Sigma] (exact mean Schur)
        eta = D/S if abs(S)>1e-12 else np.inf
        return G,rho_rest,S,D,eta,Vary,invKbar
    print(f"  {'A':>4} {'b':>4} {'K':>5} {'lam':>4} {'G''':>8} {'rho_rest':>9} {'S':>5} {'D':>6} {'eta':>6} {'G>=rho?':>8}")
    for (A,b,K,lam) in [(1.0,0.3,1.0,0.1),(1.0,0.8,1.0,0.1),(0.55,0.9,0.5,1.0),
                         (0.90,1.0,0.15,8.0),(0.5,1.0,0.5,0.3),(0.4,1.2,0.4,0.2),(1.0,1.5,1.0,0.1),(1.0,0.9,0.5,1.0)]:
        G,rr,S,D,eta,Vary,iK = toy(A,b,K,lam)
        ok = (G >= rr-2e-3)
        star = " <-- G''>0, rho_rest<0" if (G>0.05 and rr<0) else ""
        print(f"  {A:4.2f} {b:4.1f} {K:5.2f} {lam:4.1f} {G:8.4f} {rr:9.4f} {S:5.2f} {D:6.3f} {eta:6.3f} {str(ok):>8}{star}")
    print("  => G'' >= rho_rest holds in EVERY row (Brascamp-Lieb, fibre K(y)>0); rho_rest is the")
    print("     binding object.  Row (A=0.90,K=0.15,lam=8) has eta=D/S=0.26 (a SMALL defect) yet")
    print("     G''=+0.67>0 and rho_rest=-0.12<0: the BL slack b^2(E[1/K(y)]-Var(y))>0 dominates,")
    print("     so the SCHUR bound fails although G''>0 -- a small force-variance defect does NOT")
    print("     certify rho_rest>0.  The fibre gap (K(y) bounded below) is what controls the slack.")
    print()

# =====================================================================================
def part5_acid():
    print("="*92)
    print("(5) ACID TEST: free massless field must NOT acquire a positive floor  [T]")
    print("="*92)
    print("  Free massless V=1/2 k^2 x^2 (Gaussian, fibre arbitrary/decoupled): alpha=k^2, c=0,")
    print("  K=k^2, Sigma=alpha - c^2/K = k^2, D=Var(f)=0, G''=S-D=k^2.  All -> 0 as k->0.")
    print(f"  {'k':>8} {'G''=rho_rest':>14} {'chi=1/G':>12} {'(1+eta)/rho':>12}")
    for k in (1.0,0.3,0.1,0.02,0.005):
        A=k*k; G=A; chi=1.0/G
        print(f"  {k:8.3f} {G:14.6f} {chi:12.2f} {(1.0+0.0)/A:12.2f}")
    print("  => G''=rho_rest=k^2 -> 0: no manufactured gap; the hypothesis rho_rest>0 fails at k=0")
    print("     exactly as required.  (Contrast: 1/Var=1/k^2 is the same Gaussian value here, but the")
    print("     one-plaquette beta=0 Haar floor G''=0.2051 is a *compact-measure* curvature, not free.)")
    print()

def main():
    print("#"*92)
    print("# Q3 SCHUR-DEFECT witness: integrate out the unstable mode; bound the Schur complement.")
    print("# deterministic; py -3.12 -I; exact quadrature/linear algebra; no MC.  No QED.  Clay 0%.")
    print("#"*92)
    print()
    part1_haynsworth()
    part2_one_plaquette()
    part3_strip()
    part4_defect_lemma()
    part5_acid()
    print("="*92)
    print("VERDICT")
    print("  (1) In(H)=In(K)+In(Sigma) is exact: Sigma>0 iff all negative modes sit in the fibre K.")
    print("  (2) One-plaquette: theta=0 (fibre Savvidy) -> Sigma=alpha>0, excision WORKS; theta=pi")
    print("      (barrier, G''<0 at beta>4.30661) -> K>0, Sigma<0, collective mode: excision CANNOT help.")
    print("  (3) 1xL strip: constant-mode Schur curvature positive & saturating (~0.22 > 1/6); no")
    print("      Savvidy instability on the constant mode (S_W exactly flat there).")
    print("  (4) Rank-one defect lemma chi <= (1+eta)/rho_rest holds WHEN rho_rest>0 (T); sharp form")
    print("      chi<=1/rho_rest.  Binding hypothesis rho_rest=E[Sigma]>0 is NOT implied by eta=D/S.")
    print("      Missing 4D estimate: uniform fibre gap gamma>0 (K>=gamma I), i.e. E<K^-1 c,c> < S.")
    print("  (5) Acid test clean: free massless G''=k^2->0, no false floor.")
    print("="*92)

if __name__ == "__main__":
    main()
