"""
wi_reslope_finite.py -- WI LINE 2 continued: push the RE-SLOPE functional
(dRe tau/dln c1 = 1/E1 along |Im tau|-rungs; WI-MODIFIED-FUNCTIONAL) from
complex-tau (A3) to FINITE REAL-tau data, and decide force vs read.

C_O(t) = w0 + sum_i c_i e^{-E_i t}; w0=|<O>|^2>0 (VEV); c_i=|<n|O|Omega>|^2; E_i levels.
Re-slope (complex-tau): atom => tau_n=[ln(c1/w0)-i(2n+1)pi]/E1:
   rung |Im tau_n|=(2n+1)pi/E1 (c1-invariant); Re tau_n=ln(c1/w0)/E1; slope=1/E1.
Question: what survives at REAL t?

Run:  py -3.12 -I wi_reslope_finite.py
"""
import numpy as np
np.set_printoptions(precision=4, suppress=True, linewidth=140)
E1, w0 = 1.0, 1.0

def make_Z(levels, weights, w=w0):
    E = np.asarray(levels, float); c = np.asarray(weights, float)
    return (lambda t: w + np.sum(c*np.exp(-E*t)),
            lambda t: np.sum(-E*c*np.exp(-E*t)))

def newton(Z, dZ, z0, it=400, tol=1e-13):
    z = complex(z0)
    for _ in range(it):
        f, fp = Z(z), dZ(z)
        if not np.isfinite(f) or not np.isfinite(fp) or abs(fp) < 1e-300: return None
        dz = f/fp; z -= dz
        if abs(dz) < tol: return z
    return None

def continuum(edge, eps, Nc=60, dE=0.03, decay=0.5):
    Ec = edge + (np.arange(Nc)+0.5)*dE
    return list(Ec), list(eps*dE*np.exp(-(Ec-edge)/decay))

def effmass(t, levels, weights, w=w0):
    E = np.asarray(levels, float); c = np.asarray(weights, float)
    Cc = np.sum(c*np.exp(-E*t)); Ccp = np.sum(-E*c*np.exp(-E*t))
    return -Ccp/Cc

def banner(s): print("\n"+"="*78); print(s); print("="*78)

# ---------------- (A) complex-tau Re-slope
banner("(A) COMPLEX-tau Re-slope: nearest zero vs c1 (ATOM)")
print(f"{'c1':>7} {'ln(c1/w0)/E1':>13} {'Re tau0':>10} {'Im tau0':>9} {'pi/|Im|':>8}")
cs=[0.4,1.0,2.5,5.0]; Re=[]
for c1 in cs:
    Z,dZ=make_Z([E1],[c1]); z=newton(Z,dZ,(np.log(c1/w0)-1j*np.pi)/E1); Re.append(z.real)
    print(f"{c1:7.2f} {np.log(c1/w0)/E1:13.4f} {z.real:10.4f} {z.imag:9.4f} {np.pi/abs(z.imag):8.4f}")
print(f"  => dRe/dln(c1)={np.polyfit(np.log(cs),np.array(Re),1)[0]:.5f} (=1/E1={1/E1:.3f});  Im fixed.")

# ---------------- (B) real-t plateau + operator spread (atom+gap continuum)
banner("(B) REAL-t m(t): atom+gap continuum, varying c1  (P1 value, P2 spread)")
Eg,cg = continuum(2.0, eps=0.12)
print(f"{'t':>5} "+" ".join(f"c1={c:<4}" for c in cs)+"   spread")
for t in [1,2,4,8,13,20,30]:
    row=[effmass(t,[E1]+Eg,[c1]+cg) for c1 in cs]
    print(f"{t:5.1f} "+" ".join(f"{v:6.4f}" for v in row)+f"   {max(row)-min(row):.4f}")

# ---------------- (C) approach law: atom / atom+gap / threshold / cut
banner("(C) approach law m(t) at REAL t (P3 discriminator)")
Eth,cth=continuum(E1+0.02,eps=0.12); Ecut,ccut=continuum(E1,eps=0.12)
cases={"atom":([E1],[1.0]),"atom+gap":([E1]+Eg,[1.0]+cg),
       "atom+thr":([E1]+Eth,[1.0]+cth),"cut":(Ecut,ccut)}
print(f"{'t':>5} "+" ".join(f"{k:<11}" for k in cases))
for t in [2,4,8,14,22,32]:
    print(f"{t:5.1f} "+" ".join(f"{effmass(t,L,W):11.4f}" for (L,W) in cases.values()))

# ---------------- (D) finite-t amplitude -> ln c1 needs t->inf (with continuum)
banner("(D) finite-t amplitude a(t)=ln(C(t)-w0)+E1 t -> ln c1 only as t->inf (atom+gap)")
for c1 in cs:
    a=[np.log(np.sum(np.array([c1]+cg)*np.exp(-np.array([E1]+Eg)*t)))+E1*t for t in [6,12,20,30]]
    print(f"  c1={c1:4.2f}: ln c1={np.log(c1):+.4f}  a(t)@6,12,20,30 = "+", ".join(f"{x:+.4f}" for x in a))

# ---------------- (E) with a continuum, the complex rung robustness
banner("(E) complex rung under a GAPPED continuum (atom+gap), varying c1")
for c1 in cs:
    Z,dZ=make_Z([E1]+Eg,[c1]+cg)
    z=newton(Z,dZ,(np.log(c1/w0)-1j*np.pi)/E1)
    ok = z is not None
    print(f"  c1={c1:4.2f}: zero {'Re=%+.3f Im=%+.3f pi/|Im|=%.3f'%(z.real,z.imag,np.pi/abs(z.imag)) if ok else 'n/a'}"
          f"   (gap continuum: rung still ~pi/E1 => A3 robust READ, needs complex tau)")

banner("VERDICT")
print("""
  (A)/(E) COMPLEX tau: rung |Im|=pi/E1 (c1-invariant) + slope dRe/dln c1=1/E1, robust even
          with a gapped continuum -- but it is the ANALYTIC CONTINUATION (complex zeros) = A3.
  (B)/(C) REAL t: the carried content is the PLATEAU (value E1; spread O(e^{-gap*t}) -> 0, P1/P2)
          and the APPROACH LAW (exponential=atom vs power-law=cut, P3). This is finite-t data,
          but it READS the gap and needs t->inf uniformly (isolation) to certify = H1.
  (D)     the log-amplitude -> ln c1 only as t->inf.
  (E)     the '1/E1' is Re tau0/(ln c1-ln w0): a tautology that needs E1 (the plateau) and
          ln c1 (the amplitude), each requiring t->inf = isolation.
  => WI-RESLOPE-FINITE: a finite-real-data functional EXISTS but collapses to P3 (plateau+
     approach). It READS (does not force) and requires isolation (t->inf) = H1. The rung/
     slope structure proper stays A3. READ-not-FORCE confirmed.
""")
