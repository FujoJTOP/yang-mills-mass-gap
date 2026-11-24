"""
falsify_topological.py -- FALSIFY-TOPOLOGICAL witness.
Try to FORCE an atom from finite data via a topological/homotopy invariant (winding of the
Stieltjes/Herglotz function Pi on a contour). Show it reduces to "infinite data / A3".

Pi(z) = int dnu(s)/(s-z) is a Herglotz (positive-imaginary) function.
Atom = isolated pole at s1 (residue -f). Continuum = branch cut on the real axis.
"""
import numpy as np

s1 = 1.0; f = 0.5

def Pi_atom(z):            # isolated atom: simple pole at s1
    return f/(s1 - z)

def Pi_bump(z, h=0.03):    # atomless narrow bump of total mass f around s1 (NO pole)
    s = np.linspace(s1-10*h, s1+10*h, 60001)
    rho = f*np.exp(-0.5*((s-s1)/h)**2)/(h*np.sqrt(2*np.pi))
    z = np.atleast_1d(np.asarray(z, dtype=complex))
    val = np.sum(rho[:,None]/(s[:,None]-z[None,:]), axis=0)*(s[1]-s[0])
    return val if val.size > 1 else val[0]

def winding(Pi, center, rad, N=8000):
    """Winding number of the closed curve theta -> Pi(center + rad e^{i theta}) around 0."""
    th = np.linspace(0, 2*np.pi, N)
    w = Pi(center + rad*np.exp(1j*th))
    d = np.angle(w[1:]*np.conj(w[:-1]))      # incremental arg change
    return np.sum(d)/(2*np.pi)

print("="*72)
print("(1) WEAK-* DENSITY: the atomless bump matches the atom at off-axis points")
print("    (so any FINITE set of point-evals is fooled; OPE-RHO)")
print("="*72)
for z in [s1+1j*0.5, s1+0.3+1j*0.8, s1+1.0+1j*0.4, s1-0.6+1j*1.2]:
    d = abs(Pi_atom(z)-Pi_bump(z))
    print(f"  z={z!s:22s}  |Pi_atom-Pi_bump| = {d:.3e}")

print("\n"+"="*72)
print("(2) OFF-AXIS CONTOUR (all Im z >= 0.1): winding is 0 for BOTH  -> BLIND")
print("    (no pole, no cut inside; topological invariant is weak-* continuous here)")
print("="*72)
# circle centered at s1+0.6i radius 0.5  -> stays in Im z>=0.1, s1(=1+0i) OUTSIDE
cen, rad = s1+0.6j, 0.5
print(f"  off-axis circle center={cen}, radius={rad}")
print(f"  winding(Pi_atom) = {winding(Pi_atom, cen, rad):+.6f}")
print(f"  winding(Pi_bump) = {winding(Pi_bump, cen, rad):+.6f}")

print("\n"+"="*72)
print("(3) ATOMLESS MEASURE IS FEASIBLE: below the floor, nu>=0 with no atom (s2=s1)")
print("    => positivity admits isolation gap 0; no positivity/topological count excludes it.")
print("="*72)
# dnu = rho(s) ds on [s1, inf), rho>=0, no atom at s1: valid measure, atom weight 0.
print("  take dnu = rho(s)ds on [s1, 4], rho>0 smooth: nu>=0, atom at s1 = 0, isolation gap = 0.")
print("  => a matching atomless measure always exists; any invariant must be WEAK-* DISCONTINUOUS")
print("     to see the atom, hence NOT computable from finite (weak-* continuous) data.")

print("\n"+"="*72)
print("(0) BOUNDARY DENSITY: atom -> f*delta(s1) (a point), continuum -> smooth rho.")
print("    The atom is EXACTLY the weak-* discontinuity (a delta) of the boundary density;")
print("    detecting a delta requires concentrating test functions = infinite data / A3.")
print("="*72)
# Im boundary value Pi(x+i0) = -pi * rho(x); atom gives -pi f delta, bump gives -pi*rho_bump (smooth)
x = np.array([s1-0.1, s1-0.01, s1-1e-3, s1, s1+1e-3, s1+0.01, s1+0.1])
print("  rho_bump(s1±eps) stays FINITE (~peak); atom's rho is a delta (unbounded/undefined at a point).")
print("  => the topological 'pole vs cut' distinction is the delta-vs-smooth distinction = A3/infinite.")

print("\n"+"="*72)
print("(4) TO ENCLOSE the atom (at s1 ON the real axis) the contour MUST cross the real axis.")
print("="*72)
cen2, rad2 = s1+0j, 0.4            # circle through the real axis, s1 inside
print(f"  winding(Pi_atom) on real-axis-crossing circle = {winding(Pi_atom, cen2, rad2):+.4f}  (pole enclosed)")
x0 = s1+0.05
za = x0+1e-6j; zb = x0-1e-6j
jump = Pi_bump(za)-Pi_bump(zb)
print(f"  Pi_bump jump across the cut at x={x0}:  |Pi(+i0)-Pi(-i0)| = {abs(jump):.4e}   (=2*pi*rho(x))")
print("  => the continuum has a BRANCH CUT where the atom has a POLE; seeing the atom means")
print("     evaluating ACROSS the real axis (the boundary density, a delta for the atom) = A3/infinite.")
