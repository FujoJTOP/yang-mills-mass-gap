"""
star_literature.py -- witness for STAR-LITERATURE.md
Two checks, tiny (numpy only):
 (1) pure-YM 2-loop beta has NO Banks-Zaks zero: beta(g)<0 for all g>0, for SU(2),SU(3).
 (2) a single-coupling (1D) flow dg/dt = beta(g) with beta<0 everywhere is strictly monotone
     => NO limit cycle (no periodic orbit). Contrast: a 2D flow can have a limit cycle.
Run: py -3.12 -I research/unconditional/w7/star_literature.py
Labels: illustrative; the literature content is in STAR-LITERATURE.md, not here.
"""
import numpy as np

print("="*72)
print("(1) pure-YM 2-loop beta:  beta(g) = -b0 g^3 - b1 g^5 ;  b0,b1 > 0  =>  no BZ zero")
print("="*72)
# SU(N) pure gauge, 2-loop:  b0 = 11N/(3*16pi^2) = 11N/(48pi^2)
#                            b1 = 34N^2/(3*(16pi^2)^2) = 17 N^2/(384 pi^4)
for N in (2, 3):
    b0 = 11*N/(48*np.pi**2)
    b1 = 17*N*N/(384*np.pi**4)
    beta = lambda g: -b0*g**3 - b1*g**5
    g = np.linspace(1e-3, 50.0, 200000)
    bvals = beta(g)
    print(f"  SU({N}): b0={b0:.4e}  b1={b1:.4e}   max beta = {bvals.max():.3e}  (should be <0)")
    # look for any positive root of beta(g)=0  (Banks-Zaks zero)
    s = np.sign(bvals)
    sign_changes = np.where(np.diff(s) != 0)[0]
    print(f"          sign changes of beta on g>0: {len(sign_changes)}  (0 => no BZ zero)")

print()
print("="*72)
print("(2) 1D flow dg/dt=beta(g), beta<0 everywhere => monotone => NO limit cycle")
print("   (limit cycles need >=2 couplings; pure YM has one => Fortin-Grinstein-Stergiou escape closed)")
print("="*72)
b0 = 11*3/(48*np.pi**2)
beta = lambda g: -b0*g**3            # 1-loop suffices for the monotonicity point
g = 1.0                              # start
traj = [g]
dt = 0.01
for _ in range(20000):
    g = g + dt*beta(g)
    if g <= 1e-12:                   # ->0 (UV) is the only asymptote as t decreases
        g = 1e-12
    traj.append(g)
traj = np.array(traj)
mono = np.all(np.diff(traj) <= 1e-15)   # non-increasing
print(f"  trajectory strictly monotone (non-increasing)? {mono}")
print(f"  g: start=1.0  after t=200: {traj[20000]:.3e}   (flows toward 0; no return => no cycle)")

print()
print("="*72)
print("(3) contrast: a 2D system (e.g. the classic limit cycle) CAN have a closed orbit")
print("="*72)
# r' = r(1-r^2), theta' = 1  -> stable limit cycle at r=1; needs 2 DOF.
r, th = 0.2, 0.0
dt = 0.001
for _ in range(12000):
    r = r + dt*r*(1 - r*r)
    th = th + dt*1.0
print(f"  radial 1D subsystem r'=r(1-r^2): r -> {r:.4f} (->1); theta winds => closed orbit at r=1")
print("  => a limit cycle is a genuine 2D+ phenomenon; a single real coupling cannot realize it.")

print()
print("CONCLUSION (matches STAR-LITERATURE.md):")
print(" * pure YM (one coupling, beta<0 all g): no BZ zero AND no limit cycle.")
print(" * so the Fortin-Grinstein-Stergiou 'scale-without-conformal via limit cycle' escape is CLOSED for pure YM.")
print(" * (star) for pure YM reduces to: gapless one-coupling AF flow => scale-invariant (fixed point) attractor.")
print(" * this is NOT proved (Polchinski 1988 left 4D scale=>conformal open; Dymarsky et al. argued it [C]).")
print(" * Clay 0%; no QED.")
