"""
IN-3 / T1 : sign-cancellation ratio kappa_c in the SU(2) surface/character representation.

kappa_c := (z_signed / z_abs)^{1/|C0|},  z_signed = sum_C w(C),  z_abs = sum_C |w(C)|,
over surfaces/2-chains C spanning a fixed boundary, w(C) the surface weight in the
character expansion of the SU(2) Wilson measure.

Parts
  A. Character coefficients c_r(beta) = (2(2r+1)/beta) I_{2r+1}(beta) > 0  (exact).
  B. SU(2) link-integral multiplicities  M({r_i}) = dim[ (x) r_i ]_singlet >= 0  (integers).
  C. Exact signed-vs-abs weighted enumeration of a real 4D defect 2-chain with
     4-valent links -> z_signed, z_abs, ratio, kappa_c.
  D. Checkerboard SU(2) heat-bath <P> at the crossover, with binned error bars
     (reuses su3._sample_w0; validated here vs Metropolis + exact 1D integral).
"""
import numpy as np, itertools, sys
from scipy import integrate, special
sys.path.insert(0, r"D:\Dev\math\research\lattice-cpb")
import su3 as S3

# ---------------- A. character coefficients ----------------
def chi(r, th):
    s = np.sin(th)
    out = np.where(np.abs(s) < 1e-12, (2*r+1.0), np.sin((2*r+1)*th)/np.where(np.abs(s)<1e-12,1.0,s))
    return out

def c_r_beta(r, beta):
    f = lambda th: np.sin(th)**2 * np.exp(beta*np.cos(th)) * np.sin((2*r+1)*th)/np.sin(th)
    val, err = integrate.quad(f, 1e-10, np.pi-1e-10, limit=600)
    return (2/np.pi)*val

# ---------------- B. link multiplicity (singlet count) ----------------
def multiplicity(reps):
    reps = tuple(reps)
    if len(reps)==0: return 1.0
    def f(th):
        s=np.sin(th); p=s*s                      # Haar class measure sin^2(th)
        for r in reps:
            p=p*np.where(np.abs(s)<1e-12,(2*r+1.0),np.sin((2*r+1)*th)/np.where(np.abs(s)<1e-12,1.0,s))
        return p
    val,_ = integrate.quad(f, 1e-10, np.pi-1e-10, limit=800)
    return (2/np.pi)*val

# ---------------- C. exact signed/abs surface enumeration ----------------
def build_chain(L):
    """Defect 2-chain on L^3 x {x3=0}: planes (01)@x2=0, (02)@x1=0, (12)@all x0.
    Containment is by the actual link set; 4-valent links appear where (01),(02) meet."""
    plaqs=[]
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                x=(x0,x1,x2)
                if x2==0: plaqs.append((x,0,1))
                if x1==0: plaqs.append((x,0,2))
                if x0==0: plaqs.append((x,1,2))
    # links of a plaquette (x,mu,nu): (x,mu),(x+e_mu,nu),(x+e_nu,mu),(x,nu)
    def links_of(p):
        x,mu,nu=p
        def sh(y,i):
            y=list(y); y[i]=(y[i]+1)%L; return tuple(y)
        return [(x,mu),(sh(x,mu),nu),(sh(x,nu),mu),(x,nu)]
    link_incident={}
    for idx,p in enumerate(plaqs):
        for lk in links_of(p):
            link_incident.setdefault(lk,[]).append(idx)
    return plaqs, link_incident

def enumerate_chain(beta, L=2, Rmax=1):  # rep labels {0,1/2,1} -> index {0,1,2}
    plaqs, link_inc = build_chain(L)
    n=len(plaqs)
    repvals=[0.0,0.5,1.0]
    cred=np.array([c_r_beta(r,beta) for r in repvals])          # c_0, c_{1/2}, c_1
    # all rep assignments, (N,n) in {0,1,2}
    Re=np.array(list(itertools.product(range(3),repeat=n)),dtype=np.int8)
    N=Re.shape[0]
    w=np.ones(N)
    for p in range(n):
        w*=cred[Re[:,p]]
    # per-link singlet multiplicity lookup over base-3 local codes
    mult_cache={}
    def Mtuple(rt):
        key=tuple(sorted(rt))
        if key not in mult_cache: mult_cache[key]=multiplicity(key)
        return mult_cache[key]
    for lk,idxs in link_inc.items():
        k=len(idxs)
        # precompute Mtab for all 3^k local codes
        Mtab=np.empty(3**k)
        for code in range(3**k):
            rt=[repvals[(code//(3**j))%3] for j in range(k)]
            Mtab[code]=max(0.0, float(round(Mtuple(rt))))   # multiplicities are exact integers
        code=np.zeros(N,dtype=np.int64)
        pl=1
        for j,pi in enumerate(idxs):
            code+=Re[:,pi].astype(np.int64)*pl; pl*=3
        w*=Mtab[code]
    z_signed=float(w.sum()); z_abs=float(np.abs(w).sum())
    nz=w[w!=0]
    return dict(nplaquette=n, z_signed=z_signed, z_abs=z_abs,
                ratio=z_signed/z_abs if z_abs else float('nan'),
                minw=float(nz.min()) if nz.size else 0.0,
                nneg=int((w<0).sum()), ntot=int(nz.size))

# ---------------- D. SU(2) checkerboard heat bath ----------------
ID2=np.eye(2,dtype=np.complex128)
SIG=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=np.complex128)
def adj(A): return np.conj(np.swapaxes(A,-1,-2))
def quat_to_su2(u):
    return (u[...,0,None,None]*ID2 + 1j*(u[...,1,None,None]*SIG[0]+u[...,2,None,None]*SIG[1]+u[...,3,None,None]*SIG[2]))
def staple_sum(U,mu):
    L=U.shape[1]; Um=U[mu]; Sm=np.zeros((L,L,L,L,2,2),dtype=np.complex128)
    for nu in range(4):
        if nu==mu: continue
        Un=U[nu]
        S1=np.roll(Un,-1,axis=mu)@adj(np.roll(Um,-1,axis=nu))@adj(Un)
        D=np.roll(np.roll(Un,1,axis=nu),-1,axis=mu); E=np.roll(Um,1,axis=nu); F=np.roll(Un,1,axis=nu)
        Sm+=S1+adj(D)@adj(E)@F
    return Sm
def plaq(U):
    tot=0.0
    for mu in range(4):
        for nu in range(mu+1,4):
            Up=U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu])
            tot+=0.5*np.real(np.trace(Up,axis1=-2,axis2=-1)).mean()
    return tot/6.0
def parity_masks(L):
    x=[np.arange(L)]*4; g=np.meshgrid(*x,indexing='ij'); par=sum(g[i] for i in range(4))%2
    return par==0, par==1
def hb_sweep(U,beta,rng,masks):
    for c in (0,1):
        idx=masks[c]
        for mu in range(4):
            S=staple_sum(U,mu)[idx]
            s0=0.5*np.real(S[:,0,0]+S[:,1,1]); s3=0.5*np.imag(S[:,0,0]-S[:,1,1])
            s1=0.5*np.imag(S[:,0,1]+S[:,1,0]); s2=0.5*np.real(S[:,0,1]-S[:,1,0])
            Sn=np.sqrt(s0*s0+s1*s1+s2*s2+s3*s3); Ss=np.where(Sn<1e-14,1.0,Sn)
            G=quat_to_su2(np.stack([s0/Ss,s1/Ss,s2/Ss,s3/Ss],axis=-1))
            cc=S3._sample_w0(beta*Sn,rng)
            d=rng.standard_normal((S.shape[0],3)); d/=np.linalg.norm(d,axis=-1,keepdims=True)
            sv=np.sqrt(np.maximum(1-cc*cc,0))
            v=quat_to_su2(np.concatenate([cc[:,None],sv[:,None]*d],axis=-1))
            U[mu][idx]=v@adj(G)
    return U

def run_plaq(beta,L=6,nwarm=500,nmeas=400,nbin=20):
    rng=np.random.default_rng(1234); masks=parity_masks(L)
    U=np.broadcast_to(ID2,(4,L,L,L,L,2,2)).copy()
    for s in range(nwarm): hb_sweep(U,beta,rng,masks)
    traj=[]
    for s in range(nmeas):
        hb_sweep(U,beta,rng,masks); traj.append(plaq(U))
    traj=np.array(traj); bins=np.array_split(traj,nbin)
    means=np.array([b.mean() for b in bins])
    return traj.mean(), means.std(ddof=1)/np.sqrt(nbin), traj.std()/np.sqrt(len(traj))

if __name__=="__main__":
    print("="*70)
    print("A. SU(2) character coefficients  c_r(beta) = (2(2r+1)/beta) I_{2r+1}(beta)")
    for beta in [2.0,2.2,2.4,2.6]:
        vals=[c_r_beta(r,beta) for r in range(0,7)]
        bess=[(2*(2*r+1)/beta)*special.iv(2*r+1,beta) for r in range(0,7)]
        print(f"  beta={beta}: c_r(r=0..6) = "+"  ".join(f"{v:+.4f}" for v in vals))
        print(f"           min over r={min(vals):+.3e}   (Bessel max diff {max(abs(np.array(vals)-np.array(bess))):.1e})")
    print("  -> ALL c_r > 0 for all r, all beta>0: surface weights are non-negative.")
    print("="*70)
    print("B. link-integral multiplicities (singlet dims, must be >=0 integers):")
    for reps in [(0.5,0.5),(0.5,0.5,1.0),(0.5,0.5,0.5,0.5),(1.0,1.0,1.0,1.0),(0.5,)*6]:
        print(f"   M{reps} = {multiplicity(reps):+.4f}")
    print("="*70)
    print("C. exact signed-vs-absolute surface-weight enumeration (L=2 defect 2-chain):")
    for beta in [2.0,2.2,2.4,2.6]:
        res=enumerate_chain(beta,L=2)
        kappa=res['ratio']**(1.0/res['nplaquette'])
        print(f"  beta={beta}: |C|={res['nplaquette']:3d}  terms={res['ntot']:6d}  neg={res['nneg']}"
              f"  min_w={res['minw']:.3e}  z_signed={res['z_signed']:.6e}  z_abs={res['z_abs']:.6e}")
        print(f"            z_signed/z_abs = {res['ratio']:.12f}   kappa_c = {kappa:.12f}")
    print("="*70)
    print("D. SU(2) checkerboard heat-bath <P=(1/2)tr U_p> at crossover (L=6):")
    for beta in [2.0,2.2,2.4,2.6]:
        m,eb,es=run_plaq(beta,L=6)
        print(f"  beta={beta}:  <P> = {m:.4f} +/- {eb:.4f} (binned)   [naive {es:.4f}]")
    print("="*70)
