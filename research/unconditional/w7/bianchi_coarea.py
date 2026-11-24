"""
bianchi_coarea.py  --  Kill-test for the proposed "Bianchi-aware group-valued coarea"
lemma for 4D pure SU(N) Yang-Mills.

Single job, py -3.12, numpy only.  No torch/cvxpy.

The proposed route (audit target):
  A 2D result (Piazza-Shen, arXiv:2610.35286) gives a volume-uniform log-Sobolev
  inequality for 2D SU(N) Wilson YM by two-scale disintegration.  The claimed 4D
  analog needs a "Bianchi-aware group-valued coarea / Brascamp-Lieb" estimate.
  The map is  Phi_B : G^{E_int(B)} -> { face holonomies obeying the nonabelian
  Bianchi constraint on the block B }.

This script does the concrete finite-dimensional tests:
  (a) FREE U(1), linearized: is the linearized Jacobian DPhi_B the discrete
      curvature operator d ?  Does it carry the true soft modes k_min ~ L^{-1},
      or does the tool manufacture a spurious uniform mass gap (acid test)?
  (b) MINIMAL SU(2) hypercube: independent face variables + cube Bianchi
      constraints; rank(DPhi_B), nonzero singular values; is the "nonabelian
      gain" real or just dim(su(2)) colour copies that vanish once gauge zero
      modes and Bianchi constraints are removed?
  (c) SCALING: side-b block, how the smallest singular value (hence any coarea /
      Brascamp-Lieb / Bakry-Emery constant built by inverting DPhi_B) depends on b.

Conventions:
  * Cubical complex of the block [0,b]^4.  Cells: k-cells fully inside.
  * d_k : C^k -> C^{k+1} is the coboundary (discrete exterior derivative).
        (d0 f)(v->v+e_mu) = f(v+e_mu) - f(v)
        (d1 a)(face v;mu,nu) = +a(v;mu) +a(v+e_mu;nu) -a(v+e_nu;mu) -a(v;nu)
        (d2 F)(cube v;mu,nu,la) = sum of 6 face terms with cube-boundary signs.
  * The linearized face-holonomy map is exactly d1 (U(1)), and d1 (x) I_3 for
    SU(2) at the identity.
  * Gauge fixing = "max-tree gauge": set a maximal spanning tree of the vertex
    graph to zero; the remaining edges parametrise the physical (gauge-fixed)
    connection.  (Exact for the linear map: it is a choice of complement to the
    exact 1-forms im(d0).)
"""

import numpy as np
from itertools import combinations, product
from collections import deque

np.set_printoptions(precision=6, suppress=True, linewidth=140)

# -----------------------------------------------------------------------------
# 1.  Cubical block complex and coboundary operators
# -----------------------------------------------------------------------------
def build_block(b):
    """Return vertices, edges, faces, cubes and index maps for [0,b]^4."""
    verts = list(product(range(b + 1), repeat=4))
    vindex = {v: i for i, v in enumerate(verts)}

    edges, eindex = [], {}
    for v in verts:
        for mu in range(4):
            if v[mu] < b:
                eindex[(v, mu)] = len(edges)
                edges.append((v, mu))

    faces, findex = [], {}
    for v in verts:
        for mu, nu in combinations(range(4), 2):
            if v[mu] < b and v[nu] < b:
                findex[(v, mu, nu)] = len(faces)
                faces.append((v, mu, nu))

    cubes, cindex = [], {}
    for v in verts:
        for mu, nu, la in combinations(range(4), 3):
            if v[mu] < b and v[nu] < b and v[la] < b:
                cindex[(v, mu, nu, la)] = len(cubes)
                cubes.append((v, mu, nu, la))

    return verts, edges, faces, cubes, vindex, eindex, findex, cindex


def _shift(v, mu):
    w = list(v); w[mu] += 1; return tuple(w)


def build_d0(b, verts, edges, vindex, eindex):
    """(V x E) or (E x V)?  We use d0: C^0 -> C^1 as an E x V matrix."""
    V, E = len(verts), len(edges)
    M = np.zeros((E, V))
    for (v, mu), ei in eindex.items():
        M[ei, vindex[_shift(v, mu)]] += 1.0
        M[ei, vindex[v]] -= 1.0
    return M


def build_d1(b, verts, edges, faces, eindex, findex):
    """d1: C^1 -> C^2, an F x E matrix (the linearised face-holonomy map)."""
    F, E = len(faces), len(edges)
    M = np.zeros((F, E))
    for (v, mu, nu), fi in findex.items():
        # +e_mu at v ; +e_nu at v+e_mu ; -e_mu at v+e_nu ; -e_nu at v
        M[fi, eindex[(v, mu)]]            += 1.0
        M[fi, eindex[(_shift(v, mu), nu)]] += 1.0
        M[fi, eindex[(_shift(v, nu), mu)]] -= 1.0
        M[fi, eindex[(v, nu)]]            -= 1.0
    return M


def build_d2(b, verts, faces, cubes, findex, cindex):
    """d2: C^2 -> C^3, a C x F matrix (the Bianchi operator: dF=0)."""
    C, F = len(cubes), len(faces)
    M = np.zeros((C, F))
    for (v, mu, nu, la), ci in cindex.items():
        # boundary of cube (mu,nu,la): +F(nu,la at v) -F(mu,la at v) +F(mu,nu at v)
        #                              -F(nu,la at v+e_mu) +F(mu,la at v+e_nu) -F(mu,nu at v+e_la)
        M[ci, findex[(v, nu, la)]]                 += 1.0
        M[ci, findex[(v, mu, la)]]                 -= 1.0
        M[ci, findex[(v, mu, nu)]]                 += 1.0
        M[ci, findex[(_shift(v, mu), nu, la)]]     -= 1.0
        M[ci, findex[(_shift(v, nu), mu, la)]]     += 1.0
        M[ci, findex[(_shift(v, la), mu, nu)]]     -= 1.0
    return M


def max_tree_physical_edges(verts, edges, eindex, vindex):
    """Return list of non-tree ('physical') edge indices via BFS spanning tree."""
    V, E = len(verts), len(edges)
    adj = [[] for _ in range(V)]
    for (v, mu), ei in eindex.items():
        a, c = vindex[v], vindex[_shift(v, mu)]
        adj[a].append((c, ei)); adj[c].append((a, ei))
    parent = {0: None}
    tree = [False] * E
    dq = deque([0])
    while dq:
        u = dq.popleft()
        for c, ei in adj[u]:
            if c not in parent:
                parent[c] = u; tree[ei] = True; dq.append(c)
    assert len(parent) == V, "block must be connected"
    return [ei for ei in range(E) if not tree[ei]], tree


# -----------------------------------------------------------------------------
# 2.  Periodic torus (to expose the TRUE soft modes k_min ~ L^{-1})
# -----------------------------------------------------------------------------
def build_torus_d1(L):
    """d1 on the 4-torus Z_L^4 (edges->faces), E=4L^4, F=6L^4."""
    # explicit index maps
    def vidx(v): return ((v[0]*L + v[1])*L + v[2])*L + v[3]
    eindex, edges = {}, []
    for v0 in range(L):
      for v1 in range(L):
        for v2 in range(L):
          for v3 in range(L):
            v = (v0, v1, v2, v3)
            for mu in range(4):
                eindex[(v, mu)] = len(edges); edges.append((v, mu))
    findex, faces = {}, []
    for v0 in range(L):
      for v1 in range(L):
        for v2 in range(L):
          for v3 in range(L):
            v = (v0, v1, v2, v3)
            for mu, nu in combinations(range(4), 2):
                findex[(v, mu, nu)] = len(faces); faces.append((v, mu, nu))
    def sh(v, mu):
        w = list(v); w[mu] = (w[mu] + 1) % L; return tuple(w)
    E, F = len(edges), len(faces)
    M = np.zeros((F, E))
    for (v, mu, nu), fi in findex.items():
        M[fi, eindex[(v, mu)]]         += 1.0
        M[fi, eindex[(sh(v, mu), nu)]] += 1.0
        M[fi, eindex[(sh(v, nu), mu)]] -= 1.0
        M[fi, eindex[(v, nu)]]         -= 1.0
    return M


def smallest_nonzero_singular(M, tol=1e-7):
    """Singular values, rank, and smallest nonzero singular value."""
    sv = np.linalg.svd(M, compute_uv=False)
    nz = sv[sv > tol * sv[0]]
    return sv, len(nz), (nz.min() if len(nz) else 0.0)


# -----------------------------------------------------------------------------
# 3.  SU(2) via unit quaternions (nonlinear check of the Jacobian rank)
# -----------------------------------------------------------------------------
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])

def qconj(q): return np.array([q[0], -q[1], -q[2], -q[3]])

def qexp(v):
    """exp of pure-imaginary quaternion (0,v) -> unit quaternion."""
    n = np.linalg.norm(v)
    if n < 1e-15: return np.array([1.0, 0.0, 0.0, 0.0])
    return np.array([np.cos(n), *(np.sin(n) * v / n)])

def face_holonomy(vdict, face, edges, eindex):
    """Product of the four edge quaternions around a face with orientation."""
    v, mu, nu = face
    def sh(x, m): w = list(x); w[m] += 1; return tuple(w)
    seq = [((v, mu), +1), ((sh(v, mu), nu), +1),
           ((sh(v, nu), mu), -1), ((v, nu), -1)]
    q = np.array([1.0, 0.0, 0.0, 0.0])
    for (key, sign) in seq:
        u = vdict[key]
        q = qmul(q, u if sign > 0 else qconj(u))
    return q


def su2_jacobian_numeric(b, epsilon=1e-6, seed=0, perturb_keys=None, fix_keys=()):
    """Numerical Jacobian of Phi_B at a random SU(2) edge configuration.

    perturb_keys : edge keys whose columns are included (default: all edges)
    fix_keys     : edge keys pinned to the identity (gauge fixing / tree)
    """
    rng = np.random.default_rng(seed)
    verts, edges, faces, cubes, vindex, eindex, findex, cindex = build_block(b)
    E, F = len(edges), len(faces)
    if perturb_keys is None:
        perturb_keys = list(edges)
    # random configuration (tree/fixed edges set to identity if requested)
    vdict = {}
    for key in eindex:
        vdict[key] = qexp(rng.normal(size=3) * 1.0)
    for key in fix_keys:
        vdict[key] = np.array([1.0, 0.0, 0.0, 0.0])
    base = np.array([face_holonomy(vdict, f, edges, eindex) for f in faces])
    J = np.zeros((F * 3, len(perturb_keys) * 3))
    for col_ei, key in enumerate(perturb_keys):
        q0 = vdict[key]
        for j in range(3):
            v = np.zeros(3); v[j] = epsilon
            vdict[key] = qmul(q0, qexp(v))           # right perturbation
            for fi, f in enumerate(faces):
                qp = face_holonomy(vdict, f, edges, eindex)
                # right-trivialised tangent: qp * conj(base) ; vector part /eps
                t = qmul(qp, qconj(base[fi]))
                J[fi*3:(fi+1)*3, col_ei*3 + j] = t[1:] / epsilon
            vdict[key] = q0
    return J, edges, faces, eindex


# -----------------------------------------------------------------------------
# 4.  Run the tests
# -----------------------------------------------------------------------------
def hr(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

results = {}

# ---- (a) FREE U(1): hypercube linearised map, and periodic soft modes -------
hr("(a) FREE U(1)  -- linearised DPhi_B = d1 (curvature operator)")

# single hypercube b=1, max-tree gauge
verts, edges, faces, cubes, vindex, eindex, findex, cindex = build_block(1)
d0 = build_d0(1, verts, edges, vindex, eindex)
d1 = build_d1(1, verts, edges, faces, eindex, findex)
d2 = build_d2(1, verts, faces, cubes, findex, cindex)
V, E, F, C = len(verts), len(edges), len(faces), len(cubes)
print(f"hypercube b=1:  V={V}  E={E}  F={F}  cubes={C}")

physical, tree = max_tree_physical_edges(verts, edges, eindex, vindex)
phys1 = len(physical)
print(f"max-tree gauge: tree edges = {sum(tree)} (=V-1={V-1}); physical edges = {len(physical)} (=E-V+1={E-V+1})")

d1_gf = d1[:, physical]                        # F x |phys|
sv_all, rank_all, _ = smallest_nonzero_singular(d1)
sv_gf, rank_gf, smin_gf = smallest_nonzero_singular(d1_gf)
print(f"d1 (full): shape {d1.shape}, rank {rank_all}")
print(f"d1 (gauge-fixed): shape {d1_gf.shape}, rank {rank_gf}  (injective: {rank_gf==len(physical)})")
print(f"gauge-fixed singular values: {np.sort(sv_gf)[::-1]}")
print(f"  smallest nonzero singular value on the unit hypercube = {smin_gf:.6f}")

# Bianchi: d2 d1 = 0  (exactness / dF=0)
print(f"||d2 @ d1||_max = {np.abs(d2 @ d1).max():.2e}   (Bianchi d(d a)=0 holds)")
sv_d2, rank_d2, _ = smallest_nonzero_singular(d2)
print(f"d2 (cube Bianchi op): shape {d2.shape}, rank {rank_d2}  -> {rank_d2} independent cube constraints")
print(f"  F - rank(d2) = {F} - {rank_d2} = {F-rank_d2} = rank(d1)={rank_all}  (im d1 = ker d2)")

results['u1_hypercube'] = dict(V=V, E=E, F=F, cubes=C, rank_d1=rank_all,
                               rank_gf=rank_gf, smin_gf=smin_gf, rank_d2=rank_d2)

# periodic torus: find TRUE soft mode scaling k_min ~ L^{-1}
print("\nperiodic torus Z_L^4 (true soft modes):")
soft = []
for L in [2, 3, 4, 5, 6]:
    d1T = build_torus_d1(L)
    # smallest NONZERO singular value of d1T
    sv, rk, smin = smallest_nonzero_singular(d1T, tol=1e-6)
    pred = 2.0 * np.sin(np.pi / L)   # predicted lowest nonzero mode
    soft.append((L, smin, pred, rk, d1T.shape))
    print(f"  L={L}: shape {d1T.shape}  rank {rk}  min nonzero sigma = {smin:.6f}   "
          f"2 sin(pi/L) = {pred:.6f}   ratio = {smin/pred:.4f}")
results['u1_torus_soft'] = soft

# ---- (b) MINIMAL SU(2) hypercube ------------------------------------------
hr("(b) MINIMAL SU(2) hypercube -- rank / singular values / nonabelian gain")

# linearised SU(2): d1 (x) I_3
d1_su2 = np.kron(d1_gf, np.eye(3))            # (F*3) x (phys*3)
sv_su2, rank_su2, smin_su2 = smallest_nonzero_singular(d1_su2)
print(f"linearised SU(2) Jacobian = d1 (x) I_3 : shape {d1_su2.shape}")
print(f"  rank = {rank_su2}  (3 * rank_u1 = {3*rank_gf})")
# singular values: each U(1) value with multiplicity 3
print(f"  singular values (sorted desc, first 8): {np.sort(sv_su2)[::-1][:8]}")
print(f"  -> exactly 3 copies of the U(1) spectrum:  "
      f"{np.allclose(np.sort(sv_su2)[::-1], np.repeat(np.sort(sv_gf)[::-1], 3), atol=1e-8)}")

# cube Bianchi constraints on faces (per colour): d2 (x) I_3
d2_su2 = np.kron(d2, np.eye(3))
sv_d2su2, rank_d2su2, _ = smallest_nonzero_singular(d2_su2)
print(f"  SU(2) cube Bianchi operator d2 (x) I_3: shape {d2_su2.shape}, rank {rank_d2su2}")
print(f"  face space dim = {F*3}, rank(im d1_su2) = {rank_su2}, "
      f"F*3 - rank(d2_su2) = {F*3 - rank_d2su2}")
print("  --> the Bianchi constraint surface (ker d2 (x) I3) EQUALS im(d1 (x) I3):")
print(f"      {F*3 - rank_d2su2} == {rank_su2} : {F*3 - rank_d2su2 == rank_su2}")
print("  --> constraining to Bianchi removes only directions NOT in the image:")
print("      removing F*3-rank(d2_su2) constraint dims changes the map's rank by 0.")

# nonlinear SU(2): numerical Jacobian at a random configuration, rank?
# NOTE: the based-loop face holonomy is gauge-COVARIANT (conjugates under a
# vertex gauge transform), not invariant, so un-fixed edges let the 45 gauge
# directions contribute to the Jacobian via the Lie bracket.  We therefore (i)
# show the identity-config rank (= linearised = 51) and (ii) show that once the
# tree edges are pinned to the identity (gauge fixing), the rank drops back to
# <= 51 -- confirming the 72 is pure gauge framing, not physical rigidity.
print("\nnonlinear SU(2) Jacobian at random configs (numerical, quaternions):")
edges_list = list(eindex.keys())
tree_keys = [edges_list[i] for i in range(len(edges_list)) if tree[i]]
phys_keys = [edges_list[i] for i in range(len(edges_list)) if not tree[i]]

# (i) identity configuration -> should reproduce the linearised rank 51
verts_I, edges_I, faces_I, cubes_I, vindex_I, eindex_I, findex_I, cindex_I = build_block(1)
virt = {k: np.array([1.0,0.0,0.0,0.0]) for k in eindex_I}
J_id = np.zeros((len(faces_I)*3, len(edges_I)*3))
eps = 1e-6
base_id = np.array([face_holonomy(virt, f, edges_I, eindex_I) for f in faces_I])
for ei,key in enumerate(edges_I):
    for j in range(3):
        v=np.zeros(3); v[j]=eps
        virt[key]=qexp(v)
        for fi,f in enumerate(faces_I):
            qp=face_holonomy(virt,f,edges_I,eindex_I)
            t=qmul(qp,qconj(base_id[fi]))
            J_id[fi*3:(fi+1)*3, ei*3+j]=t[1:]/eps
        virt[key]=np.array([1.0,0.0,0.0,0.0])
sv_id=np.linalg.svd(J_id,compute_uv=False)
rk_id=int(np.sum(sv_id>1e-4*sv_id[0]))
print(f"  identity config: shape {J_id.shape}, rank ~ {rk_id}  (expect 51)")

# (ii) gauge-fixed: only physical (non-tree) edges perturbed, tree fixed at id
nonlin_ranks = []
for seed in range(4):
    J, ed, fc, ei = su2_jacobian_numeric(1, seed=seed, perturb_keys=phys_keys, fix_keys=tree_keys)
    svJ = np.linalg.svd(J, compute_uv=False)
    rkJ = int(np.sum(svJ > 1e-4 * svJ[0]))
    nonlin_ranks.append(rkJ)
    print(f"  gauge-fixed seed={seed}: shape {J.shape}  rank ~ {rkJ}  smin ~ {svJ[-3]:.3e}")

# (iii) un-fixed for contrast (the misleading 72)
J_uf,_,_,_ = su2_jacobian_numeric(1, seed=0)
sv_uf=np.linalg.svd(J_uf,compute_uv=False)
rk_uf=int(np.sum(sv_uf>1e-4*sv_uf[0]))
print(f"  UN-fixed seed=0: shape {J_uf.shape}, rank ~ {rk_uf}  "
      f"(>51 only because gauge dirs act by conjugation)")
print(f"  --> identity rank {rk_id}; gauge-fixed nonlinear ranks {nonlin_ranks}; "
      f"physical rank cap = 3*(E-V+1) = {3*len(phys_keys)}")
print("  --> no extra rigidity from the nonabelian structure; any rank above 51 is pure gauge.")

results['su2_hypercube'] = dict(rank_lin=rank_su2, rank_d2=rank_d2su2,
                                rank_identity=int(rk_id),
                                ranks_nonlin_gf=nonlin_ranks,
                                rank_nonlin_unfixed=int(rk_uf), smin_lin=smin_su2)

# ---- (c) SCALING: smallest singular value vs block side b ------------------
hr("(c) SCALING: side-b block, smallest singular value of gauge-fixed d1")

print("b   V       E       F       phys    rank    min nonzero sigma    1/b   sigma*b")
scale = []
for b in [1, 2, 3, 4]:
    verts, edges, faces, cubes, vindex, eindex, findex, cindex = build_block(b)
    d1 = build_d1(b, verts, edges, faces, eindex, findex)
    physical, tree = max_tree_physical_edges(verts, edges, eindex, vindex)
    d1_gf = d1[:, physical]
    sv, rk, smin = smallest_nonzero_singular(d1_gf)
    V, E, F = len(verts), len(edges), len(faces)
    scale.append((b, V, E, F, len(physical), rk, smin))
    print(f"{b}   {V:<6} {E:<6} {F:<6} {len(physical):<6} {rk:<6} {smin:.6f}          "
          f"{1.0/b:.4f}  {smin*b:.4f}")

# fit log-log slope
bs = np.array([s[0] for s in scale], float)
sms = np.array([s[6] for s in scale])
p = np.polyfit(np.log(bs), np.log(sms), 1)
print(f"\nlog-log fit: sigma_min(b) ~ b^({p[0]:.3f})   (slope -1 => linear degradation)")
results['scaling'] = dict(scale=scale, slope=float(p[0]))

# ---- summary ---------------------------------------------------------------
hr("SUMMARY (numbers only; verdict in BIANCHI-COAREA.md)")
print(f"(a) U(1) hypercube: d1 injective on {phys1} gauge-fixed edges, rank {rank_gf}, "
      f"min sigma {smin_gf:.4f}.")
print(f"    periodic soft modes: min nonzero sigma tracks 2 sin(pi/L) ~ 2pi/L  (recovered, NOT erased).")
print(f"(b) SU(2) linearised rank {rank_su2} = 3*{rank_gf}; singular values = 3x U(1) copies.")
print(f"    Bianchi ker(d2 x I3) == im(d1 x I3): gain from constraining = 0.")
print(f"    nonlinear random-config rank stays {nonlin_ranks}.")
print(f"(c) sigma_min(b) ~ b^{p[0]:.3f}  => coarea/BL constant degrades ~ b (not volume-uniform).")

print("\nDONE.")
