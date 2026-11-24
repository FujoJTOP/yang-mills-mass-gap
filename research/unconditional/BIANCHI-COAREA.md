# BIANCHI-COAREA — audit of the proposed "Bianchi-aware group-valued coarea" route

> 2026-11-24. Audit + finite-dimensional kill-test of the proposed 4D lifting of the
> 2D Piazza–Shen log-Sobolev mechanism via a "Bianchi-aware group-valued coarea /
> Brascamp–Lieb" estimate.
> Witness: `research/unconditional/w7/bianchi_coarea.py` (`py -3.12`, single job).
> Labels: **[T]** true theorem · **[C]** conditional theorem · **[O]** exact missing step ·
> **[NG]** no-go for a named class · **[V-here]** new computation this round ·
> **[V-doc]** original source read this session · **[R]** external, as reported.
> **No QED is claimed.** Clay remains open.

---

## VERDICT (up front)

The proposed "Bianchi-aware coarea/Brascamp–Lieb" lemma is **not a new handle**. It reduces to
**H1** (the boundary-/volume-uniform upgrade). The concrete finite-dimensional test settles the
at-issue statement and it fails — but **not** by the route one fears (it does *not* manufacture a
spurious mass gap on the free field). It fails the other way:

- **Free U(1), linearized:** the Jacobian `DΦ_B` *is* the discrete curvature operator `d`. It
  **faithfully carries the free-field soft modes** `σ_min = 2 sin(π/L) ≈ 2π/L → 0`
  ([V-here], periodic `Z_L^4`, exact match to `2 sin(π/L)` for `L=2..6`). **No spurious uniform
  gap** — the acid sub-test is *passed*. But that is precisely why any coarea/BL constant built by
  inverting `DΦ_B` scales as `1/σ_min ~ L`: it **cannot be volume-uniform**.
- **Minimal SU(2) hypercube:** the claimed "nonabelian gain" is **not real**. At the linearized
  (free/weak-field) order `DΦ_B = d ⊗ I₃`: rank `= 51 = 3 × 17`, singular values are **exactly three
  copies** of the U(1) spectrum, and the cube-Bianchi constraint surface **coincides with the image**
  (`ker d₂ = im d₁`), so imposing Bianchi removes **zero** rank ([V-here]). At generic nonlinear
  points the only extra rank is **pure gauge framing** (the based-loop face holonomy is gauge-
  *covariant*); after gauge-fixing the nonlinear rank returns to `51` ([V-here]).
- **Brascamp–Lieb is not even instantiable here.** The only published compact-Lie-group BL theorem
  (CLL, arXiv:2602.10647) **requires the maps to be group homomorphisms** ([V-doc]); the
  plaquette-holonomy map is not a homomorphism, so that theorem gives no BL constant for it at all.

**Net:** the route is **(a) reduces to H1**, with the additional structural finding that the
nonabelian/Bianchi "gain" is empty at the only order where the soft modes live. The exact missing
step is a scale-`b`-uniform lower bound on `DΦ_B` (equivalently `σ_min ≥ c > 0`), which is literally
the boundary-uniform DLR log-Sobolev inequality — the same wall as H1. See §4.

---

## §1 Citation verification (charter §2 — nothing fabricated)

### (a) Piazza–Shen 2026 — **REAL** [V-doc]
- **J. Piazza** (Univ. Wisconsin–Madison) & **H. Shen** (Univ. Chicago), *"Log-Sobolev inequality for
  2D SU(N) lattice gauge theory at all couplings,"* **arXiv:2610.35286** [math.PR; cross-list math-ph],
  26 Oct 2026, 13 pp. DOI **10.48550/arXiv.2610.35286**.
- Content confirmed: the 2D Wilson lattice YM measure with arbitrary fixed `β` satisfies an LSI whose
  constant depends on `β,N` but is **uniform in the lattice volume**; proof by **two-scale
  disintegration**, coarse step via **uniformly positive Ricci curvature**.
- **Important honest caveat the proposal omitted:** by their own Remark 1.2 the Ricci step covers
  `SU(N)` (and `SO(N)`, `N≥3`) but **not `U(1)`** or groups with nontrivial Lie-algebra centre. So the
  cited 2D result **does not itself cover the free `U(1)` limit** — the very regime our acid test
  probes. A Nov-2026 alphaXiv follow-up extends to arbitrary compact connected Lie groups by
  replacing Ricci with a Caputo–Menz–Tetali entropy-tensorisation criterion [R]; that follow-up was
  not needed for this audit.

### (b) Cowling–Li–Liang 2026 — **REAL**, and it does **NOT** cover the plaquette map [V-doc]
- **M. G. Cowling** (UNSW), **J. Li** (Macquarie), **C.-W. Liang** (NTU), *"The Brascamp–Lieb
  inequality on compact Lie groups and its extinction on homogeneous Lie groups,"* **arXiv:2602.10647**
  [math.GR; cross-list math.CA], 11 Feb 2026, 23 pp. DOI **10.48550/arXiv.2602.10647**.
- **Homomorphism requirement — CONFIRMED.** Definition 1.1: the datum is "a `J`-tuple of
  **homomorphisms** `σ_j : G → G_j`" (continuous assumed). Theorem C ("canonical datum of compact
  Lie groups") gives finiteness iff the **codimension condition on Lie ideals** holds:
  `Σ_j (1/p_j)(dim dσ_j(𝔤) − dim dσ_j(𝔥)) ≤ dim 𝔤 − dim 𝔥` for all ideals `𝔥 ⊲ 𝔤` (eq. 1.5 / 1.6);
  when everything is connected with Haar mass 1, the constant equals `1`.
- The **only** non-homomorphic discussion is a background remark: a *nonlinear* variant with smooth
  maps "gives only a **local** inequality" under a BCCT-type derivative condition at a point; it is
  **not part of Theorems B–E**. The Heisenberg projections example is also explicitly "not
  homomorphisms" (Remark 4.8).
- **Consequence for the route:** the map `Φ_B` (face holonomy = product of four edge group elements,
  with the cube Bianchi relation) is a **polynomial/non-homomorphic** map `G^{E} → G^{F}`. CLL's
  finiteness criterion is therefore **not applicable**. Invoking "Brascamp–Lieb for compact Lie
  groups" as the engine of the proposed coarea lemma is a **misattribution** — the cited theorem does
  not apply to the object it is being asked to bound.

### (c) The claimed 2026 4D weak-coupling lattice LSI/mass-gap scheme — **REAL but NOT a proof** [V-doc]
- **ai.viXra.org** (AI-assisted e-print archive, explicitly **"not a journal"**, items "may not yet
  have been verified by peer-review" — i.e. **unrefereed**). Series under reported sole author
  **Lluis Eriksson**: *"Uniform Log-Sobolev Inequality and Mass Gap for Lattice Yang-Mills Theory: A
  Conditional Reduction,"* **ai.viXra:2602.0041** (v1 2026-02-09, v4 2026-07-31); and *"From Uniform
  LSI to Mass Gap … a Conditional and Windowed Assembly,"* **ai.viXra:2602.0054** v2.
- The pasted summary is **accurate**. Directly verified:
  - The mass-gap conclusion is **conditional** on `(H-XSD)` (a cross-scale derivative bound) and
    `(H-DOB)` / `(H-DOB-blk)` (a block-level **Dobrushin** contraction), plus `(H-P0)`. The v4 abstract
    **explicitly disclaims** any unconditional mass gap, continuum construction, Osterwalder–Schrader
    reconstruction, or Clay result.
  - **Volume windowed:** statements are restricted to `L ≤ e^{C/g²+O(1)}`.
  - **Holley–Stroock fiber cost corrected:** v1's "`O(1)` regardless of `β`" is **withdrawn**; the
    corrected form is `osc = 2β · n_plaq + C_poly`, i.e. **grows linearly in `β`** (bad at weak
    coupling). v2 walks back v1's "unconditional / bypasses Dobrushin" claims.
- **Conclusion:** it is a **conditional-and-windowed** assembly, **unrefereed**, explicitly **not** a
  proof. Its own retreat (DLR step needs an unverified block-Dobrushin condition; windowed volume;
  `β`-linear cost) is independent evidence that the last step is exactly the boundary-uniform DLR LSI
  — see §3.

---

## §2 The finite-dimensional kill-test (single 4D hypercube, max-tree gauge) [V-here]

All numbers from `w7/bianchi_coarea.py` (`py -3.12`, one job). Cubical block `[0,b]⁴`, cells fully
inside; `d_k` = coboundary. Max-tree gauge = set a maximal spanning tree of the vertex graph to zero.

### 2.0 The object
`b=1` hypercube: `V=16, E=32, F=24, cubes=8`; `d₀ (16→32)`, `d₁ (32→24)` = linearized face-holonomy
map, `d₂ (24→8)` = cube Bianchi operator. Exactness `‖d₂ d₁‖_max = 0` (Bianchi `d(da)=0` holds
identically). `rank d₁ = 17`, `rank d₂ = 7` (the 8 cube constraints have one linear dependency — the
4-cell `d₃`), and `F − rank d₂ = 24 − 7 = 17 = rank d₁`, i.e. **`im d₁ = ker d₂`**: every curl is
Bianchi-closed, and Bianchi is exactly the image condition.

### 2.1 (a) Free U(1), linearized — PASSES the "no spurious gap" acid sub-test
- Max-tree gauge leaves `E − (V−1) = 32 − 15 = 17` physical edges; **`d₁` restricted to these is
  injective, rank 17**, singular values
  `{2.734, 2.655, 2.449, 2.284, 2.255, 2.170, 1.825, 1.576, 1.558, 1.482, 1, 1, 1, 0.904, 0.697,
  0.680, 0.509}`; smallest on the unit hypercube `= 0.509`.
- **True soft modes (the real test):** on the periodic torus `Z_L⁴` the smallest **nonzero**
  singular value of `d₁` matches `2 sin(π/L)` to 4 digits for every `L`:

  | L | shape | rank | min nonzero σ | `2 sin(π/L)` | ratio |
  |---|-------|------|---------------|--------------|-------|
  | 2 | 96×64     | 45   | 2.000000 | 2.000000 | 1.0000 |
  | 3 | 486×324   | 240  | 1.732051 | 1.732051 | 1.0000 |
  | 4 | 1536×1024 | 765  | 1.414214 | 1.414214 | 1.0000 |
  | 5 | 3750×2500 | 1872 | 1.175571 | 1.175571 | 1.0000 |
  | 6 | 7776×5184 | 3885 | 1.000000 | 1.000000 | 1.0000 |

  So `σ_min(d₁) = 2 sin(π/L) ≃ 2π/L → 0`. The **`k_min ~ L⁻¹` soft modes are preserved, not erased**,
  and there is **no spurious uniform mass gap**. The tool does **not** manufacture a gap on the free
  field — **[V-here]**, acid sub-test passed.
- **The moral is the opposite of reassuring for the route:** because `DΦ_B` faithfully carries the
  `1/L` soft modes, any coarea/BL constant `∝ 1/σ_min` **blows up as `L`** — it is *structurally*
  incapable of being volume-uniform. The only escape would be extra rigidity from the nonabelian
  Bianchi structure; §2.2 shows there is none.

### 2.2 (b) Minimal SU(2) hypercube — the "nonabelian gain" is **not real**
- **Linearized** (weak-field/free limit, the regime that hosts the soft modes): `DΦ_B = d₁ ⊗ I₃`,
  shape `72 × 51`, **rank `51 = 3 × 17`**, and the singular-value multiset is **exactly three copies**
  of the U(1) spectrum (`np.allclose(...) = True`). The cube-Bianchi operator `d₂ ⊗ I₃` has rank `21`;
  `F·3 − 21 = 72 − 21 = 51 = rank(DΦ_B)`. **`ker(d₂ ⊗ I₃) = im(d₁ ⊗ I₃)`**: the Bianchi constraint
  surface *is* the image, so **"imposing Bianchi" and "removing gauge zero-modes" delete exactly the
  directions that were never in the image — the rank and every singular value are unchanged.**
  The gain is therefore **only `dim su(2) = 3` colour multiplicity** — no rigidity.
- **Nonlinear, generic configs:** identity-config rank `51` (reproduces linearization). With the
  **tree edges pinned to identity** (gauge fixed) the Jacobian on the `51` physical edge-directions has
  rank `51` at every random seed `{51,51,51,51}` = the physical cap `3(E−V+1)`. The misleading "`72`"
  appears **only** for the *un-fixed* based-loop map, because the face holonomy is gauge-*covariant*
  (conjugates by `g_v`), so the 45 gauge directions act by the Lie bracket — **pure gauge framing, not
  physical rigidity**.
- **Verdict on (b):** the "nonabelian gain" **vanishes** once constraints and gauge zero-modes are
  removed. Rank and singular values are literal colour copies of the abelian problem. The route's
  hoped-for 4D advantage (the 2D convolution structure being replaced by a "richer" 6-plaquette
  projector) does **not** produce any linear-order rigidity.

### 2.3 (c) Scaling with block side `b` — degrades, not volume-uniform
Smallest nonzero singular value of the gauge-fixed `d₁` on side-`b` blocks:

| b | V | E | F | phys | rank | `σ_min` | `σ_min·b` |
|---|---|---|----|------|------|---------|-----------|
| 1 | 16  | 32   | 24   | 17   | 17   | 0.508865 | 0.509 |
| 2 | 81  | 216  | 216  | 136  | 136  | 0.215086 | 0.430 |
| 3 | 256 | 768  | 864  | 513  | 513  | 0.118896 | 0.357 |
| 4 | 625 | 2000 | 2400 | 1376 | 1376 | 0.075535 | 0.302 |

Log-log fit: **`σ_min(b) ~ b^{−1.368}`**; the clean periodic case gives the true asymptotic `2π/L`
(exact), so the continuum block behaviour is `σ_min ≃ π/b` (smallest Dirichlet 1-form mode of a
side-`b` cube). Hence `1/σ_min ~ b` — **linear degradation** (slightly steeper at small `b` from
discreteness). **Not volume-uniform.**

**β-dependence (folding in the measure).** The linear map is `β`-blind; `β` enters through the
curvature weight `e^{−β tr F²}`. The natural "control" combination is `β σ_min² ~ β/b²`, which
**shrinks at weak coupling** (`β = 2N/g²` large), and the weak-coupling **Holley–Stroock fiber cost
`2β·n_plaq`** (the vixra author's own corrected Lemma, §1c) grows **linearly in `β`**. Combined with
the asymptotic-freedom **volume window `L ≤ e^{C/g²+O(1)}`**, this reproduces exactly the
windowed/conditional status the vixra series concedes. So the task's hypothesized "`b²` or
exponentially at weak coupling" degradation is borne out in the combined `β/b²` + `β`-linear-cost
form — but the *dominant, irreducible* obstruction is the bare `1/σ_min ~ b` soft-mode denominator
already present in the free theory.

---

## §3 Does the LSI → gap chain require boundary-uniform (DLR) LSI, and is that the same wall?

**Yes to both.** [C]/[O]

- The standard chain is `LSI ⇒ hypercontractivity ⇒ exponential clustering ⇒ transfer-matrix gap`.
  A **bulk (volume-uniform)** LSI controls clustering *inside a finite box*. To pass to an
  **infinite-volume transfer-matrix gap uniform in the transverse/gluing boundary condition** one needs
  the **boundary-uniform (DLR) LSI** — the *uniform-in-boundary-condition* log-Sobolev inequality in
  the sense of Stroock–Zegarlinski. This is a **stronger** object than the bulk LSI, and it is the
  step the vixra assembly explicitly cannot close (its `(H-DOB-blk)` block-Dobrushin hypothesis + the
  `L ≤ e^{C/g²}` window are precisely the DLR gap).
- The boundary-uniform LSI is, by definition, a **uniform lower bound on the Dirichlet form / the
  conditional variance** across the boundary — i.e. a uniform control of the **inverse susceptibility
  / soft-mode gap**. That is **exactly H1** (the scale-uniform inverse-susceptibility/soft-mode
  lower bound; equivalently `f₁² > 0` transported uniformly). So:
  **boundary-uniform (DLR) LSI ≡ H1.** The LSI→gap chain does not avoid the wall; it *lands on it*.
- The free-U(1) computation in §2.1 makes this concrete: the Dirichlet form's lowest mode is the
  `k_min ~ 1/L` soft mode (`σ_min ≃ 2π/L`); a `b`-uniform lower bound on it is exactly what is
  missing, and it is the same object whether phrased as "coarea/BL constant uniform in `b`",
  "boundary-uniform DLR LSI", or "H1".

---

## §4 The exact missing step

The proposed lemma would need:

> **[O] (the missing step)** a lower bound `σ_min(DΦ_B) ≥ c > 0` **uniform in the block side `b`**
> (hence uniform in volume), for the face-holonomy map on side-`b` blocks — equivalently a `b`-uniform
> Poincaré/log-Sobolev/coarea constant for `Φ_B`, equivalently the boundary-uniform DLR LSI.

The finite-dimensional test shows this **cannot** be obtained from the linearized/nonabelian structure:
`σ_min(b)` carries the genuine free-field soft modes `≃ π/b` (exact `2 sin(π/L)` on the torus), and the
nonabelian/Bianchi structure adds **zero** rigidity at that order (rank `= 3 ×` abelian, singular
values identical, Bianchi surface `=` image). Therefore:

- The route is **(a): it reduces to H1.** Any volume-uniform coarea/BL bound here is *literally* the
  boundary-uniform LSI = H1, and the only gap between the proposed lemma and H1 is the `b`-uniform
  lower bound — which is H1 restated, not a new tool.
- It is **not** (b) a failure by spurious gap — the free-U(1) soft modes are faithfully preserved; the
  tool does not manufacture a gap.
- It is **not** merely (c) an unproven conjecture at this scale — the finite-dimensional statement
  (`rank`, singular values, scaling) is **settled here**: no gain, linear degradation.
- Additional structural **[NG] for this route**: the cited compact-Lie-group Brascamp–Lieb theorem
  (CLL) **requires homomorphisms** and therefore **cannot be applied** to the non-homomorphic
  plaquette-holonomy map; the "Brascamp–Lieb" engine the proposal names is not available for the
  object it must bound.

**One-line honest statement.** *The Bianchi-aware coarea carries the free-field soft modes correctly
(no fake gap), but the nonabelian/Bianchi structure supplies no rigidity beyond colour multiplicity;
hence its constant degrades like the block size and the lemma is H1 (the boundary-uniform DLR LSI) in
disguise. The exact missing step is a `b`-uniform lower bound on `σ_min(DΦ_B)` — i.e. H1 itself.*

---

### Refs & provenance
- **[V-doc]** J. Piazza, H. Shen, *Log-Sobolev inequality for 2D SU(N) lattice gauge theory at all
  couplings*, arXiv:2610.35286 (2026), DOI 10.48550/arXiv.2610.35286.
- **[V-doc]** M. G. Cowling, J. Li, C.-W. Liang, *The Brascamp–Lieb inequality on compact Lie groups
  and its extinction on homogeneous Lie groups*, arXiv:2602.10647 (2026), DOI 10.48550/arXiv.2602.10647.
- **[V-doc]** ai.viXra:2602.0041 (v4, 2026-07-31) and ai.viXra:2602.0054 (v2), ai.viXra.org
  (unrefereed, AI-assisted e-print archive; not a journal). Reported sole author L. Eriksson.
- **[V-here]** `research/unconditional/w7/bianchi_coarea.py` — all ranks, singular values, scaling,
  and the SU(2) nonlinear Jacobian checks in §2.
- Context: `docs/BRIEF-FOR-MODEL.md`, `PROOF-SKELETON.md`, `research/unconditional/
  TRANSPORT-FIVE-WAY-VERDICT.md`, `CLAUDE.md §5`.
