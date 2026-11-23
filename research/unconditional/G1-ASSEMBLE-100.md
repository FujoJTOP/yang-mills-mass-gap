# G1-ASSEMBLE-100 — the strongest complete G1 statement, the exact modulo list, and the 100% ledger

> 2026-11-23 · Assembles G1 (4D pure SU(N) YM existence) from G1-EXISTENCE / G1-L3 / G1-UNIFORM-ZA / G1-CONTINUUM-OS / G1-CONTINUUM-SCALING / G1-UNIQUENESS-VACUUM / PROOF-G1-BPHZ / A2 / PROOF-L3 / L3-scale-uniform / L3-band-decay / IN-2.
>
> **判词 (verdict).** The G1 **implication graph closes 100%** — every arrow is a proven theorem, no implication is missing. Every input is either **(a) a theorem**, **(b) reducible to L3**, or **(c) = the gap (G2)**. There is **no fourth category** for the full OS0/1/2/3/4 statement, with **one** exception: **L3 itself**, the sole genuinely-open hard-analysis node — and L3 *is* the single named species (the **one-sided variance upper bound**). At the IR end the same species reappears as **H1 = G2**. Hence:
> **`4D pure SU(N) YM exists (OS0/1/2/3/4)` modulo `{L3, G2}` — and `L3 = H1 = G2` are one species at two scales.** Fixed-volume existence (OS0/1/3) is modulo **L3 alone**; OS2 and OS4 add **no new species**, only the same L3(UV)+G2(IR) pair. The **classification** is 100% complete; the **Clay statement** stays 0% because the one species is unconstructed. This doc is honest bookkeeping, not a proof of the gap.

---

## 0 · One numbering caveat (read first)

The project uses its **own** OS numbering, fixed by G1-CONTINUUM-OS and the os2/os3 reports; it is **not** the textbook (Glimm–Jaffe) numbering, so do not reserialize it as if it were universal. Project convention:
`OS0` = **regularity/temperedness** (a-uniform moment bounds + Minlos factorial); `OS1` = **invariance** (torus translations/point group, then E(4)); `OS2` = **full continuous O(4)** rotation invariance (**an `a→0` irrelevance statement**, not an `L→∞` one — see §3); `OS3` = **reflection positivity**; `OS4` = **clustering / ergodicity / unique vacuum**. Reconstruction (OS73/75) needs OS0, OS1, OS3 (+OS4 to name "the" vacuum). This is the object whose mass gap is G2. *(OS = Osterwalder–Schrader.)*

---

## 1 · The exact master statement

> **MASTER (G1 existence, full OS).** Fix compact simple simply-connected `G=SU(N)`, `N≥2`. Along the two-loop renormalized trajectory `β(a)=2N/g²(a)`, `g²(a)∼1/(2b₀log(1/aΛ))`, in scheme C (L1). **There exists** a probability measure `μ` on `𝒮′(ℝ⁴; gauge-invariant local observables)` such that:
> - **(OS0)** `μ` is regular/tempered: for every `n` and compact `K`, `|S^{(n)}(f₁,…,f_n)| ≤ C_n(L)∏‖f_i‖_{C^k}`, **uniform in `a`**, plus the Minlos factorial bound (⇔ `sup_a Z_a⁻¹‖H_a‖_D<∞` on the common core `D`);
> - **(OS1)** `μ` is invariant under `E(4)=ℝ⁴⋊O(4)`;
> - **(OS2)** `μ` is invariant under the **full continuous `O(4)`** (no preferred frame);
> - **(OS3)** `μ` is reflection positive: `∫\overline{θA}A dμ ≥ 0` on the forward cone `C_RP`;
> - **(OS4)** `μ` is **clustering / ergodic** (unique vacuum);
>
> **modulo [LIST] = `{ L3 , G2 }`.**

Here `L3` = the a-uniform scale-uniformity across the `J∼log(1/a)` shell tower (`= (★)∧(★dec)`, §3), and `G2` = the mass gap `Δ>0` (`⇔ f₁²>0`, `= H1`). The rest of this file shows the list cannot be shortened *except* by the unification `L3 = H1 = G2` (§6), and that no other axiom contributes a third entry.

**The two nested theorems:**

| # | Statement | Modulo (precise) | Status |
|---|---|---|---|
| **G1-fix** | At **fixed** `T⁴_L`, along the RT/scheme C, ∃ subsequence `a_n→0` and `μ_L` satisfying **OS0, OS1, OS3** | **L3 alone** | conditional theorem (§2) |
| **G1-full** | On `ℝ⁴` (`L→∞` after `a→0`), ∃ `μ` satisfying **OS0, OS1, OS2, OS3, OS4** | **`{L3, G2}`** (OS2→{L3,G2}, OS4→G2) | conditional theorem (§3) |
| **Clay** | G1-full **∧** mass gap `Δ>0` | **`{L3, H1}`** = **one species** | **0%** — the species is unbuilt |

---

## 2 · G1-fix, restated (the proven modulo-L3 core)

This is G1-EXISTENCE verbatim: **L1** (scheme C, definitional) + **L2** (all-orders lattice↔RT match, a theorem) give the `a⁻⁴` subtraction in C; a **single-block** `G″≥δ_W` (Lemma B Haar curvature `W≥δ_W` + free-box gap `4π²/L²` a-uniform + weak coupling `g_j≤g*`) **together with L3** give `UVS(L)`; **GL3** turns `UVS(L)` into relative compactness + Minlos ⇒ `μ_L`; **GL4** passes reflection positivity (OS3) to the weak limit unconditionally. The single-block bound is the `j*`-block; **L3 is exactly the extension of it across the remaining `J−j*∼log(1/a)` shells**. OS0 from the `a`-uniform moment bounds; OS1 from lattice translation/rotation invariance. The gap never enters; the single open input is L3.

The retargeting that makes this a theorem (not the false A1-(c)): `G″=E[H]−Var+W`. The literal `E[H]−Var≥δ` is **false** (one-loop equipartition ⇒ `0`; W7: SU(2) β=2.5 `E[H]−Var=+0.02`, SU(3) β=5.7 `+0.04`). The positive curvature is the **Haar/Jacobian term** `W≥δ_W`, which `(c)` omitted. That repair is what supplies the `δ` **single-block**; extending it across the whole `J∼log(1/a)` tower *is* L3.

---

## 3 · G1-full: how OS2 and OS4 enter, and what each costs

**OS1 (full E(4)).** Translations: proven on the torus (lattice symmetry), promoted through the `L→∞` limit once the limit exists. Rotations split into the torus point group `W=(ℤ/2)⁴⋊S₄` (proven) and the continuous `O(4)` (**OS2**, below). No new input beyond OS2 + the `L→∞` existence.

**OS2 (full continuous `O(4)`) — the limit-ordering correction.** The loose phrase "OS2 at `L→∞`" is wrong and rests on a limit-ordering error (G1-OS2-INFVOL). The `O(4)→W` breaking (`W=(ℤ/2)⁴⋊S₄`, `|W|=384`, a **fixed finite closed subgroup**) is a **UV (`a`) effect, not an IR (`L`) effect**; in the ordered scheme (`a→0` at fixed `L` ⇒ `μ_L`; then `L→∞`) rotations are decided when `μ_L` is built. Weak limits preserve symmetry **downward only** (`compactness ⇒ W`, never `compactness ⇒ O(4)`), and the `W`-invariant anisotropic Gaussian `C=f(|x|²)+g(Σx_i⁴)` (RP, Ward-complete, translation-invariant) survives `L→∞` still anisotropic. ⇒ **`L→∞` buys full `ℝ⁴` translations (OS1) + volume-tightness, and zero rotation.** OS2 is the **fixed-`L`, `a→0`** statement "`μ_L` is `O(4)`-invariant," and it splits into exactly the two species:
- **GL5-KIN = KIN-2 (proven).** No `O(4)`-breaking gauge-invariant operator of dimension `≤4`: `dim<4` is constants; `dim=4` is `{tr F², tr F F̃}`, both `O(4)`-invariant; the only candidate `Σc_{μν}tr F_{μν}²` has `c` a `W`-invariant symmetric tensor ⇒ `c∝δ` (`invariant_sym2_eq_smul_one`, Lean, 0 sorry) ⇒ itself invariant. Residual breaking operators start at `dim=6`, the `ℓ=4` harmonic. *(Sharparrow: had the point group been `(ℤ/2)⁴⋊S₃` — time singled out — a dim-4 anisotropic direction would exist; `S₄` spacetime permutations kill it, so OS2 is non-vacuous.)*
- **SYM-4D = the actual OS2 content (the `a→0`, fixed-`L` lift).** The leading `W`-but-not-`O(4)` operator is dim-6 `ℓ=4`; `SYM-4D` demands its `a²` coefficient be a **nonperturbative, volume-uniform** bound (the `a²` coefficient is an all-space connected correlator). This is a `(a,L)`-uniform renormalized-moment/irrelevance bound — **the L3 species at `(a,L)`-uniformity**; its box-scale/IR component is G2. **No published nonperturbative Symanzik theorem in 4D** (Balaban RG stops at `λ∼1`). ⇒ **OS2 is *not free given L3+OS0/1/3+compactness***; it is a **separate open input, logically gap-free, but MG-adjacent in cost** (its volume-uniformity = a clustering input). → **(b) reducible to L3-species (+ G2 at the IR edge).**
- **Q-CONS = `NoSpacetimeGoldstone` (route D, only the `L→∞` decay face).** The angular charge must annihilate the vacuum, `Q_gΩ=0`. Local Ward gives only `[H,Q_g]=iB_g` with `B_g` a boundary angular flux. **F-flip:** `Q-CONS ⟸ [H,Q_g]=0` (structural, preserved by the RP continuum limit) `+ ker H=ℂΩ` (unique vacuum) `+ clustering (MG)`. ⇒ **downstream of the gap** (MG-grade). → **(c) = G2.**

So **OS2 = {proven KIN-2} + {L3-species (SYM-4D)} + {G2 (Q-CONS)}**: a **separate open input** (not free, not the gap), but it introduces **no third species** — it is the same `L3`(UV) + `G2`(IR) pair. Honest status: **open, gap-free in logic, MG-adjacent in cost.**

**OS4 (clustering / unique vacuum).** `OS4 = IN-6`. In the implication poset `MassGap ⇒ OS4 ⇒ UNI/VAC` is **strict** (`OS4 ⇏ MassGap`: critical 2D Ising is a unique, clustering, power-law Gibbs state; the generalized free field satisfies OS0/1/3/4 with `Δ>0` yet no atom, Theorem B′). So OS4 is **strictly below** the gap in logic. **Per the brief, OS4 is granted "given the gap":** and the split is exact (G1-OS4-GIVEN-GAP). **(C) clustering** is free and unconditional from *any* gap: `T=P_0+R`, `‖R‖≤e^{−aΔ}`, so `|ω(Aτ_sB)−ω(A)ω(B)|≤‖A‖‖B‖e^{−Δs}`. **(U) unique vacuum** is free iff the gap is **volume- and sector-uniform** (`(U-gap)`: `λ₁(T_Λ)≤e^{−aΔ}` independent of `Λ` and bc `ξ`) — which the project's G2 *is* (SZZ/Osterwalder–Seiler "uniform in `L,a`"; CONDITIONAL-LATTICE-GAP `Δ_Λ≥m`). Then Dobrushin uniqueness ⇒ every bc gives the same limit. If the gap is *not* uniform (fixed-`L` only), coexistence = tunneling splitting `~e^{−cL}→0`, uniqueness fails (2D Ising), and the per-sector screening `W1₀∧W1₁` must be re-bought (MG-grade, MG-sub-strength). The `VAC`/center half is otherwise half-free (center `Z_N` "vacua" die free under OS4 + center-invariant periodic bc at `T=0`). The `θ`-axis is orthogonal and RP already pins `θ∈{0,π}`. → **(c) = G2** (gap-downstream; free given the *uniform* gap).

---

## 4 · The classification ledger — every input, one category each

`(a)` theorem · `(b)` reducible to L3 · `(c)` = gap (G2) · `(d)` genuinely open hard analysis.

| # | Input / node | Category | Why / source |
|---|---|---|---|
| 1 | **L1 scheme C** | **(a)** theorem | finite reparametrisation `g'²=g²+cg⁴`, `Λ'=Λe^{c/2b₀}`; invariants `m/Λ,√σ/Λ,a/ξ`. Definitional. (A2 §L1) |
| 2 | **L2 all-orders lattice↔RT match** | **(a)** theorem | Reisz lattice BPHZ + Symanzik; same counterterm structure order-by-order. (A2 §L2, PROOF-G1-BPHZ §2) |
| 3 | **dim-4 operator class `{TrF²,TrFF̃}`, θ drops** | **(a)** theorem | BRST cohomology (ghost-#0, dim 4); θ total derivative. (W3(0), G1-UNIFORM-ZA) |
| 4 | **scale-sum `Σ_j g_j⁴=ζ(4)/(16β₀⁴)<∞`** | **(a)** theorem | `g_j²=1/(2β₀j)+O(j⁻²)`; the `Σg_j²∼loglog(1/a)` "wall" *is* the running coupling, absorbed. (BREAK-IN-1, G1-UNIFORM-ZA (A)) |
| 5 | **Lemma B `W≥δ_W`** | **(a)** theorem | Haar volume `∏(sin x/x)` log-concave on `B_{π/2}`; SU2 ≥½g, SU3 ≥¾g, scale/β-independent. (PROVE-Gpp §2) |
| 6 | **free-box gap `c_L=4π²/L²` a-uniform** | **(a)** theorem | lattice Poincaré on `ker(d)^⊥`; zero modes pure gauge; sharp. (P1-L0) |
| 7 | **single-block `G″≥δ_W`** | **(a)** theorem | one-loop equipartition + cone-tail `e^{−c/g_j²}` + Lemma B + free-box gap. **Single block only.** (A1 §4–5, G1-EXISTENCE) |
| 8 | **GL3: `UVS(L)⇒` relative compactness + diagonal subsequence + Minlos ⇒ `μ_L`** | **(a)** theorem (implication) | Lean `exists_subseq_tendsto_pi`, `gl3_imp_with_minlos`. Hypothesis = UVS(L) = node 9. |
| 9 | **OS0 = `UVS(L)` = a-uniform moments + Minlos factorial** | **(b) = L3** | proven **given** L3; the substance is L3. (G1-EXISTENCE proof) |
| 10 | **GL4: RP weakly closed (OS3)** | **(a)** theorem (unconditional) | `rp_transfer`, `rp_l2_closure_trunc`; OS78/MP87 fine-lattice RP, descends under block-spin. |
| 11 | **OS1 torus translation/point group** | **(a)** theorem | lattice E(4)-symmetry (+OS2 for continuous `O(4)`, node 12–14). |
| 12 | **OS2 GL5-KIN / KIN-2** | **(a)** theorem | no relevant `O(4)`-breaking op; `dim-6` irrelevant `O(a²)`. Lean `SYMDim4.lean` 0 sorry. (IN-2 §3) |
| 13 | **OS2 SYM-4D (nonperturbative `a²`)** | **(b) = L3 (+G2 at IR)** | `(a,L)`-uniform moment bound = L3 species; box-scale part = G2. No 4D nonperturbative Symanzik. (IN-2 §4.1) |
| 14 | **OS2 Q-CONS / NoSpacetimeGoldstone** | **(c) = G2** | F-flip: `[H,Q_g]=0 + kerH=ℂΩ + clustering(MG) ⇒ Q_gΩ=0`. (IN-2 §4.2, §5) |
| 15 | **`L→∞` limit existence / IN-6 (tightness, `UVS-L`)** | **(c) = G2 (adjacent)** | `L`-uniform tightness = L3-species at volume level; IR/box-scale obstruction = G2. Or given by gap (clustering ⇒ unique limit). (G1-fixed-volume §3; P1-L0 obstruction) |
| 16 | **OS4 clustering / unique vacuum** | **(c) = G2** | `MassGap ⇒ OS4` (Dobrushin–Föllmer); `OS4⊉ gap` but cost `W1₀∧W1₁` is MG-grade. Given the gap per brief. (G1-UNIQUENESS-VACUUM) |
| **17** | **L3 = `sup_a Z_a⁻¹‖H_a‖_D<∞` = `(★)∧(★dec)`** | **(d) — THE species** | one-sided **upper** bound on force variance per shell, uniform `j,U′,a`; §5. The **only** (d). |
| **18** | **G2 = mass gap `Δ>0` ⇔ `f₁²>0`** | **(d) — THE species (IR)** | `= H1`; `Var_μ(a_{k→0})≤1/λ`; same variance-upper-bound species as L3, IR face. |

**Reading of the table.** Rows 1–12 are **theorems** `(a)`. Rows 9, 13 reduce to **L3** `(b)`. Rows 14–16 reduce to **G2** `(c)`. Only rows **17 (L3)** and **18 (G2)** are genuinely-open hard analysis `(d)` — and by §6 they are **one functional form**. **No row is a `(d)` distinct from L3/G2.** This *is* the 100% result: the modulo list of the full OS0/1/2/3/4 statement has exactly two entries, and they are the same named species.

---

## 5 · L3, the exact remaining input (UV face)

**What L3 is.** Not reflection positivity — OS3 is free (node 10). L3 is the **OS0 regularity axiom in `a`-uniform form**, equivalently the one-derivative generator bound `sup_a Z_a⁻¹‖H_a‖_D<∞` on the smeared gauge-invariant core `D`. Its substance is **one inequality per shell**, uniform in `j,U′,a`:
```
(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E_{μ(·|U′)}⟨PᵀH_jP v,v⟩ + W_j(v) − δ_W ,   ∀|v|=1,
```
plus band decay `(★dec)`: `|Cov_μ(s_B,s_{B′})|≤C e^{−κ|B−B′|/a}`, `κ>0` uniform.

**Why it resists (the precise missing step).** `G″_j = E[H] − Var + W`, and `Var≥0` **subtracts**: `G″≤E[H]`. Coercivity `∫TrF²≥0` and cone-emptiness bound the *value* / *typical sign* of `H` (a **lower** edge), never the variance from **above**. The only bridge from the generic `Var(G_v)≤E[G_v²]∼β_j²` down to the needed `Var(G_v)≲E[H]∼β_j` is **one-loop equipartition**, exact only as `g_j→0` — and the `a⁻⁴` coefficient is an **asymptotic IR renormalon** (Borel ambiguity = `⟨TrF²⟩`, Stokes constant undetermined on `ℝ⁴`), so the bridge is *asymptotically*, not *uniformly*, true. The proven `Λ_UV` (Holley–Stroock, PROVE-H2) is `e^{−2β_jN_p}` — at `β_j∼2.6` this is `10^{−107}`, catastrophically non-uniform.

**Why it is believed true.** The deep-UV shell is *ordered* (`β_j→∞`), `|F_j|=O(g_j)`, the Savvidy cone `√(2g_j|F_j|)∼g_j` sits below `Λ_{j+1}`, so the physical `G″_j≈W≥δ_W`. L3 is *true* but not *proven*; the bridge is strict one-loop exactness, blocked by the renormalon. Balaban is the right multiscale frame but gauge-fixed and short of the large-field two-sided bound; Magnen–Rivasseau–Sénéor is UV-complete but gauge-fixed (no gauge-invariant OS, no RP). **No published constructive-QFT result reaches L3** — it is the frontier, not a cited lemma. *(Correction per charter: the "Chandra–Fröhlich" joint paper does **not exist**; the live references are Chandra–Chevyrev–Hairer–Shen (Invent. Math. 237, 2024) and Brydges–Fröhlich–Seiler — both 2D/3D, not the 4D marginal log-critical case. Shen–Zhu–Zhu (CMP 400, 2023) is strong-coupling `β<1/48` only, silent on the weak-coupling UV.)*

**What closing L3 buys.** Fixed-volume existence (OS0/1/3, subsequential `μ_L`) becomes **unconditional**. It removes **none** of G2, `L→∞`, uniqueness, or the gap.

---

## 6 · The unification: `L3 = H1 = G2`, one species at two scales

The deepest structural fact of the campaign (PROOF-SKELETON §0, §3): **L3 and H1 have the same functional form — an *upper* bound on a variance:**
- **L3 (UV, one derivative, weak coupling):** the variance upper bound `Λ_UV` uniform across the shell tower, `Var(G_v)≤E[H]+W−δ`.
- **H1 = G2 (IR, two-sided Hessian, physical coupling):** `Var_μ(a_{k→0})≤1/λ`, via the Ward identity `G″·Var(a)=1` making sign and variance bound one number.

Both fail for the *same* reason: in `G″=E[H]−Var+W` the term `Var` is **subtracted**, and `Var≥0` only ever yields the *lower* edge — no positivity axiom, no cone-emptiness, no coercivity, no BPHZ/BRST identity supplies the *upper* edge. The difference is only the **scale** (UV vs IR) and the **subspace** (`∇_⊥` fiber vs `∇_∥` `k→0`). The `a²`-scaling bridge (`G″_a/a²=Δ₀²`, `λ(μ)=Z_O²λ(μ₀)`) is exact and free (multiplicative transport through the trace anomaly, acid-safe: `0↦0`), but it transports only the **value**, never the **positivity** of the seed — **the wall is the nonzero seed, not the transport.**

**Consequence for the modulo list.** Rows 17 and 18 are the **same species**, so the full Clay modulo list `{L3, H1}` = **one** unknown. Building one tool — nonperturbative, sign-controlling, non-circular, non-SUSY (free-field-failing, `f₁²=0`), `(L,a)`-uniform, giving a *deterministic upper bound on a variance of an anti-correlated non-abelian measure* — would close **both** G1's L3 and G2's H1. That tool **does not exist** (C5/C6: all nine dependency-inequality frameworks fail — no lattice order, no attraction, no FKG, no log-submodularity for frustrated non-abelian SU(N)).

---

## 7 · The exact completion % (precise, not rhetorical)

Two different "100%"s; never conflate them (charter + PROOF-SKELETON honest rule: node-% ≠ Clay-%). 

**(i) Implication-graph / assembly completion = 100%.** Every arrow (`L1+L2 ⇒ subtraction`, `single-block + L3 ⇒ UVS(L)`, `UVS(L) ⇒ μ_L`, `RP ⇒ OS3`, `KIN-2 ⇒ OS2-restoration`, `MG ⇒ OS4`, `MG ⇒ Q-CONS`) is a proven theorem or an elementary spectral-theorem step. **No missing implication.** The chain from `{L3, G2}` to the master statement is complete.

**(ii) Non-species analytic-lemma completion = 100%.** Of the nodes that are *not* the species, **all are theorems**: rows 1–8, 10, 11, 12 (11 analytic/combinatorial lemmas), plus rows 9, 13–16 which are *reductions* (proven "given L3/G2"). **Every non-species input is proven or explicitly reduced.**

**(iii) G1-fix (OS0/1/3, fixed volume) scalar completion ≈ 92%.** 12 of 13 sub-lemmas proven; the single open node is **L3**, the substantive residue (the skeleton rates the L3 node itself ~40% progressed). The honest phrasing is **not** "92% done" but "**a theorem modulo one named lemma, L3**."

**(iv) G1-full (OS0/1/2/3/4) modulo list = `{L3, G2}`, exactly.** OS2 supplies `{proven KIN-2} + {L3:SYM-4D} + {G2:Q-CONS}`; OS4 supplies `{G2}`. So G1-full is modulo two entries, and by §6 they are **one species**. Its "completion" is therefore: **100% of the existence-vs-gap structure is classified; 0% unconditionally complete**, because the one species (`L3=H1`) is unconstructed.

**(v) Clay (existence ∧ gap) = 0%.** Millennium open problem; the named species is not built; gap `= H1 = 0%`. Per the honest rule, **no finite % claims otherwise.**

---

## 8 · The precise remaining input(s) — one line each

1. **L3** — the UV face: a per-shell a-uniform **upper** bound on the force variance `Var(Pᵀ∇S_j)≤E[PᵀH_jP]+W_j−δ_W` `(★)`, plus band decay `(★dec)`, across the `J∼log(1/a)` tower. *Genuinely open hard analysis; the modest face of the Savvidy wall (one derivative, UV, fixed volume, `⊥` gap). No published route; Balaban/MRS stop here; the renormalon blocks the one-loop bridge.*
2. **G2 = H1** — the IR face: `Var_μ(a_{k→0})≤1/λ`, equivalently `f₁²>0` and `Δ>0` at physical coupling `β≈2.5`. *The gap itself (= Clay). Its strong-coupling end is proven (P3 unconditional at `β<1/48`); the crossover `[1/48, 2.5]` is the open wall.*

**And these two are one species.** So the honest final line:

> **`4D pure SU(N) YM exists (OS0/1/2/3/4)` is a theorem modulo `{L3, G2}`; `L3 = G2` is a single named unknown — the nonperturbative upper bound on a variance of a frustrated, non-attractive, non-abelian measure, which does not currently exist. Fixed-volume existence (OS0/1/3) needs only L3. Everything else in G1 — the entire OS0/1/2/3/4 and `L→∞` apparatus — is proven or reduces to that one species. The classification is 100% complete; the Clay theorem is 0%.**

---

## 9 · Honest boundary — what this assembly is and is not

- **It is** the cleanest possible statement: every G1 input is a theorem, reduces to L3, or reduces to G2; OS2 and OS4 add no new species; the modulo list collapses to one named species at two scales.
- **It is not** a proof of existence or of the gap. L3 is *true-but-unproven* (deep-UV ordering) and G2 is *the open problem*. The one-loop bridge to L3 is blocked by the asymptotic renormalon; the variance-upper-bound tool for G2 is provably absent from all nine known dependency-inequality frameworks.
- **No exaggeration of L3's payoff:** closing L3 alone gives fixed-volume existence (OS0/1/3, subsequential) — **not** `L→∞`, **not** OS4, **not** the gap, **not** uniqueness. OS2-at-`L∞` and OS4 are honest separate rungs, but both reduce to the same species (SYM-4D→L3; Q-CONS, IN-6, OS4→G2).
- **No claim that solving L3 solves G2:** they are the *same functional form*, not literally the same statement. The unification says one *new tool* might close both; it does not say one closure implies the other. That is a structural conjecture, stated as such.

*Refs (external, as checked in the source files).* Osterwalder–Schrader, CMP **31** (1973) 83, **42** (1975) 281; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Menotti–Pelissetto, CMP **113** (1987) 369; Lüscher, CMP **54** (1977) 283; Balaban, CMP **109** (1987) 249, **116** (1988) 1, **122** (1989) 175/355; Magnen–Rivasseau–Sénéor, CMP **155** (1993) 59/325 (DOI 10.1007/BF02097397); Bauerschmidt–Brydges–Slade, LNM **2242** (2019); Aizenman–Duminil-Copin, arXiv:1912.07973; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; Chandra–Chevyrev–Hairer–Shen, Invent. Math. **237** (2024); Reisz (lattice BPHZ); Dobrushin, TPA **13** (1968) 197; Elitzur, PRD **12** (1975) 3978. *(Corrected per charter: no "Chandra–Fröhlich" joint paper exists.)*
*Internal.* G1-EXISTENCE-THEOREM; G1-L3; G1-UNIFORM-ZA; G1-CONTINUUM-OS; G1-CONTINUUM-SCALING; G1-UNIQUENESS-VACUUM; G1-fixed-volume; PROOF-G1-BPHZ; A1-twosided-UV; A2-scheme-match; P1-L0-proof; PROOF-L3; L3-scale-uniform; L3-band-decay; PROOF-P3-UNCONDITIONAL; PROOF-ISOLATION; IN-2-attack; CONTINUUM-CP; PROOF-COARSE-POINCARE; docs/PROOF-SKELETON.md §0–3.
