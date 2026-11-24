# SIGNED-FOURIER-CORE — the representation-ring / Bochner species: its benign face *is* RP=KP, its signed face lives in the strong-coupling surface sum; both acid tests fail, no crossover lever

> 2026-11-23 · Core adjudication of the candidate 6th species "representation-ring / signed-Fourier (Bochner-type)": *a LOWER bound on a binding representation coefficient of the connected `tr F²` correlator, certifying coherence rather than disorder.* **Verdict: DEAD — score 1/5.** The species has exactly two faces, and each is a strong-coupling structure under a new name. (A) The **benign face** — the character coefficients `c_R(β)≥0` — is not merely *like* positivity: by the SU(N) Bochner–Godement theorem it **is** the statement that the Wilson weight is of positive type, i.e. the PSD/reflection-positivity kernel, i.e. verbatim the Kotecký–Preiss hypothesis that the polymer weights are non-negative. It is the same positivity as three of the nine dead families at once, and — decisively — it holds **strictly** for the free heat-kernel measure (`c_R=d_R e^{−tC₂(R)}>0`), where `f₁²=0`; so it is **vacuous** (free-field acid test FAILS). (B) The **signed face** — Weingarten sign-alternation — is exact but lives in the area-law **surface** sum at `β<β₁`, controls the mass gap (an upper bound on the correlator), **not** the per-tube residue, and has no theorem at `β≈2.5` (T5). There is **no** signed-Fourier lower bound on a binding coefficient off strong coupling that is not (A) or (B). The exact stuck step is unchanged: `H*` = a volume- and `(L,a)`-uniform lower bound `f₁²≥r₀>0` at `β≈2.5` that does not pass through a PSD/Bochner inequality and vanishes at `g=0`.

---

## §1 · The species, stated precisely — and exactly where the sign lives

**Setup.** `T⁴_L`, spacing `a`, Wilson action `S=βΣ_p(1−N⁻¹Re tr U_p)`, Haar measure `dU_e` per link. The observable is the gauge-invariant `0⁺⁺` interpolator `O(x)=tr F²(x)` (`=Σ_{A1g} o_p` a finite sum of one-cell plaquette operators, real character content `o_R` = fixed *group* data: singlet from the "`1−·`", fundamental/antifundamental and higher from `Re tr U_p` and its powers). The connected two-point function is the Källén–Lehmann Laplace transform
`C(t):=⟨O(t)O(0)⟩_c = ∫e^{−mt}\,dν(m²)`, `ν≥0` (OS reflection positivity), with the atom `f₁²=ν({m₁²})=|⟨Ω|O|G⟩|²`.

**Character expansion.** Expand each plaquette weight and each `O` in characters:
`e^{(β/N)Re tr U_p}=Σ_R c_R(β)χ_R(U_p)`, `c_R(β)≥0`, `c_R=O(β^{|R|})` (Drouffe–Itzykson; PROOF-P3 §1), so `u(β):=c_f(β)≍β→0` at strong coupling. Expanding every plaquette and the operator and integrating Haar with Weingarten/recoupling bookkeeping gives
`C(t)=Σ_{P} w(P)`, over **link-connected polymer clusters** `P` meeting both the `t=0` and `t=t` slices, with

  `w(P) = [∏_{p∈P} c_{R_p}(β)] · 𝕀(P) · [∏ o_{R_{p'}}]`,

`𝕀(P)` = the Haar weight of the polymer's dual surface (product of link recouplings / Weingarten functions `Wg(σ;N)`).

**The single flux tube (the binding polymer).** A connected `P` spanning `[0,t]` must occupy ≥ one plaquette per time slice (PROOF-P3 §3). At equality it is **one fundamental plaquette per slice**, consecutive plaquettes sharing a link with the singlet recoupling `f⊗f̄→1`, coefficient `1` (the normalized contraction `∫dU\,χ_fχ_f̄/\,=1>0`). Hence

  `w(P_tube) = o_f² · 1 · u(β)^{t/a} = a₀(β)\,u^{t/a}`, `a₀=o_f²∏(\text{positive recouplings})>0` (**PROOF-P1P2 §2–3, PROOF-P3 §3**).

So the **leading connected (binding) amplitude is a strictly POSITIVE product**: positive character coefficients `c_f>0`, positive recoupling `+1`, fixed real `o_f²`. **There is no sign in the term that carries the atom.**

**Exactly where the sign lives.** Two places, neither of them the per-tube residue:
1. **Operator coefficients `o_R`,** real but *signed* across `R` (the singlet "`1`" vs the `−Re tr U_p` content); connectedness subtracts the singlet, and the surviving leading tube content is positive (§ above). A *relative* sign, not an oscillatory one.
2. **The Weingarten function `Wg(σ;N)` inside `𝕀(P)` for non-tube polymers** (sheets, branched/merged tubes, multi-plaquette surfaces): `Wg` **alternates** in the permutation cycle structure and in `N`. This is the *only* genuine sign-alternation. It enters `w(P)` only at cost `≥ u^{t/a+1}` (decoration/wider tube) and, in the full surface sum, is the area-law sign structure. **PROOF-P1P2, sharpest line:** *"SU(N) Weingarten sign-alternation lives in the area-law SURFACE sum, not in the per-tube character expansion (positive for β>0)."*

That single placement kills the species' premise: the "signed-Fourier structure" is either **positive** (the per-tube/Cartesian-character face, `c_R≥0`) or **in the surface sum** (the Weingarten face). The hope — a *signed* lower bound on a *binding* coefficient — has no home: the binding coefficient is positive, and the signed coefficients are at higher cost, on the surface, where they do not bound the tube.

---

## §2 · Attempting the coefficient lower bound off strong coupling

The species needs `∃ε>0:⟨Ω|O|G⟩≥ε` (equivalently `f₁²≥ε²`) at `β≈2.5`, derived *without* passing through measure positivity. What does the representation ring give?

**(a) The recoupling sign gives only `a₀>0`, i.e. `⟨Ω|O|G⟩≥0`, not `≥ε`.** The positive-type/recoupling structure is a statement about the *sign* (the glueball is coherent, not destructive), and a sign lemma is not a magnitude lemma. The magnitude `a₀` is the residue of a flux tube = a *dynamical* (interacting) quantity: it is `e^{−(\text{string self-energy})}` times the endpoint overlaps, controlled by the confining area-law coefficient, which is **not** a group datum. Representation theory fixes the sign of `⟨Ω|O|G⟩` — it does not lower-bound its modulus. Lower-bounding the modulus *is* the nonperturbative dynamics.

**(b) The atom is a Schur-concave functional; every representation-space total is an upper bound (atom ≤ total).** The intertwiner inner product is `ν≥0`; `f₁²` is majorization-**minimal**. Any quadratic/Cauchy–Schwarz/Bochner functional `Φ(ν)≥0` is **maximized by spreading** (moving weight off the atom into a continuum) and cannot force or forbid the atom. This is the Schur-convex sign wall (SCHUR-CONVEX; EDGE-f₁²-attack §"unifying obstruction"; H1-T9 §"lower-bound f₁²"). The needed object is the **reverse**: a *total bounded above by the atom* (`∫φ\,dν ≤ C(φ)ν({m₁²})`), which no representation-space positivity supplies. The representation ring gives selection rules (it can force `f₁²=0` in the wrong channel) — never `f₁²≥r₀`.

**(c) Even the time-Bochner (Laplace) route cannot isolate the atom.** `C(t)=∫e^{−mt}dν` makes `C` a Laplace transform of a positive measure, i.e. **completely monotone** (Bernstein). But an atom is a *jump* in `ν`; a purely absolutely-continuous `ν_ε` (smooth bump of width `ε` at `m₁`, mass `f₁²`) reproduces `C(t)` to `O(e^{−m₁t}·\text{width})` with **no atom**. Complete monotonicity (`C>0`, `(−1)^nC^{(n)}≥0`) is satisfied by both the atomless and the atomic case; Laplace inversion is exponentially ill-conditioned and **cannot** yield a discrete-component lower bound. Positivity (Bochner) is blind to singular-vs-absolutely-continuous — exactly the free-field (`ν` a two-photon cut) vs interacting distinction that constitutes the whole problem.

**Outcome of §2:** off strong coupling there is **no** coefficient lower bound `≥ε`. The only positive lower bound the character structure yields at any accessible regime is `f₁²≥0` — the free-field value. The bound is not missing a constant; it is missing its direction (a Schur-convex coercivity) and its magnitude (the dynamical tube residue).

---

## §3 · The two acid tests — both fail, decisively

### Acid test 1 — is it free of measure-positivity? **NO (benign face is literally RP=KP).**
The claim of the species is that the sign structure is "NOT a measure-positivity statement." For the **benign face it is false.** By the compact-group **Bochner–Godement theorem** (Bochner 1933; Godement 1948; Folland): a continuous central function `W` on `SU(N)` is **of positive type** (the PSD kernel `Σ_i c_i\bar c_j W(g_i^{-1}g_j)≥0`) **iff** its character coefficients are all `≥0` in the Peter–Weyl expansion `W=Σ_R a_Rχ_R`. The Wilson weight is exactly such a function, and its `c_R(β)≥0` (PROOF-P3 §1; Drouffe–Itzykson; PROOF-HEAT-KERNEL). So

  `c_R(β)≥0` ⟺ `e^{(β/N)Re tr U}` **positive-definite (Bochner)** ⟺ polymer weights `w(P)≥0` (**KP hypothesis**) ⟺ PSD transfer kernel (**RP**).

These four are *one statement*. The "benign face" is simultaneously measure positivity, reflection positivity, and the Kotecký–Preiss positivity hypothesis — **three of the nine dead families, identified**. It degenerates; **Acid Test 1 FAILS.**

### Acid test 2 — does it vanish on the free field? **NO (the benign face is strictly positive for the free field, where f₁²=0).**
This is the clean kill, and it is specific to this species. The **free** (gauge) field, in the heat-kernel regularization, is
`K_t(U)=Σ_R d_R\,e^{−tC₂(R)}χ_R(U)`, with `d_R>0`, `C₂(R)≥0` ⇒ **strictly positive character coefficients** `c_R^{free}=d_R e^{−tC₂(R)}>0` (PROOF-HEAT-KERNEL §"exact single-plaquette expansion"). So the free measure is **strictly Bochner-positive / positive-type** — *more* positive than the interacting one at small `t`. Yet the free Maxwell field has `ν` a two-photon **cut on `[0,∞)`**, purely absolutely continuous ⇒ **`f₁²=0`** (SCHUR-CONVEX "master obstruction"; PROOF-F1SQ). Therefore

  **Bochner positivity ⇒ `f₁²>0` is FALSE** (the free heat-kernel measure is a strict counterexample).

Any "lower bound" derived from the positivity of the character coefficients returns `>0` on the free field, where the true answer is `0` ⇒ **VACUOUS**. The free-field discriminator — `f₁²(free)=0` — is exactly what the benign face fails. **Acid Test 2 FAILS.**

The **signed face** (Weingarten) *does* pass the letter of Acid Test 2 (its effect — the area law — vanishes at `g=0`), but it fails on the *object*: it controls the surface/area-law (an upper bound on the correlator ⇒ a *mass-gap* lower bound), **not** the per-tube residue `f₁²`, and only at `β<β₁` (§4). Passing Acid Test 2 on the wrong object is not survival.

---

## §4 · The crossover attempt (Savvidy region) — why it reduces to KP

At `β≈2.5` the activity `u(β)` reaches `O(0.1–0.4)` and `u·K≳3–40≫1` (`K≤CN`), so the character/polymer expansion is **outside its absolute-convergence radius** (W2; T5). Two faces again:

**(A) Benign face.** The character coefficients are **non-negative for all `β>0`**, so the divergence is a **real-activity (Pringsheim) barrier**: the singularity sits on the positive real axis and there is **no sign cancellation to exploit**. The convergence radius cannot be crossed by any reorganization of the *positive* series (W2, "double-edged, no cancellation loophole"). In other words, the benign face at the crossover **is** the Kotecký–Preiss expansion, and it has *exactly* the strong-coupling wall. There is no "boost" from representation-ring structure — the positive-coefficient face is KP verbatim. **Reduces to KP; KP does not cross.**

**(B) Signed face.** The only route past the real barrier is sign-changing weights, i.e. the **blocked/Savvidy effective action** (negative Hessian, W4/W6). But: (i) the signed measure is an artifact of the RG block, **not the physical real-β measure** (T5 "Reframing": real-β pure YM has no phase to complexify); (ii) no theorem controls the blocked resummation — thimbles/CL are heuristics, resurgence is formal, signed KP still closes only at `β<β₀` (T5 verdict); (iii) the Weingarten signs themselves remain a **strong-coupling surface** datum (large-`χ`/area-law series), so they too live on the `β<β₁` side. The spin-foam (`6j`) signed structure is, by the master's own registry, **one of the five faces** of the single non-abelian-sign no-go (Savvidy / `6j` / Weingarten / random-cluster / `H_ff`). **Signed face = the same wall in surface/spin-foam clothes.**

**Three forge attempts (rule 1 — actively build, don't stop at "hard").** Each lands on an existing no-go:
- *Loop-space / Migdal duality total-positivity:* 4D has no closed duality; the 3D `6j` cascade with controlled signs does not lift to 4D; the 4D spin-foam (`Barrett–Crane`/EPRL) amplitude is oscillatory with an uncontrolled "flatness" phase ⇒ the `6j` face of the non-abelian-sign no-go.
- *Transfer-matrix Perron–Frobenius coherence lower bound:* PF is **proved** (PROOF-GAUSS-INVARIANT-PF) but pins only the **vacuum** (top eigenvector); the glueball is the *second* eigenvector in the *same* `0⁺⁺` sector, which PF (a "top-of-sector" statement) does not pin ⇒ the existing PF-limitation no-go.
- *Strict Bochner-interior gap (`c_R≥ε`):* the free field sits on the **cone boundary** (`f₁²=0`) yet is Bochner-positive; distinguishing "interior by `ε`" from "boundary" in a family of positive-type functions is a gap no positivity argument sees ⇒ the Schur-convex master obstruction again.

No forge survives; all three are catalogued no-gos.

---

## §5 · Scorecard, sharpest line, exact stuck step

**Conditions (of 5):**
1. **Nonperturbative — PARTIAL (½).** The recoupling sign `f⊗f̄→1=+1` and Weingarten alternation are exact, nonperturbative group facts (acid-test-clean *as signs*). No nonperturbative *bound*.
2. **Sign-controlling LOWER bound `m_eff²(k→0)≥λ` — FAIL.** The benign face gives `≥0` (free-field value); the signed face gives a *mass* (upper-bound-on-correlator) result at `β<β₁`, and *never* a residue lower bound at `k→0`, `β≈2.5`.
3. **Non-circular — FAIL.** To turn either face into `f₁²≥ε` at the crossover needs the margin/geometry = `H1`/(CP) (the activity must be made small = the gap), or presupposes `f₁²>0` (sharp cone, INVENT-SECOND-CONE). Circular.
4. **Non-SUSY + fails free field — FAIL.** Benign face does **not** fail the free field (free heat kernel is *strictly* Bochner-positive, `f₁²=0`) ⇒ vacuous (Acid Test 2). Signed face passes the letter but on the wrong object/regime.
5. **`(L,a)`-uniform — FAIL.** KP is `L`-uniform only at `β<1/48`; no `a→0`; the atom `f₁²` is usc and may drop to 0 at the crossing (PROOF-TRANSPORT-ANALYTIC K1).

**Score: 1/5** — *one point* for the genuine, exact, acid-test-clean **sign** structure (the positive recoupling and the real Weingarten alternation, analogous to BS-KERNEL-SIGN's surviving color-Casimir lemma, 2/5); everything that could convert that sign into a **bound** is either the free-field-shared Bochner positivity (vacuous) or the strong-coupling surface sum (wrong regime and wrong object). Strictly weaker than BS-KERNEL-SIGN because its usable face fails the free-field discriminator outright.

**Sharpest line.** *The representation-ring species is the strong-coupling residue theorem wearing a Fourier hat.* Its benign face — `c_R(β)≥0` — is, by Bochner–Godement, **the same single positivity** as reflection positivity and the Kotecký–Preiss hypothesis (one statement, three dead families); and that positivity is held **strictly** by the free heat kernel `Σd_Re^{−tC₂(R)}χ_R` whose `f₁²=0`, so it cannot be the lever. Its signed face — Weingarten alternation — is real but lives in the area-law **surface** sum at `β<β₁`, bounding the *mass*, not the per-tube *residue*, and dies at the crossover (real Pringsheim barrier; no theorem past `β₀`, T5). "Coherence not disorder" is precisely the **positive recoupling `f⊗f̄→1=+1`** — which is already inside P3's positive polymer sum. The species adds no sign the positive-measure route does not already have, and its only genuine sign (Weingarten) sits where it cannot lower-bound the atom.

**Exact stuck step.** `H*`: a volume- and `(L,a)`-uniform **lower** bound `f₁² = lim_{t→∞}e^{m₁t}C(t) ≥ r₀(β)>0` at `β≈2.5`, equivalent to a Schur-**convex** coercivity `∫φ\,dν ≤ C(φ)\,ν({m₁²})` (total ≤ atom), that does **not** pass through a PSD/Bochner inequality (or it is vacuous on the free field) and that **vanishes at `g=0`** (or it returns `>0` on the free two-photon cut). The representation ring supplies neither the coercivity direction (it is Schur-concave: atom ≤ total only) nor the magnitude (the flux-tube residue is dynamical, not group-theoretic). Transfer to `β≈2.5` at the crossover is the same analyticity/real-barrier wall the whole campaign has hit (`H1`/`(CP)`). Species DEAD; the wall is unchanged and now has one more precisely-registered face.

*Refs.* S. Bochner, *Monotone Funktionen, Stieltjes Integrale und harmonische Analyse*, Math. Ann. **108** (1933) 378; R. Godement, *Les fonctions de type positif et la théorie des groupes*, Trans. AMS **63** (1948) 1; Folland, *A Course in Abstract Harmonic Analysis* (Bochner–Godement for compact groups); D. Weingarten, J. Math. Phys. **19** (1978) 999; Drouffe–Itzykson, Phys. Rep. **38C** (1978) 133 (positive character coefficients); Drouffe–Zuber, Phys. Rep. **102** (1983) 1; Bernstein (complete monotonicity / Laplace positivity); Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Kotecký–Preiss, CMP **103** (1986) 19; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737. Internal: PROOF-P1P2 (per-tube positive vs Weingarten surface); PROOF-P3-UNCONDITIONAL (P3 residue I1, `β<β₁`); PROOF-HEAT-KERNEL (free heat kernel `c_R>0`; character gap = fiber, not `0⁺⁺`); H1-T9-spin-network (representation reformulation = character expansion; atom ≤ total); SCHUR-CONVEX (master free-field obstruction; `f₁²` Schur-concave); EDGE-f1²-attack (Schur-convex coercivity missing; sign obstruction); BS-KERNEL-SIGN (color-Casimir sign lemma, 2/5, wrong direction); W2-positive-scheme (`c_R≥0`, Pringsheim real barrier, no cancellation); T5-signed-weight (no theorem at `β≈2.47`); PROOF-F1SQ (free `f₁²=0`); PROOF-GAUSS-INVARIANT-PF (PF pins vacuum, not same-sector 2nd state); INVENT-SECOND-CONE (sharp cone = `f₁²>0`); PROOF-TRANSPORT-ANALYTIC (K1, `f₁²` usc); INEQ-T2 (residue-dual is an infimum); master §4 (non-abelian-sign five faces: Savvidy/`6j`/Weingarten/random-cluster/`H_ff`).
