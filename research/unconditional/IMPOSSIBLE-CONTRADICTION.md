# IMPOSSIBLE-CONTRADICTION — assume pure YM is GAPLESS and derive a contradiction from exact properties: **no contradiction is obtained; every route pins the SAME exact gap** — "gapless interacting 4D IR ⇒ conformal" (route i, the sharpest) / the a-theorem excludes neither branch (ii) / positivity is lsc (iii) / the anomaly is at θ=π not θ=0 (iv) / the renormalon IS H1 (v). The contradiction-side re-derives H1.

> **USER DIRECTIVE.** Jump out of the fixed mindset (KL ν, atom f₁², χ, weak-*, H1). Take a base that looks
> **completely impossible** and take it seriously: **assume 4D pure YM is GAPLESS, derive a contradiction.**
>
> **VERDICT (no QED, Clay 0%).** Five independent contradiction routes were pushed as hard as they go.
> **None yields a contradiction.** Each pins the *exact* failing step, and the steps are not independent:
> the sharpest (i) reduces the whole thing to one folklore step — **"a gapless interacting 4D Wilsonian IR
> attractor is a conformal field theory"** — after which the exact β-function of pure YM (no Banks–Zaks zero
> at `g>0`) would force a contradiction. That folklore step is **not a theorem**, and closing it is exactly as
> hard as the gap. The contradiction-side therefore **re-derives `H1`** (the question "gapped vs gapless" is
> literally "is `⟨TrF²⟩` `=0` or `>0` in the IR"). This is the honest outcome: the impossible base does not
> break the wall — it relocates the wall onto the one step no tool closes.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** witness. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The five exact properties and the exact gap in each

**Setup.** "Gapless pure YM" = the `T=0`, infinite-volume IR has a continuous spectrum down to zero: an
interacting phase with no mass. In 4D, Lorentz-invariant, unitary, with a stress tensor, such a phase is
*believed* to be a CFT (scale- then conformal-invariant). The contradiction attempt tests whether the **exact**
YM properties forbid it.

### (i) Trace anomaly + β-function — the sharpest route; gap = "gapless ⇒ CFT"
- **Exact:** pure SU(N) YM is asymptotically free; at 2-loop `β(g) = −β₀g³ − β₁g⁵` with `β₀, β₁ > 0`, so
  `β(g) < 0` for **every** `g>0` — **no Banks–Zaks zero** at finite coupling (only the UV fixed point `g=0`).
  Witness R1: `β<0` on all of `g∈(0,8]` for `N=2,3,4`; zero-at-`g>0` = `False`.
- **If** the IR is a CFT, it has an exactly-marginal coupling `g*` with `β(g*) = 0`, `g*>0`. Pure YM has none
  ⇒ **contradiction** — the IR cannot be a CFT at finite coupling; and it cannot flow to `g=0` either
  (that would make `g` *decrease* toward the IR, reversing asymptotic freedom).
- **Exact gap `[O]`:** the premise "gapless interacting ⇒ CFT" is **not** established. Gaplessness (no mass)
  does **not** imply scale invariance — a gapless phase can, in principle, be non-conformal (a continuously
  running coupling that never hits a fixed point). The step "gapless ⇒ conformal" is **folklore, not a
  theorem** (the proven statement is the converse direction: scale invariance + unitarity ⇒ conformal,
  Polchinski 1988; and that needs *scale* invariance, not just gaplessness).

### (ii) a-theorem / monotone effective central charge — gap = it excludes NEITHER branch
- **Exact:** `⟨T^μ_μ⟩ = a E₄ − c W² + (β/2g)Tr F²`; Komargodski–Schwimmer prove `a_UV > a_IR` along any
  non-trivial 4D flow. Pure YM UV = free Maxwell, `a_UV = 62(N²−1)/360 > 0`.
- **Gapped IR:** `a_IR = 0`. **Gapless CFT IR:** `0 < a_IR < a_UV`. **Both** satisfy `a_IR < a_UV` strictly.
  Witness R2: both branches pass for `N=2,3`. The theorem **excludes neither**. (This reproduces
  `A-THEOREM.md` §S1 — the a-anomaly is a one-sided **upper** bound on a *dimensionless* number, the wrong
  direction for a mass.)
- **Exact gap `[O]`:** forcing "gapless" out would need `a_IR ∉ (0, a_UV)`, i.e. the IR is *not* a CFT — which
  is the gap itself, assumed.

### (iii) Unitarity + Källén–Lehmann positivity — gap = positivity is lower-semicontinuous (blind)
- **Exact:** the spectral measure is positive, `ν ≥ 0`; any gapless continuum is a perfectly valid positive
  measure. Positivity certifies a state *exists*, never that it is *isolated above zero*.
- **Exact gap `[O]`:** positivity is lsc (`PNT-OPE-RHO`): a bump near the floor matches every finite datum.
  A gapless measure is as positive as a gapped one. Unitarity forbids negative norm, not gaplessness.

### (iv) 't Hooft anomaly matching — gap = the anomaly is at θ=π, trivial at θ=0
- **Exact:** the 4D `Z_N` **1-form center symmetry** has a mixed anomaly with the θ-angle, but it is a
  `θ=π` phenomenon; at `θ=0` the anomaly vanishes, so a trivially-gapped symmetric phase is **allowed at
  θ=0** (and so is a gapless one). `ROOF-CONFINEMENT-GAP`/`PROOF-1FORM-ANOMALY` established this.
- **Exact gap `[O]`:** at `θ=0` anomaly matching gives **nothing**. (At `θ=π` the anomaly forces
  *deconfinement*, not gaplessness — a different phase.)

### (v) OPE / renormalon / scale structure — gap = circular (the renormalon IS H1)
- **Exact:** the `TrF²` OPE has the `u=2` IR renormalon; a scale-invariant IR would have a scale-free OPE
  (no condensate). The gluon condensate `⟨TrF²⟩` is the one dimension-4 scale of YM; a CFT has `⟨TrF²⟩=0`
  (dimension-4 vev, no scale), the confined phase has `⟨TrF²⟩>0`. Witness R3: the exact trace-anomaly Ward
  identity in flat space reduces to `(β/2g)⟨TrF²⟩ = 0`; in a CFT `β=0` makes this hold for **any**
  `⟨TrF²⟩`, including 0 — **no exact relation forbids `⟨TrF²⟩=0`.**
- **Exact gap `[O]`:** "`⟨TrF²⟩=0` vs `>0`" **is** "gapless vs gapped". The OPE/renormalon is the `H1`
  coordinate verbatim (`ROOF-CONDEMSATE-FLOOR`, `ROOF-RESURGENCE`: `σ` free). Circular.

---

## §2 · The sharpened crux (what the impossible base actually produces)

The five routes collapse to **one** possibility: if the IR is gapless and *interacting*, its Wilsonian
attractor is a fixed point — a CFT — and route (i) would then give a genuine contradiction (no `BZ` zero). So
the entire "impossible base" reduces to the single sharp target:

> **(★)** **A gapless interacting 4D Wilsonian IR attractor of an asymptotically-free gauge theory is
> scale-invariant (hence, by unitarity, conformal).**

If (★) is a theorem, **pure YM is gapped** — route (i) closes it. (★) is the exact step that is currently
*folklore*, not proof. It is a statement about the **space of 4D RG attractors**, not about any measure or
observable, so it sits genuinely *outside* every failed base — this is the one place the contradiction
direction is not immediately circular.

Two honest defeats of (★) that must be checked before trusting it: a gapless **non-conformal** 4D attractor
would have to exist (none is known; all known 4D gapless phases are CFTs or free), and the argument assumes
the RG flow reaches *an* attractor (no limit cycles). Closing (★) rigorously would be a Clay-scale theorem
in its own right (it is, essentially, "the 4D IR of YM is not a non-conformal interacting phase" = `¬(iii)`).

---

## §3 · Additional "impossible" bases tried (and their gaps)

- **Triviality / UV Landau-type:** YM is asymptotically free, not trivial; no IR Landau pole. The triviality
  argument uses the *wrong* end (UV). `[NG]`.
- **S-matrix unitarity + Froissart + asymptotic freedom:** AF sets the UV S-matrix → 1; a gapless CFT IR
  S-matrix is unitary and consistent; nothing in the bounded S-matrix forbids the CFT. `[O]`.
- **Compact gauge group + finite-volume discreteness:** the finite-volume spectrum is discrete; gaplessness
  is an infinite-volume phenomenon — consistent, no contradiction. `[O]` (this is the `a→0, L→∞` limit = H1).
- **Vacuum-energy / pressure sign:** a gapless CFT has `p = ρ/3`; the confined phase has `ρ<0<...`; both are
  thermodynamically allowed. `[NG]`.

None closes.

---

## §4 · Honest bottom line

- **The impossible base does not break the wall.** "Gapless pure YM" is **not** in contradiction with the
  exact trace anomaly, the β-function, the a-theorem, unitarity/KL positivity, θ=0 anomaly matching, or the
  OPE structure — because each either **excludes neither branch** (ii, iii), **is trivial at θ=0** (iv),
  **is the gap itself** (v), or **hinges on one unproven folklore step** (i): `gapless interacting 4D IR ⇒
  conformal`.
- **What the contradiction direction *does* produce:** a **single sharp target (★)** — prove the 4D gapless
  IR attractor is conformal. This is genuinely *outside* the failed bases (it is a statement about RG
  attractors, not a measure/observable), and if (★) holds, route (i) closes pure YM as gapped. It is the one
  non-circular residue of the impossible base.
- **Honest scope:** (★) is not proved here and is not obviously easier than the gap; it is `¬(iii)` — "the IR
  is not a non-conformal interacting phase" — i.e. a re-labelling of the wall, not a demolition of it.
- **Clay = 0%.** No QED. The contradiction attempt re-derives the same invariant from the other side.

*Refs:* `A-THEOREM.md`, `PROOF-1FORM-ANOMALY.md`, `ROOF-CONFINEMENT-GAP.md`, `PNT-OPE-RHO.md`,
`ROOF-CONDEMSATE-FLOOR.md`, `ROOF-RESURGENCE.md`; Komargodski–Schwimmer JHEP 12 (2011) 099 [1107.3987];
Polchinski, *Scale and conformal invariance in QFT*, Nucl. Phys. B303 (1988) 226; witness
`w7/impossible_contradiction.py` (runs clean, `py -3.12 -I`).
