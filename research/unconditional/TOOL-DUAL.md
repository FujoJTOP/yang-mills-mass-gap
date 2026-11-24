# TOOL-DUAL — convex / Legendre–Fenchel / Nehari duality CANNOT force the KL atom: strong duality is EXACT and the atomless measure is feasible, so primal inf = dual optimum = 0. The blind spot is CONVEX-STRUCTURAL, not a tool-choice artifact

> **VERDICT (no QED).** The "tool blind spot" is that no positivity/integral bound forces
> `f₁² = ν({m₁²}) > 0`. This fork asks whether **convex duality** — a genuinely non-integral
> (Lagrange / Legendre–Fenchel / Nehari) route — can force it instead. **Answer: no, and for a sharp
> reason.** The feasible set (positive measures matching `M, χ, m₁²`, and any finite moments) is **convex and
> contains the atomless measure** (all mass at `s = M/χ`), so the primal **infimum of the atom is 0**;
> by **exact LP duality** (strong duality holds, gap = 0) the **dual optimum is also 0**, and the dual
> multipliers vanish. **Duality does not create information — it is an equivalence — so it inherits the
> blindness exactly.** The Nehari/AAK analytic dual only *reads* the atom (the Stieltjes pole = A3,
> read-not-force). Forcing requires a **nonlinear / structural** constraint (the actual YM measure) or the
> `J → ∞` (full Stieltjes) limit — i.e. `H1`. **The blind spot is convex-structural; no dual trick removes it.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric.
Acid test: a forcing bound must be `0` on free massless *and discriminatingly* `0` on `:φ²:` (vacuously `0` everywhere is **vacuous**).

---

## §1 · The primal (measure LP) — infimum is 0 `[T]`

`ν ≥ 0` on `[m₁², ∞)`, atom `a := ν({m₁²})`, data `M = ∫dν`, `χ = ∫s⁻¹dν`, floor `m₁²`. Primal: `inf a`.

**Atomless point is feasible.** Put all mass at `s₀ = M/χ`: then `∫dν = M`, `∫s⁻¹dν = M/s₀ = χ`, and
`a = 0`. Feasible whenever `s₀ = M/χ ≥ m₁²`, i.e. `χ ≤ M/m₁²` (the data-feasibility of Lemma G). Hence
`inf a = 0`. The witness solves the discretised LP for `M=1, χ=0.5, m₁²=1`: **primal `= 0`**.
Adding the trace-anomaly moment `∫s⁻²dν = A` keeps it `0` (a spread-out atomless bump matches all three
constraints). **Any finite set of linear moments keeps it `0`** — this is the `Δ_min = 0 ∀K` theorem
(`SLOT-BOOTSTRAP`, `PNT-OPE-RHO`) in LP language.

---

## §2 · The dual (Legendre–Fenchel) — optimum is 0, because duality is exact `[T]`

The dual of `min cᵀw  s.t. Aw = b` is `max bᵀy  s.t. Aᵀy ≤ c`. **Strong LP duality holds** (feasible,
bounded), so `dual optimum = primal optimum = 0`. The witness reads the equality marginals: **`y ≈ (0,0)`**
— the optimal dual multipliers vanish, so **no finite dual witness pushes the atom up**. This is the crux:
**the duality gap is 0, and that is exactly why no dual trick can work.** A duality "forcing" bound would
need a *positive* duality gap — but for this convex measure program there is none. Duality is an
**equivalence**, so it cannot manufacture the atom any more than the primal integral functionals can.

---

## §3 · Acid test — the dual bound is vacuous `[NG]`

The dual lower bound is `0` **everywhere** (free massless, `:φ²:`, interacting YM alike). It is therefore
`0` on the acid adversaries only **vacuously** — it is not *discriminatingly* acid-clean. A genuine forcing
tool must be `> 0` on YM and `0` on both adversaries; the dual bound cannot do this, because `0` is the
primal infimum over the whole feasible class.

---

## §4 · Nehari / AAK dual — reads (A3), does not force `[T]`

The non-integral analytic object is the **Stieltjes/Herglotz transform** `Π(z) = ∫ dν/(z−s)`, and
`f₁² = −Res_{z=−m₁²} Π`. In **Nehari / AAK** terms the moment Hankel operator `H` (`H_{ij} = ∫s^{i+j}dν`)
is PSD (positivity of `ν`), and each **atom of `ν` is a rank-one Hankel perturbation** — a spike in the
Schur parameters (reflection coefficients), i.e. a pole of `Π`. The Nehari distance / AAK singular value is
`0` for an **atomless** `ν` and equals the atom magnitude when an isolated atom is present. So the analytic
dual **reads** the atom (this is precisely **A3, read-not-force**) — it gives **no lower bound** from the
integral data, only a detector whose input is the full `Π` (the continuum eigenstates). Nehari's `L^∞`-distance
duality is exact and `0` on the atomless side, like §2.

---

## §5 · What forcing would actually require `[O]`

The only ways to make duality force `f₁² > 0`:

1. **A nonlinear / structural constraint** that the atomless measure violates — i.e. "`ν` is the KL measure of
   4D pure YM", not an arbitrary positive measure. That constraint *is* the nonperturbative theory = `H1`.
2. **The `J → ∞` limit** (Legendre-Fenchel of the free energy w.r.t. a source `J`): the max-entropy /
   exponential-family solution is **continuous for every finite `J`**; an atom is recovered only in the
   `J = ∞` limit, i.e. from the **full Stieltjes transform** — again A3, again `H1`.

So the dual escapes only by becoming the exact (nonconvex, full-`Π`) object — which is the wall.

---

## §6 · Bottom line

- **[NG]** Convex / Legendre–Fenchel / Nehari duality **cannot force** the KL atom: primal inf = dual
  optimum = 0, multipliers vanish, the atomless measure is feasible.
- **[T]** **Strong duality is exact (gap = 0)** — *this is the reason*: an equivalence cannot create the atom.
  The blind spot is **convex-structural**, not a defect of the particular positivity/integral tool chosen.
- **[T]** Nehari/AAK only *reads* the atom (Stieltjes pole = A3, read-not-force); it is `0` on the atomless side.
- **[O]** Forcing needs a nonlinear/structural constraint (the YM measure) or the `J→∞` full-`Π` limit = `H1`.
- **The tool blind spot is not fillable by duality.** It is a property of the convex class of positive
  measures; the atom lives in the *boundary/nonconvex* structure that only the actual YM measure (or the
  full analytic `Π`) carries. **Clay 0%; no QED.**

*Files:* `research/unconditional/w7/tool_dual.py` (witness, `py -3.12 -I`, clean).
*Refs:* `PNT-OPE-RHO.md` (integral blind spot), `SLOT-BOOTSTRAP.md` (`Δ_min=0 ∀K`), `PROOF-F1SQ.md`/A3,
`WI-ESCAPE.md` (A3 coordinate), `CONSTRUCT-SYNTHESIS.md`.
