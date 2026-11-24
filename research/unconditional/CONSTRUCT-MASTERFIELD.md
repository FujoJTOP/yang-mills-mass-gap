# CONSTRUCT-MASTERFIELD — building the 4D master field: the Gaussian/Hartree closure reproduces the DSE wall, but a **rigid (extrinsic-curvature) string functional stabilises the neutral 0⁺⁺** that pure Nambu–Goto cannot, giving an explicit `m_{0⁺⁺}/√σ = 2π√(2κ)` with one missing input `κ`

> **VERDICT (construction, no QED).** Two closures of the large-N 4D Makeenko–Migdal loop equation `W(C)` were built.
> **(i) Gaussian/Hartree** master field → reduces to the Dyson–Schwinger gap equation: the gapless solution
> `m=0` **persists** at the crossover (the known `SLOT-DSE` wall). **(iv) Rigid surface-sum** (a string
> functional with an extrinsic-curvature term) → the **neutral, contractible 0⁺⁺ closed loop acquires an
> interior minimum** — which **pure Nambu–Goto provably lacks** (`PROOF-EFFECTIVE-STRING`,
> `TRANSMUTE-STRING-SPECTRUM`) — and the mass is explicit:
> **`m_{0⁺⁺}/√σ = 2π√(2κ)`**, with rigidity `κ`, at loop radius `ρ*=√(κ/2σ)` (`L*≈1.9/√σ`, the glueball
> scale). Fitting the lattice band `m/√σ∈[3.3,3.9]` gives `κ∈[0.138,0.193]`. **The one missing input is
> `κ` (the flux-tube rigidity / core thickness) = `H1` in string coordinates.** The rigidity turns the
> corpus's "no interior minimum, missing = core thickness, no model" into an **explicit one-parameter
> relation** — a genuine, if modest, advance on the construction line. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric witness. Units: `σ=1` (so `√σ=1`); `κ` is the dimensionless coefficient of `∫K²ds`.

---

## §1 · The two closures of the loop equation

The large-N MM equation `∂_μ(δW/δσ_{μν}(x)) = λ∮_C dy_ν δ(x−y) W(C_{xy})W(C_{yx})` is **closed in `W(C)`**
(no n-point tower, `PNT-LARGEN-LOOP` `[T]`). Two ways to *solve* it:

**Closure (i) Gaussian / Hartree.** `W(C)` from a self-consistent propagator `D(k)=1/(k²+m²)` with `m²` from
its own tadpole — i.e. the large-N Hartree master field. Its gap equation is the rainbow DSE.
**Result: the gapless solution `m=0` is always present; a positive branch exists only if `F'(0)>1`** — the
`k→0` tachyon **absent at β=2.5**. Reproduces the `SLOT-DSE` wall (witness §3). **This closure does not
stabilise a mass** `[NG]`.

**Closure (iv) Rigid surface-sum.** `W(C) = Σ_{surfaces ∂=C} exp(−σ·Area − κ∮K² ds)`, i.e. the
Makeenko–Migdal/master-field weight **completed by an extrinsic-curvature (rigidity) term** — the natural
reflection of a flux tube with finite thickness. A **circular neutral loop** (radius `ρ`, `L=2πρ`, N-ality 0,
`0⁺⁺` by symmetry) has energy
`E(ρ) = 2πσρ + πκ/ρ` (tension + bending). **This is the object built below.**

---

## §2 · Pure NG has no interior minimum; rigidity creates one `[T]` + `[V-num]`

- **Pure NG** (`κ=0`): `E₀(L)=σL√(1−π(D−2)/(3σL²))`, `D=4`. Witness: minimum sits at the **left edge**
  `L→L_c⁺=√(2π/3)=1.447`, monotone beyond — **no interior minimum** (confirms `PROOF-EFFECTIVE-STRING` /
  `TRANSMUTE-STRING-SPECTRUM`: the neutral closed string is the shrink modulus, `√σ` is a density, no mass).
- **Rigid** (`κ>0`): `E(ρ)=2πσρ + πκ/ρ`. Interior minimum at `ρ*=√(κ/(2σ))`,
  `E_min = 2π√(2σκ)` ⇒ **`m_{0⁺⁺} = 2π√(2σκ)`**. The bending term `πκ/ρ` **resists the shrink** that pure
  NG cannot. The stabilised loop is **neutral and `0⁺⁺`** (a contractible circular loop) — so this sidesteps
  the corpus's objection that the *√σ-stabilised* string is the *charged* torelon: here the stabiliser is
  the tube's own rigidity, not `σL`, and the loop stays neutral.

`[V-num]` table (σ=1):

| `κ` | `ρ*` | `L*=2πρ*` | `m_{0⁺⁺}/√σ` |
|---|---|---|---|
| 0.100 | 0.224 | 1.405 | 2.81 |
| 0.138 | 0.262 | 1.65 | **3.30** |
| 0.193 | 0.310 | 1.95 | **3.90** |
| 0.300 | 0.387 | 2.433 | 4.87 |

The glueball scale `L*≈1.7–1.95/√σ` falls out at the physical band — a consistency check, not an input.

---

## §3 · The construction's content: an explicit one-parameter relation

**Constructed object (best available approximation to `W(C)`):**

> `W(C) = Σ_{surfaces} exp(−σ·Area − κ∮K² ds)` — the **rigid-string master field**.
> Its neutral-`0⁺⁺` closed-loop saddle has **`m_{0⁺⁺} = 2π√(2σκ)`**, i.e. **`m/√σ = 2π√(2κ)`**.

This **is an advance over pure NG**: the corpus proved pure NG has *no interior minimum and no neutral
`0⁺⁺` mass*; the rigid completion *does* produce a **stabilised neutral `0⁺⁺`** with a closed-form mass.
The relation `m/√σ = 2π√(2κ)` is **falsifiable**: `κ∈[0.138,0.193]` reproduces the lattice band
`3.3–3.9` (Morningstar–Peardon `≈3.9`), and it predicts the **shape** (linear in `√κ`).

**Cross-check (closure i):** the Gaussian/Hartree closure gives the DSE gap equation (witness §3),
`m²=λm²/(1+m²)`: `m=0` is **always** a solution; the positive branch needs `λ>1` (the absent tachyon). So
the *Gaussian* closure hits the wall, while the *rigid* closure produces a mass — the two closures differ
exactly by the rigidity, consistent with the picture.

---

## §4 · Honest: what is missing and how far this is from a proof

1. **`κ` is not derived — it is the missing input `[O]`.** The rigidity `κ` (core thickness) must come from
   the YM measure; deriving it is `H1` in string coordinates. The construction **localises** the missing
   input to one dimensionless number and gives its exact relation to the mass, but does not compute it.
   (This is the same "core thickness" that `TRANSMUTE-STRING-SPECTRUM` named; here it has a formula.)
2. **The rigid string is an EFT with its own caveats.** The Polyakov–Kleinert rigidity runs and has known
   UV subtleties in `D=4`; the simple circular ansatz + quadratic bending `πκ/ρ` is a **model**, not the
   exact surface sum. The result `[C]`, not `[T]`.
3. **Neutral vs charged — resolved only within the model.** The stabiliser is the tube's rigidity (a
   neutral object), not `σL`; but the identification of the physical flux-tube rigidity with `κ` is an
   assumption. If the true YM flux tube is *not* a rigid string (e.g. its rigidity is fluctuation-generated
   with the wrong sign), the construction fails — that is the real open question.
4. **Not a proof.** `m_{0⁺⁺}/√σ = 2π√(2κ)` is a one-parameter fit to a model; the missing step is `κ=H1`.

---

## §5 · What it predicts, and the distance to a proof

- **Prediction:** `m_{0⁺⁺}/√σ = 2π√(2κ)`, `κ∈[0.138,0.193]` for the physical band, `L*≈1.9/√σ`. A
  one-parameter relation tying the *neutral* `0⁺⁺` mass to the flux-tube rigidity.
- **Distance to proof:** one named parameter (`κ`), plus the rigid-string EFT validity. The construction
  moves the target from "pure NG has no `0⁺⁺` mass and no model" to "the rigid completion has one, with
  `κ` = the missing `H1` input."

**Sharpest line.** The Gaussian closure reproduces the DSE wall (`m=0` persists); the **rigid** closure
stabilises the neutral `0⁺⁺` that pure NG cannot, giving `m_{0⁺⁺}=2π√(2σκ)` — a real structural step that
makes the missing core thickness an explicit, falsifiable parameter. But the rigidity itself is the wall:
`κ=H1`, and no derivation of `κ` from the YM measure is given. **Clay 0%; no QED.**

*Refs (corpus):* `PROOF-EFFECTIVE-STRING` (pure-NG no interior min, `L_c`); `TRANSMUTE-STRING-SPECTRUM`
(neutral vs charged torelon, core thickness = H1); `PNT-LARGEN-LOOP` (MM closed in `W`, 1-form channel);
`ROOF-CONFINEMENT-GAP` (`m≥4√σ` strong-coupling); `SLOT-DSE-FIXPOINT` (Gapless `m=0` persists).
*Refs (external):* Lüscher–Symanzik–Weisz, NPB 173 (1980) 365; Lüscher–Weisz, JHEP 07 (2002) 049; Aharony–Komargodski, JHEP 05 (2013) 118; Polyakov, NPB 268 (1986) 406 / Kleinert, Phys. Lett. A 174 (1993) 335 (rigid string); Morningstar–Peardon, PRD 60 (1999) 034509 (`m_{0⁺⁺}/√σ≈3.9`).
*Witness:* `research/unconditional/w7/construct_masterfield.py`.
