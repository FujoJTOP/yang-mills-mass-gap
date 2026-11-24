# P6-PERIODIC-TO-H1 — from a positive periodic constant-mode curvature floor to H1 and to the mass gap: the exact logical chain, and why the floor is a susceptibility *ceiling* (H1-shaped) but never an *atom*

> **2026-11-24 · Wave P6 — the logic question.** Assume the physical periodic (untwisted) constant
> mode has a positive, `L`-independent curvature floor `G″_{k=0}(a₀*) ≥ c₀ > 0`, as the exact `β=0`
> Haar anchor `0.2688` suggests and as the strip/twist fits (`c=0.248…0.41`) support. **What does that
> floor imply for `H1` and for the mass gap, and with what logical strength?** This note resolves the
> chain: it separates the curvature floor, the amplitude-variance bound, and the atom; proves what the
> floor *does* imply (a `k=0` susceptibility ceiling, conditionally); states precisely where the
> implications stop (the non-constant soft tower, the 4D uniformity, the atom); and names the one
> missing input. **No QED. Clay 0%.**
>
> **Reads.** `D-WAVE-SYNTHESIS.md`, `SM-WAVE-SYNTHESIS.md`, `STAR2-CLASSIFICATION.md`,
> `D6-DECISION.md` (the dichotomy `χ<∞⟺saturation`, `G″(a₀*)=1/χ`), `D4-HAAR-WILSON.md` (exact
> `G″=G″_Haar+G″_Wilson`, Haar floor `0.2688`), `SM4-CONSTANT-MODE.md` (constant-mode floor, scope),
> `H6-HEAD-SPECTRAL.md` (typical point vs essential infimum), `H1-RG-DESCENT.md` (descent breaks on
> the soft mode; curvature≠atom), `docs/paper-v3/main.tex` §6.7–§6.8, `docs/paper-v2/ch02-squeeze.tex`
> (the squeeze, H1 necessary-not-sufficient, thin-atom proposition, Lean `H1_does_not_imply_atom`).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness (cited, not re-run here). No QED.
> **Terminology fix (charter §2).** Throughout, `H1` is the *paper's* object, `Var_μ(ξ_{k→0}) = χ ≤ 1/λ`
> (ch02 `eq:h1`) — an **upper** bound on the susceptibility, satisfied by *saturation* `G″→1/χ>0`. The
> `D2-ANALYTIC-LIMIT` phrase "H1 is false for the constant mode" refers to the *decay-to-zero*
> question and must be read as "the `χ=∞` (decay) scenario does not occur"; it does **not** contradict
> the paper's `H1` ceiling. We do not propagate that wording as a theorem.

---

## 0 · Verdict (no QED; Clay 0%)

**A positive periodic constant-mode curvature floor is, at full strength, the `k=0` face of `H1` — a
susceptibility ceiling `χ_{k=0} ≤ 1/c₀` — and nothing more.** It is:

1. a **curvature floor at `k=0`** `[T]`/`[V]` — not an atom;
2. a **bound on `Var(a)`** (the `k=0` zero-momentum susceptibility) only **conditionally** `[C]`
   (Gaussian marginal / CLT; the rigorous Brascamp–Lieb route is `[NG]` at the Savvidy flank);
3. **never an atom** `[T]` — at `β=0`, `G″=0.2688>0` while `f₁²=0` exactly.

It gives **no** bound on the non-constant soft tower (`[NG]`/`[O]`), and it is **H1's shape in one
channel, conditionally** — the 4D uniform statement *is* the open step. Crucially, **even a fully
proven `H1` does not give the gap** `[NG]`: `H1` is an upper bound on a harmonic mean, and an upper
bound on a moment can never force the *bottom* of support (thin-atom proposition, machine-checked).
The floor constrains `χ` — a **moment** — and the smooth squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` shows a χ-ceiling
cannot pin the bottom of support: it is consistent with *any* bottom (thin atom), and it only makes the
upper rung `M₀/χ` *larger* (`≥ M₀ c₀`), so it gives no upper bound either. It touches the wrong
functional. The gap needs the *other* half — the residue/atom floor `f₁² ≥ r₀ > 0` with isolation `I₂`
— which the curvature floor neither supplies nor points to.

**Logical skeleton (the whole answer in one chain; `╳` = does not imply).**

```
 [P]  G″_{k=0}(a₀*) ≥ c₀ > 0 , L-independent         premise; [T]/[V] in 1D & β=0; [O] in 4D periodic
        │  [C]  Gaussian marginal / CLT  (G″·Var=1; measured 2–5% defect; BL route [NG] at Savvidy flank)
        ▼
 [Q]  Var(a) = χ_{k=0,L} ≤ 1/c₀   (k=0 susceptibility ceiling)
        │  [T] Var(a)=χ_L=Σ_d(1−|d|/L)C(d)→χ ;  [O] 4D: is χ_L bounded uniformly? (= the decompactification theorem)
        ▼
 [H1] Var_μ(ξ_{k→0}) = χ ≤ 1/λ , uniform in (L,a)          necessary for the gap
        ╳  H1 ⇏ atom  ╳   ([T] thin-atom; machine-checked H1_does_not_imply_atom)
 [atom] f₁² ≥ r₀ > 0   with isolation  Λ_★ − m₁² ≥ δ > 0    ← THE MISSING HALF  [O]
        │  [T] squeeze: f₁²/χ ≤ m₁² ≤ M₀/χ
        ▼
 [gap]  m₁² ≥ f₁²/χ ≥ r₀ λ > 0                              needs BOTH halves  (H* = I1 ∧ I2)
```

The floor only ever touches the upper arrow (the ceiling). The lower half (the atom) is untouched.

---

## 1 · What a positive constant-mode floor actually is

### 1.1 The object

The campaign's (correct) coarse curvature is the **coarse-field-conditioned marginal curvature**

```
   G″(v;a₀) = −d²/da² log ρ_v(a)|_{a₀} = E[∂_v²V|a₀] − Var(∂_vV|a₀) ,   V = S_W + Φ ,   Φ = −log J ,   [T]
```

`ρ_v` the exact constant-mode marginal, `a = ⟨v, ξ⟩`, `v_l = ê/√N` the unit-norm *uniform* (`k=0`)
mode (`‖v‖=1`). For the constant mode `a₀*` is the *interior* maximiser of `ρ_v` (the flat-chart Haar
edge `J(2π)=0` forces it interior), so `ρ_v'(a₀*)=0`, `ρ_v''(a₀*)<0`, and

```
   G″(v;a₀*) = −ρ_v''(a₀*)/ρ_v(a₀*) ≥ 0 ,   strict >0 for a non-degenerate maximum.       [T]
```

This is `SM4` Thm (constant-mode positive floor); the sign is unconditional for every constant `v`,
every `L`, boundary/twist, `β`; the **value** is `[V]`.

**The premise `[P]`.** That this floor is positive **and `L`-independent in the physical periodic
(untwisted) model**. Proven: at `β=0` exactly (`G″→9/(4π²−6)=0.2688`, `L`-independent; `D1`,`D2`,
`D4`), and on the `1×L` strip (1D exponential mixing ⇒ `χ<∞`); the twist fits give positive intercepts
`c=0.248` (`θ=0`), `0.41` (`θ=π`). `[O]` in 4D periodic: at `b=0` the variance grows and `G″` decays,
and in `d≥2` the constant-mode susceptibility can diverge (KT/Goldstone). **So `[P]` is itself the
open question; we treat it as a hypothesis and ask what it buys.** This is the entire point of the note.

### 1.2 Three candidate readings — and which one `[P]` is

| reading | statement | status |
|---|---|---|
| (i) **curvature floor at `k=0`** | `G″_{k=0}(a₀*) ≥ c₀ > 0`, `L`-independent | **[T]** sign (`SM4` Thm) + **[V]** value; **[O]** 4D uniformity |
| (ii) **bound on `Var(a)`** | `Var(a) = χ_{k=0,L} ≤ 1/c₀` | **[C]** — Gaussian marginal / CLT only; see §1.3 |
| (iii) **an atom** | `f₁² = ν({m₁²}) ≥ r₀ > 0` | **NOT implied, ever** — **[T]** counterexample at `β=0` (§1.4) |

So a positive constant-mode floor **is** (i); it **conditionally gives** (ii); it **is not and cannot
give** (iii).

### 1.3 Why (i) ⇒ (ii) is only `[C]`, not `[T]`

The bridge is the exact **Gaussian-marginal** identity

```
   marginal of a = ⟨v,ξ⟩ is Gaussian   ⟹   G″(a₀*) · Var(a) = 1 ,   i.e.  G″ = 1/χ_{k=0} .     [T]
```

The converse fails: for a non-Gaussian marginal the product `G″·Var` can sit above or below `1`. At
`β=2.5` it is measured `G″≈0.281` vs `1/Var≈0.272` (a `2–5%` gap; `ch02 §G'' vs 1/χ`, `FORGE5`). Two
consequences:

- **The clean implication is conditional.** `G″(a₀*) ≥ c₀` bounds the *local* curvature at the mode,
  not the global variance: a density that is sharp at `a₀*` but has heavy tails can have large `G″`
  **and** large `Var`. The `k=0` mode is a delocalised sum of links, so the CLT *makes* it near-Gaussian
  (Lindeberg / local CLT; `D6` §1, `D4` §2.1), which is why `(i)⇒(ii)` holds in practice — but this is
  **[C]**, not a theorem.

- **The rigorous one-way route is `[NG]`.** The only *proved* inequality `curvature ⇒ variance` is
  Brascamp–Lieb: `V_eff'' ≥ c₀` **a.e.** `⇒ Var(a) ≤ 1/c₀`. It needs the **essential infimum**
  `ess inf G″`, which is **negative** on the Savvidy flank (`D4` §4.2: `−(log w)'' = −0.37` on ~half the
  support; `H6` §3: `ess inf G″ = m² < 0` at the double-well barrier). So the BL route to (ii) fails
  exactly where the constant-mode floor is evaluated. `[NG]` for the "global-log-concavity ⇒ variance"
  class.

Hence: **`[P]` gives `Var(a) ≤ 1/c₀` only through the Gaussian/CLT reading `[C]`; the unconditional
inequality does not exist.**

### 1.4 Why curvature ≠ atom — the `β=0` acid anchor `[T]`

At `β=0` the links are i.i.d. Haar, the constant-mode marginal is the `L`-fold convolution of the
one-link Haar marginal — a **bounded continuous density, no point mass** — with exact limit

```
   χ = Var_Haar(ξ₃) = (4π²−6)/9 = 3.7198 ,     G″ → 1/χ = 9/(4π²−6) = 0.26883   (L-independent).   [T]+[V]
```

Exact FFT (`D1`): `G″ = 0.205, 0.240, 0.249, …, 0.265 → 0.2688` for `L=1…16`; `Var(a)=3.7198` at every
`L`; kurtosis `2.57→2.97→3`. **Meanwhile `f₁² = ν({m₁²}) = 0`** at `β=0` — the Haar/free anchor is
atomless (`D-WAVE` §2; `H1-RG-DESCENT` §5; `V2` §7). The free massive composite `:φ²:` likewise has
`G″ = 1/χ ≠ m₁²` and `f₁² = 0` at its bottom (`ch02 §G'' vs 1/χ`).

**Why these are different objects (the exact reason a curvature floor cannot be an atom).**

- `G″(a₀*) = −(log ρ_v)''(a₀*)` is a **second-moment / local-concavity** functional of the *coarse
  amplitude density* `ρ_v` on the chart. In the Gaussian case it equals `1/Var(a) = 1/χ`, and
  `χ = ∫ s^{-1} dν(s)` is the **residue-weighted harmonic mean** of the Källén–Lehmann masses.
- The **atom** `f₁² = ν({m₁²})` is the *point mass* of the Källén–Lehmann measure `ν` (a measure on
  `s = m²`) at the **bottom of its support** `m₁² = inf supp ν`.
- The two are functionals of different spaces (`ρ_v` on the chart vs `ν` on mass-squared) bridged only
  by `χ = ∫ s^{-1} dν`. The atom contributes `f₁²/m₁²` to `χ`, and is **invisible** to any floor/ceiling
  on `χ` unless `f₁²` is bounded below. A tiny residue at tiny mass (`f₁² ~ m₁²`) keeps `χ` finite while
  driving the bottom to zero — the **thin-atom measure** `ν = wδ_ε + (1−w)δ_{s₂}`, `w = ε/(2λ)`, has
  `χ ≤ 1/λ` yet `inf supp ν = ε → 0`, `f₁² = w → 0` (`ch02` Prop. thin atom; machine-checked
  `H1_does_not_imply_atom`, model `(f,m,M,χ,λ)=(0,1,1,1,1)`).
- **Structural reading** (`ch02 §schur`): every functional built from finitely many moments or from
  positivity is **convex-order / Schur-concave**, hence gives an **upper** bound on the atom — never a
  lower bound. Curvature (a moment) and the atom (the edge of support) lie on opposite sides of this
  barrier. `β=0` is the exact witness: `G″ = 0.2688 > 0` is a *moment* fact, `f₁² = 0` is the *edge*
  fact, and they coexist.

---

## 2 · The implications, with their exact logical strength

### 2.1 No bound on the non-constant soft tower — `[NG]`/`[O]`

The constant mode is `k=0` **exactly**. The soft tower is the `k≠0`, `k→0` modes (long-wavelength
profiles). A floor at one direction `v = e₃/√N` gives nothing at another direction without a
**correlation inequality / monotonicity**. The two written mechanisms both fail here:

- **Fast-fibre silence.** The per-block (fibre) transfer gap is positive and volume-uniform (`~Λ_UV`)
  but **provably silent** on the coarse curvature: the marginal `ν` and the fibre conditionals are
  independent coordinates (two-coordinate independence; `B3-ISOPERIMETRIC`, `H6` §2 reading (iii), §6·D).
  Sending the fast gap to `∞` leaves `G″_marg` unchanged.
- **No monotonicity `χ(k) ≤ χ(0)`.** The tool that *would* transport a `k=0` bound to `k≠0` is a
  reflection-positivity / ferromagnetic ordering `Var(a)_k ≤ Var(a)_0`. For 4D Yang–Mills this is **not
  established** (`[O]`). Absent it, **no implication**.

Equivalently (`SM4` §4 scope): the constant-mode floor **does not transport** to the non-constant
`k→0` modes; that transport is exactly the `D`-defect / uniform-Fisher bound, which is `H1`. `[NG]` for
the claim "constant-mode floor ⇒ tower floor"; `[O]` for the transport itself.

### 2.2 A `k=0` susceptibility ceiling, conditionally — `[C]`

Grant (i)⇒(ii) (§1.3, Gaussian/CLT). With the unit-norm mode, the amplitude variance **is** the
finite-volume zero-momentum susceptibility:

```
   Var(a) = χ_L := Σ_d (1−|d|/L) C(d)  →  χ := Σ_d C(d) .       C(d)=⟨ξ_l ξ_{l'}⟩_c connected .   [T]
```

So `[P]` gives `χ_{k=0,L} = Var(a) ≤ 1/c₀` — an **upper bound** on the `k=0` susceptibility, in the
**H1 direction**. This is a genuine `[C]` implication (`[T]` at `β=0`: `χ=3.72`, trivially bounded;
`[T]` in the 1D strip: `χ<∞` by exponential mixing; `[C]` in 4D at `β=2.5`).

**But note the direction.** `χ ≤ 1/c₀` is a *ceiling*. H1 is exactly such a ceiling. So §2.2 is
"the floor gives the H1-shaped ceiling in the constant channel." It is **not** a lower bound on any
mass.

### 2.3 This is `H1`'s shape in one channel; the 4D uniform version is `[O]`, not implied

The dichotomy is exact (`D6` §1, `[T]` Gaussian / `[C]` CLT):

```
   χ < ∞  ⟺  χ_L bounded  ⟺  G″ → 1/χ > 0     (saturation)
   χ = ∞  ⟺  χ_L ≍ L^p   ⟺  G″ ≍ L^{−p} → 0    (decay: p=2 free-box, p=1 marginal)
```

The paper's `H1` is `χ = Var_μ(ξ_{k→0}) ≤ 1/λ` **uniform in `(L,a)`**. §2.2 supplies this **at the mode
`a₀*`, in the `k=0` channel, conditionally**, and **on the strips / at `β=0`**. It does **not** supply:

- **uniformity in `a`** — the floor is at the mode; the essential infimum over the support is negative
  (Savvidy tail, `H6` §3), so "typical point `>0`" does not lift to "`ess inf > 0`";
- **4D uniformity in `L`** — whether `χ_L` stays bounded in the periodic 4D model is the
  Kosterlitz–Thouless/Goldstone question, which **is** the open `H1` (`D-WAVE` §2, `D6` §5).

**Therefore `[P]` does not *imply* `H1`.** At best it *is* `H1`-in-the-`k=0`-channel conditionally, and
the 4D uniform statement it would have to become **is** the open step. Asserting `[P]` for all `L` in
4D is asserting (a chunk of) `H1` itself — no free lunch. `[C]` conditional, `[O]` in 4D.

### 2.4 Even full `H1` does not give the gap — `[T]`/`[NG]`

This is the decisive, conservative point. By the exact squeeze (`ch02` `eq:squeeze`):

```
   f₁²/χ  ≤  inf supp ν  ≤  M₀/χ ,        χ = ∫ s^{-1} dν ,   M₀ = ∫ dν .            [T]
```

- **Necessity** `[T]`: a gap `m₁² ≥ λ` forces `χ ≤ M₀/λ` (`gap ⇒ χ ≤ M₀/λ`).
- **Non-sufficiency** `[T]`/`[NG]`: the thin-atom measure has `χ ≤ 1/λ` (H1 holds) yet `inf supp ν → 0`
  (gap fails). Machine-checked `H1_does_not_imply_atom` (`H1Sharpening.lean`; axioms only
  `[propext, Classical.choice, Quot.sound]`, zero `sorry`). Structural reason: H1 is an upper bound on
  a moment; the atom is the edge of support; **an upper bound on a moment can never force a bottom**
  (`ch02 §schur`, the "upper-bound machine / lower-bound target" mismatch).

**The gap is the wrong functional for a moment ceiling to reach.** With the ceiling `χ ≤ 1/c₀` from
§2.2 (so `1/χ ≥ c₀`), the squeeze's **lower** rung is `m₁² ≥ f₁²/χ ≥ f₁² c₀`; this is a genuine
**lower** bound on `m₁²` **iff `f₁² > 0`**, and with `f₁² = 0` it degrades to the trivial `m₁² ≥ 0`.
The **upper** rung `m₁² ≤ M₀/χ` is only made *larger* by the ceiling (`M₀/χ ≥ M₀ c₀`), so the ceiling
supplies no upper bound on the gap either. A χ-ceiling is consistent with *any* bottom — the thin atom
makes this exact. The lower bound requires `f₁² ≥ r₀ > 0`. So:

```
   gap needs both:   (f₁² ≥ r₀ > 0, isolation I₂)  ∧  (χ ≤ 1/λ)  ⟹  m₁² ≥ f₁²/χ ≥ r₀ λ > 0 .   [T]
   the floor supplies, at best and conditionally:  the SECOND factor only.                       [C]
```

This is the campaign's `H*` (`ch02` `eq:hstar`): `I1 ∧ I2` with `I2` the isolation
`Λ_★ − m₁² ≥ δ > 0`. The floor addresses `χ`; it cannot manufacture `f₁²`.

---

## 3 · What additional input would turn the floor into a gap bound

Ordered from "closes the gap" to "closes one sub-step":

### 3.1 The residue/atom floor `f₁² ≥ r₀ > 0`, with isolation `I₂` — **the missing half `[O]`**

This is the single decisive input. With it and the susceptibility ceiling (which `[P]` can at best,
conditionally, supply), the squeeze gives the gap:

```
   (f₁² ≥ r₀)  ∧  (χ ≤ 1/λ)   ⟹   m₁² ≥ f₁²/χ ≥ r₀ λ > 0 .          [T] given the pair
```

`I2` (isolation `Λ_★ − m₁² ≥ δ > 0`) is what upgrades "an atom somewhere at the bottom" to "an isolated
pole = a definite mass"; the residue–level-gap lemma (`ch02`) turns `I1` into `I2` and back under a
self-certifying condition. **No known tool supplies `f₁²`.** By §1.4, moments / positivity / convex
order are Schur-concave and give **upper** bounds on `f₁²` only; the free field is a feasible atomless
point; finite moment data cannot distinguish an atom from a smoothed atomless measure. A **new object
outside the convex-order class** is required — the charter's "6th species" (a Hardy-/total-positivity-/
new-cone-type certificate; the campaign's attempted Hardy escape). `[O]`, and `[NG]` for the entire
moment/positivity/convex-order method class.

### 3.2 Tower transport — a sum rule or a monotonicity `[O]`

To reach the **non-constant** `k→0` tower from the `k=0` floor one of:

- **monotonicity / RP ordering** `Var(a)_k ≤ Var(a)_0` (ferromagnetic-type correlation inequality), so a
  `k=0` ceiling bounds all `k≠0` modes; not established for 4D YM;
- **spectral sum rule** expressing the tower's contribution to the constant-mode susceptibility with
  explicit **positive** weights, so a `k=0` floor forces the lightest tower shell up. The corpus has
  sum-rule attempts (`A1-SUMRULE-A`, `A2-SUMRULE-CONDENSATE`) but none closes the light-shell bound; the
  head is a **spectral no-go** for the tool class (`H6` §4, `main.tex` §6.5 "The spectral no-go").

Either way: `[O]`.

### 3.3 Ward identities — what they do and do not give `[T]`/`[O]`

The exact Ward identity is `E[ξ f] = 1` with `f = ⟨v,∇S⟩`. It is exact for every measure and `β`, and
says `ξ` and `f` have **unit correlation** — *not* that the conditional variance is `G″` (`ch02
§G'' vs 1/χ`). So Ward is **not** a bridge from the floor to the atom by itself; it fixes the
Gaussian-limit normalization, not the residue. A genuine **Ward sum rule** that pins the atom's residue
against the longitudinal channel would be a `[O]` input — related to, but not identical with, the
spectral sum rule of §3.2.

### 3.4 4D uniformity — `χ_L` bounded = the decompactification theorem itself `[O]`

Even the *ceiling* `χ ≤ 1/c₀` needs `χ_L` bounded uniformly in `L` in 4D, which **is** `H1`/the mass gap
(`D6` §5). In `d≥2` the constant-mode susceptibility can diverge (KT/Goldstone); the periodic model at
`b=0` is not settled. This is not "additional" in the sense of a separate tool — it is the target
re-stated. `[O]`.

**Summary of §3.** The floor can, at most and conditionally, contribute the susceptibility **ceiling**
(one factor of the sufficient pair). The gap is gated by the **atom** `f₁² ≥ r₀` (and its isolation),
which no moment/positivity argument can reach; the tower transport and the 4D uniformity are separate
`[O]` steps that the constant-mode floor does not touch.

---

## 4 · Logical-strength ledger

| item | statement | label |
|---|---|---|
| constant-mode sign | `G″(v;a₀*) = −ρ_v''/ρ_v ≥ 0`, interior maximiser, all `L`/twist/`β`/constant `v` | **[T]** (`SM4` Thm) |
| constant-mode value | `+0.362` (off-diag)…`+0.694` (Cartan) at `β=2.5`, `L=1`; `0.2688` exact at `β=0` | **[V]** / **[T]** (β=0) |
| premise `[P]` (4D periodic `L`-uniform floor) | `G″_{k=0}(a₀*) ≥ c₀ > 0` | **[T]/[V]** 1D & β=0; **[O]** 4D |
| `G″·Var = 1` (Gaussian marginal) | inverse-susceptibility identity | **[T]**; **only-if-Gaussian** |
| `[P] ⇒ Var(a) ≤ 1/c₀` | `k=0` susceptibility ceiling | **[C]** (Gaussian/CLT; measured 2–5% defect) |
| BL route `V_eff''≥c₀ a.e. ⇒ Var≤1/c₀` at the floor | needs `ess inf G″`; negative on Savvidy flank (`D4`, `H6`) | **[NG]** for the class |
| `Var(a) = χ_L = Σ_d(1−|d|/L)C(d) → χ` | amplitude variance = zero-momentum susceptibility | **[T]** |
| `χ<∞ ⟺ saturation G″→1/χ>0`; `χ=∞ ⟺ decay` | the dichotomy | **[T]** (Gaussian) / **[C]** (CLT) |
| `[P] ⇒ H1` (full, uniform in `(L,a)`) | `χ_{k→0} ≤ 1/λ` everywhere, all `L`, 4D | **[O]** — 4D uniformity = the open step |
| curvature floor ⇒ atom `f₁²>0` | **false**: `β=0` has `G″=0.2688`, `f₁²=0` | **[T]** (counterexample) |
| `H1 ⇒ gap` | **false**: thin atom `ν=wδ_ε+(1−w)δ_{s₂}`, `w=ε/2λ` | **[NG]**; machine-checked `H1_does_not_imply_atom` |
| `gap ⇒ χ ≤ M₀/λ` | necessity | **[T]** |
| `(f₁²≥r₀) ∧ (χ≤1/λ) ⇒ m₁² ≥ r₀λ` | the sufficient pair `H*` | **[T]** (given the pair) |
| floor ⇒ non-constant tower bound | **false**: fast-fibre silent; no `χ(k)≤χ(0)` | **[NG]** / **[O]** (transport) |
| the missing input | `f₁² ≥ r₀ > 0` + isolation `I₂` (outside convex-order class) | **[O]** (=`H*` residue half) |
| 4D `χ_L` bounded | the decompactification theorem (=`H1`) | **[O]** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 5 · Sharpest line

**A positive periodic constant-mode curvature floor is a `k=0` curvature floor `[T]`/`[V]` — at most a
susceptibility *ceiling* `χ_{k=0} ≤ 1/c₀`, and that only conditionally `[C]` (Gaussian/CLT; the rigorous
Brascamp–Lieb route is `[NG]` at the Savvidy flank). It is never an atom: at `β=0`, `G″=0.2688>0` while
`f₁²=0` exactly, because curvature is a moment (`1/χ`, a residue-weighted harmonic mean) and the atom is
the edge of the Källén–Lehmann support — opposite sides of the convex-order barrier. The floor gives no
bound on the non-constant soft tower (`[NG]`, fast-fibre silent; `[O]`, no `χ(k)≤χ(0)`), and it is
`H1`'s *shape* in one channel only — the 4D uniform statement it would have to become **is** the open
`H1`, so it implies nothing beyond what it already assumes. Most decisively, even a full `H1` does not
give the gap: `H1` is an upper bound on `χ` (a moment), the gap is the bottom of support, and the
squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` shows a χ-ceiling reaches the edge only through `f₁²/χ` — vacuous at
`f₁²=0`, and it only enlarges the upper rung `M₀/χ`. The mass gap is gated by the other half — the
residue floor `f₁² ≥ r₀ > 0` with isolation `I₂` — which the curvature floor neither supplies nor
approaches; it requires a certificate outside the moment/positivity/convex-order class. So
the chain is: **floor ⇒ (conditional) `k=0` ceiling, `╳` tower, `╳` atom, `╳` gap.** Do not upgrade it.
**No QED. Clay 0%.**

---

*Refs.* `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`, `eq:h1`, `eq:squeeze`, thin-atom Prop.,
`H1_does_not_imply_atom`, `eq:sufficient`, `eq:hstar`, `§schur`); `docs/paper-v3/main.tex` §6.5
(spectral no-go), §6.7 (decompactification), §6.8 (pincer), items 8–9 of "what this paper establishes";
`D-WAVE-SYNTHESIS.md` §2 (`G″>0` coexists with `f₁²=0`; `χ<∞` is `H1` in 4D); `SM-WAVE-SYNTHESIS.md`
§2–§3; `STAR2-CLASSIFICATION.md` §3 (head = `H1`); `D6-DECISION.md` §1, §5 (dichotomy,
`G″(a₀*)=1/χ`, `χ<∞`=H1); `D2-ANALYTIC-LIMIT.md` ($\beta$-dependent plateau; *terminology note* — "H1
is false" there means "no decay to 0"); `D4-HAAR-WILSON.md` §4–§5 (Haar floor `0.2688`, Wilson sign at
large `L` = H1, `(LC-W)` false pointwise); `D5-UNIVERSALITY.md` (floor value non-universal,
boundary-dependent); `SM4-CONSTANT-MODE.md` §4 (scope: constant mode only, not the `k→0` tower, not
decompactification); `H6-HEAD-SPECTRAL.md` §3 (typical point vs essential infimum; double well);
`H1-RG-DESCENT.md` §4.2, §5 (descent breaks on soft mode; `N/12` curvature ≠ atom). Standard inputs
(CLT, Brascamp–Lieb, convex order/Schur, Källén–Lehmann) used as in the corpus. **No third-party sources
downloaded.** No witness run for this note.
