# A1 · Attack surface A — the monotone-to-residue sum rule `Δa` vs `f₁²`

> **VERDICT (no QED).** Both proposed forms are **`[NG]` — this route is dead.**
>
> - **`Δa ≥ c·f₁²`** cannot be the right inequality *in either direction*. As written it re-reads as an
>   **upper** bound `f₁² ≤ Δa/c ≤ a_UV/c` on the residue — the opposite of the residue **floor**
>   `f₁² ≥ r₀ > 0` the campaign needs. And the exact Komargodski–Schwimmer (KS) sum rule carries the
>   residue scale `f` (⇔ `f₁²`) as a **prefactor `f⁴`** that cancels from the identity; the only
>   inequality in it is `Im A ≥ 0`. So `a`'s spectral representation is *homogeneous in `f`* and pins
>   neither an upper nor a lower bound on `f₁²`.
> - **`a_IR ≥ F(f₁², m₁)`** is false at the only endpoint pure YM actually has: a gapped IR is the
>   trivial theory, **`a_IR = 0`** for *any* residue. Hence `F` would have to be `≤ 0`, i.e. vacuous.
> - **Structural reason (the sharpest statement).** `a` is a **three-point** ⟨TTT⟩ anomaly coefficient
>   (in a CFT the *two*-point ⟨TT⟩ is fixed by `c` alone; `a` does not enter it). `f₁²` is the atom of a
>   **two-point** Källén–Lehmann density — the ⟨ΘΘ⟩ trace-trace correlator, whose spectral weight is
>   `∝ (β/2g)² f₁²` and is governed by the β-function and `f₁²`, **not** by `a`. The two objects live in
>   different multiplicity channels; there is **no kernel** built from anomaly data that connects them.
>
> The bridge, *as posed*, is dead. The genuinely useful statement that survives is a **`[T]` operator
> identity** (§2.1) that contains `f₁²` and `β` but no `a`; and the **exact missing step `[O]`** is an
> *IR-controlled, UV-insensitive upper bound on the dilaton spectral integral* — the a-theorem integral
> is deliberately UV-sensitive and is the wrong object for a residue bound.

Labels: `[T]` theorem · `[C]` conditional · `[NG]` no-go for this route · `[O]` exact missing step · `[V]` numerics.

---

## 1. The attempted statements (verbatim)

Let `a` be the Euler-anomaly coefficient in the trace anomaly and `f₁²` the atom weight of the
Källén–Lehmann measure `ν` of `tr F²` in the `0^{++}` channel (project normalization, `ch02-squeeze.tex`,
eqs. `chi`, `squeeze`).

- **(S1)** `Δa := a_UV − a_IR  ≥  c·f₁²` with a universal constant `c > 0`.
- **(S2)** `a_IR  ≥  F(f₁², m₁)` for a function `F` with `F(0,·)=0` and `F(s,·)>0` for `s>0`.

Known inputs: KS monotonicity gives `0 ≤ Δa ≤ a_UV`, with `a_IR = 0` at a gapped endpoint and
`a_UV = (N²−1)·31/180` for pure SU(N) (`ch10-betarg.tex` §"Central charges"; `LIT-ERA-10.md` §4).

---

## 2. The spectral representation that `a` actually satisfies

### 2.1 Operator identity — `[T]`

The trace anomaly of pure YM reads
```
T^μ_μ = (β(g)/2g) tr F²  +  a E₄  −  c W²                (+ renormalization-scheme contact terms).
```
The three terms are three *different* objects: the first is the (classically non-scale-invariant)
β-function term, the second/third are the genuine anomalies.

Squaring and taking the connected vacuum correlator at flat metric:
```
⟨Θ(x) Θ(0)⟩_conn  =  (β/2g)² ⟨trF²(x) trF²(0)⟩_conn
                  =  (β/2g)² ∫ dν(s) Δ(x;s) ,      Θ := T^μ_μ .
```
The **spectral (`s>0`) part of `⟨ΘΘ⟩` is exactly `(β/2g)² ν`** — i.e. exactly the `tr F²` measure whose
atom is `f₁²`. The `aE₄ − cW²` pieces are **local curvature terms**: in flat space they contribute only
`δ⁴(x)`-and-derivative contact polynomials, *zero spectral weight at `s>0`*. Hence:

> **`f₁²` is the atom of the two-point function `⟨ΘΘ⟩`, and that function is controlled by `β` and `f₁²`;
> the anomaly coefficient `a` is invisible in its Källén–Lehmann density.** `[T]`

This is the honest positive sum rule, and it does **not** involve `a`. It is the campaign's own
`(β-function) × (residue)` pair (coordinates 2 and 1), not `a`.

### 2.2 Why `a` has *no* two-point spectral representation — `[T]`

In 4D, the stress-tensor **two**-point function of any CFT is fixed by conformal + Poincaré symmetry
to a single tensor structure with coefficient **`c` alone**:
```
⟨T_{μν}(x) T_{ρσ}(0)⟩  =  c · I_{μνρσ}(x) / (x²)⁴ ,       a  does not appear.
```
The Euler coefficient `a` is one of the coefficients of the **three**-point function `⟨TTT⟩` (the
"a-structure"/anomaly-pole piece; `a` and `c` are two of its three constants). This is standard and is
re-confirmed by the modern `⟨TTT⟩` literature (see §7). Consequently:

> **No Källén–Lehmann sum rule over any two-point spectral density can produce `a`.** `a` is
> intrinsically a 3-point (⟨TTT⟩) datum. `[T]`

### 2.3 The *only* exact sum rule for `a`: the KS dilaton relation — and `f₁²` is a prefactor — `[T]`

The KS proof of the a-theorem represents the flow by a dilaton `τ` (Goldstone of spontaneously broken
conformal symmetry) with Wess–Zumino term
```
δS_WZ = ∫√−g σ [ (c_UV−c_IR) W² − (a_UV−a_IR) E₄ ] .
```
The flat-space low-energy `2→2` dilaton amplitude is `A(s) = 2 Δa s²/f⁴ + …` (forward limit), where
`f` is the **dilaton decay constant**. Unitarity (optical theorem `Im A = sσ ≥ 0`) plus a dispersion
relation gives the **exact** sum rule
```
                        f⁴       ⌠∞  ds′
        Δa  =  a_UV − a_IR  =  ────   ⎮   ───  Im A(s′)  >  0 .        (KS)
                        π        ⌡0  s′³
```
**Crucial observation.** `f` is precisely the object `f₁²` measures. By the dilaton PCAC relation,
`⟨0|Θ|dilaton⟩ = f m₁²`, hence `⟨0|tr F²|m₁⟩ = f m₁² /(β/2g)` and
```
        f₁²  =  |⟨0| tr F² |m₁⟩|²  ∝  f² m₁⁴ / (β/2g)² .                (f↔f₁²)
```
Substituting the low-energy `A(s') = 2 Δa s′²/f⁴` into (KS) makes the sum rule an **identity** in which
`f⁴` and `Δa` cancel against the amplitude normalization; the surviving content is only `Im A ≥ 0`,
i.e. `Δa > 0`. Therefore:

> **(KS) is homogeneous in `f`**: it fixes the **sign** of the dilaton spectral integral
> `∫ ds′ Im A/s′³`, and treats `f` (⇔ `f₁²`) as a *free normalization*. It bounds neither
> `f₁² ≤ …` nor `f₁² ≥ …` given `Δa`. `[T]`

This is the exact technical statement of the synthesis's "`a` is phase-blind": the one exact spectral
representation `a` possesses has the residue scale sitting in a prefactor that drops out of the only
inequality (`Im A ≥ 0`) it yields.

---

## 3. Why no universal bridge exists (kernel, hypotheses, direction)

### 3.1 Kernel: none within the anomaly data — `[NG]`
`a` lives in the 4-point/dilaton (equivalently 3-point ⟨TTT⟩) sector; `f₁²` lives in the 2-point
`⟨ΘΘ⟩`/`tr F²` sector. A sum rule is an integral transform between *one* channel's density and a
coefficient; here the coefficient and the density are in **different channels of different
multiplicity**. No anomaly identity supplies a kernel. `[NG]` for "derive `a` as `∫ K(s) dν(s)`."

### 3.2 Dilaton-EFT parameter count — explicit counterexample class — `[NG]`
The universal low-energy description of a conformal flow is
```
S_eff = ∫ d⁴x [ ½ f² e^{2τ}(∂τ)² + … ] + Δa·τ E₄ − Δc·τ W² + …
```
`Δa` (WZ anomaly coupling) and `f` (kinetic decay constant) are **independent Wilson coefficients**.
Unitarity requires `f² > 0`; the a-theorem requires `Δa > 0`. For **every** pair `(Δa, f)` with both
positive there is a consistent dilaton EFT (subject only to UV completion, which is exactly what is in
question). Hence no universal function `F` and no constant `c>0` can relate them within this data class.
A proof of (S1) or (S2) would therefore require **more than** (i) the trace-anomaly operator identity,
(ii) the a-theorem, (iii) the dilaton EFT — i.e. genuinely new YM-specific input (§6). `[NG]`.

### 3.3 The direction/sign obstruction (EFT-independent, decisive) — `[NG]`
The a-theorem is a **decreasing** monotone, `0 ≤ Δa ≤ a_UV`. Read the algebra:
```
(S1)  Δa ≥ c f₁²   ⇔   f₁² ≤ Δa/c ≤ a_UV/c    ← an UPPER bound on the residue.
Reverse (what the campaign actually needs):  f₁² ≥ something,  i.e.  Δa ≤ c′ f₁² would give,
      with Δa = a_UV > 0 fixed for gapped pure YM,   f₁² ≥ a_UV/c′ > 0   ← a FLOOR.
```
So the helpful inequality is the **reverse** of (S1). But (KS) is homogeneous in `f` (§2.3), so it
yields **neither** direction: to get `Δa ≤ c′ f₁²` one would need an **upper** bound on the dilaton
spectral integral `∫ ds′ Im A/s′³` in terms of the *IR* residue `f₁², m₁`. That integral is
**UV-dominated** by construction (it integrates the whole flow up to the conformal regime) and is not
controlled by `f₁²`. Hence (S1) is an upper-bound machine aimed at a lower-bound target — the exact
"upper-bound machine, lower-bound target" mismatch of `ch02-squeeze.tex` §"Why H1 is blind", now
transposed from `χ` to `a`. `[NG]`.

### 3.4 `a_IR` is phase-blind at a gapped endpoint — `[NG]`
A mass gap means the IR is the **trivial** theory, which has **no degrees of freedom and `a_IR = 0`**,
*for any value of `f₁²`*. Thus (S2) `a_IR ≥ F(f₁², m₁)` with `F>0` at `f₁²>0` is **false** at the only
endpoint pure YM actually has. Only `F ≡ 0` survives — which carries no information. This is the
`LIT-ERA-10.md` §2.4 statement ("`a`-theorem is vacuous at a gapped endpoint") made quantitative:
`a_IR` is an IR *count*, and a gap makes it zero independent of the residue. `[NG]`.

---

## 4. Dimensions and normalization (explicit)

| object | dimension | value / definition |
|---|---|---|
| `a` | `0` (dimensionless) | `a_UV^{(SU(N))} = (N²−1)·31/180` (free-vector count, `ch10`, `LIT-ERA-10` §4) |
| `c` (Weyl) | `0` | `c_UV = (N²−1)·1/10`; `a/c = 31/18` at the UV edge |
| `f₁²` (ch02 normalization) | `0` | `ν` normalized so `χ=∫s⁻¹dν` is the propagator and `f₁²/χ ≤ m₁²` ⇒ `f₁² ∈ [0, M₀]`, dimensionless atom fraction (`ch02` `chi`,`squeeze`) |
| `f` (dilaton decay constant) | `+1` (mass) | `⟨0|Θ|dilaton⟩=f m₁²`; `f₁² ∝ f²m₁⁴/(β/2g)²` |
| `s = m²` | `+2`; `dν` | measure on `s`; `dν` dimensionless in ch02 normalization |

**Consequence.** `a` and `f₁²` are **both dimensionless** in the project normalization, so (S1) has a
dimensionless `c` with **no unit obstruction**. The failure of (S1) is *structural and directional*
(§3), not dimensional. In the alternative *unnormalized* convention `f₁² = |⟨0|trF²|m₁⟩|²` has
`[mass]⁶`; then `c` must carry `[mass]⁻⁶ ~ 1/m₁⁶`, the inequality is again `f₁² ≤ Δa·m₁⁶/(…)` (an upper
bound), and matching `1/m₁⁶` at the cutoff is divergent — an independent reason the naive form cannot be
a clean continuum statement. Either way: no floor.

---

## 5. Acid tests (mandatory)

**(a) Free massless field.** `a_UV = a_IR = a_free` (UV=IR, no flow) ⇒ `Δa = 0`; `f₁² = 0` (gapless, no
atom). (S1) reads `0 ≥ c·0 = 0`: **satisfied, but saturates at zero.** By the charter ("a tool that
also 'works' on the free theory proves nothing") a relation that returns `0 ≥ 0` at the free point is
**vacuous there** — it cannot certify a positive residue from `Δa`. (S2) reads `a_free ≥ F(0)=0`:
trivially true, no information. → **both forms pass trivially but prove nothing; no strict positivity.**

**(b) Gapped composite `:φ²:` (gap, no atom).** The operator `:φ²:`, dim 2, has `Z₂` symmetry
`φ→−φ`, so `⟨0|:φ²:|φ⟩=0`: **no one-particle atom** — the KL density is a two-particle continuum from
`4m²`, `f₁² = 0`. The theory has a mass scale (gap) but no atom. Meanwhile it flows free-UV → trivial
IR, so `Δa = a_free > 0`. (S1) reads `a_free ≥ 0`: **satisfied with maximal slack**, giving *zero*
information about the (absent) atom. (S2) reads `0 ≥ F(0)=0`: trivial. → **vacuous on the
gap-but-no-atom test.** A relation true on *both* the free field *and* a gapped no-atom composite, and
returning no strict positivity on either, **cannot** be an atom-certifying tool. **This is the acid-test
failure** — not because (S1) is false, but because it is one-sided/vacuous in exactly the direction that
matters. `[NG]`.

**(c) KS sum rule at the gapped endpoint.** There is no *asymptotic* dilaton in a gapped theory (the
would-be dilaton is the massive lightest glueball; `f` is not a Goldstone decay constant and the
effective δa description degenerates). The (KS) integral `∫ ds′ Im A/s′³` is then saturated by the
massive glueball continuum and constrains the *sign* of `Δa` only — nothing about the isolated atom
weight. Consistent with §2.3. `[T]`.

**(d) Secondary (heuristic, not load-bearing).** Near a conformal fixed point `Δa → 0` at an
`O(1)` light-scalar residue, so any *universal* `c>0` in (S1) is violated in that limit. Pure YM itself
is not near-conformal (gapped, `Δa = a_UV` fixed), so this does not by itself falsify (S1) for YM; it
only shows (S1) is not a universal anomaly identity. `[O]`.

---

## 6. Verdict and what would be needed to close it

| item | status | content |
|---|---|---|
| operator identity `T^μ_μ = (β/2g)trF² + aE₄ − cW²` | **`[T]`** | standard trace anomaly |
| `⟨ΘΘ⟩_conn = (β/2g)²∫dν Δ`, spectral weight `∝ (β/2g)² f₁²`, **no `a`** | **`[T]`** | §2.1 — the positive sum rule, but it is `β × f₁²`, not `a` |
| `a ∉ ⟨TT⟩` (2-pt fixed by `c`); `a ∈ ⟨TTT⟩` (3-pt) | **`[T]`** | §2.2 — `a` blind to 2-pt spectra |
| KS sum rule, `f⁴` a prefactor, only `Im A ≥ 0` survives | **`[T]`** | §2.3 — homogeneity in `f` |
| (S1) `Δa ≥ c f₁²` | **`[NG]`** | upper bound / wrong direction (§3.3); 2-pt vs 3-pt mismatch (§3.1,§2.2) |
| (S2) `a_IR ≥ F(f₁²,m₁)` | **`[NG]`** | `a_IR = 0` at any gapped endpoint (§3.4) |
| **exact missing step** | **`[O]`** | an **IR-controlled, UV-insensitive upper bound** on the dilaton spectral integral `∫ ds′ Im A/s′³` in terms of `f₁², m₁` (to get the useful reverse inequality `Δa ≤ c′ f₁²`, §3.3) |

**What would be needed to close (if one insists).**
1. A **YM-specific low-energy theorem** relating the `0^{++}` glueball's WZ/E₄ coupling (`Δa`) to its
   decay constant `f` — i.e. a relation between two *independent* dilaton-EFT Wilson coefficients
   (§3.2). No such theorem exists; it would be a genuinely new dynamical input, in effect a shadow of
   the missing transport machine. `[O]`.
2. Or replace the UV-sensitive KS integral by an **IR-localized** functional of `ν` that is monotone
   *and* 2-point. Candidate: the entanglement (Casini–Huerta) `c`-function, which is built from the
   stress-tensor two-point function and is monotone — but it is still **phase-blind** (`c_IR = 0` at a
   gap) and so reproduces §3.4. No rescue by that path either. `[NG]` for that sub-route.

**Bottom line.** `a` is the wrong monotone to bridge to `f₁²`: it is a dimensionless, UV-finite
**three-point** anomaly coefficient whose only sum rule carries `f₁²` as a cancelling prefactor, and it
is *decreasing* from the UV with a *zero* IR value at any gap. It can only ever *upper*-bound the
residue, exactly when the campaign needs a *lower* bound. The bridge as posed is **dead**; the honest
residue-carrier is `⟨ΘΘ⟩ ∝ (β/2g)² f₁²` (§2.1), which has no monotone.

---

## 7. Verification list (external citations)

| claim used | source | status |
|---|---|---|
| `⟨TT⟩` fixed by `c` alone; `a,c` are two of the three constants of `⟨TTT⟩` | Chalabi–Herzog–Iliesiu–Jiang, "Weyl anomalies of 4d conformal boundaries and defects," arXiv:2111.14713; Corianò–Maglio, "TTT in CFT: trace identities and the conformal anomaly effective action," Nucl. Phys. B **942** (2019) 303, arXiv:1703.08860; Bzowski–McFadden–Skenderis, "The general 3-graviton vertex (TTT) …, d=4," Nucl. Phys. B **937** (2018) 56, arXiv:1808.10221 | **VERIFIED** this session (WebSearch, all three) |
| KS a-theorem: WZ term `δS_WZ = ∫√−g σ[(c_UV−c_IR)W² − (a_UV−a_IR)E₄]`; forward dilaton amplitude `A(s)=2Δa s²/f⁴`; sum rule `Δa=(f⁴/π)∫ds′ ImA(s′)/s′³>0` | Komargodski–Schwimmer, "On RG flows in four dimensions," JHEP **12** (2011) 099, arXiv:1107.3987 (`LIT-ERA-10` E7, VERIFIED) | **VERIFIED** (previously + this session) |
| `a` coefficients: free vector `31/180`, `a/c=31/18`, `a_UV=(N²−1)·31/180` | `LIT-ERA-10.md` §4 (free-field block, XC-verified); `ch10-betarg.tex` §"Central charges" | **verified by cross-check** (N=4 sum rule) |
| `a_IR = 0` at a gapped endpoint; a-theorem vacuous there | `LIT-ERA-10.md` §2.4 (E7 discussion) | standard, **VERIFIED** (`LIT-ERA-10` E7) |
| trace-trace correlator captures `tr F²` with `(β/2g)²`, `aE₄−cW²` contact | standard trace-anomaly Ward identity (Duff 1977, `LIT-ERA-10` E23); operator algebra | standard; E23 UNVERIFIED in era ledger but relation is textbook |
| `a`-theorem reformulated via stress-tensor 2-pt / entanglement `c`-function | Casini–Huerta line, e.g. arXiv:1304.4411 (appeared in WebSearch) | **seen** this session; not load-bearing |
| dilaton EFT has independent WZ coefficient and kinetic/decay constant | Schwimmer–Theisen, Nucl. Phys. B **847** (2011) 590, arXiv:1011.0696 (from WebSearch); structure of `S_eff` in KS | **VERIFIED** (WebSearch) |

*Note on citation hygiene:* the ⟨TT⟩/⟨TTT⟩ statement is the load-bearing external fact; it is confirmed
by three independent `⟨TTT⟩` papers. The KS sum rule is quoted in the exact form returned by the search
of arXiv:1107.3987 (WZ variation and the `f⁴` prefactor); re-check eq. numbers against the published
JHEP version before any v3 submission. No bibliographic datum here is invented.

*Numerical witness:* `research/v3/w7/a1_sumrule_witness.py` (numpy; single job) — confirms (i) the
`a_UV` table and that for `N≥3` `a_UV > 1`, so the (S1)-implied dimensionless cap on `f₁²` is *weaker*
than the trivial cap `M₀=1` (for `SU(2)`, `a_UV=0.517<1`, so the cap can bite for `c>a_UV`, but it
remains an upper bound and still allows `f₁²=0` — verdict unchanged); (ii) the KS homogeneity: with `Δa`
fixed and the toy positive density `Im A(s)=A₀ s³ e^{−s/M²}`, the solved `f` varies over orders of
magnitude — `Δa` does not bound `f` (⇔ `f₁²`); (iii) that `f₁²=0` is compatible with `Δa=a_UV>0`, so no
`c>0` in (S1) can certify an atom.
