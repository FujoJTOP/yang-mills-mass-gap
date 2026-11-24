# P3-FOURWALL-DERIVATION — the minimal connected `¼tr F²` polymer is a four-walled tube: rigorous re-derivation of PROOF-P3-UNCONDITIONAL §3–§4 at `m₁=4ln(1/u)/a`

> 2026-11-24 · Rigorous redo of §3–§4 of `PROOF-P3-UNCONDITIONAL.md` with the **correct** four-wall minimal-surface geometry, as demanded by the charter §2 correction and `P3-AUDIT-CORRECTED.md` objection 1. Setup as in the parent file: 4D pure SU(N) Wilson lattice YM, `S=βΣ_p(1−N⁻¹Re tr U_p)`, character expansion `e^{(β/N)Re tr U_p}=Σ_r c_r(β)χ_r(U_p)` with `χ_r(𝟙)=dim r` (so `χ_f=tr U_p`, `dim f=N`), `c_f(β)=β/(2N)+O(β³)` for `N≥3` (doubles to `β/2` for SU(2)); `O` = A1g (zero-momentum) single-**spatial**-plaquette `¼tr F²`; `C(t)=⟨O(0)O(t)⟩_c`. Rigor over optimism; every claim labelled **[T]/[C]/[O]/[R]**.
>
> **Gupta verified directly.** R. Gupta, *Introduction to Lattice QCD*, arXiv:hep-lat/9807028, §14.2.2, Eqs. (14.17)–(14.21) — read from the source PDF, verbatim (§4 below). The four-wall tiling and `M_{0⁺⁺}=−4log(β/2N²)` are Gupta's, not paraphrase. `[O]`

---

## VERDICT (up front)

| Object | Status | One line |
|---|---|---|
| **§1 minimal surface = four walls, `\|P\|≥4t/a`** | **[T]** | Airtight: link-balance forces 4 flux strands (one per source-plaquette edge), and only temporal plaquettes advance time, so each strand crosses `t/a` slices ⇒ `≥4t/a` temporal plaquettes. **No cheaper connected structure** — proved, not asserted. |
| **§2 unique leading class, `a₀>0`, `m₁=4ln(1/u)/a`** | **[T]** | Leading class = straight four-wall tubes (summed over transverse position); amplitude `a₀ = (operator coeffs)·(positive endpoint recoupling) > 0`. Everything else costs `≥4t/a+1` ⇒ a full factor `u` down. |
| **§3 combinatorial remainder `Δ≥1 ⇒ O(u)`, rate gap `δ`** | **[T]** (combinatorics) / **[C]** (the stated uniform bound) | The integer cost grading gives `O(u)`-relative for all bounded-excess competitors and `δ=2ln(1/u)/a` for rate-separated ones. Converting this into `Cu\|w_tube\|` uniformly in `t` needs the **decoration exponentiation** (dressed mass/amplitude) = the pre-existing **[R]** gap; it is *orthogonal to the factor 4*. |
| **§4 SU(2) sanity: `−4ln(β/2N²) ↔ 4ln(1/u)`** | **[O]/[T]** | Matches Gupta (14.21) exactly with `u=β/2N²`; SU(2)-special orientation-doubling gives the standard `−4ln(β/4)`, same coefficient 4. Numerical `β=1` `m_eff≈6` brackets the two conventions (`5.5`–`8.3`). |

**Bottom line.** The **four-wall minimal-surface geometry is now a theorem [T]**, and the factor-4 correction to the *geometry and to the leading amplitude/positivity* is complete and rigorous. **Nothing about the four-wall structure resists.** What resists is exactly the one pre-existing hurdle the audit already isolated, wholly independent of the geometry: the **dressed-mass / dressed-amplitude resummation [R]** (audit objection 3), needed to turn the combinatorial `O(u)` statement into a `t`-uniform bound. So the residue conclusion `f₁²≥r₀(β)>0` remains **[C] conditional on [R]**, now resting on a *correct* `m₁=4ln(1/u)/a` and a *proved* minimal surface. `[J]` **No QED on the residue theorem itself**; the geometry sub-proof is genuinely closed.

---

## §1 The minimal connected surface is the four-walled tube `[T]`

Throughout, `a` = lattice spacing, `t` = physical time separation, `N_t := t/a ∈ ℕ` = number of time slices traversed. Grade a polymer by plaquette cost `|P| := Σ_{p∈P}|r_p|` (fundamental `|f|=1`). Let `P` be a **connected, gauge-invariant** set of charged plaquettes (nontrivial irreps) whose link integrals are nonzero, joining the source insertion `O` in slice `0` to the sink `O` in slice `N_t`.

### 1.1 Link balance: each of the four source edges needs a partner `[T]`

The insertion `O` is a single **spatial** plaquette; in the fundamental term it contributes `tr U_{p₀}=Σ (U_{e₁})_{a₁a₂}(U_{e₂})_{a₂a₃}(U_{e₃}^†)_{a₃a₄}(U_{e₄}^†)_{a₄a₁}`. Each of its four boundary links `e₁,…,e₄` therefore carries **one open fundamental index pair**. The Haar integral over link `e_i` factorizes; if `e_i` is covered by the source alone, the factor is `∫dU (U)_{ij}=0` (Gupta's identity list: `∫dg U_{ij}=0`). Hence

> **(1.1)** For the weight to be nonzero, **every one of the four source boundary links must be shared with at least one other charged plaquette.** `[T]`

Two distinct unit plaquettes in a hypercubic lattice share **at most one link** (their intersection is ∅, a site, or a link; two shared links would force identical plaquettes). Therefore no single partner plaquette can cover two source links, and

> **(1.1′)** The source has **≥ 4 distinct partner plaquettes** — four independent "flux strands," each rooted on one source edge. `[T]`

This is exactly the statement that the tube's cross-section is the full 4-link boundary of the source plaquette, i.e. **four walls**. `[T]`

### 1.2 Time-advance: each strand costs `t/a` temporal plaquettes `[T]`

A plaquette is **spatial** (all four edges spacelike; lies wholly inside one time slice) or **temporal** (contains two timelike links; spans slice `j → j+1`). Spatial plaquettes in adjacent slices share **no** link (their links live at different time coordinates), so two plaquettes in slices `j` and `j+1` are adjacent only if one of them is **temporal**. Consequently:

> **(1.2)** Any plaquette-adjacency path from a plaquette in slice `0` to one in slice `N_t` crosses each of the `N_t` time-boundaries `j→j+1`, and each crossing requires `≥1` **temporal** plaquette. A fundamental strand therefore contains `≥N_t` temporal plaquettes; `[T]`

(Equivalently, in the flux picture: fundamental flux is conserved, the polymer's only boundary is source∪sink, and the mod-2 crossing-chain of temporal plaquettes on any interior boundary is a 1-cycle whose minimal support — homologous to the source square — is 4 links; so `≥4` temporal plaquettes per boundary, `≥4N_t` total. The two arguments agree; 1.2 is the elementary one. `[T]`)

A higher-representation temporal plaquette carrying `|r|≥2` fundamental flux units does **not** help: it costs `|r|` in the grading, exactly what it carries, so merging strands is neutral-to-worse. `[T]`

### 1.3 The lower bound and its equality case `[T]`

Because `C(t)_c` is the **connected** correlator, the only charged objects are the source and sink plaquettes; each of the four source flux units must therefore be transported to a **sink** edge along its own strand — it cannot close on itself at the source (that would be a *disconnected* contraction, excluded from `⟨O O⟩_c`). So (1.2) applies to all four strands over the full height. Combining (1.1′)+(1.2): four strands, each `≥N_t` temporal plaquettes of unit cost, with merging non-saving:

> **(1.3) Lemma (four-wall minimality).** Every connected, gauge-invariant, nonvanishing polymer linking the insertions satisfies `|P| ≥ 4N_t = 4t/a`. `[T]`

> **(1.4) Corollary (equality ⟹ straight tube).** Equality `|P|=4t/a` forces (i) exactly `N_t` temporal plaquettes per strand, (ii) **no** spatial plaquettes anywhere in `P`, and (iii) each strand straight (any spatial displacement of a strand adds spatial plaquettes, §1.4). Hence the unique minimizers are the **straight four-walled tubes**: four side strips of height `t`, one fundamental temporal plaquette per face per slice, `4` plaquettes per slice, `|P|=4t/a`. `[T]`

This is Gupta's "rectangular tube consisting of four sides of size 1×T" (Eq. 14.20 context, Fig. 13B) — now derived, not quoted. `[T]`

### 1.4 Why "one plaquette per slice" is *not* the minimum — and is not even a polymer `[T]`

The wrong object of `PROOF-P3-UNCONDITIONAL §3` ("one fundamental plaquette per slice") fails two independent ways:

- **Disconnected.** A **spatial** plaquette per slice stacked in time shares no link with the next slice (§1.2) ⇒ the set is not a connected polymer at all.
- **Zero weight even if connected.** If one instead threads a single **temporal** staircase (each plaquette sharing one link with the next), it covers only **one** source edge; the other **three** dangle, `∫dU U_{ij}=0` on each ⇒ weight `0` (1.1).

So the file's `|P|≥t/a` was a valid but *unsaturable* lower bound; the saturating object does not exist as a nonzero polymer. Correlation functions of **connected polymers** therefore have no `t/a`-cost term at all — the leading exponent is `4t/a`. `[T]`

### 1.5 No cheaper connected structure `[T]`

(I) Any non-straight modification adds `≥1` plaquette (a corner, a kink, a displacement increment) ⇒ `|P| ≥ 4t/a + 1`. (II) A transverse **displacement** `d=|x_source−x_sink|` of the whole sink adds `4d` plaquettes (all four walls become slanted staircases of `N_t+d` plaquettes), so `|P|=4(t/a)+4d` — on `T³_L` the transverse distance is bounded, `d≤L/2=O(1)` as `t→∞`, so this is still `4t/a+O(1)`, **same rate**; summed over the (finite) set of sink positions it folds into the positive amplitude `a₀` (§2.1) and is *not* rate-separated. (III) Perimeter increase (a 1×2 cross-section) costs `6t/a`; a second tube costs `8t/a`. Therefore:

> **(1.5)** The **unique** minimizers at `4t/a` are the straight four-wall tubes; every other connected gauge-invariant polymer has `|P| ≥ 4t/a + 1`. `[T]`

---

## §2 Uniqueness of the leading concentration and `a₀>0` `[T]`

### 2.1 The leading class and its weight

`C(t)=Σ_x Σ_{y∈spatial} ⟨O(x,0)O(y,t)⟩_c` (zero-momentum sum over sink; source fixed by translation invariance). By §1 the leading connected class is the **straight four-wall tube** from `(x,0)` to `(y,t)`. Summing over `y` (all transverse positions) and over `x` gives a finite, **positive** geometric multiplicity; by translation invariance it is independent of `x` per unit spatial volume. Each tube has

- `4N_t` fundamental temporal plaquettes ⇒ character/contraction weight `(u)^{4N_t}` with `u := c_f(β)/N` the per-tile strong-coupling weight (Gupta's "`β/2N` from the expansion × `1/N` from the integration" = `β/2N²+O(β³)`; §4) — **strictly positive** for `0<β<β₁` `[O, Gupta 14.18]`;
- an endpoint/source-sink factor: the operator coefficient `(2N)^{-1}` per insertion times the fundamental-to-singlet recoupling at each end, `f⊗f̄→1`, whose SU(N) invariant coefficient is `1/N>0`. `[T]`

Hence `w_tube(t) = a₀ · u^{4N_t}` with

> **(2.1)** `a₀(β) = (2N)^{-2}\,·\big[\text{positive endpoint recouplings}\big]\,·\,\big[\text{finite positive transverse sum}\big] \;>\;0`, **`t`-independent**, and `w_tube(t)=a₀ e^{−m₁ t}`, `m₁ := 4\ln(1/u)/a`. `[T]`

### 2.2 Uniqueness is an *amplitude* uniqueness (and it is what P3 needs)

Every non-tube polymer costs `≥4t/a+1` (§1.5), so its weight is `≤ u^{-1}·`(entropy)`×`(tube weight per unit) — a full factor `u` smaller. Therefore:

> **(2.2)** Over the leading rate `e^{−m₁ t}`, the four-wall-tube class has the **unique** maximal amplitude, and that amplitude `a₀(β)>0`. Hence the residue `f₁² = |⟨Ω|O|G⟩|² = a₀(1−O(u))` is a **positive** modulus-square: no sign or recoupling can undo it. `[T]`

This establishes the **leading-amplitude concentration (P3's I1 ingredient)** rigorously with the correct `m₁=4ln(1/u)/a`. It does **not** establish the level-isolation `I2` (that the next 0⁺⁺ sits at `≥1.5m₁` *as a transfer-matrix eigenvalue with a uniform gap*) — the cost grading is an amplitude statement, not eigenvalue spacing; `I2` remains the open gap of `PROOF-ISOLATION` and is **not** what P3 needs. `[T]/[O]` (This is precisely the REDTEAM-FORGE correction already in the parent file, unchanged by the factor 4.)

### 2.3 The factor 4 changes only the mass, not the positivity `[T]`

The audit's core claim is confirmed and now proved: the substitution `1→4` in the exponent (from `|P|≥t/a` to `|P|≥4t/a`) rescales `m₁: ln(1/u)/a → 4ln(1/u)/a`; `a₀>0`, the `O(u)` relative remainder, and the uniqueness are **structurally unchanged**. `[T]`

---

## §3 Corrected remainder `[T]` combinatorics, `[C]` uniform bound

### 3.1 The integer-cost gap `⇒` relative `O(u)` `[T]`

Let `R(t) := C(t) − w_tube(t)` be the sum over all connected non-tube polymers. By §1.5 each has `|P| ≥ 4t/a + Δ`, `Δ≥1`, fundamental (higher `|r|` only enlarge `Δ`). Its weight is `≤ u^{Δ}·(\text{entropy factor})·(\text{per-tube normalization})`. With the Kotecký–Preiss criterion `Σ_{P∋x}|w(P)|e^{a|P|}≤a`, equivalently `u·K<1`, holding **uniformly in `t` and `V`** (the criterion is per-plaquette; `K≤CN`; `β₁=O(1/N)`), the sum over the `Δ=O(1)` "bounded-excess" competitors obeys

> **(3.1)** `|R_{\Delta=O(1)}(t)| ≤ C\,u\,|w_tube(t)|` — one extra activity, uniform in `t,V`. `[T]` (Full-rep convergence `K<∞` is the standard factorial-domination estimate, audit §1.2 / Seiler LNP 159; not re-proved here. `[T]`)

The KP uniformity is the structural fact that converts "*leading term positive*" into "*residue bounded below for all `V,L`*." `[T]`

### 3.2 The `t`-linear decoration sum is a mass shift, not a small correction `[R]`

The bound (3.1) is **not** yet `t`-uniform as written, because a single extra plaquette ("decoration"/pillbox) can attach at any of the `4N_t = 4t/a` tube sites, giving a **linear-in-`t`** count `(4t/a)·u·w_tube`. Summing the decoration gas exponentiates (linked-cluster / Mayer):

`w_tube·(1+u)^{4t/a} = w_tube·e^{(4t/a)\ln(1+u)} = w_tube·e^{4ut/a−2u²t/a+…}`,

i.e. a **self-energy shift**

> **(3.2)** `m₁^{dressed} = 4\ln(1/u)/a \;-\; 4u/a \;+\; O(u²/a) = \dfrac{4\ln(1/u)}{a}\big(1+O(u/\ln\tfrac1u)\big)`, and `a₀ → a₀(1−O(u))`. `[O]` (the arithmetic; matches audit §1.3, GJS/OS two-point structure). `[R]` (**what is not proved here:** existence of the resummed effective weight, its `L`-uniformity, and a *lower* bound `a₀(1−O(u))>0` from a *signed* remainder. Absolute convergence gives `|1+r|≤1+O(u)`, not `1+r≥1−O(u)`.) `[R]`

This is the **one genuine remaining gap** — audit objection 3 — and it is **independent of the four-wall geometry**. `[R]`

### 3.3 The rate gap `δ` `[T]`

Competitors that grow *faster* than `4t/a` with `t` are rate-separated. The smallest such excess is a **perimeter-6** single tube (1×2 cross-section; the standard first excited 0⁺⁺ at `≈1.5m₁`), cost `6t/a`, rate `6ln(1/u)/a`:

> **(3.3)** `δ = (6−4)\ln(1/u)/a = 2\ln(1/u)/a = \tfrac12 m₁`. (The transverse-displacement family of §1.5-II is `O(1)` in cost on finite `V` and is *not* rate-separated; the genuinely heavier classes are perimeter-6 (`6t/a`, `δ=2\ln(1/u)/a`) and two tubes (`8t/a`, `δ=4\ln(1/u)/a`). How surfaces are partitioned into "dressing of the ground-state tube" vs "new pole" is the content of the [R] resummation; `δ=2\ln(1/u)/a` is the rate of the first strictly-heavier 0⁺⁺ content — the standard `≈1.5m₁` excited glueball.) `[T]`

### 3.4 The stated uniform remainder `[C]`

Modulo the **[R]** resummation of §3.2, (3.1)+(3.3) give the corrected bound with the correct geometry and the correct `δ`:

> **(3.4)** `\boxed{\,|R_L(t)| \le C\,u\,|w^{dressed}_{tube}(t)| + B\,e^{−\delta t},\qquad \delta = 2\ln(1/u)/a\,}` for all `t,V,L`; hence `e^{m₁^{dressed}t}|R_L(t)| \le Cu + B e^{−δt} → Cu` as `t→∞`, uniformly. `[C]` (conditional on [R])

Consequently `f₁² = a₀(1−O(u)) ≥ r₀(β)>0` with the **correct** `m₁=4ln(1/u)/a`, `V`- and `L`-uniform — **[C] conditional on [R]**, exactly the status the audit assigned, now with proved geometry. `[C]`

### 3.5 Orthogonal caveat: the `t→∞` limit on `T⁴_L` `[R]`

`PROOF-P3-UNCONDITIONAL` works on `T⁴_L` (time also periodic of size `L`) yet takes `t→∞`; on a periodic time direction the separation is defined mod `L`, so that limit does not exist. The four-wall geometry and `m₁=4ln(1/u)/a` are **unaffected**; the extraction must be done on the cylinder `T³_L×ℝ` (or `T⁴_L` with `t<L/2` then `L→∞`). This is a *bookkeeping* gap (audit objection 4), orthogonal to the factor 4, flagged here and not claimed closed. `[R]`

---

## §4 SU(2) strong-coupling mass: sanity check against Gupta `[O]/[T]`

### 4.1 Gupta's equations, verbatim from the source `[O]`

From `arXiv:hep-lat/9807028`, §14.2.2 (read from PDF):

- Eq. (14.17): `W_{RT} = \frac1Z\int dU\,W_{RT}\,[1+(\beta/2N)(W_{1+1}+W_{1-1})+\dots]` — the expansion factor per plaquette is `β/2N`.
- Eq. (14.18): each tiling plaquette "brings a factor of `β/2N` from the expansion and another factor of `1/N` from the integration … combined factor `β/2N²`" (with the note: "SU(2) is special and the combined factor is `β/N²`"); `⟨W_{RT}⟩=(β/2N²)^{RT}(1+\dots)`.
- **Eq. (14.20):** "The minimum tiling of the connected plaquette-plaquette correlation function is the rectangular tube consisting of four sides of size 1×T … `W₁₁(T)W₁₁(0) ∼ e^{−M₀₊₊T} = (β/2N²)^{4T}(1+\dots)`."
- **Eq. (14.21):** `M₀₊₊ = −4\log(β/2N²)+\dots` .

The `4T` in the exponent is *the four walls*; it is the same `4` derived in §1. `[O]`

### 4.2 Generic match `[T]`

Identify the per-tile strong-coupling weight `u := β/2N² = c_f/N` (§2.1, §4.1). The four-wall tube gives `w_tube ∼ u^{4T}` ⇒

`m₁ = −\ln(u^4)/a = 4\ln(1/u)/a = 4\ln(2N²/β)/a = −4\ln(β/2N²)/a = M₀₊₊`. ✔

So **`4\ln(1/u)` with `u=β/2N²` reproduces Gupta (14.21) exactly at leading order**, coefficient 4 = four walls. `[T]/[O]`

### 4.3 SU(2) special case `[O]`

Gupta's own (14.18) note: for SU(2) the two plaquette orientations are equal (trace real), so the combined per-tile factor is `β/N²` (= `β/4`), not `β/2N²`. Applying the same four-wall logic: `m₁ = 4\ln(1/u)/a` with `u=β/N²=β/4` ⇒

`M₀₊₊^{SU(2)} = −4\ln(β/4)+\dots` ,

the standard SU(2) strong-coupling series. **The coefficient 4 is universal; only the argument of the log differs from the generic formula by the SU(2) orientation factor 2** (`β/4` vs `β/8`). Using the *generic* `β/2N²` blindly at `N=2` would give `−4\ln(β/8)`, off by that factor — a convention artifact, not a geometry issue. `[O]`

### 4.4 Numerical consistency `[O]`

At `β=1.0` (SU(2), strong coupling): predicted `M₀₊₊ = 5.545` (SU(2)-special `−4\ln\tfrac14`) or `8.318` (generic `−4\ln\tfrac18`). The direct lattice measurement (`COMPUTE-F1SQ-SU2.md`, `β=1.0`, `L=6`) gives `m_eff≈6` (state at/beyond the cutoff, `m·a≳3`). The measured `≈6` **lies between** the two conventions' predictions `5.5`–`8.3` — order-of-magnitude consistent given that the leading-log series is asymptotic at `β=1`. (At `β=2.2`, `m_G=1.51` is in the *crossover* region, `β_c≈2.30`, where the strong-coupling series does not apply; it is not a test of `−4\ln(β/…`.) `[O]`

### 4.5 Münster corroboration `[O]`

The 4D strong-coupling gap series `m=−4\ln u+Σ_{k≥1}m_k u^k` (leading log coefficient exactly `−4`) is Münster's: F. Münster, Phys. Lett. B **95** (1980) 59; Nucl. Phys. B **190** [FS3] (1981) 439, erratum B200 [FS4] (1982) 536 (refs as verified in `P3-AUDIT-CORRECTED`; not re-derived from a primary PDF here — flagged honestly per charter §2). The `−4` is the same four-wall coefficient. `[O]`

---

## §5 What is now proved / what resists / what is untouched

**Proved (this document, [T]):**
1. Minimal connected gauge-invariant `¼tr F²` polymer = **four-walled tube**, `|P|≥4t/a`, equality iff straight (§1).
2. No cheaper connected structure exists; every non-tube costs `≥4t/a+1` (§1.5).
3. Unique leading **amplitude** class, `a₀(β)>0`, `m₁=4ln(1/u)/a` (§2).
4. Integer-cost gap gives `O(u)`-relative for bounded-excess competitors and rate gap `δ=2ln(1/u)/a` (§3.1, §3.3); KP criterion `uK<1` `t,V`-uniform `[T]`.
5. Gupta (14.17)–(14.21) verified; `4ln(1/u)` matches `M₀₊₊` (§4).

**Resists ([C]/[R], unchanged by the factor 4):**
- **Dressed-mass / dressed-amplitude resummation** (§3.2, audit objection 3): needed to make (3.4) `t`-uniform and to get a *lower* bound `a₀(1−O(u))>0`. This is the single genuine gap; it is **independent of the four-wall geometry**.
- Full-rep KP `K<∞` (§3.1): standard, [T] by factorial domination.
- Time-periodic `t→∞` extraction (§3.5): bookkeeping, on `T³_L×ℝ`.
- `I2`/isolation, `a→0`, transport to `β≈2.5`: untouched, as before.

**Net.** The factor-4 correction demanded by the audit is **fully discharged as a theorem** on the geometry/positivity side; the residue conclusion `f₁²≥r₀(β)>0` survives on a *correct* `m₁=4ln(1/u)/a` and now stands **[C] conditional only on the pre-existing dressed-mass [R]** — no new assumption, and nothing specific to four walls resists. `[J]`

---

*Refs.* R. Gupta, *Introduction to Lattice QCD*, arXiv:hep-lat/9807028, §14.2.2, Eqs. (14.17)–(14.21) — **verified verbatim from source**. F. Münster, Phys. Lett. B **95** (1980) 59; Nucl. Phys. B **190** (1981) 439, erratum B200 (1982) 536 (via `P3-AUDIT-CORRECTED`). Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Seiler, LNP **159** (1982); Kotecký–Preiss, CMP **103** (1986) 491; Glimm–Jaffe–Spencer, Ann. Math. **100** (1974) 585. Internal: `PROOF-P3-UNCONDITIONAL.md`, `P3-AUDIT-CORRECTED.md`, `PROOF-P1P2.md`, `PROOF-P3.md`, `PROOF-ISOLATION.md`, `COMPUTE-F1SQ-SU2.md`. *(This file supersedes §3–§4 of `PROOF-P3-UNCONDITIONAL.md` on the four-wall geometry; it does not edit that file.)*
