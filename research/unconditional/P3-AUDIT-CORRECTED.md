# P3-AUDIT-CORRECTED — adversarial audit of `PROOF-P3-UNCONDITIONAL.md`; the factor-4 tube, the rep-sum, and the dressed mass

> 2026-11-24 · Independent audit of our own claimed strong-coupling residue theorem (`f₁² ≥ r₀(β)>0` uniformly in `L`). Every external citation re-checked against source. **Labels:** `[T]` theorem (proved here or in cited standard work) · `[C]` conditional theorem · `[O]` observation / factual · `[NG]` no-go · `[R]` required-but-missing step / roadmap · `[J]` judgment. Charter §2: no `QED` on anything unproven.

---

## 裁决 / VERDICT (up front)

| # | Objection | Status | One-line reason |
|---|---|---|---|
| **1** | Minimal-cost tube is **4 walls**, `m₁=4ln(1/u)/a`, not 1-per-slice `ln(1/u)/a` | **VALID — decisive** | Gupta §14.2.2 Eq. (14.20)–(14.21): minimal tiling = "rectangular tube of four sides of size 1×T", `M₀₊₊=−4 log(β/2N²)`; Münster `m=−4 ln u+Σmₖuᵏ`. The "1 per slice" object is **not gauge-invariant** (weight 0), not merely unphysical. |
| **2** | KP hypothesis proved only on the fundamental sector, not the full rep sum | **VALID as written; CLOSABLE** | File bounds `C_g≤CN` only on `fundamental`; full-rep activity is controlled by factorial decay `c_r∼β^{|r|}/|r|!`×poly(`dim r`) — standard (Seiler LNP 159). SZZ does **not** help (it is not a rep expansion). |
| **3** | Dressed-mass / amplitude resummation asserted, not proved | **VALID — the real remaining gap** | KP absolute convergence does **not** produce the effective mass, its `L`-uniformity, or a positive amplitude floor. Standard-in-principle (GJS/OS two-point analysis) but unwritten. |
| **4** | 4D plaquette adjacency degree is 20, not 12 | **VALID — minor** | Each link lies in `2(d−1)=6` plaquettes; each plaquette has 4 links ⇒ degree `4·(6−1)=20`. File's "≤3 plaquettes meet / degree 12" is wrong by the same factor. Immaterial (absorbed in `K`). |

**Does `f₁²>0` survive the corrected (4-tube) derivation?** **YES.** `[J]` The factor 4 only rescales the mass, `m₁^{(0)}=4ln(1/u)/a` instead of `ln(1/u)/a`; the four-wall tube is still the **unique leading amplitude class**, still **positive** (`c_f>0` × singlet recoupling `f⊗f̄→1`), the remainder is still `O(u)` relative, and KP is still `t,V`-uniform. Hence `f₁²=a₀u^p(1−O(u))≥r₀(β)>0` survives with the same `r₀`-structure, **provided** objections 2 and 3 are discharged. `[C]`

**Is the residue conclusion currently a closed theorem?** **NO.** `[J]` As written, `PROOF-P3-UNCONDITIONAL.md` is **not** a complete unconditional proof: §3 uses the wrong geometry, §2 restricts to the fundamental sector, §4 asserts the resummation. The conclusion is very likely true and is a **genuine strong-coupling residue theorem once completed**, but the present document must carry **[条件 / 待重证]**.

**Is P3 new vs OS'78 / SZZ'23?** `[O]` **Not subsumed.** SZZ'23 (a stochastic/Langevin proof of LSI + Poincaré + clustering) gives the **mass gap** (a lower bound on the decay *rate*, i.e. an *upper* bound on correlations) and cannot give a **residue floor** — it is the wrong direction (a functional inequality never lower-bounds a single spectral weight). OS'78 give the uniform-in-volume strong-coupling expansion and the mass gap but never state `f₁²>0`. So `f₁²>0` is a **residue/atom refinement of OS'78's expansion** (derivable from it with essentially the same machinery) and is **not** derivable from SZZ. Net: **genuinely new as a stated theorem, not a new mechanism, and no simpler derivation exists.** `[J]`

---

## §1 Objection-by-objection adjudication

### 1.1 Objection 1 — the factor-4 tube geometry. **VALID, decisive.** `[T]`

**The file's claim** (§3): *"A connected `P` linking the time-`0` and time-`t` cells must meet every intermediate slice, so `|P|≥t/a`; equality forces exactly one fundamental plaquette per slice … `m₁=ln(1/u)/a`."*

**The correct minimum.** Verified against two independent standard sources read directly:

- **R. Gupta, *Introduction to Lattice QCD*, arXiv:hep-lat/9807028, §14.2.2, pp. 90–91, Eqs. (14.20)–(14.21).** Verbatim: *"The minimum tiling of the connected plaquette-plaquette correlation function is the rectangular tube consisting of four sides of size 1 × T"*; then `W₁₁(T)W₁₁(0) ∼ e^{−M₀₊₊T} = (β/2N²)^{4T}(1+…)`, and **`M₀₊₊ = −4 log(β/2N²) + …`**.
- **F. Münster, *Strong coupling expansions for the mass gap in lattice gauge theories*, Phys. Lett. B **95** (1980) 59; full paper Nucl. Phys. B **190** [FS3] (1981) 439–453 (erratum B200 [FS4] (1982) 536).** The 4-D gap series is `m = −4 ln u + Σ_{k≥1} m_k u^k` (leading log coefficient exactly `−4`).

Both agree: **the minimal connected polymer is a rectangular *tube* of four walls, four plaquettes per time slice, `|P|=4·(t/a)+O(1)`, hence `m₁^{(0)} = 4 ln(1/u)/a`** (with `u=u(β)≍β` the fundamental character coefficient). This is the `0⁺⁺` (A1g) plaquette–plaquette channel; the **two-glueball / wider-tube** states sit at `≥ 1.5 m₁` (perimeter-6 cross-section) or `2 m₁` (two tubes), i.e. at strictly higher rate. `[T]`

**Why the file's "one per slice" is *zero*, not merely loose.** `[T]` The operator insertion `O=Σ_{A1g}(2N)^{-1}(χ_f+χ_f̄)` is a single **spatial** plaquette; in the character expansion every one of its four links must be contracted (a lone `χ_f` integrates to zero). A chain with one plaquette per slice attaches to **one** source link and leaves the other **three** dangling ⇒ its link integral vanishes ⇒ **weight zero**. To contract all four source links, the flux tube's cross-section must be the full 4-link boundary loop of the source plaquette, extruded in time ⇒ **four walls**. So `|P|≥t/a` is true but far from saturated; the saturating object does not exist as a nonzero polymer. `[T]`

**Consequence for the file.** §3's identification of the leading polymer, its weight `a₀u^{t/a}` and its mass `m₁=ln(1/u)/a` are all wrong by the substitution `1→4` in the exponent. **The downstream files that repeat `m₁=ln(1/u)/a` are all off by this factor** (`E2-FLIP-MINORANT.md` line 72/101, `PROOF-IN11.md` "`a_c=1`", etc.); see §4.

### 1.2 Objection 2 — the KP representation sum. **VALID as written; closable (not fatal).** `[R]→[T]`

The file (§2) proves `Σ_{P∋x}|w(P)|e^{a|P|}≤a` using `C_g≤CN` *"on the fundamental sector"* only, while the polymer sum runs over **all** nontrivial irreps. As written this is incomplete: `c_r=O(β^{|r|})` alone does not bound the recoupling multiplicity of large reps. **However the gap is standard and closable**, and the conclusion is unaffected:

- At a link, `≤2(d−1)=6` plaquettes meet; the number of invariant contractions of reps `r₁,…,r₆` is `≤∏_i dim(r_i)`.
- The character coefficients decay **factorially**: `|c_r(β)| ≤ K' β^{|r|}/|r|!` (Schur/Bessel structure of the Wilson weight), and for fixed `N`, `dim(r)` is **polynomial** in `|r|` (degree `N(N−1)/2`).
- Hence `Σ_r dim(r)|c_r(β)|e^{a|r|} = Σ_r poly(|r|)β^{|r|}/|r|!·e^{a|r|} < ∞` for every `β`, giving a finite `C_g(β,N)` and the KP hypothesis for the **full** rep sum. `[T]`

This is the standard content of **Seiler, *Gauge Theories as a Problem of CQFT and Statistical Mechanics*, LNP **159** (1982)** (character expansion of the Wilson action, convergence below a small-`β` radius). The only correction needed in the file: replace "bounded by `CN` on the fundamental sector" by the full-rep factorial-domination estimate above; `K` grows but stays finite, `β₁` is unaffected in structure. `[T]`

**SZZ does not discharge this.** Shen–Zhu–Zhu (arXiv:2204.12737; CMP **400** (2023) 805) is a **stochastic/Langevin** proof of uniqueness, LSI, Poincaré, and exponential clustering for `SU(N)/SO(N)` at `|β|<1/(16(d−1))` (in `d=4`: `β<1/48`). It uses **no** character/polymer sum and proves **no** residue — it cannot be cited for the KP full-rep estimate. `[O]`

### 1.3 Objection 3 — the dressed-mass / amplitude resummation. **VALID — the genuine remaining gap.** `[R]`

The file (§4): *"resumming decorations of a tube … renormalizes the dressed mass `m₁` … the finite part multiplies `a₀` by `1+O(u)` … because the KP criterion is independent of backbone length, the decoration gas converges."* This is **asserted, not proved**. KP absolute convergence bounds the *sum of decorations*; it does **not** by itself give:
1. existence of a resummed/exponentiated weight (need the linked-cluster/Mayer expansion of the decoration gas along the tube),
2. an effective mass `m₁^{dressed}=4ln(1/u)/a+O(u/a)` with a **positive** gap to the next level,
3. **volume-uniformity** of that mass and amplitude floor,
4. a **positive** amplitude lower bound `a₀(1−O(u))>0` (bounds in absolute value give `|1+r|≤1+O(u)`, but a *lower* bound needs `r` real small, `1+r≥1−O(u)`).

The mechanism is exactly the **Glimm–Jaffe–Spencer two-point / particle-structure** analysis (Ann. Math. **100** (1974) 585) as adapted to strong-coupling gauge theory in OS'78: the self-energy shift is genuine — decorating each of the `4t/a` slices with a "pillbox" adds one activity `u`, giving `(1+u)^{4t/a}=e^{4ut/a}` ⇒ `m₁=4ln(1/u)/a − 4u/a + … = 4ln(1/u)/a·(1+O(u/ln(1/u)))`. `[O]` **Proving** this uniformly (items 1–4) is the load-bearing missing step; the file's §4 does not contain it. `[R]`

### 1.4 Objection 4 — plaquette adjacency degree. **VALID, minor.** `[T]`

A link is contained in `2(d−1)` plaquettes; in `d=4`, that is `6`. A plaquette has 4 links, each shared with `6−1=5` others ⇒ degree `4·5=20`. The file's "at a link `≤d−1=3` plaquettes meet" (should be `2(d−1)=6`) and "degree `≤12`" (should be `20`) are both wrong by a fixed factor. This enters only the animal-count `λ_a≤eΔ`, shifting `β₁` by a constant (`20/12`); **no structural consequence.** `[T]` (Factual, immaterial.)

---

## §2 The corrected 4-tube derivation and remainder `[C]`

Write `κ:=4ln(1/u)/a`, `u=u(β)=c_f(β)≍β` for `β<β₁`.

**Lemma 2.1 (four-wall minimality).** *Let `P` be a gauge-invariant link-connected polymer of plaquettes carrying nontrivial reps, linking a spatial-plaquette insertion at `t=0` to one at `t=T`. Then `|P| ≥ 4(T/a) + O(1)`; equality is attained only by a tube whose cross-section is the 4-link boundary of the source plaquette (four walls of `T/a` plaquettes each).*

*Proof.* Nonvanishing of the link integrals forces an even (contractible) covering of every link. The source plaquette's four links each require a partner; the minimal partner set is one temporal-plaquette wall per source link, and the walls must close through the sink (flux conservation of the boundary loop) ⇒ four walls, `4(T/a)` temporal plaquettes, `+O(1)` for source/sink. Any `k<4` walls leave `4−k` source links uncovered ⇒ zero weight (§1.1). ∎ `[T]`

**Corollary 2.2 (leading weight; positivity).** The minimal class is the set of four-wall tubes of arbitrary spatial path and transverse orientation. Each tube has `w_tube(t)=a₀ u^{4(T/a)}=a₀ e^{−κT}` (with `4T/a` temporal plaquettes, `T≡t`), where
`a₀(β) = (2N)^{−p}\,[∏_{wall plaquettes}c_f(β)]\,[f⊗f̄→1 recoupling] > 0`,
a product of positive fundamentals `c_f>0` (for `0<β<β₁`) and the positive singlet recoupling. `[T]`

**Corollary 2.3 (uniqueness is an *amplitude* uniqueness, not a rate uniqueness).** Every polymer that is not a four-wall tube costs `≥4T/a+Δ` with `Δ≥1`; it carries an extra factor `u^{Δ}` (same or higher rate, `O(u)`-smaller amplitude). Hence the tube is the **unique leading amplitude class**, with all other sectors `O(u)` relative — but *not* uniquely separated in rate (the `O(u)` amplitude sits at the same rate). `[T]` (This is precisely the point of `PROOF-P1P2 §3`, restated correctly.)

**Lemma 2.4 (remainder).** With the KP criterion of §1.2 holding uniformly in `t,V`, and after the decoration resummation of §1.3, there is `C` with
`|R_L(t)| ≤ C u\,|w_tube(t)| + B e^{−δ t}`, `δ=(ρ−1)κ`, `ρ≥3/2` (the perimeter-6 tube or two tubes),
so `e^{m₁^{dressed} t}|R_L(t)| ≤ Cu + O(e^{−δt})` uniformly in `V,L`. `[C]` (conditional on `[R]` §1.3).

**Theorem 2.5 (corrected residue theorem). `[C]`** Fix `N`, `a>0`. There is `β₁=β₁(N)>0` (KP radius, same order as `1/48`) and, for each `β<β₁`, a constant `r₀(β)=½c₁(β)u^p>0` (`c₁>0`, `p` the operator-content exponent, `p=0` for the normalized A1g plaquette) such that for every `T⁴_L`,
`f₁²_{1,L}=lim_{t→∞}e^{m₁(L)t}C_L(t) ≥ r₀(β)`, **uniformly in `L`**,
where now `m₁(L)=4ln(1/u)/a·(1+O(u/L))` is the **dressed** mass. `[C]`

**What the corrected proof still needs, precisely. `[R]`**
1. **Lemma 2.1** (four walls) — *done here* `[T]`.
2. **Full-rep KP estimate** (factorial domination) — *standard, Seiler LNP 159* `[T]`.
3. **Dressed-mass + amplitude resummation** (existence, `L`-uniformity, positive floor `1−O(u)`) — *genuine gap; GJS/OS two-point analysis* `[R]`.
4. **KL extraction** at fixed `L` then `t→∞`, then `L`-uniformity — bookkeeping, standard `[T]`.

Item 3 is the one that must be written out to promote `[C]` to `[T]`.

---

## §3 Novelty vs OS'78 / SZZ'23 `[J]`

- **OS'78** (Osterwalder–Seiler, Ann. Phys. **110** (1978) 440–471; ref verified) prove OS positivity, the **uniform-in-volume strong-coupling cluster/character expansion**, the infinite-volume limit, Wilson's confinement bound, and the **mass gap** (exponential clustering). They do **not** state `f₁²>0`: a mass gap is a lower bound on the decay **rate** (= an upper bound on correlations); the **residue floor** is a lower bound on a single spectral weight, strictly finer. `[O]`
- **SZZ'23** (Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; ref verified) prove LSI + Poincaré + exponential clustering (⇒ mass gap), `SU(N)/SO(N)`, any `d`, `β<1/48`, by **stochastic/Langevin** analysis — **not** a polymer/rep expansion, **no** residue. A functional inequality is the **wrong direction** to lower-bound a spectral atom. `[O]`
- **Therefore:** `f₁²>0` is **not** a corollary of SZZ, and is **not** stated by OS'78; it is a **residue refinement derivable from OS'78's expansion** (leading positive tube + reflection positivity + KL). It is **genuinely new as a stated theorem** and non-vacuous — the required **acid test** holds: the free massless field (`f₁²=0`, no binding, no tube) fails it, while the strong-coupling four-wall tube is a genuine YM positive. `[J]` But the mechanism (positive leading connected polymer) is **old** (OS'78), and the residue extraction costs the same as the P3 proof — **no simpler route exists**. `[J]`

> Honest framing: P3 is a **modest but real** addition — the atom/residue half that gap-methods (OS as usually quoted, SZZ) leave implicit or one-directional. It is **not** a new tool, and it does **not** touch the analyticity wall (`a→0`, `β→2.5`). `[J]`

---

## §4 Precise status `[J]`

**P3 (the residue floor `inf_L f₁²_{1,L}>0`):** **TRUE and survives the factor-4 correction**, class **(a) with the extra `[R]` §1.3** — not (b) (no simpler OS/SZZ route), not (c) (the geometry does **not** threaten the conclusion).
- **Proof as written (`PROOF-P3-UNCONDITIONAL.md`): NOT valid.** §3 uses the wrong (1-per-slice) tube; §2 handles only the fundamental sector; §4 asserts the resummation. Label **[条件 / 待重证]** (matching the 2026-11-24 header already added). `[J]`
- **Proof after correction:** genuine strong-coupling residue theorem `[T]` **iff** item 3 (§1.3, dressed-mass/amplitude resummation) is supplied; otherwise `[C]`.

**Citation audit (charter §2).**
- **Gupta** hep-lat/9807028 §14.2.2 pp.90–91, Eqs. (14.20)–(14.21) — **verified correct.** `[O]`
- **Münster** Phys. Lett. B **95** (1980) 59 / Nucl. Phys. B **190** (1981) 439, erratum B200 (1982) 536 — `m=−4 ln u+…` **verified.** The file cites only "Münster–Weisz" without a precise ref; give the above. `[O]`
- **Osterwalder–Seiler** Ann. Phys. **110** (1978) 440–471 — **verified.** `[O]`
- **SZZ** CMP **400** (2023) 805 / arXiv:2204.12737 — **verified** (stochastic, not rep expansion). `[O]`
- **Kotecký–Preiss** — file says "CMP **103** (1986) **19**"; correct is **CMP 103(3) (1986) 491–498** (DOI 10.1007/BF01211762). **Citation error to fix.** `[O]`
- **Brydges–Giffen–Durhuus–Fröhlich** NPB **275** (1986) 459–487 ("Surface representations of Wilson loop expectations") — **verified.** (Weingarten is a *different* author; do not attribute signs to him in this ref.) `[O]`
- **"`c_r(β)≥0` for all r" (§1 of the file)** — not needed and **not established** for `N≥3`; only `c_f(β)>0` (true for `0<β<β₁`) and `|c_r|` factorial decay are used. Drop the blanket claim. `[O]`

**Downstream corrections required.**
- `docs/YANG-MILLS-MASTER.md` §2.1 line 36–37 says *"P2 代价分层给 I2"* — this **contradicts** `PROOF-P3-UNCONDITIONAL.md §5`'s own REDTEAM-FORGE correction (the cost grading gives **I1 only**, not I2). Fix the master doc. `[O]`
- The mass value `m₁=ln(1/u)/a` is propagated (`E2-FLIP-MINORANT.md` §, `PROOF-IN11.md` `a_c=1`, etc.). **All must become `4ln(1/u)/a`.** In `PROOF-IN11.md` the qualitative conclusions (membership `{0⁺⁺,2⁺⁺}`, leading-order degeneracy, `m_{2⁺⁺}/m_{0⁺⁺}→1`) are **unchanged** because the factor 4 is common to both channels; only `a_c=1→4`. `[O]`
- `PROOF-P3.md` / `PROOF-P1P2.md` repeat the 1-per-slice geometry in their "sharpest line"; both need the 4-wall correction. `[O]`

**Bottom line.** The audit **found a real error in our own proof** (the factor-4 tube) — the highest-value outcome per charter §2 — but the error is *repairable*: the residue conclusion `f₁²≥r₀(β)>0`, uniform in `L`, **survives**. What does **not** survive is the label "unconditional/closed": the correct status is **[C] conditional on the dressed-mass resummation (§1.3) + the standard full-rep KP (§1.2)**. The result remains a **regime-limited strong-coupling residue theorem (I1 only)**; `I2`/isolation, `a→0`, and transport to `β≈2.5` are untouched, as before. `[J]`

*Refs.* Gupta, arXiv:hep-lat/9807028 §14.2.2 (Eqs. 14.20–14.21); Münster, Phys. Lett. B 95 (1980) 59, Nucl. Phys. B 190 (1981) 439, erratum B200 (1982) 536; Osterwalder–Seiler, Ann. Phys. 110 (1978) 440; Seiler, LNP 159 (1982); Kotecký–Preiss, CMP 103 (1986) 491; Shen–Zhu–Zhu, CMP 400 (2023) 805 / arXiv:2204.12737; Glimm–Jaffe–Spencer, Ann. Math. 100 (1974) 585; Brydges–Giffen–Durhuus–Fröhlich, NPB 275 (1986) 459. Internal: PROOF-P3-UNCONDITIONAL, PROOF-P3, PROOF-P1P2, PROOF-ISOLATION, PROOF-F1SQ, PROOF-IN11, E2-FLIP-MINORANT.
