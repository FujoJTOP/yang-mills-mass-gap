# P3-FULLREP-KP — the Kotecký–Preiss hypothesis for the **full** representation sum, with explicit `K ≤ CN`

> 2026-11-24 · Closes audit objection 2 of `P3-AUDIT-CORRECTED.md` ("KP proved only on the fundamental sector"). Target: prove `Σ_{P∋x}|w(P)|e^{a|P|} ≤ a` over **all** nontrivial irreps, with an explicit `K ≤ CN`.
> **Labels:** `[T]` theorem (proved here or in cited standard work) · `[C]` conditional · `[O]` observation · `[R]` required-but-missing step. No `QED` on anything unproven (charter §2).

---

## 裁决 / VERDICT (up front)

| Question | Answer | Status |
|---|---|---|
| Is `Σ_{P∋x}|w(P)|e^{a|P|} ≤ a` provable over the **full** rep sum? | **YES.** | `[T]` |
| Explicit `K` in `|w(P)| ≤ C_g^{|P|}u^{|P|}`, `K = 2λ_aC_g e`? | **`K ≤ 80\,e^{2+β}N ≤ 6.1×10²\,N`** (with `Δ=20`). | `[T]` |
| What makes the rep sum converge? | `c_r ≤ \frac{β^{|r|}}{|r|!}\frac{e^{β}}{\dim r}` ⇒ `c_r\dim r ≤ \frac{β^{|r|}}{|r|!}e^{β}` ⇒ `Σ_r c_r\dim r\,e^{a|r|} < ∞`. | `[T]` |
| Does anything resist? | **Nothing essential.** Only bookkeeping: the constant is explicit but not optimized; and one normalization remark (§4.3, does not affect `K`). | `[O]` |

**One-line reason the audit gap is closed.** The factor that the audit worried about — "`c_r=O(β^{|r|})` alone does not control large-rep recoupling multiplicity" — is supplied by the **factorial** `1/|r|!` in the character coefficient, **together with** the pointwise bound `|χ_r(U)| ≤ \dim r` (valid because Haar measure is a *probability* measure). These combine to `|w(P)| ≤ ∏_p c_{r_p}\dim r_p ≤ ∏_p \frac{β^{|r_p|}}{|r_p|!}e^{β}`, where the offending `\dim r` from the contraction is **cancelled** by the `1/\dim r` inside `c_r`. The remaining `1/|r|!` dominates the (sub-exponential) entropy of representation labels. No Littlewood–Richardson/Weingarten refinement is needed — the crudest global contraction bound already suffices, and is in fact *better* than the per-link bound (which over-counts by `dim^4`, §3).

---

## §1 The standard character-coefficient bound `[T]`

### 1.1 Setup and claimed bound

Wilson plaquette weight, `U_p` the plaquette holonomy:
`e^{(\beta/N)\operatorname{Re}\operatorname{tr} U_p} = e^{(\beta/2N)(\operatorname{tr}U_p + \operatorname{tr}U_p^\dagger)} = Σ_{r} c_r(\beta)\,χ_r(U_p)`, `c_r(β) := ` coefficient.

**Lemma 1.1 (factorial decay, rep-uniform).** For every SU(N) irrep `r` and every `β ≥ 0`,
`\boxed{\;0 \le c_r(β) \le \frac{β^{|r|}}{|r|!}\,\frac{e^{β}}{\dim r}\;}`
where `|r|` = number of boxes of `r`, `\dim r` its dimension. Consequently
`c_r(β)\,\dim r \le \frac{β^{|r|}}{|r|!}\,e^{β}`. `[T]`

### 1.2 Proof (self-contained; the multiplicity route)

Write `F(U) := e^{(\beta/2N)(χ_f(U)+χ_{f̄}(U))}` (with `χ_f = \operatorname{tr}U`, `χ_{f̄} = \operatorname{tr}U^\dagger = \overline{χ_f}`). Taylor-expand:
`F = Σ_{m≥0}\frac{(β/2N)^m}{m!}(χ_f+χ_{f̄})^m = Σ_{m≥0}\frac{(β/2N)^m}{m!}\sum_{a+b=m}\binom{m}{a}χ_f^{\,a}χ_{f̄}^{\,b}.`

Each `χ_f^{\,a}χ_{f̄}^{\,b} = χ_{f^{⊗a}⊗f̄^{⊗b}} = Σ_s m_s(a,b)\,χ_s`, where `m_s(a,b) ≥ 0` is the multiplicity of the irrep `s` in the tensor product `f^{⊗a}⊗f̄^{⊗b}`. Hence, by orthonormality `∫χ_s\overline{χ_r}\,dU = δ_{sr}`,

`c_r(β) = Σ_{m≥0}\frac{(β/2N)^m}{m!}\sum_{a+b=m}\binom{m}{a}\,m_r(a,b).` `(★)`

Three elementary facts:

1. **Nonnegativity:** every term on the RHS of `(★)` is `≥0`, so `c_r(β) ≥ 0`. *(This also settles the small loose end flagged in `P3-AUDIT-CORRECTED.md` §4 — `c_r≥0` for all r is in fact a by-product of `(★)`; it is not otherwise needed.)* `[T]`
2. **Factorial cutoff:** `m_r(a,b) ≠ 0 ⟹ |r| ≤ a+b`. Reason: `f^{⊗a}⊗f̄^{⊗b}` contains only irreps with `≤ a+b` boxes (each fundamental/fundamental-bar factor carries one box; a tensor product cannot manufacture more boxes than it has factors). Thus the sum in `(★)` starts at `m = a+b ≥ |r|`.
3. **Dimension cap:** `m_r(a,b) \le N^{\,a+b}/\dim r`. Reason: `\dim(f^{⊗a}⊗f̄^{⊗b}) = N^{a+b}`; if `V = ⊕_s m_s V_s` then `m_r·\dim r \le Σ_s m_s\dim s = N^{a+b}`. `[T]`

Combining (2)+(3) into `(★)` and using `Σ_{a+b=m}\binom{m}{a} = 2^m`:
`c_r(β) \le \frac{1}{\dim r}\sum_{m\ge |r|}\frac{(β/2N)^m}{m!}\,2^m N^m = \frac{1}{\dim r}\sum_{m\ge|r|}\frac{β^m}{m!} \le \frac{β^{|r|}}{|r|!}\,\frac{e^{β}}{\dim r},`
since `\sum_{m\ge n}β^m/m! \le \frac{β^n}{n!}e^{β}`. This is Lemma 1.1. ∎

### 1.3 Relation to the literature `[O]`

- **Seiler, *Gauge Theories as a Problem of CQFT and Statistical Mechanics*, LNP **159** (1982), §II** establishes convergence of the character (high-temperature) expansion below a small-`β` radius via exactly the mechanism above: the coefficients decay like `β^{|r|}/|r|!` up to `\dim r`-factors, so that sums `Σ_r(\dim r)^k\,β^{|r|}/|r|!` converge for every `k`. Lemma 1.1 is a sharp, self-contained instance (`k=1`, with the `1/\dim r` from the multiplicity cap); the displayed constants are proved in §1.2 and are not quoted verbatim from a numbered formula.
- **Osterwalder–Seiler, Ann. Phys. **110** (1978) 440** use the same factorial structure for the uniform-in-volume cluster expansion.
- **Münster, Phys. Lett. B **95** (1980) 59 / Nucl. Phys. B **190** [FS3] (1981) 439 (erratum B **200** (1982) 536)** writes the 4-D gap series as `m = -4\ln u + Σ_{k≥1}m_k u^k` with `u = c_f` the fundamental coefficient; consistent with Lemma 1.1 at `|r|=1`: `c_f \le \frac{β}{\dim f}e^{β} = \frac{β}{N}e^{β}` (roughly `2e^β×` the exact `c_f \approx β/2N`; immaterial).

*Normalization caveat, see §4.3.* `[O]`

---

## §2 The full representation-sum estimate `[T]`

### 2.1 Polymer activity bound

For a polymer `P` (link-connected plaquette set, all reps `r_p` nontrivial), `|P| = Σ_{p}|r_p|`,
`w(P) = \Big(∏_{p}c_{r_p}\Big)\cdot I(P)\cdot(\text{operator coeffs}), \qquad I(P) := \int\prod_{e\in P}dU_e\prod_{p\in P}χ_{r_p}(U_p),`
the "link recoupling" value. The **key pointwise bound** is

**Lemma 2.1 (global contraction bound).** `|I(P)| \le ∏_{p\in P}\dim(r_p)`. `[T]`

*Proof.* `∏_e dU_e` is a **probability** measure (normalized Haar on each link). For `U \in \times_e SU(N)`, `|∏_p χ_{r_p}(U_p)| \le ∏_p|\operatorname{tr}D^{r_p}(U_p)| \le ∏_p\dim(r_p)` because `D^{r_p}(U_p)` is unitary (spectrum on the unit circle ⇒ `|\operatorname{tr}| ≤ \dim`). Hence `|I(P)| ≤ ∫∏_e dU_e\,∏_p\dim r_p = ∏_p\dim r_p`. ∎

(Operator coefficients are `O(N^{-p})`, bounded; they only improve the estimate and are absorbed in `C_g` below.)

**Corollary 2.2.** `|w(P)| \le ∏_{p}\big(c_{r_p}\dim r_p\big) \le ∏_{p}\dfrac{β^{|r_p|}}{|r_p|!}e^{β} = e^{β\,|P|_{\mathrm{pl}}}\dfrac{β^{|P|}}{∏_p|r_p|!} \le \big(e^{β}β\big)^{|P|}\Big/∏_p|r_p|!`, `[T]`
where `|P|_{\mathrm{pl}} = \#\{\text{plaquettes in }P\} \le |P|` was used. If `u := c_f \asymp β/(2N)` (literal normalization, §4.3), then `β \le 2Nu(1+O(β))`, giving the file's form
`|w(P)| \le \big(C_g\,u\big)^{|P|},\qquad \boxed{C_g \le 2N\,e^{β}} \;\;(=O(N)).` `[T]`

### 2.2 The KP sum

Sum over all polymers through a fixed plaquette `x`, graded by `|P|`, with grading parameter `a>0`:
`S_x(a) := Σ_{P\ni x}|w(P)|\,e^{a|P|} \le Σ_{m\ge1}\Big(\#\{\text{connected } m\text{-sets}\ni x\}\Big)\cdot Φ(a)^m,`
where we grouped the `m = |P|_{\mathrm{pl}}` plaquette positions and summed the reps per plaquette, defining the **per-plaquette rep activity**
`Φ(a) := Σ_{r\ne 0} c_r\,\dim r\;e^{a|r|} \le e^{β}\sum_{r\ne 0}\frac{β^{|r|}}{|r|!}e^{a|r|} = e^{β}\big(e^{βe^{a}}-1\big).` `(★★)`

(The factorization `Σ_{\{r_p\}}\prod_p(\cdots) = \prod_p\big(Σ_{r_p}\cdots\big)` is valid because the connection condition is on the plaquette *set*, independent of the labels.)

**Animal count.** `#\{connected m-sets through x\} ≤ λ_a^{\,m}` for a finite constant `λ_a = λ_a(Δ)` depending **only on the adjacency degree** `Δ = 20` (§4.1) — the standard animal-growth/tree-graph bound; one may take `λ_a = eΔ` up to a harmless `e^{m}` factor (absorbable into `λ_a`). Crucially `λ_a = O(1)` in `N`. Hence
`S_x(a) \le \frac{λ_aΦ(a)}{1-λ_aΦ(a)} \qquad\text{whenever } λ_aΦ(a) < 1.` `[T]`

### 2.3 The factorial beats the rep-label entropy

The concern "the sum over all irreps, with their dimensions and recoupling growth, might diverge" is quantified: at grading `|r|` there are at most `\sim e^{c\sqrt{|r|}}` SU(N) irreps (sub-exponential; `p(|r|)`-type partition growth for fixed `N`) each with `\dim r \le \binom{|r|+N-1}{N-1} = \mathrm{poly}_N(|r|)`. In `(★★)` all of this is already inside the summand `c_r\dim r \le β^{|r|}/|r|!e^β`, which carries **no residual `\dim` or per-rep constant**. Since `|r|!` dominates every sub-exponential factor (`\mathrm{poly}(|r|)e^{c\sqrt{|r|}}`), the series converges for every `a,β`. Concretely `Φ(a) = O(β)` as `β\to0` (the `|r|=1` levels `f,\bar f` dominate: `Φ(a) \approx 2βe^{a}e^{β}`), so `λ_aΦ(a) \to 0`. `[T]`

**Lemma 2.3 (full-rep KP hypothesis).** For `β ≤ β₁` with `λ_aΦ(a) ≤ \tfrac12`,
`Σ_{P\ni x}|w(P)|e^{a|P|} \le \frac{Φ(a)}{1-λ_aΦ(a)} \le 1 = a` (taking `a=1`). `[T]`

---

## §3 Recoupling / invariant-multiplicity bounds `[T]`

The task asks to bound the per-link invariant-tensor contraction (Littlewood–Richardson / Weingarten / `6j`) growth across all reps. Two routes; the second is what closes the proof.

**(a) Pointwise route (used, closing).** Lemma 2.1 bounds the *entire* contraction — every link integral and every recoupling multiplicity at once — by `∏_p\dim r_p`. This is the crudest possible bound, but it is **exactly matched** to the `1/\dim r` inside `c_r` (Lemma 1.1), so it loses nothing: `c_r\dim r \le β^{|r|}/|r|!e^β` has no `\dim`. This bypasses LR/Weingarten entirely. `[T]`

**(b) Per-link route (over-refined, strictly worse).** One can also bound link-by-link. At a link `e` with incident irreps `s_1,\dots,s_{k_e}` (`k_e ≤ 2(d-1)=6`), the number of singlets in `s_1⊗\dots⊗s_{k_e}` satisfies the standard LR/Racah bound `#\mathrm{singlets} \le \prod_{j\ne i}\dim s_j/\dim s_i` for **any** chosen `i` (multiplicity of `\bar s_i` in `⊗_{j\ne i}s_j \le \dim(⊗_{j\ne i}s_j)/\dim s_i`). Multiplying over the `≤4` links of each plaquette gives a per-plaquette factor between `\dim(r_p)^0` (if `p` is the designated `i` at two of its links — a balanced assignment exists for *surface* polymers since `\#links = 2\,\#plaquettes + L_b/2`) and `\dim(r_p)^4` (crudest, if none designated). Dropped into `c_r\dim^4` and summed, the factorial still wins, but only after paying `\dim(r)^3 \le |r|^{3(N-1)}`, i.e. a **worse** `C_g = O(N^3)`. The pointwise route (a) gives `C_g = O(N)` and is both simpler and sharper; route (b) is unnecessary. `[T]`

**Verdict on §3:** recoupling multiplicity is *controlled by the factorial decay* — via `c_r\dim r \le β^{|r|}/|r|!e^β`. No residual assumption. `[T]`

---

## §4 Conclusion `K ≤ CN`, adjacency, and what (if anything) resists

### 4.1 Adjacency degree `[T]`
Each link lies in `2(d-1) = 6` plaquettes; each plaquette has 4 links, each shared with `6-1=5` others ⇒ plaquette-adjacency degree `Δ = 4\,(6-1) = 20` (not 12). So `λ_a = eΔ = 20e ≈ 54.4`. (The file's original "degree ≤12, ≤3 plaquettes per link" was wrong; corrected here to `Δ=20`, `6` per link. It only shifts `K` by `20/12`, no structural effect.)

### 4.2 Explicit `K`
Assemble §2: `C_g \le 2N\,e^{β}`, `λ_a = O(Δ) = O(1)` (in `N`), `K := 2λ_aC_g e`. Then
`\boxed{K \le C_0\,e^{β}\,N\quad\text{with } C_0 = 4\lambda_a e \text{ explicit},\, \lambda_a=O(20);\ \text{in particular } K \le 10^{3}\,N\ \ (\beta\le 1/48).}` `[T]`
(With the generous `λ_a = 20e` one gets `K\le 80\,e^{2+β}N\le 6.1\times10^{2}N`; the precise animal constant is immaterial and only rescales `β₁`.) KP criterion `u(β)\,K < 1` ⟺ `u < 1/(2λ_aC_g e)`; with `K=O(N)` this is a **finite, `N`-dependent-small** condition identical on every `T^4_L` and every strip length, hence uniform in `V` and `L`. This is exactly the audit's target `K ≤ CN`. `[T]`

### 4.3 Normalization remark (does **not** affect `K`) `[O]`
In the *literal* action `S = βΣ_p(1 - N^{-1}\operatorname{Re}\operatorname{tr}U_p)` the expansion coefficient is `c_f = β/(2N)+O(β²)`, i.e. `u \asymp β/N`, **not** `u\asympβ`. If one takes the bare coefficient as `u`, then `uK \asymp (β/N)(CN) = Cβ`, and the KP radius is `β₁ = O(1)` — **N-independent**, consistent with Shen–Zhu–Zhu's `β<1/48`. If instead `β` is Hooft-rescaled so that `u\asympβ` (the file's convention), the radius reads `β₁ = O(1/N)`. Either way `K\le CN` and the theorem hold; only the stated numerical `β₁` shifts. Flagged so the master doc's `β₁=O(1/N)` is understood as convention-dependent (the literal-coordinate radius is N-independent).

### 4.4 What resists `[O]/[R]`
- **The KP full-rep bound itself: nothing resists.** `[T]` The audit's objection 2 is fully closed: `C_g=O(N)`, `K≤6.1×10²N`, unconditional, with the mechanism being `1/|r|!` × `1/\dim r` ∘ `\dim r`.
- **Not optimized:** the constant `6.1×10²` is explicit but crude (`a=1`, unoptimized animal constant). Lowering it only enlarges `β₁`; existence with `K≤CN` is what was needed. `[O]`
- **Genuinely open elsewhere (unchanged, out of scope here):** the dressed-mass/amplitude resummation (`P3-AUDIT-CORRECTED.md` §1.3, class `[R]`), the transfer-matrix residue bridge, `I2`/isolation, `a\to0`, transport to `β≈2.5`. This note closes **objection 2 only**.

---

## Sharpest line
The "large-representation recoupling" obstruction evaporates on contact with two standard facts: `|χ_r(U)|\le\dim r` (Haar is a probability measure) and `c_r\le β^{|r|}/(|r|!\,\dim r)` (multiplicity cap). They multiply to `c_r\dim r\le β^{|r|}/|r|!e^β`, where the `\dim r` that a naive recoupling count would blow up is exactly cancelled, leaving a factorial that dominates the sub-exponential rep-label entropy. Full-rep KP holds with `K\le 6.1\times10^{2}\,N` on every `T^4_L`; the audit gap is closed with elementary tools.

*Refs.* Kotecký–Preiss, CMP **103** (1986) 491–498 (DOI 10.1007/BF01211762); Seiler, LNP **159** (1982) §II; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440–471; Münster, Phys. Lett. B **95** (1980) 59, Nucl. Phys. B **190** [FS3] (1981) 439, erratum B **200** (1982) 536; Gupta, arXiv:hep-lat/9807028 §14.2.2 (factor-4 tube, used only for context); Shen–Zhu–Zhu, CMP **400** (2023) 805 / arXiv:2204.12737 (N-independent radius `β<1/48`). Internal: `PROOF-P3-UNCONDITIONAL.md` §1–§2, `P3-AUDIT-CORRECTED.md` §1.2 (objection 2).
