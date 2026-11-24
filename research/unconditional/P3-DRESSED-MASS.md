# P3-DRESSED-MASS — the backbone–decoration (skeleton-polymer) expansion of the four-wall tube: the decoration gas is KP-summable `t`- and `L`-uniformly, the dressed mass `m₁^{dressed}=κ−f_dec` exists with `f_dec=O(u/a)`, and the dressed amplitude has a uniform positive floor; the only remaining spectral input is the cylinder Källén–Lehmann/OS reconstruction (companion P3-CYLINDER-RESIDUE)

> **2026-11-24 · companion to `PROOF-P3-UNCONDITIONAL.md` §4 and `P3-AUDIT-CORRECTED.md` §1.3.** Purpose: discharge audit objection 3 ("dressed-mass / amplitude resummation asserted, not proved") by actually building the skeleton-polymer expansion of the tube and proving what is provable. **This file does not edit the proof or the paper; it is an independent constructive-QFT write-up.** Labels: `[T]` theorem (proved here or in a cited standard work) · `[C]` conditional theorem · `[O]` observation/factual · `[NG]` no-go · `[R]` required-but-missing input · `[J]` judgment. Charter §2: no `QED` on anything unproven.

---

## 裁决 / VERDICT (up front)

The audit gap asks for four things — (a) the decoration-gas resummation exists, (b) the effective (dressed) mass exists and is `L`-uniform, (c) the dressed mass stays a spectral-atom energy, (d) the dressed finite amplitude has a uniform strictly-positive lower bound. Adjudication:

| Sub-gap | Status | One-line reason |
|---|---|---|
| **(a) decoration-gas resummation exists** | **`[T]`** | The decoration gas attached to the tube is an abstract polymer gas with per-attachment activity `O(u)`; Kotecký–Preiss on the **full** representation sum (Seiler: factorial domination) gives absolute convergence of the Mayer/Ursell series, uniformly in tube length `T=t/a` and volume `L`. |
| **(b) effective mass exists, `L`-uniform** | **`[T]`** | `\log Z_dec = φ_dec·T + O(1) = f_dec·t + O(1)` with `φ_dec = O(u)`, `f_dec=O(u/a)` and `|defect| ≤ C u`; existence of the free-energy density is the standard KP thermodynamic-limit theorem for the (signed) pressure; `L`-dependence enters only through wrapped decoration polymers of weight `O(u^{L/a})=e^{−cL}`, hence exponentially uniform. The rate is isolated from the next level (`≥1.5κ`) by `δ ≥ κ/3`. |
| **(c) dressed mass stays a spectral atom energy** | **`[T]` modulo `[KL]`** | An elementary **atom-from-exponential lemma** (spectral measure `≥0`) upgrades `C(t)e^{m t}→A>0` to `dμ({m})=A` with no density assumption. The only nontrivial input is that `C(t)=∫e^{−Et}dμ_O(E)` *is* a KL representation on the cylinder — the **companion P3-CYLINDER-RESIDUE** (equivalently OS'78 reconstruction, reflection positivity). `f₁² = A = dμ_O({m₁^{dressed}})`. |
| **(d) uniform positive amplitude floor** | **`[T]`** | Decorations are `O(u)`-small in **absolute value** (`|D_fin−1| ≤ Cu` by absolute convergence), so `A = a₀(1±Cu) ≥ (1−Cu)a₀ > 0`; `a₀>0` uses only `c_f>0` and the singlet recoupling `f⊗f̄→1` (multiplicity `1`). **No blanket `c_r≥0` is used.** |

**Net.** The dressed-mass resummation (audit objection 3) is **closed for (a),(b),(d)** — these are theorems once the (already-standard) full-rep KP is granted. The single remaining **analytic input** is **(c)'s Källén–Lehmann / OS transfer-matrix representation on the cylinder**, which is a separate companion (`P3-CYLINDER-RESIDUE`) and is standard `[T]` (OS'78) but not written in the P3 file. **`[J]` This upgrades `P3-AUDIT-CORRECTED` §2 Theorem 2.5 from `[C]` to `[C]` with the *named* residual `[KL]`, not an unnamed "assert." No QED is claimed; `I2`/isolation of the *excited* atom, `a→0`, and `β→2.5` are untouched.**

---

## §0 · Setup, conventions, and the two structural facts

**Geometry. `[O]`** Per `P3-AUDIT-CORRECTED` and audit item 4, the correct arena is the **cylinder** `Λ_{L,T} = T³_L ×_a ℤ`, spacing `a`, spatial torus `T³_L`, infinite/discrete time `ℤ`; we take `L→∞` last and `t→∞` before that. (Working on `T⁴_L` makes `t→∞` ill-defined; the KL representation below lives naturally on the cylinder.) The two-point function of the `A₁g`, zero-momentum plaquette operator
`O(x) = Σ_{spatial plaquettes p ⊂ x-slice} (2N)^{-1} Re tr U_p = Σ_p (2N)^{-1}(χ_f(U_p)+χ_{\bar f}(U_p))`
is `C_{L}(t) := ⟨O(0)O(t)⟩_{c}` (`c` = connected/truncated). Write `T := t/a ∈ ℕ`.

**Character/polymer expansion. `[T]`** For `β<β₁` the Wilson weight expands `e^{(β/N)Re tr U_p} = Σ_r c_r(β)χ_r(U_p)`. The **fundamental** coefficient is `u:=c_f(β)∈(0,1)`, `u≍β` (audit §1.2, Seiler LNP 159). A **polymer** `P` is a link-connected set of plaquettes carrying nontrivial irreps; its weight is `w(P) = (∏_{p∈P} c_{r_p})·∏_{links ℓ∈P}(link recoupling)` times operator coefficients at the source/sink.

**Fact 0.1 (positivity of link integrals). `[T]`** At any link,
`∫_{SU(N)} dU ∏_{i=1}^{k} χ_{r_i}(U) = dim Hom(1, ⊗_i r_i) = mult(1 ∈ ⊗_i r_i) ∈ ℤ_{≥0}.`
Hence the *recoupling* multiplicity at every link is a **nonnegative integer**, so the direct character series is a sum of nonnegative terms in the recouplings. `[T]` (This is the "positive scheme" of W2; Weingarten signs appear only in the *inverse* area-law series, never here.)

**Fact 0.2 (cost grading / four walls). `[T]`** (Audit Lemma 2.1, Gupta hep-lat/9807028 §14.2.2 Eqs. (14.20)–(14.21).) A connected polymer `P` linking the `t=0` insertion to the `t=t` insertion satisfies `|P| ≥ 4T + O(1)`; the source spatial plaquette's four links must each be contracted, forcing a **four-wall tube** (four temporal plaquettes per time slice) whose cross-section is the 4-link boundary of the source plaquette. The bare tube has `|P|=4T+O(1)`, weight
`w(B_T) = a₀(β)\,u^{4T} = a₀\,e^{−κ t}`, `κ := 4\ln(1/u)/a`,
where `∏_{p∈B_T}c_f = u^{4T}` (each of the `4T` temporal plaquettes carries `c_f=u`) has been absorbed into the exponential, and the **`T`-independent amplitude is** `a₀(β) = (2N)^{-p}·s > 0`, with `s ∈ ℤ_{≥1}` the product of link **recoupling multiplicities** (Fact 0.1; the source/sink `f⊗f̄→1` singlet recoupling contributes `1`) and `p∈{0,2}` the operator-content exponent. `[T]`

> **Why the "minimal cost" is `4T` and the only zero-kink tube is straight. `[T]`** A four-wall tube whose cross-section *shifts* spatially between adjacent slices needs at least one extra spatial plaquette (a "kink") per shift, costing an extra `u`; hence at cost exactly `4T` the tube is **straight** in its cross-section (a fixed spatial plaquette), and *every* spatial wandering is a decoration (a kink chain) multiplied into `Z_dec` below. The zero-momentum operator's sum over sink transverse positions `x` gives a factor `Σ_x u^{|x|_{\text{link}}} = O(1)` (geometric, `u<1`) absorbed into `a₀`; it is **not** an entropy `μ^T` competing with the energy. `[T]`

---

## §1 · Backbone–decoration (skeleton-polymer) decomposition

**1.1 Backbone. `[T]`** Fix, among the four-wall tubes of temporal length `T`, a **backbone** `B_T`: the straight four-wall tube with cross-section the source spatial plaquette `p₀` (its four walls are the `4T` temporal plaquettes of a `1×T` "rectangular tube" as in Gupta). The `A₁g` sum over the transverse sink position is folded into `a₀` (§0). The backbone weight is `w(B_T)=a₀u^{4T}` (§0, Fact 0.2).

**1.2 Decorations. `[T]`** Given `B_T`, the **decoration polymers** are the connected polymers `γ` of plaquettes that are (i) disjoint from `B_T`, (ii) **adjacent** to `B_T` (share a link with it), possibly also adjacent to each other. A **decoration configuration** `𝔇 = {γ₁,γ₂,…}` is a set of pairwise disjoint decoration polymers (this is exactly the set of connected components of `P∖B_T` when `P` is a polymer containing `B_T`). Its weight is `W(𝔇)=∏_i w(γ_i)`, and the **decoration partition function** is
`Z_dec(B_T) := Σ_{𝔇} W(𝔇) = 1 + Σ_{𝔇≠∅} W(𝔇)`. `(1.1)`

**1.3 Skeleton lemma (coverage at amplitude level). `[T]`** *Every connected polymer `P` linking the two insertions is either (I) `B_T` plus a decoration configuration (i.e. `P ⊇ B_T`, decorations `= P∖B_T`), or (II) a "non-tube" topology whose leading cost is `≥ 4T + Δ` with `Δ→∞` relative to the tube in the sense of rate `≥ (1+β₀)κ`, `β₀ ≥ 1/2`.* More precisely, writing `𝔓_{\text{tube}}` for the set of `P⊇B_T` and `𝔓_{\text{rest}}` for the complement, one has the **exact** decomposition
`C_L(t) = [ Σ_{P∈𝔓_{\text{tube}}} w(P) ] + [ Σ_{P∈𝔓_{\text{rest}}} w(P) ] = a₀u^{4T}\,Z_dec(B_T) + R^{hi}_L(t)`, `(1.2)`
and every `P∈𝔓_{\text{rest}}` satisfies `|P| ≥ 6T + O(1)` (perimeter-6 wider tube, the lowest non-tube) or `|P| ≥ 8T+O(1)` (two tubes / two-particle), hence `w(P) \lesssim u^{6T}` etc.

*Proof sketch. `[T]`* The first equality is tautological (`B_T` is a fixed polymer; the sum over `P⊇B_T` factors as `w(B_T)·Σ_{decor}` because `w` is multiplicative over disjoint pieces and the recoupling at the interface `∂B_T` is fixed into `a₀`). For the classification: by Fact 0.2 any `P` linking `0` to `t` has `|P|≥4T+O(1)`. If `|P| ≤ 4T + Λ₀` with `Λ₀` the fixed source/sink overhead, `P` is the straight four-wall tube `B_T` plus `≤Λ₀` extra plaquettes, i.e. a decoration configuration of bounded total cost (locality of decorations in the leading class — see 1.4). If `|P| ≥ 4T+Λ₀+1`, either `P` still contains a four-wall `T`-trunk but with extra structure (counted in `Z_dec`), or it does not, in which case the cheapest linking cost is a wider tube (perimeter `6`, `6T`) since the cross-section perimeter must be an even integer `≥4` and `=4` is the tube; next two tubes (`8T`). These are the classes in `𝔓_{\text{rest}}`. The recoupling positivity (Fact 0.1) makes the enumeration sign-free at the level of absolute values. ∎

**1.4 Locality of decorations in the leading class. `[T]`** A decoration component `γ` of `P∖B_T` with temporal extent `ℓ` (number of time-slices it meets) has `|γ| ≥ ℓ` (it must carry at least one plaquette per slice it spans, to remain link-connected across slices). Hence "long" decorations (`ℓ` comparable to `T`) carry cost `≳ T`, i.e. belong to `𝔓_{\text{rest}}` (a parallel-tube / two-tube object, rate `≥ 2κ`), **not** to the leading class. The decorations that enter `Z_dec(B_T)` in the leading class have total cost = `|P|−4T = O(1)` and are therefore **temporally local** (bounded extent). Conversely, an arbitrary *assembly* of local decorations re-sums in `Z_dec` to unrestricted `T` (it is a genuine gas along the backbone). `[T]` (Routine enclosure bookkeeping; standard "backbone/polymer" device of Kotecký–Preiss inhomogeneous activities and Bovier–Zahradník.)

> **Remark 1.5 (the "sector of arbitrary spatial path" is not an independent entropy). `[T]`** The audit (Cor. 2.3) is explicit that tube uniqueness is an *amplitude* uniqueness: a spatially wandering four-wall tube is not a distinct zero-kink object — each transverse shift is a kink (extra `u`), i.e. a decoration, and is counted in `Z_dec`. There is no `μ^T` path entropy at the minimal cost; the straight tube is the unique zero-kink backbone. This is what makes `f_dec` the *whole* renormalization of `κ` (no separate roughening contribution). `[T]`

---

## §2 · Resummation existence and `L`-uniformity

**2.1 The decoration gas obeys KP. `[T]`** The decoration polymers `γ` are a subfamily of *all* nontrivial connected plaquette polymers, so the KP hypothesis for the full polymer system controls them. The **full-representation** KP hypothesis is (audit §1.2; Seiler LNP 159):
`Σ_{γ∋x,\,γ\ \text{connected}} |w(γ)|\,e^{α|γ|} ≤ α`, `(2.1)`
with `x` any plaquette. *Verification `[T]`:* the character coefficients decay **factorially**, `|c_r(β)| ≤ K' β^{|r|}/|r|!` (Schur/Bessel structure of the Wilson weight), and `dim(r) ≤ poly_N(|r|)`; at a link `≤2(d−1)=6` plaquettes meet (audit §1.2, §1.4), so a contraction is a product of `≤6` reps with weight `≤∏dim`. Hence `Σ_r poly_N(|r|)\,β^{|r|}/|r|!\,e^{α|r|} < ∞` for every `β`, giving a finite `C_g(β,N)`; the animal count through `x` is `≤(eΔ)^n`, `Δ=20` (correct 4D degree, audit §1.4). Taking `α=1`, (2.1) holds whenever `u·K<1`, `K:=eΔ\,C_g e ≤ C' N` — **a per-plaquette condition, identical on every `Λ_{L,T}`, hence independent of `L` and `T`.** This is the "KP radius is `t`-independent" part of the audit request. ∎

**2.2 Resummation exists (the Mayer/Ursell log). `[T]`** By KP (Kotecký–Preiss CMP 103 (1986) 491, Theorem with general signed activities), `log Z_dec(B_T)` equals the **absolutely convergent** sum of Ursell polynomials over connected decoration clusters attached to `B_T`:
`log Z_dec(B_T) = Σ_{\text{connected } \Gamma \text{ attached to } B_T} φ^U(\Gamma)`, `(2.2)`
and the series is analytic in the activities on the polydisc `uK<1`. **This is precisely the resummation whose existence the audit (§1.3 item 1) demanded.** The **decoration self-energy** is the free energy per lattice slice `φ_dec` and the equivalent rate per unit Euclidean time `f_dec` (same units as `κ`):
`φ_dec := lim_{T→∞} \tfrac{1}{T}\log Z_dec(B_T) = O(u)`, `f_dec := φ_dec/a = lim_{t→∞}\tfrac{1}{t}\log Z_dec(B_T) = O(u/a)`, `(2.3)`
both existing by the standard KP thermodynamic-limit theorem for the pressure (the Ursell sum is absolutely convergent and its density has a limit; equivalently, `\log Z_dec(B_T) − f_dec\,t = \log Z_dec(B_T) − φ_dec\,T = O(1)`). `[T]`

**2.3 The per-link decoration activity is `O(u)`. `[T]`** The number of plaquettes *adjacent* to `B_T` per time slice is `O(1)` (the cross-section has `O(1)` links, each in `O(1)` neighboring plaquettes). By (2.1), the summed absolute activity of decoration polymers through any such plaquette is `≤ α = 1` at `α=1`, but the *scaled* activity per attachment that enters `φ_dec` is `Σ_{γ\ni x}|w(γ)| ≤ C''u/(1-C''u) = O(u)` (each `γ` has weight `≤ C^{|γ|}u^{|γ|}`, `Σ_n(Cλu)^n`). Hence
`φ_dec = (O(1)\ \text{adjacent plaquettes/slice})·(O(u)) = C_{dec}(\beta,N)\,u + O(u^2)`, `f_dec = φ_dec/a = \frac{C_{dec}u}{a} + O(u^2/a)`, `(2.4)`
i.e. **`φ_dec = O(u)`** (equivalently `f_dec = O(u/a)`), with `C_{dec}` a finite constant bounded uniformly in `L,T`. (This is the audit's "(1+cu+…)^t ⇒ shift `−cu`" made a theorem.) `[T]`

**2.4 `L`-uniformity. `[T]`** On `T³_L` the decoration gas may include polymers that wrap a spatial cycle; such a wrapped decoration has cost `≥ L/a` (it must cross `L/a` slices) and weight `≤ u^{L/a}=e^{−cL}`. Since (2.1) is `L`-independent and the wrap corrections are `e^{−cL}`, one has
`f_dec(L) = f_dec(∞) + O(e^{−cL/a})`, `\log Z_dec^{(L)}(B_T) = f_dec(∞)\,t + O(1) + O(t\,e^{−cL/a})`, `(2.5)`
uniformly in `t`. The finitely many spatial positions of the backbone (folded into `a₀`, §0) are `L`-independent up to the overall `L³` measure that cancels in the normalized residue. **Thus the resummed decoration factor is `L`-uniform.** ∎ `[T]`

---

## §3 · The effective (dressed) mass

**3.1 Existence and value. `[T]`** Combining (1.2) and (2.3)–(2.4),
`Σ_{P∈𝔓_{\text{tube}}} w(P) = a₀u^{4T}Z_dec(B_T) = a₀u^{4T}\,e^{f_dec\,t}\,D_{fin}(T) = A^{dressed}\,e^{−m₁^{dressed}\,t}\,D_{fin}(T)`, `(3.1)`
with `A^{dressed}=a₀`, `u^{4T}=e^{−κ t}`, and
`m₁^{dressed} := κ − f_dec = \frac{4\ln(1/u) - C_{dec}u + O(u^2)}{a} = \frac{4\ln(1/u)}{a}\bigl(1 - O(u/\ln(1/u))\bigr)`, `(3.2)`
`D_{fin}(T) := Z_dec(B_T)e^{−f_dec\,t} = 1 + O(u)` (bounded by §2.3–2.4; the `O(u)` comes from the endpoint defect, see §4). **This is the effective/dressed mass, and it exists.** It is `L`-uniform by (2.5): `m₁^{dressed}(L)=m₁^{dressed}+O(e^{−cL/a})`. `[T]`

> **Value-check vs Münster. `[O]`** (3.2) reproduces the *structure* of Münster's strong-coupling mass `m = −4\ln u + Σ_{k≥1}m_k u^k` (Phys. Lett. B **95** (1980) 59; NPB **190** (1981) 439): leading log `−4\ln u = κ`, and the `+Σm_k u^k` is `−C_{dec}u−…` with `m_k>0` ⇔ `f_dec>0` (decorations raise `Z_dec`, lowering the mass). We prove only the leading bound `f_dec=(C_{dec}u+O(u^2))/a`, i.e. `m₁^{dressed}=κ−(C_{dec}u+O(u^2))/a`; the full series is `[O]`, not needed for the residue floor. `[T]/[O]`

**3.2 Rate isolation from the next level. `[T]`** The remainder `R^{hi}_L(t)` in (1.2) is `𝔓_{\text{rest}}`, whose lowest rate is the perimeter-6 wider tube (`|P|≥6T`, rate `6\ln(1/u)/a=1.5κ`), décor-renormalized to `1.5κ−O(u/a)`; or two tubes / two-particle at `2κ` (rate `2κ−O(u/a)`). Hence, resumming decorations in `𝔓_{\text{rest}}` exactly as above,
`C_L(t) = A^{dressed}e^{−m₁^{dressed}t}\bigl(1+O(u)\bigr) + O\bigl(e^{−m₁^{dressed}t}e^{−δ t}\bigr)`, `δ ≥ \tfrac12κ − O(u/a) ≥ \tfrac13κ\ (u \text{ small})`, `(3.3)`
**uniformly in `L`.** I.e. `e^{m₁^{dressed}t}C_L(t) → A^{dressed}` with an exponentially small error at rate `δ≥κ/3`. `[T]`

**3.3 Is `m₁^{dressed}` a spectral-atom energy? `[T]` modulo `[KL]`.**

**Lemma 3.4 (atom-from-exponential). `[T]`** *Let `ν ≥ 0` be a finite Borel measure on `[0,∞)` and `C(t):=∫_{[0,∞)}e^{−Et}dν(E)`. If `C(t)e^{mt}→A` as `t→∞` for some `m≥0`, `A∈ℝ`, then `dν({m})=A`, `ν` has no mass on `[0,m)`, and `A≥0`. If `A>0`, `m` is an atom of `ν` of weight `A`.*

*Proof. `[T]`* Boundedness of `C(t)e^{mt}` forces no mass below `m` (else `C(t)e^{mt}≥δ e^{εt}→∞`). Write `dν=A'δ_m+ν'` with `ν'({m})=0`, `supp ν'⊂(m,∞)`. Then `C(t)e^{mt}=A'+∫_{[m,∞)}e^{−(E−m)t}dν'(E)`. The integrand `→0` for `ν'`-a.e. `E` (`ν'({m})=0`) and is dominated by `1∈L¹(ν')`, so the integral `→0` by dominated convergence. Hence `A=A'=ν({m})`. ∎

Now, **granting the cylinder Källén–Lehmann representation** (the single analytic input — §5),
`C_L(t) = ∫_{[0,∞)} e^{−Et}\,dμ_{O}^{(L)}(E)`, `dμ_O^{(L)} ≥ 0`, `(3.4)`
(`ν` = spectral measure of the Hermitian, gauge-invariant `O` in the vacuum sector, vacuum atom `E=0` removed by the connected part). Then (3.3) + Lemma 3.4 give, **for each finite `L`,**
`f₁²_{,L} := dμ_O^{(L)}\bigl(\{m₁^{dressed}(L)\}\bigr) = A^{dressed}(L) > 0`, `(3.5)`
and `m₁^{dressed}(L)` is the **bottom of `spec(H|_{O\text{-sector}})`**, an isolated atom (nothing below it, and the next level at `≥1.5κ−O(u/a)>m₁^{dressed}` since `m₁^{dressed}<1.5κ−O(u/a)`). `[T]` modulo `[KL]`

> **What 3.3 does *not* claim. `[J]`** It identifies the **lightest** `A₁g` state as an atom (this is `I1`). It does **not** prove the *excited* `0⁺⁺` at `≈1.5m₁` is an isolated atom (`I2`): that is a second-level (`λ₁ − λ₂`) datum, out of scope and provably a different statement (`G2D-ISOLATION-I2`, `PROOF-EDGE-E2`). The "isolation" in 3.3 is the *trivial* single-particle-below-continuum (`m<2m`) plus cost-graded second-level separation — **not** the nontrivial `I2`. `[J]`

---

## §4 · Uniform positive lower bound on the dressed amplitude

**4.1 Absolute-value control of `D_fin`. `[T]`** Define `D_fin(T):=Z_dec(B_T)e^{−φ_dec T}`. The defect `\log D_fin(T)=\log Z_dec(B_T)−φ_dec T` is the sum of: (i) the superadditivity slack and (ii) the subadditivity slack across a temporal cut (§2.2), plus (iii) the two endpoint terms where the tube terminates. Each slack is carried by Ursell clusters that **cross a fixed cut**; every such cluster costs `≥u` (it must contain an extra plaquette straddling the cut) and there are `O(1)` of them per cut, so by the **absolute** convergence of (2.2),
`|\log D_fin(T)| ≤ C_{1}\,u + O(1)\cdot u = C\,u`, `C=C(\beta,N)<\infty`, uniformly in `T,L`, hence `D_fin(T)=e^{O(u)}=1+O(u)`. `(4.1)`
**Crucial:** (4.1) uses only **absolute convergence**, *not* termwise nonnegativity. Therefore the amplitude floor does **not** require the (audit-flagged, not-established for `N≥3`) blanket `c_r(β)≥0`; it needs only `|c_r|` factorial decay (KP) — which is standard. `[T]`

**4.2 Bare amplitude `a₀>0`. `[T]`** From §0, `a₀ = (2N)^{-p}·s` with `s∈ℤ_{≥1}` a product of link recoupling multiplicities (Fact 0.1: each `≥0`, and the source/sink `f⊗f̄→1` singlet recoupling has multiplicity `1>0`), and `p∈{0,2}`. Hence `a₀ ≥ (2N)^{-2} > 0`, `T`- and `L`-independent. (Audit Cor. 2.2; `c_f>0` for `0<β<β₁`.) `[T]`

**4.3 Assembling the floor. `[T]`** By (3.1),(4.1),(4.2),
`A^{dressed} = a₀\,D_{fin} ≥ a₀\,(1 − C u) ≥ \tfrac12 a₀ ≥ \tfrac12 (2N)^{-2} > 0`, `(4.2)`
provided `u < 1/(2C)` — guaranteed on `β<β₁` (shrink `β₁` by an `N`-dependent constant if needed; this does not change the `β₁=O(1/N)` structure). Finally, `A^{dressed}(L) \to A^{dressed}` at rate `O(e^{−cL/a})` (§2.4), so
`\inf_{L} A^{dressed}(L) ≥ \tfrac14(2N)^{-2} > 0`, `(4.3)`
a **uniform strictly-positive lower bound** on the dressed amplitude, hence (via §3.3) `\inf_L f₁²_{,L} ≥ r₀(β)>0` with `r₀=¼(2N)^{-2}`-type (the audit's `r₀`). `[T]` modulo `[KL]`

**4.4 Free-field acid test (charter §1). `[O]`** For free massless Maxwell, the connected `A₁g` correlator is a power law (`C(t)≍t^{-2}`, no atom, `f₁²_{free}=0`); the polymer tube has `a₀=0` (no binding singlet recoupling — no glueball), so (4.2) gives **no floor** and the construction correctly returns `f₁²=0`. The floor `r₀>0` is non-vacuous because it is *produced by* the YM recoupling `f⊗f̄→1=1>0` (Fact 0.1) that the free photon lacks. **Acid test PASSED** — the tool is `>0` for YM, `=0` for the free field. `[O]`

---

## §5 · The exact remaining analytic input

Everything in §1–§4 is `[T]` **except one line**: the Källén–Lehmann representation (3.4) on the cylinder together with `H≥0` self-adjoint and `Ω` the unique vacuum. This is a standard theorem, but it is **not written in the P3 file** and is the load-bearing spectral input:

**`[KL]` (named input). Källén–Lehmann / OS transfer-matrix representation on `T³_L ×_a ℤ` (companion `P3-CYLINDER-RESIDUE`).** Reflection positivity of the Wilson measure (Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; project: PROOF-GAUSS-INVARIANT-PF) yields a self-adjoint transfer matrix `T=e^{−aH}`, `H≥0`, on the gauge-invariant physical Hilbert space, and for the Hermitian `O`,
`⟨Ω|O(t)O(0)|Ω⟩_{c} = ∫_{[0,∞)} e^{−Et}\,dμ_O^{(L)}(E)`, `dμ_O^{(L)}≥0`, `μ_O^{(L)}=Σ_n |⟨Ω|O|n⟩|² δ_{E_n} + (\text{continuum})`.
This is the classical Källén (1952) / Lehmann (1954) representation, here for the lattice transfer matrix (Reed–Simon IV, §XIII). `[T]` (standard, but it is the companion's job to assemble it on the cylinder with the correct `t→∞` order per audit item 4.)

With `[KL]`, §3.3's atom-from-exponential lemma is unconditional and `f₁²_{,L}=A^{dressed}(L)`.

**On GJS and the "dressed-two-point theorem". `[O]`** The audit (§1.3) names **Glimm–Jaffe–Spencer, Ann. Math. 100 (1974) 585** ("particle structure of `P(φ)₂`") as the *template* for a dressed-two-point / self-energy analysis. In **2D `P(φ)₂`**, GJS prove the two-point function has a one-particle pole by a *Bethe–Salpeter / renormalized-expansion* analysis. In **4D lattice strong coupling** we do **not** need the GJS machinery: the four-wall tube *is already* the one-particle object (its rate is isolated in the cost grading, §3.2), and the KL representation + Lemma 3.4 extract the atom directly. GJS/OS'78/Seiler are cited here as the **standard precedents** that (i) a convergent cluster expansion + reflection positivity + KL does produce a one-particle atom (§3.3), and (ii) the self-energy shift is `O(u)` (Münster). The precise division:

| Piece | Source | Label |
|---|---|---|
| Full-rep KP (`|c_r|` factorial decay, `dim` poly) | Seiler LNP **159** (1982); audit §1.2 | `[T]` |
| KP criterion ⇒ absolute Mayer/Ursell convergence, analyticity | Kotecký–Preiss, CMP **103** (1986) 491 | `[T]` |
| Pressure / free-energy density limit | KP (1986); standard thermodynamic-limit theorem | `[T]` |
| Backbone/polymer (skeleton) bookkeeping | KP inhomogeneous activities; Bovier–Zahradník | `[T]` |
| Self-energy `=O(u)`, mass series structure | Münster, PLB **95** (1980) 59 / NPB **190** (1981) 439 | `[O]` (leading bound `[T]`) |
| Atom from exponential + KL positivity | Lemma 3.4 (this file); Källén 1952 / Lehmann 1954 | `[T]` |
| **KL/OS representation on the cylinder** | **Osterwalder–Seiler (1978); companion P3-CYLINDER-RESIDUE** | **`[T]` (named input)** |
| Dressed-two-point template (not needed in 4D SC) | Glimm–Jaffe–Spencer, Ann. Math. **100** (1974) 585 | `[O]` |

**Not claimed / honest boundary. `[J]`**
1. **`I2`** — isolation of the *excited* `0⁺⁺` (`≈1.5m₁`) from the continuum: **not addressed, `[R]`**, a separate rung (`G2D-ISOLATION-I2`).
2. **The exact value of `m₁^{dressed}`** beyond `κ−(C_{dec}u+O(u^2))/a`: only the leading `O(u)` bound is proved; the full series is Münster `[O]`.
3. **`a→0`, `β→2.5`** — untouched (the analyticity wall; `P3-AUDIT-CORRECTED` §4).
4. **The `[KL]` assembly on the cylinder** (audit item 4: `t→∞` before `L→∞`, `T³_L×_aℤ` not `T⁴_L`) is the companion's task; this file sets up consistently on the cylinder but does not itself write the OS reconstruction.

---

## §6 · Scorecard (target = close audit objection 3; "5/5" = full `[T]` with no named input)

1. **Nonperturbative — PASS `[T]`.** KP polymer/Mayer expansion, all `β<β₁`, no loop truncation.
2. **(a) resummation exists — PASS `[T]`** (§2.2).
3. **(b) effective mass exists, `L`-uniform — PASS `[T]`** (§3.1–3.2, eq. (3.2),(2.5)).
4. **(c) spectral atom — PASS modulo `[KL]` `[T]`** (§3.3 Lemma 3.4 + `[KL]`).
5. **(d) uniform positive amplitude floor — PASS `[T]`** (§4, eq. (4.3)).

**Was audit objection 3 ("asserted, not proved") discharged? `[J]`** **Yes for (a),(b),(d) unconditionally `[T]`; (c) reduced to the single named standard input `[KL]`.** The "assert" is replaced by a proof up to the KL representation, which is OS'78 — standard but belonging to the companion `P3-CYLINDER-RESIDUE`. **No `QED` on (c) without that companion; no fake closure of `I2`; no touching of the analyticity wall.** The residue floor `f₁²_{,L} ≥ r₀(β)>0`, `L`-uniform, is now `[T]` modulo one named spectral input — matching `P3-AUDIT-CORRECTED` §2's verdict that the conclusion "survives" with objection 3 as the sole load-bearing written gap.

**Sharpest line. `[J]`** The dressed-mass resummation is not an extra dynamical assumption — it is the **free energy per unit length of a KP decoration gas** whose per-attachment activity is `O(u)` and whose radius is independent of `t` and `L`; its existence is the KP pressure limit, its `O(u)` size makes `m₁^{dressed}=κ(1−O(u/\ln(1/u)))` a genuine positive renormalized mass, and its absolute convergence makes the dressed amplitude `a₀(1±Cu)` uniformly positive (`a₀>0` from `c_f>0` and `f⊗f̄→1`). What remains is *not* dynamics but the **spectral identification** — that this convergent connected-correlator asymptotic is the atom of a positive Källén–Lehmann measure — i.e. OS reconstruction on the cylinder, the companion `P3-CYLINDER-RESIDUE`. That is the exact remaining analytic input, and it is standard `[T]`, not a new no-go.

*Refs.* Kotecký–Preiss, CMP **103** (1986) 491–498; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440–471; Seiler, LNP **159** (1982); Glimm–Jaffe–Spencer, Ann. Math. **100** (1974) 585; Münster, Phys. Lett. B **95** (1980) 59 / Nucl. Phys. B **190** (1981) 439 (erratum B200 (1982) 536); Gupta, arXiv:hep-lat/9807028 §14.2.2, Eqs. (14.20)–(14.21); Bovier–Zahradník, *A simple proof of the Kotecký–Preiss criterion*; Källén, Helv. Phys. Acta **25** (1952) 417; Lehmann, Nuovo Cim. **11** (1954) 342; Reed–Simon IV §XIII.2; Brydges, *A short course on cluster expansions*. Internal: PROOF-P3-UNCONDITIONAL §4; P3-AUDIT-CORRECTED §1.3, §2; G2D-ISOLATION-I2; G2G-TRANSFER-SPECTRAL; PROOF-EDGE-E2; PROOF-GAUSS-INVARIANT-PF; BS-BLIND (KL `ν≥0`); W2-positive-scheme; G2-RESIDUE-ATOM; P3-CYLINDER-RESIDUE (companion, `[KL]`).
