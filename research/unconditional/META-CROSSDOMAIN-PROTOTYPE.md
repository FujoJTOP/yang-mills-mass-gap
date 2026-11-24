# META-CROSSDOMAIN-PROTOTYPE — where in mathematics does a **non-positivity, interaction-driven, parameter-free lower bound on a spectral atom / bound state** already exist, and does any import to 4D YM?

> 2026-11-24 · META-METHOD. The charter (`CLAUDE.md §1`) asks for the missing species: a lower bound on the spectral atom `f₁²=ν({m₁²})` that is **(i) non-positivity** (not convexity/FKG/log-Sobolev — those are Schur-concave ⇒ atom-blind, `BRIEF §3`), **(ii) uses the interaction**, **(iii) has no small parameter**. This note does a **cross-domain prototype search**: find such an object in *other* mathematics, state its exact structure, run the free-field acid test, and give the honest import verdict. It does **not** re-tread the death map (positivity, Bakry–Émery, Wasserstein, replica-swap, tensor-RG, stochastic-quantization G1; `BRIEF §3`, `TRANSPORT-FIVE-WAY-VERDICT §1`).
>
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 / **[C]** 条件定理 / **[O]** 精确缺失步 / **[NG]** no-go for a *named class* / **[R]** external standard, as-reported (not re-verified this session) / **[J]** 判断. **No QED on anything unproven. This note does NOT prove Clay and does NOT prove no tool can exist (`∀class ≠ ∀possible`).**
>
> **House fact carried in from the corpus (load-bearing for §5):** the atom `f₁²` already has an **exact Fredholm-determinant identity** in the registry — `A3`: `f₁² = −λ₁ Res_{z=z₁}G(z) = −(1/z₁)·N_O(z₁)/D_O′(z₁)`, `D_O(z)=\det(1−zT_O)`, `z₁=e^{am₁}` — which is *structurally a Birman–Schwinger determinant* (`docs/new-tools-registry.md §A3`; `INVENT-NEW-IDENTITY.md`; `INVERSE-SCATTERING-IST.md §5`: "`det(1+K_GLM)=1/T(k)`, zeros=atoms ≡ A3"). So the cross-domain "best prototype" below is **not a new mechanism** — it is the corpus's own A3, recognized as Birman–Schwinger. That is the central, honest META result.

---

## VERDICT (up front)

> **In native mathematics there ARE production-ready, non-positivity, interaction-driven, parameter-free lower bounds on a spectral atom / bound state — three of the four prototypes are genuine. But NONE imports to 4D pure YM at the crossover `α_s~O(1)` as a *lower bound on `f₁²` / on `m₁`*. This is a clean "none imports" with a precise, per-prototype obstruction, and it is valuable.**
>
> - **Single best prototype = Birman–Schwinger / bound-state criterion (§1).** It has exactly the sought character: an **exact** operator-norm criterion `‖V^{1/2}(H₀+E)^{-1}V^{1/2}‖≥1 ⟺ an eigenvalue ≤ −E`, **non-positivity** (uses the attractive part `V_-` / a sign-indefinite kernel), **interaction-driven** (`V=0⇒` no atom), **parameter-free** (valid for all `V`, no weak-coupling), **acid-clean** (`V=0⇒K=0⇒f₁²=0`). It is *literally* A3/IST in another alphabet. **But its native topology is inverted for the glueball** (§1.4, §5): BS certifies a bound state **below a fixed continuum edge**; the glueball atom sits **above a positive vacuum, separated by confinement** (`BS-THRESHOLD §1–2`). Applied to the YM transfer matrix, the BS norm-`≥1` criterion is directionally locked onto the **Savvidy negative mode** (the Nielsen–Olesen tachyon, `λ_min=−g|B|<0`), not the `0⁺⁺` atom. **[NG] for the import; the residue form is A3, which stalls at the `a→0` zero-finding = marginal BPHZ wall.**
> - **Semiclassical / instanton double-well (§3)** is the only prototype that yields a **direction-correct *positive* lower bound on a separation** (`m₁ ≥ ω√S e^{−S}`, `[T-in-model]`, `E2-FLIP-TRANSPORT §2`), acid-clean non-vacuously. Its import is the **`e^{−S_inst}` transport machine**, which is `[C]/[T]` in the rigid-saddle/dilute class and **provably stops at the crossover** (size modulus `ρ*=O(1)/Λ`, packing `nρ⁴=O(1)`, `[V-here]`). **[O] = non-abelian liquid control = H1.** It bounds the *gap*, not the *residue* `f₁²`.
> - **Fractional-moment / Anderson localization (§2): [NG].** Native output is **upper** (Wegner DOS bound, modulus bound on `Var`) and localization, not a floor; the coherent Savvidy "disorder" is a Landau field that is **maximally delocalized** (`ANDERSON-LOCALIZATION §2`).
> - **Spectral-shift / Lieb–Thirring / Kato–Temple min–max (§4): [NG].** Lieb–Thirring bounds the *number/moment* of bound states from **above** (direction wrong); Krein's spectral shift is the **BS determinant** relabeled (⇒ §1's obstruction); Kato/Temple/Weinstein produce a lower bound only **given `E₂`** (separation oracle `ρ≤E₂`, `E2-FLIP-MINORANT Theorem ENCL`), and `E₂=ρ_r·m₁` by Single-Scale Covariance ⇒ circular.
>
> **The common obstruction (§5):** every prototype needs a **controller** — a *fixed continuum edge* (BS), a *quenched disorder ensemble* (FMM), a *rigid isolated saddle* (semiclassical), or a *separation oracle `E₂`* (Kato/Temple). In 4D YM at the crossover that controller is either **(a)** the gap itself, **(b)** the delocalized sign-indefinite Savvidy tachyon (wrong sign + wrong sector), or **(c)** an uncontrolled `O(1)`-coupling object. Each of (a),(b),(c) is `H1` under another name. **None crosses; all relocate the wall to the crossover.**

---

## §0 · The atom object and the four bars each prototype must clear

**Target (`BRIEF §1`).** `m₁=inf supp ν`; atom `f₁²=ν({m₁²})=|⟨0|¼trF²|G⟩|²>0`; Clay `H*` = atom at `m₁` **and** isolation. The KL measure `ν` is the spectral measure of `O_S=(β/2g)trF²` in the physical Hilbert space `H_phys=ker Gauss` (`INVENT-GAUGE-INVARIANT-MASS`). Direction bar: any **variational / Rayleigh / moment** instrument bounds `m₁` **from above** (`E2-FLIP-MINORANT §0`, `GAP-G2-edge`); a genuine tool must give `Φ≤m₁` or `f₁²≥r₀`.

**Bars (charter §1 + `BRIEF §2`).** A usable tool needs: **(B1)** non-perturbative; **(B2)** sign-controlling **lower** bound; **(B3)** non-circular (no presupposed gap/condensate); **(B4)** fails on the free *massless* field (`f₁²=0`, `G″=0`); **(B5)** `(L,a)`-uniform. The acid test **(B4)** is the sharp filter: a mechanism must be **true for the interacting theory and false (or empty) for the free massless field**.

I score each prototype on *native* structure first (does such an object exist in its home field?), then on **import** (does it transfer to 4D YM at the crossover?).

---

## §1 · Birman–Schwinger / bound-state criterion (the best prototype)

### 1.1 Exact structure [R: Birman 1961, Schwinger 1961, Simon; Reed–Simon IV §XIII]
For `H₀=−Δ` on `L²(ℝ^d)` and `H=H₀−V`, `V≥0` form-bounded with relative bound `0`, define the **Birman–Schwinger operator**
```
      K(E) := V^{1/2}(H₀+E)^{-1}V^{1/2},      E>0.
```
`K(E)` is compact, positive, and **decreasing in `E`** (`E↗⇒(H₀+E)^{-1}↓`). The **Birman–Schwinger principle**:
```
   ‖K(E)‖ ≥ 1   ⟺   spec(H) ∩ (−∞,−E] ≠ ∅        (an eigenvalue at −κ², κ²≥E);
   # eigenvalues of K(E) > 1   =   # bound states of H with energy < −E;
   H at −κ_j²  ⟺  K(κ_j²) has eigenvalue 1 ;   ground energy  −E₁ = −sup{E : ‖K(E)‖≥1}.
```
The **residue/normalization** of the bound state is `∝ 1/|∂_E(eigenvalue of K)|` at the crossing: with `K′(E)=−V^{1/2}(H₀+E)^{-2}V^{1/2}`, a simple eigenvalue `1` at `E=κ²`, eigenvector `g`, gives the spectral-measure atom weight
```
      f_{bound}²  =  |⟨·|ψ⟩|²  ∝  1 / |⟨g, K′(κ²) g⟩| .       (positive, but magnitude = kernel derivative)
```

### 1.2 Non-positivity mechanism [T]
The criterion is driven by the **negative (attractive) part `V_-`** and a **sign-indefinite** kernel; it is **not** a convexity/FKG/Γ₂ object — no measure is log-concave, `V` may change sign, and the conclusion is existence of a *negative* spectral point. It is **interaction-driven** (`V=0⇒K=0⇒` no atom) and **parameter-free** (`‖K(E)‖≥1` is a fixpoint, no small-`V` expansion; but see 1.4: it is a *threshold*, not automatic for every `V`).

### 1.3 Acid test in native form [T]
`V≡0` ⇒ `K(E)≡0` ⇒ `‖K(E)‖=0<1` for all `E` ⇒ **no bound state, `f₁²=0`**. The mechanism correctly returns *empty* on the free field. **Acid test PASSES non-vacuously** — this is the rare prototype that clears **(B4)** by construction, because the interaction literally *creates* the kernel.

### 1.4 Import to 4D YM — the inversion [NG for the import; A3 residue form = [O]]
Three facts, all already banked, turn the import into a no-go for the *literal* criterion:

1. **Topology is inverted.** BS certifies a state **below a fixed continuum edge** (energy `< −E`, essential spectrum `[0,∞)`). The glueball is an atom **above a positive vacuum** (`m₁²>0`), and in confinement there is **no** free two-gluon continuum below it — the "edge" `(2m_g)²` / multi-gluon threshold `θ_c≥4m₁²` *is itself the gap* (`BS-THRESHOLD §1–2`: "the two-gluon threshold is not a physical object; the only physical threshold is multi-glueball, `θ_c≥(2m₁)²`"). So the only *fixed* controller edge (the free continuum at `s=0`) lies **below** the atom — BS would look for a state *under* 0, i.e. an **instability**, not the `0⁺⁺` mass.
2. **It locks onto the Savvidy tachyon.** On the YM transfer matrix, `V` = the fluctuation operator `2adF_B` whose transverse Hessian has `λ_min=−g|B|<0` (`PROOF-SAVVIDY-EXACT-SIGN`). The BS norm-`≥1` criterion becomes a detector of the **Nielsen–Olesen negative mode**, whose dominant spectral weight is the **maximally delocalized Landau level** (`ANDERSON-LOCALIZATION §2`, `ξ=∞`), a **1-form / UV** object — `[T]` wrong sign and wrong sector (`D-BBS-LOWRANK-KATO §a`; `FORGE-KREIN-SPECTRAL B5`: "`J`-positive bottom `=Θ(g|B)`, UV/Landau; the gap is the `k→0` IR `min`"). **The non-positivity criterion is exactly right for detecting negative weight — and the glueball is positive mass above a positive vacuum.**
3. **The residue form is A3 and stalls at `a→0`.** Read `z=e^{-aE}`, `T_O=e^{-aH}`: BS's `det(1−K)` *is* `D_O(z)=det(1−zT_O)`, its **zero is the atom**, its **derivative is the residue** (`docs/new-tools-registry.md §A3`; `INVERSE-SCATTERING-IST §5`). To **lower**-bound `f₁²` one must (i) locate `z₁=e^{am₁}` and (ii) lower-bound `N_O(z₁)/D_O′(z₁)` — i.e. bracket the zero (= know `m₁`) **and** upper-bound the kernel derivative at the crossing. Both are controls on the exact BS kernel at `α_s~O(1)`, whose negative part is the sign-indefinite tachyon; making them `(L,a)`-uniform as `a→0` is the **marginal BPHZ analyticity / zero-finding wall** (`new-tools-registry §A3 残留`, `B6`). **This is precisely why the corpus's best A3 tool "sees the atom" yet keeps the wall on the continuum side.**

> **§1 verdict.** Best prototype by *character* — exact, non-positivity, interaction-driven, parameter-free, acid-clean, and *already the corpus's A3/IST*. But the import is **[NG]**: native topology inverted (bound-below-edge vs atom-above-vacuum), the norm criterion detects the delocalized Savvidy tachyon (wrong sign/sector), and the residue form is A3, which stalls at the `a→0` zero-finding `= [O]` = marginal BPHZ / H1. **It relocates the wall; it does not cross it.**

---

## §2 · Fractional-moment method / Anderson localization (Aizenman–Molchanov, Aizenman–Warzel)

### 2.1 Exact structure [R: Aizenman–Molchanov CMP 157 (1993) 245; Aizenman–Warzel, *Random Operators* (AMS GSM 168, 2015); Fröhlich–Spencer; Germinet–Klein]
For a random Schrödinger operator `H_ω=−Δ+V_ω`, at energies `E` in the localized regime there is `s∈(0,1)` with
```
      E_ω |G_ω(x,y;E)|^s  ≤  C_λ e^{−c|x−y|}        (fractional-moment bound, from a Wegner estimate),
```
giving **pure-point spectrum**, exponential eigenfunction decay, and **dynamical localization**; the **Wegner estimate** `E[N_ω(I)] ≤ C|I|` bounds the **local density of states from above**.

### 2.2 Mechanism & direction [T]
The mechanism is **quenched disorder averaging** plus the non-linearity of `s<1` (subharmonicity making resolvent moments summable). It is **not** interaction-driven in our sense (it is a random external potential), and its native outputs are **upper** (Wegner DOS, finite `Var`) and about the **nature** of spectral weight (point vs continuous), **not** a lower bound on the spectrum's bottom.

### 2.3 Acid test
`V_ω≡0` ⇒ free resolvent **delocalized**, FMM gives nothing. But this is a **letter pass only**: "not localized" is not "no atom," and FMM never distinguishes a *gapped* from a *gapless* theory — a **massive free field is extended yet gapped**, so "localized ⇒ gap" is a false inference. **Acid test does not discriminate.**

### 2.4 Import to YM — [NG]
Already adjudicated in `ANDERSON-LOCALIZATION`: the coherent Savvidy background is the **Landau problem** (maximally extended, `ξ=∞`, bottom the delocalized tachyon `−g|B|`); genuine localization would need the **fluctuating `F` as quenched disorder**, which is **structurally absent** (the disorder is the mode's own dynamical field = G1). Even granting localization, a modulus bound on `Var` cannot **sign** the spectral bottom (`G″=E[H]−Var+W`, sign-blind; `ANDERSON-LOCALIZATION §3`), and at the mobility edge `ξ(λ)→∞` the modulus is **vacuous**. **No floor, wrong direction, no quenched ensemble: [NG].**

---

## §3 · Semiclassical / WKB / instanton double-well spectroscopy (Coleman)

### 3.1 Exact structure [R: Coleman, *Aspects of Symmetry* Ch. 7; Coleman, PRD 15 (1977) 2929; Callan–Coleman]
Symmetric double well `V=(g/4)(x²−v²)²`: the tunneling **splitting** (the gap between the two lowest states) is
```
      m₁  =  E₁−E₀  ≈  (2ω/√π)·√(S_inst)·e^{−S_inst},   S_inst = (4v³/3)√(g/2),  ω=√(V″(well)),
```
with controlled relative error `O(1/S_inst)` in the **rigid-saddle / `S_inst≫1`** regime. The corpus's `[V-here]` witness (`w7/e2_flip_transport.py`, `E2-FLIP-TRANSPORT §2`) gives `m₁/Φ_tun = 1.48±0.02` essentially constant over `g∈[4,10]`, `corr(m₁,e^{−S})=+0.973`.

### 3.2 Mechanism & direction [T-in-model]
**Non-positivity:** the splitting comes from the **instanton/WKB** (the potential's barriers `S_inst`), not from convexity. **Interaction-driven:** `S_inst` from the potential. **Parameter-free** in the asymptotic sense (`S_inst≫1` is a scale, not a coupling expansion — but it *is* a regime). **Direction:** a genuine **lower bound on a separation** `m₁ ≥ ω√S e^{−S}` — the only prototype delivering a *positive* bound of the right sign.

### 3.3 Acid test [T-in-model; V-here]
Scale-invariant/free limit `S_inst→∞` ⇒ `e^{−S_inst}→0` ⇒ splitting `→0`; the free massless size-integral is scale-invariant (no running ⇒ `∫dρ/ρ⁵` diverges ⇒ no preferred size ⇒ `m=0`, `E2-FLIP-TRANSPORT §2b`). **Acid test PASSES non-vacuously.** This is the one genuinely direction-correct, acid-clean positive prototype.

### 3.4 Import to YM — [C]/[T] in class, [O] at the crossover
Import is the **`e^{−S_inst}` transport machine** (`E2-FLIP-TRANSPORT`), and its status is banked and precise:
- **[T-in-model]** rigid-saddle: `m₁ ≥ ω√(S_inst)e^{−S_inst}` in any theory with an isolated finite-action saddle (double well; `ℝ³×S¹` fixed-size bion; 3D compact `U(1)`).
- **[T]** exact multiplicative residue transport `f₁²(μ)=Z_O²f₁²(μ₀)` (`0→0`, acid-clean) — but this **moves** a residue, it does not **create** one; it is `Z_O²` times an existing sign.
- **[NG] at the crossover:** the 4D instanton has a **size modulus**; the one-loop distribution peaks at `ρ*=e^{−2N/(b−5)}/Λ` (SU(3): `2.72/Λ`, **large = IR**), where the packing `nρ⁴=O(1)` — a **liquid, not a dilute gas** (`[V-here]`, `TRANSPORT-PACKING-AUDIT` corrects the earlier `O(10)` to `0.01–0.06` **at `ρ̄`** but the peak-scale control is the issue). The cluster/KP "theorem" needs `S_inst≫1` **and** `nρ⁴≪1` **at the mass-generating scale**: control ⟺ (UV / small-`S¹` / Coulomb) = **wrong scale/phase/channel** (`E2-FLIP-TRANSPORT §4.2`). Two further locks: **channel** (the controlled condensates give a 1-form `√σ`, not the `0⁺⁺` edge) and **one-scale circularity** (`e^{−S}=e^{−8π²/g²}=(Λ/μ)^{power}` is transmutation; `m₁=c·Λ` with `c` unknown = the gap).
- **Crucially for the *atom*:** the semiclassical splitting bounds the **level separation `m₁`**, not the **residue `f₁²`**. The visibility `f₁²` is a *different* overlap (`~e^{−S}`), and transporting it to `2.5` is again the `a→0` marginal wall.

> **§3 verdict.** The strongest **positive-in-spirit** prototype: right-signed lower bound on a separation, acid-clean, `[T]`/`[C]` in the rigid-saddle/dilute class. But it **provably stops at the crossover** (size-modulus/diluteness ⟺ wrong-scale dichotomy `[V-here]`), and it does not by itself lower-bound the residue. **Its import is exactly the missing machine, which `= H1`.**

---

## §4 · Spectral-shift function / eigenvalue counting / Kato–Temple min–max / Lieb–Thirring

### 4.1 Lieb–Thirring — direction is UPPER [R: Lieb–Thirring PRL 35 (1975) 687; Stud. Math. Phys. (1976)]
```
      Σ_j |E_j|^γ  ≤  L_{γ,d} ∫ |V_-|^{γ+d/2}        (γ>0);  #bound states ≤ C_d ∫|V_-|^{d/2}.
```
This bounds the **number / moment of bound states from above** — it **caps** the count, it does **not floor** any single atom. *Direction wrong for our target* (we need `f₁²≥r₀`, a **lower** bound on one atom). The complementary "at least one bound state" statement is **exactly the BS criterion** (§1) — the same object, so §4 cannot add a lower bound that §1 lacked.

### 4.2 Krein spectral shift = BS determinant relabeled [R: Krein 1953; Birman–Krein; Lifshits]
`ξ(λ;H,H₀)` satisfies `∫f dξ = Tr[f(H)−f(H₀)]` and the **Birman–Krein formula** `ξ(λ)=−(1/π)arg det(1+K(λ+i0))`. Its jumps at bound states *are* the eigenvalues; its determinant is the **same `det(1+K)`** as §1 / A3. So the "spectral shift" route is **not a fourth mechanism** — it is Birman–Schwinger in another alphabet, inheriting §1's inverted-topology obstruction.

### 4.3 Kato / Temple / Weinstein min–max with a sign-changing test function — needs `E₂` [T: Kato 1949, Temple 1928, Weinstein 1932; corpus `E2-FLIP-MINORANT Theorem ENCL`]
```
      Temple:  m₁ ≥ R − Var/(ρ−R)   sound  ⟺  ρ ≤ E₂ ;   sharp and strictly increasing in ρ.
```
Min–max with *any* test function (sign-changing included) gives `R(ψ) ≥ E₁` = **upper**. The only *lower* bounds (Kato/Temple/Weinstein) **consume a certified `ρ ≤ E₂`** — the separation oracle. By **Single-Scale Covariance** `E₂=ρ_r·m₁` (`E2-REGRESS §0.4`, `G2D-ISOLATION-I2 §2`), the oracle is a fixed multiple of the very `m₁` being bounded ⇒ **circular** (`E2-FLIP-MINORANT Theorem ENCL`). With a sign-changing test function the sign controls nothing about `f₁²`; it only enters the (upper) Rayleigh quotient.

### 4.4 Acid test & import — [NG]
Lieb–Thirring with `V=0` is `0≤0`: satisfied by the free field trivially ⇒ **does not distinguish** interacting from free. Temple/Kato "pass" the free-field test only by **borrowing `m_g>0` / `E₂>0`**, which is the gap (`BS-VARIATIONAL-LOWER §6`: "manufactures `m₁>0` only by importing a mass gap through `m_g`"). **All four land on H1 or are direction-wrong: [NG].**

---

## §5 · The best candidate and the precise obstruction

### 5.1 Summary table

| prototype | native object | mechanism | direction | acid test | import to YM @ crossover |
|---|---|---|---|---|---|
| **BS / bound-state criterion** (§1) | `‖K(E)‖≥1 ⟺ e.v. ≤−E`; residue `∝1/|K′(κ²)|` | non-positivity (`V_-`, indefinite kernel) | **lower** on binding *below a fixed edge* | **PASS** (`V=0⇒K=0`) | **[NG]** — topology inverted; detects Savvidy tachyon; residue form = A3 `[O]` |
| **FMM / Anderson** (§2) | `E\|G\|^s≤Ce^{−c\|x−y\|}`; Wegner DOS | quenched-disorder averaging | upper (DOS, `Var`) / localization | no discrimination | **[NG]** — delocalized (Landau), no quenched ensemble (G1), sign-blind |
| **Semiclassical / WKB / instanton** (§3) | `m₁≈(2ω/√π)√S e^{−S}` | instanton / tunneling | **lower** on a *separation* | **PASS** (non-vacuous) | **[C]/[T] in rigid-saddle/dilute; [NG]/[O] at crossover** (size modulus) |
| **Lieb–Thirring** (§4.1) | `Σ\|E_j\|^γ≤L∫\|V_-\|^{γ+d/2}` | kinetic-vs-potential | **upper** on count/moment | trivially satisfied by free | **[NG]** direction wrong |
| **Spectral shift** (§4.2) | `ξ=−(1/π)arg det(1+K)` | = BS determinant | count, jumps=atoms | — | **[NG]** = §1 relabeled |
| **Kato/Temple min–max** (§4.3) | `m₁≥R−Var/(ρ−R)` iff `ρ≤E₂` | separation oracle | **lower given `E₂`** | borrows gap | **[NG]** `E₂=ρ_r m₁` (SSC) circular |

### 5.2 The precise obstruction (common to all four) [T]
> **Controller-necessity.** Each prototype certifies an atom against a **controller**: BS against a **fixed continuum edge** (the `0` of `spec(H₀)`); FMM against a **quenched random ensemble**; semiclassics against a **rigid isolated saddle** (fixed `S_inst`, no modulus); Kato/Temple against a **separation oracle `E₂`**. In 4D pure YM at `α_s~O(1)` the controller is, in every case, one of three things:
> **(a) the gap itself** — the multi-gluon threshold `θ_c≥(2m₁)²` (`BS-THRESHOLD §1`), or `E₂=ρ_r m₁` (`G2D-ISOLATION-I2 §2`) ⇒ circular;
> **(b) the delocalized sign-indefinite Savvidy tachyon** — `λ_min=−g|B|<0`, dominant weight on the Landau level `ξ=∞` (`ANDERSON-LOCALIZATION §2`), a **1-form/UV** object, so any BS/spectral-shift **norm/lower** criterion triggered by the negative part lands on the **instability**, not the `0⁺⁺` atom (`D-BBS-LOWRANK-KATO`, `FORGE-KREIN-SPECTRAL B5`);
> **(c) an uncontrolled `O(1)`-coupling object** — the semiclassical saddle's **size modulus** (`ρ*=O(1)/Λ`, `nρ⁴=O(1)`) `[V-here]`, or the exact BS kernel's derivative at the crossing.
> **(a) = (b) = (c) = `H1`** read in three coordinates. The prototypes are genuine; the **controller is the wall**.

### 5.3 Why "none imports" is the right verdict, and is valuable
- The prototypes are **not vacuous**: three of four are real theorems in their home fields, and BS + semiclassics pass the free-field acid test **non-vacuously**. The cross-domain search **succeeds in its literal task** — it locates *where* a non-positivity, interaction-driven, parameter-free atom-lower-bound lives (Schrödinger spectral theory: BS; semiclassical spectroscopy: WKB/instanton).
- But the import fails **structurally**, not for lack of effort: the sought object `f₁²≥r₀` is an atom **above a positive vacuum separated by confinement** — a *spectral-edge* phenomenon, whereas BS/localization are theories of states **below a fixed edge / in disorder**, and Kato/Temple are theories of **separation given a second level**. **YM's atom is none of these topologies**; it is the removal of the low continuum (confinement) creating an isolated weight — which no "bind-below-continuum" or "count-below-edge" criterion reaches.
- The one prototype whose **direction and acid test are exactly right** (semiclassical instanton splitting) is **already the corpus's `e^{−S_inst}` machine**, is `[T]`/`[C]` only in the rigid-saddle/dilute class, and stops at the crossover with a `[V-here]`-computable obstruction. **No new mechanism is hiding in this cross-domain search; the best prototype is the corpus's own A3 (BS) and its own instanton splitting, and both already sit on the banked wall.**

### 5.4 Exact missing step (unchanged, sharpened)
```
   [O]  a lower-bound instrument Φ with 0<c ≤ Φ ≤ m₁  (or f₁² ≥ r₀),
        non-positivity, interaction-driven, parameter-free, (L,a)-uniform,
        whose "controller" is NEITHER E₂ (circular) NOR the Savvidy tachyon
        (wrong sign/sector) NOR an uncontrolled O(1) liquid (H1).
   Equivalently: a stable, sign-definite, gauge-invariant quadratic minorant of
   the YM Dirichlet form at k→0  —  the "6th species" (E2-FLIP-MINORANT §5). None of
   {BS, FMM, WKB/instanton, Lieb–Thirring, spectral-shift, Kato–Temple} supplies it.
```
**Flip condition.** This verdict flips to "a prototype imports" iff someone exhibits, for pure `SU(N)` YM on `ℝ⁴` and uniform in `(L,a)`: **(BS)** either a fixed continuum edge *above* the `0⁺⁺` atom that is **not** a multiple of `m₁` (i.e. a kinematic, non-circular controller — `BS-VARIATIONAL-LOWER §2–6` shows the multigluon/binding threshold is dynamical `=I2`), **or** a Birman–Schwinger **norm-`≥1`** proof whose excess comes from the **positive-mass** `0⁺⁺` sector rather than the negative Savvidy mode; **(WKB)** or a non-abelian liquid closure `nρ⁴≤𝔡₀<1` at `ρ~1/Λ` (`E2-FLIP-TRANSPORT §5`). None is known; each is `[O]` `= H1`.

---

## §6 · References

**External — [R] standard / as-reported (not re-verified this session):**
- **Birman–Schwinger:** M. S. Birman, "On the spectrum of singular boundary-value problems," *Mat. Sb. (N.S.)* **55(97)** (1961) 125–174; J. Schwinger, "On the bound states of a given potential," *Proc. Nat. Acad. Sci. USA* **47** (1961) 122–129; B. Simon, *Quantum Mechanics for Hamiltonians Defined as Quadratic Forms* (Princeton, 1971); M. Reed, B. Simon, *Methods of Modern Mathematical Physics IV*, §XIII.3, §XIII.14 (Academic Press, 1978).
- **Anderson / fractional moments:** M. Aizenman, S. Molchanov, "Localization at large disorder and at extreme energies: an elementary derivation," *Comm. Math. Phys.* **157** (1993) 245–278; M. Aizenman, S. Warzel, *Random Operators: Disorder Effects on Quantum Spectra and Dynamics*, AMS GSM **168** (2015); J. Fröhlich, T. Spencer, *CMP* **88** (1983) 151; F. Germinet, A. Klein, *Ann. Henri Poincaré* **2** (2001) 489; P. Stollmann, *Caught by Disorder* (Birkhäuser, 2001).
- **Semiclassical / instanton:** S. Coleman, "Fate of the false vacuum: semiclassical theory," *Phys. Rev. D* **15** (1977) 2929; S. Coleman, *Aspects of Symmetry* (Cambridge, 1985), Ch. 7 (double-well splitting `~ω√S e^{−S}`); C. Callan, S. Coleman, *Phys. Rev. D* **16** (1977) 1762.
- **Lieb–Thirring:** E. H. Lieb, W. Thirring, "Bound for the kinetic energy of fermions which proves the stability of matter," *Phys. Rev. Lett.* **35** (1975) 687; "Inequalities for the moments of the eigenvalues of the Schrödinger Hamiltonian…," in *Studies in Mathematical Physics* (Princeton, 1976); M. Cwikel; G. V. Rosenbljum; B. Simon.
- **Kato / Temple / Weinstein / spectral shift:** T. Kato, "On the upper and lower bounds of eigenvalues," *J. Phys. Soc. Japan* **4** (1949) 334; T. Kato, *Perturbation Theory for Linear Operators* (Springer, 1966); G. Temple, *Proc. London Math. Soc.* (1928); A. Weinstein, *Phys. Rev.* **40** (1932) 797; M. G. Krein, "On the trace formula in perturbation theory," *Mat. Sb. (N.S.)* **33** (1953) 597; M. Sh. Birman, M. G. Krein; I. M. Lifshits.

**Internal (absolute paths; the banked results this META note classifies against):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop + acid test; §2 labels; §5 wall/key); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§1 atom `f₁²`; §2 five bars; §3 death map; §4 transport≡confinement).
- `D:\Dev\math\docs\new-tools-registry.md` (§A3 Fredholm determinant identity = the BS structure; §B6 completely-monotone no-go; §D residue dictionary); `D:\Dev\math\research\unconditional\TRANSPORT-FIVE-WAY-VERDICT.md`; `E2-FLIP-TRANSPORT.md` (semiclassical `[T-in-model]`, size-modulus `[V-here]`, `[NG]` for condensate class); `ANDERSON-LOCALIZATION.md` (FMM, Landau delocalization); `INVERSE-SCATTERING-IST.md` (GLM = A3, Temple needs `E₂`, Levinson integer frozen); `BS-VARIATIONAL-LOWER.md` / `BS-THRESHOLD.md` (BS kernel sign-indefinite; threshold dynamical); `D-BBS-LOWRANK-KATO.md` / `FORGE-KREIN-SPECTRAL.md` (Kato UV tower, count≠bottom, UV≠IR); `WKB-BORN-OPPENHEIMER.md` (WKB fast gap ≠ slow `G″`); `E2-FLIP-MINORANT.md` (Theorems ENCL/DET/CMP/POS; five lower-bound methods ⇒ `H1`); `G2-RESIDUE-ATOM.md`; `G2G-TRANSFER-SPECTRAL.md` (Bushell lower bound = `G″` = `H1`); `G2D-ISOLATION-I2.md` (SSC `E₂=ρ_r m₁`); `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_min=−g|B|`); `INVENT-NEW-IDENTITY.md` (A3); `INVENT-VARIANCE-INEQ.md`; `TRANSPORT-PACKING-AUDIT.md` (`𝔡=0.01–0.06` at `ρ̄`).

---

*Status: **META-CROSSDOMAIN VERDICT = "none imports," cleanly.** Three of four prototypes are genuine non-positivity, interaction-driven, parameter-free atom/bound-state lower bounds in native mathematics (Birman–Schwinger `[T]`, Aizenman–Molchanov `[T]`, Coleman/WKB instanton splitting `[T]`/`[R]`); **none imports** to 4D pure YM at the crossover as a lower bound on `f₁²` or `m₁`. **Best prototype = Birman–Schwinger**, which is the corpus's own **A3** Fredholm determinant; its native topology (bound state **below a fixed continuum edge**) is **inverted** for the glueball (**atom above a positive vacuum, separated by confinement**), so on the YM transfer matrix the norm-`≥1` criterion detects the **delocalized Savvidy tachyon** (wrong sign + wrong sector), and the residue form stalls at the `a→0` zero-finding = marginal BPHZ wall. The single direction-correct, acid-clean positive prototype (semiclassical instanton splitting) is **already the `e^{−S_inst}` machine**, `[T]`/`[C]` only in the rigid-saddle/dilute class, and `[NG]`/`[O]` at the crossover (size modulus). **Common obstruction = controller-necessity:** every prototype needs a controller (`E₂` / fixed edge / quenched ensemble / rigid saddle) that in YM is the gap itself, the Savvidy tachyon, or an uncontrolled `O(1)` liquid — all `= H1`. **NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay stays open; the missing species is a stable sign-definite gauge-invariant quadratic minorant at `k→0 = H1`.***
