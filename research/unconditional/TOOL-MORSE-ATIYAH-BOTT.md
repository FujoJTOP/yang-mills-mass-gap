# TOOL-MORSE-ATIYAH-BOTT — the Witten deformation / Morse complex of `S_W` on `A/G`: does its **low spectrum** (not its index) force a 0⁺⁺ floor?

> 2026-11-24 · Area-scout for the 4D pure SU(N) Yang–Mills mass gap (Clay). Scope: algebraic topology / Morse theory of the Yang–Mills functional `S_W=(1/2g²)∫Tr F²` on `A/G` — Atiyah–Bott equivariant cohomology, Morse homology, and specifically **Witten's deformation / Witten complex low-spectrum structure** (small eigenvalues governed by critical points + gradient flow). This document deliberately targets the sub-question the prior files left open: *the Witten **spectral** structure at the critical locus, as opposed to the index/invariant level* (`PROOF-MORSE-FLOER` = invariants; `H1-T7` = localization/Q-cohomology; `INVENT-NEW-INDEX` = the frozen index).
> **Labels (charter §2):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go *for a named class* / **[O]** 精确缺失步 / **[J]** 判断. External facts marked **[V]** (verified this session) / **[R]** (standard/as-reported). **No QED on anything unproven; Clay stays open; this is NOT a global impossibility theorem (`∀class ≠ ∀possible`).**

---

## Verdict up front

> ### **WRONG-OBJECT + INDEX-RIGIDITY. No genuine leverage. The one thing the spectral-flow no-go *missed* — the Witten low-eigenvalue *splitting* is continuous, not a frozen integer — escapes rigidity but lands on the instanton/condensate (θ / 0⁻⁺) channel, is an *upper* bound not a lower bound, and its prefactor is the *indefinite Savvidy Hessian itself*. The "gap above" is exactly the harmonic Hessian at the critical point, which at the YM vacuum is the *free* operator (`λ_min(k→0)=0`) for every `g`.**

**Why, in one paragraph.** The Witten deformation `Δ_{f,h}=e^{-f/h}\,d\,e^{2f/h}\,d^*\,e^{-f/h}…` of the de Rham complex on a finite-dimensional `M`, with `f=S_W`, has a low spectrum with **three strata** (Helffer–Sjöstrand **[R/V-method]**): *(I)* the bottom `b_k` eigenvalues are exponentially small — the **index/Morse** stratum, coupling-free and deformation-rigid; *(II)* the precise small eigenvalues are **Eyring–Kramers tunneling rates** `λ≈(h/π)|λ̂₁(saddle)|√(…)·e^{−2Δf/h}`, `Δf=f(saddle)−f(min)`; *(III)* the "gap above" scales like `h·spec(Hess f at the critical point)`. On `A/G` these become: *(I)* rigid → free value `0` (`INVENT-NEW-INDEX`); *(II)* = the **instanton/sphaleron barrier** `8π²/g²` (the θ/condensate/0⁻⁺ relaxation, orthogonal to the neutral `0⁺⁺` `tr F²`, and a bound **on the near-zero eigenvalues**, i.e. an *upper* statement — the wrong direction for a support floor); *(III)* = the Hessian at the vacuum `−Δ` (= free Maxwell, `λ_min(k→0)=0`, fails the acid test: identical on free and interacting). The 0⁺⁺ glueball is a **within-well vibrational** excitation, i.e. stratum (III) *plus* the non-harmonic measure correction; the Witten/Morse complex retains only the harmonic Hessian and the inter-well tunneling, i.e. it **discards exactly the measure/variance interaction `G″=E[H]−Var` that is H1.** Structurally there is not even a finite-dimensional Morse–Smale complex on 4D `A/G` to run the theorem on (Uhlenbeck bubbling; instanton number continuous on ℝ⁴). **This is the spectral-flow/index wall re-derived in the *spectral* (not just invariant) coordinate, plus the one new fact that the splitting is the wrong channel/upper-bound.**

---

## §0 What is genuinely new here vs. the existing corpus

The area is densely pre-worked. To be honest about the delta:

| prior file | what it established | what it did **not** isolate (the gap this file fills) |
|---|---|---|
| `PROOF-MORSE-FLOER` | `S_W` on `A/G` is **Morse–Bott**; invariants count critical points (topology), coupling-free ⇒ not the gap. | Only the **invariant/count** level. Did not analyze the **Witten-deformed Laplacian's small eigenvalues** (the Eyring–Kramers stratum) nor the "gap above." |
| `H1-T7-localization` | Equivariant/BV **localization** computes `Q`-cohomology / exact `Z`; the gap is not a `Q`-closed datum. | Localization = exact partition function to saddles. Did not treat the **spectral** small-eigenvalue theorem. |
| `INVENT-NEW-INDEX` | Any discrete, locally-constant invariant is **rigid** → free value `0`; the acid-test-passing index counts **gauge copies**. | Rigidity is about the **index**. It did not ask whether the **continuous small-eigenvalue splitting** escapes rigidity (it does) and whether that helps (it doesn't — wrong channel/upper bound). |
| `FORGE-SPECTRAL-FLOW` | The Hessian's negative index / APS flow is `0` (0⁺⁺ block PSD) or `=N_neg` (UV Savvidy): frozen to a sign, blind to the support edge. | Again the **index**; not the Witten low-**spectrum** strata. |
| `BPS-SUSY-PROTECTION`, `PROOF-WITTEN-SPINOR` | SUSY/BPS protects only the **charge** sector (⊥ 0⁺⁺); the Witten **spinor** identity consumes `ℛ_A≳0` = Savvidy. | The positive-mass **spinor** identity, and BPS **states** — not the Witten **Laplacian** small-eigenvalue asymptotics on `A/G`. |

**The delta:** a precise statement that the *Witten spectral trichotomy* (index / Eyring–Kramers splitting / harmonic gap-above) maps on `A/G` to *(rigid 0)*(instanton/θ splitting, upper-bound, Savvidy-prefactored)*(free/Savvidy Hessian)* — and that the glueball is the **one object none of the three strata is** (a measure-generated within-well vibrational mass). Plus two live audit targets surfaced this session (§7).

---

## §1 The object precisely: `S_W` on `A/G`, its critical set, and why Atiyah–Bott is a *surface* theorem

**(a) Critical set.** `S_W=(1/2g²)∫_M Tr F_A∧*F_A` on the space of connections `A` (modulo gauge `G`). Critical points of `S_W` on `A/G` solve the **Yang–Mills equation** `D_A^*F_A=0`; the finite-action critical set is `{vacuum F=0} ∪ {instantons F=±*F}` (`BPST`; `Atiyah–Hitchin–Singer` **[R]**). This is the same content as `PROOF-MORSE-FLOER §1`; I restate only the pieces I use.

**(b) Why "Morse" is false — it is Morse–Bott, non-compact, and (on ℝ⁴) the label is continuous.** On ℝ⁴ the instanton moduli are **non-compact** (dilatation = an exact `Hess S_W=0` direction, scale-flat), so critical points are **manifolds**, not isolated points; and the instanton number `Q=(1/8π²)∫Tr F∧F̃` is **not** a discrete label (it is continuous before one fixes a compactification/topology). This is `PROOF-MORSE-FLOER §1` **[R]**.

**(c) Atiyah–Bott is the *surface* theorem that does not lift.** Atiyah–Bott, *Phil. Trans. R. Soc. A* **308** (1983) 523 **[R]**, computes the **equivariant cohomology** `H^*_G(A)` and the perfect-Morse theory of `S_W` **over a Riemann surface** — where the action *is* perfect Morse and the moduli are compact. The 4D analogue fails: Uhlenbeck **bubbling** (`CMP` **83** (1982) 11, 31 **[R]**) destroys compactness, so there is **no perfect-Morse, no finite-dimensional Morse–Smale complex** on 4D `A/G`. The rigorous 4D relatives are **Floer homology** (3-manifolds; `CMP` **118** (1988) 215) and **Donaldson theory** (compact 4-manifolds; `JDG` **18** (1983) 279) — both **topological field theories**, whose invariants are `g`-independent and whose "Hamiltonian" carries **no propagating mass**. This structural fact alone means the literal Witten "low spectrum ⇔ critical points" theorem has **no rigorous 4D instance** (§5).

**Key structural point (reused everywhere below): critically, `S_W`'s Hessian at the vacuum `A=0` is `−Δ` — the free Maxwell operator — for every `g`.** Because `F=0` at the vacuum, the `2[*,F]` curvature term of `Hess S_W=−D²+2\,\mathrm{ad}F_B` (`PROOF-SAVVIDY-EXACT-SIGN`) **vanishes there**; the Hessian is `−∂²`, `λ_min(k→0)=0`. This is the acid-test seed (`PROOF-MORSE-FLOER §3(a)`).

---

## §2 The Witten deformation and its three spectral strata (the core of this file)

Witten (*J. Diff. Geom.* **17** (1982) 661 **[R]**) proved the Morse inequalities by deforming the de Rham Laplacian on a **finite-dimensional** compact `M`: set `d_{f,h}=e^{-f/h}\,h\,d\,e^{f/h}` and `Δ_{f,h}=(d_{f,h}+d_{f,h}^*)^2`, `h>0`, `f` a Morse function (so `Δ_{f,h}=h^2Δ+|\nabla f|^2+h\,(\text{first-order})`).
The rigorous low-spectrum theory (**Helffer–Sjöstrand**, *Comm. PDE* **10** (1985) 245 **[R/V-method]**; tunneling: **Helffer–Klein–Nier / Helffer–Nier** **[R]**) gives **three distinct strata**, and it is the *separation of these three* that the prior files did not draw:

### Stratum (I) — the index strata are exponentially small and **rigid**
There are exactly `b_k(f)` (Betti number) eigenvalues of `Δ_{f,h}` on `k`-forms with `λ=o(e^{-C/h})` as `h→0`, and their count is the Morse number; `Σ(−1)^k c_k=χ`. This is the topological stratum. **On `A/G` with `f=S_W` the critical equation `F=±*F` carries no coupling `g`.** Hence this stratum is `g`-independent, deformation-invariant, and by `INVENT-NEW-INDEX §3` (rigidity: a discrete datum locally constant on the connected `β`-line equals its free value) evaluates to the **free value `0`**. It cannot force a *coupling-dependent* mass floor. — **[NG] (rigidity), already banked.**

### Stratum (II) — the precise small eigenvalues are **Eyring–Kramers tunneling rates** (continuous — escapes rigidity)
The exponentially small eigenvalues are **not zero**; for the 0-form Witten Laplacian the Eyring–Kramers asymptotics is **[R]** (Helffer–Nier):
> `λ_k(h) ≈ (h/π)\,|λ̂₁(saddle)|·√( det Hess f(min) / det Hess f(saddle) )·(1+O(h))·exp(−2[f(saddle)−f(min)]/h)`,
where `saddle` runs over **index-1** points and the exponent is `−2·(Agmon distance)` = the barrier `Δf`.
`λ_k` is a **continuous, positive, splitting** quantity — unlike an integer index, it is **not** frozen by rigidity. **This is the single item the spectral-flow no-go genuinely missed**, and it deserves a fair hearing.

**On `A/G`:** the "wells" are the **topological vacua** (Chern sectors `Q∈ℤ`, or `θ`-sectors), and `saddle` = the sphaleron / instanton-type saddle; the barrier `Δf` is the **instanton action** `≈8π²/g²·|ΔQ|`. After the stochastic-quantization identification (§4) `λ_k` is a **relaxation rate**, and this slowest rate is the **inter-sector (topological/θ) tunneling** — the `0⁻⁺`/`θ`/gluon-condensate channel, **orthogonal** to the neutral `0⁺⁺` `tr F²` residue `f₁²` (`BPS-SUSY-PROTECTION` Q4: the only charge in `F` is `tr F F̃`, which lives in `0⁻⁺`; `tr F²` is neutral). So Stratum (II) produces a positive quantity in the **wrong channel**.

**Three further reasons Stratum (II) is not a lever** (each independently fatal):
1. **Direction.** `λ_k` *is* the small eigenvalue — the formula **upper-bounds** how small the near-zero eigenvalues are. A support **floor** (`m₁²≥λ`, or `f₁²≥r₀`) needs a **lower** bound on the gap *above* these modes, which Stratum (II) does not supply. (This is the same direction failure as `HISTORICAL-TOOL-SWEEP` Family-I.)
2. **The formula may not even apply.** Eyring–Kramers needs a **non-degenerate index-1 Morse–Smale saddle**. On `A/G` the critical points are **Morse–Bott manifolds**, and the relevant saddles at the crossover are **Savvidy** points with **extensive negative index** `N_neg=c_N(gB)²V/8π²` (`PROOF-SAVVIDY-EXACT-SIGN`, `W6-BE`), not index-1. So the clean small-eigenvalue formula is not instantiated.
3. **The wall is in the prefactor.** Even granting an effective formula, `det Hess f(saddle)` involves the **indefinite Savvidy Hessian** — the wall reappears inside the prefactor. "The Witten identity does not see around Savvidy — it sees exactly Savvidy" (`PROOF-WITTEN-SPINOR`).

— Stratum (II): **[NG] (wrong-object: instanton/θ/condensate; wrong direction: upper not lower; needs index-1 it lacks).**

### Stratum (III) — the "gap above" is the harmonic Hessian at the critical point = **free/Savvidy**
Above the exponentially small group, the next eigenvalues scale like `h·spec(Hess f(critical point))` (harmonic/WKB quasimodes; Helffer–Sjöstrand **[R]**). **At the YM vacuum the Hessian is `−Δ`, the free Maxwell operator** (§1c). Therefore:
> the "Witten spectral gap" read at the vacuum is the gap of the **free Laplacian at the vacuum** — `λ_min(k→0)=0` — **for the free theory and for interacting YM alike.** Acid test: **returns the free gapless spectrum on both sides** (`PROOF-MORSE-FLOER §3(a)`); it does **not** fabricate a floor (good) but supplies **none** for YM (useless) — the charter's "vacuous" outcome; **fails the sign-controlling-lower-bound requirement (criterion 2).** At Savvidy backgrounds it is *indefinite* (no gap at all).
— Stratum (III): **[NG] (vacuous on the acid test; free Hessian).**

### The decisive synthesis: the glueball is **none** of the three strata
The 0⁺⁺ glueball is a **within-well vibrational** excitation of the measure — Stratum (III) is its *harmonic skeleton* (`−Δ`, gap 0 at the vacuum), and the **actual mass is generated by the non-harmonic measure correction** `G″=E[H]−Var=S−D` (`INVENT-VARIANCE`, `G″=S−D`, H1). The Witten/Morse complex retains **only** the harmonic Hessian (Stratum III) and the **inter-well tunneling** (Stratum II); it **discards precisely the measure/variance interaction** that *is* the gap. Statement worth quoting:
> **The Morse/Witten complex computes the two things that are provably *not* the 0⁺⁺ glueball — the free harmonic Hessian at the vacuum (Stratum III) and the inter-well instanton tunneling in the θ/0⁻⁺ channel (Stratum II). The glueball's mass arises in the "beyond-Morse" remainder `G″=E[H]−Var`, which is H1 in the critical-point coordinate.**

---

## §3 Equivariant cohomology (Atiyah–Bott) — the direct answer to Q2

**Q2: is there a cohomological/equivariant invariant that FORCES a spectral bottom, or does it collapse to index-rigidity?**

`H^*_G(A)` (Atiyah–Bott; Berline–Getzler–Vergne; Mathai–Quillen) is computed from the **self-duality equation `F=±*F`, which contains no coupling `g`**. It is therefore **locally constant on the connected coupling family** — the exact hypothesis of `INVENT-NEW-INDEX`'s rigidity lemma — hence equals its **free value**. Its outputs are **counts** (Morse numbers `≥` Betti numbers, `Σ(−1)^k c_k=χ`) and **exact partition functions** (localization: `Z` localizes to the instanton fixed locus). A **count** cannot lower-bound a **continuous support edge** `m₁²` — the same "counts vs. the `k→0` min" mismatch as `FORGE-KREIN-SPECTRAL`. **Answer: it collapses to index-rigidity; it does NOT force a spectral bottom.** The only equivariant/Atiyah–Singer object that *does* survive the free-vs-interacting test is the **gauge-dressed Dirac index** `n(A)=Ind D₄[A]`, which (i) is **fermionic** (pure YM has no spinor, `⟨n⟩=0` by CP), and (ii) yields only the θ/`0⁻⁺` susceptibility `χ_t` as an **upper** bound on `m_{0⁻⁺}` — **wrong channel, wrong direction** (`INVENT-NEW-INDEX §1`; `GAP-IN11`). — **[NG].**

---

## §4 Is the identification legitimate? (stochastic quantization) — yes, and it makes the verdict sharper

A fair objection: "the Witten Laplacian on forms on `A/G` is not the physical transfer operator." **Reality is subtler and *favours* the question:** via **stochastic quantization** (Parisi–Wu; rigorous: Jona-Lasinio–Mitter, Zwanziger), the **Fokker–Planck / Witten–Laplacian on the *infinite-dimensional field space*** with `f=S_W` **is** the physical object, and its low eigenvalues **are** the Euclidean **KL masses** (the relaxation rates = the imaginary-time poles). So Q1 is legitimate: *could* the Witten/Morse complex on field space read off the masses?

**The trichotomy then answers it negatively and precisely:** the physical masses split into (i) the **inter-well tunneling** `λ_k` (§2-II) — which is the topological/instanton/condensate relaxation, *not* the 0⁺⁺; and (ii) the **within-well** excitations (§2-III) — whose harmonic part is the **free** spectrum (`−Δ`, gap 0) and whose *mass-generating* non-harmonic part is the **measure variance** `E[H]−Var`, the "beyond-Morse" remainder. The Witten complex captures (i) fully and (ii) only harmonically; the physical 0⁺⁺ mass is the part of (ii) the complex **throws away**. So even with the legitimate identification, the Morse/Witten complex does not pinch the glueball. — **[J]/[NG].**

---

## §5 Structural no-instance: there is no 4D Morse–Smale complex to run the theorem on

Beyond the object-level verdicts, the *machinery* does not exist in 4D:
- Witten's spectral theorem (Helffer–Sjöstrand) is a **finite-dimensional** statement. On 4D `A/G`, Uhlenbeck **bubbling** (`CMP` **83** (1982) **[R]**) and the **non-compact, scale-flat Morse–Bott** critical manifolds (`PROOF-MORSE-FLOER §1`) mean there is **no Morse–Smale complex** and hence **no instance** of the small-eigenvalue theorem. 
- The rigorous 4D "instanton-as-gradient-flow" theories (**Floer** 3-manifolds; **Donaldson** 4-manifolds; Kapustin–Witten / Haydys–Witten 5d–4d Morse–Smale–Witten complex **[R]**) are **topological field theories**: their data are `g`- and metric-independent, their "Hamiltonian" is `{Q,·}`-exact (`T_{μν}={Q,G_{μν}}`), so **"mass gap" is not even defined** in them. Relating a "Floer spectrum" to the physical Hamiltonian requires a **topological twist / SUSY** that pure non-SUSY SU(N) does not have (`PROOF-MORSE-FLOER §2`; `PROOF-NONINVERTIBLE-SYM`).
- So the precise **missing object** if one insisted on this route ≈ *a finite-dimensional, coupling-dependent Morse–Smale complex whose small-eigenvalue splitting equals the 0⁺⁺ glueball mass.* Building it requires the **twist/SUSY** (§6) and the **measure interaction** (§2) — the same wall.

---

## §6 Q3 — does it need SUSY, and does non-SUSY pure YM lose it?

**Yes to the first; the second loses it into the wrong channel.** The *spectral* interpretation of the Morse complex (the Witten index, `{Q,Q†}=H`, zero modes = critical points) is **supersymmetric quantum mechanics** (`Witten 1982` **[R]**). Pure YM has only the nilpotent **BRST** `Q` (cohomology = observables, no energy positivity) — no second scalar supercharge, no fermion. The SUSY completions that *do* supply it all fail to reach the 4D pure-YM IR:
- the **BPS/protected** atoms live in the **charge** sector (`Z=Witten–Olive area integral`), orthogonal to the neutral `tr F²`; `m=|Z|` degenerates to `m=0` on `tr F²|0⟩` (`BPS-SUSY-PROTECTION`, full argument);
- **N=2** on the Coulomb branch is itself **gapless**; its lightest BPS states are monopoles/dyons (charged), killed by charge conservation on `0⁺⁺`; **N=4** is conformal (gapless — the *opposite* of a gap); **N=1** confinement is an EFT/gaugino-condensation argument, **not** localization/Morse;
- **soft breaking** N=2→0 destroys the superalgebra (`Z` ceases to exist, multiplets lengthen) — protection gone; and the monopole condensation crossing is the dual-Meissner/√σ (cross-sector) wall.
— **[NG] (BPS/SUSY = wrong channel; rigidity at the index level).** The non-SUSY theory **loses the spectral protection entirely**, and the SUSY version protects the **wrong channel** — the "SUSY is the only protected positive-mass atom, wrong channel" no-go, confirmed from the Morse/Witten side.

---

## §7 Audit targets surfaced this session (the brief's "(A) audit" — highest-value mode)

The searches returned several 2024–2026 attempts in exactly this area. None are peer-reviewed; I flag the **structural** reason each cannot cross the wall (no claim about the authors' rigor beyond the specific logical point):

1. **Feehan, *Morse theory for the Yang–Mills energy function near flat connections*, arXiv:1906.03954** **[R]** — a serious monograph proving a **"Yang–Mills energy gap": small L²-curvature `⇒` flat connection** (via Łojasiewicz + infinite-dim Morse theory). **Do not confuse this with the Clay mass gap.** Feehan's "energy gap" is a statement about the **values of the functional `S_W`** (`S=0` vs `S≥ε` unless flat) — the **Morse–Bott/topological, `g`-free** structure (`PROOF-MORSE-FLOER`). It is **not** a spectral gap of the quantum theory; it is gauge-topological and holds identically for the free theory. Citing it as the mass gap is a category error (Stratum III ≡ free Hessian; §2).
2. **"Agawa proof" + 2025 Cambridge working-paper addendum** (`coe.cambridge.org`) **[R, unrefereed]** — claims "completion" via **Morse theory on the gauge-fixing functional** (Gribov uniqueness) + RG continuum limit. The gauge-fixing-functional Morse theory computes **Gribov copies / orbit-level** data — precisely `INVENT-NEW-INDEX §2`'s "the only acid-test-passing index counts **gauge copies, not physical states**." A count of Gribov horizons is gauge-dependent and not a cohomological invariant of the physical Hilbert space; it cannot supply the neutral `0⁺⁺` residue. **Wrong object.**
3. **Unrefereed preprints applying Witten–Helffer–Sjöstrand small-eigenvalue theory directly to a lattice YM "mass gap"** (e.g. a Zenodo "Master4.pdf", `records/18800920`) **[R, unrefereed]**. This is *exactly* the route analyzed in §2: the small eigenvalues are the **Eyring–Kramers tunneling rates** (instanton/θ/condensate, wrong channel, **upper** bound), and the "gap above" is the **harmonic Hessian at the vacuum = free** (`λ_min=0`). The standard HS theorem, applied faithfully, yields the **trichotomy**, not a 0⁺⁺ lower bound. **Treat with caution; the structural obstruction is §2/§5.** I did not execute or read that file as trusted code/data.

Honest calibration **[J]:** the existence of these attempts is itself evidence the area is attractive and the wall is *located* (people land on Stratum II/III), not that the route works.

---

## §8 Scorecard (charter's 5 hard conditions for a qualifying tool)

| # | Condition | Result | Why |
|---|---|---|---|
| 1 | **Non-perturbative** | **partial PASS** | an integer index is all-`β`; the Eyring–Kramers rate is non-perturbative. |
| 2 | **sign-controlling LOWER bound** (`m₁²≥λ` or `f₁²≥r₀`) | **FAIL** | Stratum II is an *upper* bound on the near-zero eigenvalues; Stratum III gives `λ_min(k→0)=0` (no floor). Rigidity kills Stratum I. |
| 3 | **non-circular** | **FAIL** | "critical-point-splitting ⟺ glueball mass" would need the measure interaction `E[H]−Var` = H1; the only escape (a `β*` jump) is built from the gap (`INVENT-NEW-INDEX §3`). |
| 4 | **non-SUSY + fails free field** | **FAIL (vacuous)** | at the vacuum the Hessian is the **free** `−Δ` for all `g`: identical free/interacting ⇒ vacuous; the YM-specific datum (instanton) is the wrong channel & an upper bound. |
| 5 | **uniform `(L,a)` floor, `a→0`** | **FAIL** | no floor; Morse–Bott non-compact, bubbling, no uniform small-eigenvalue control. |
| **score** | | **1/5** | = the known wall in the spectral coordinate. |

---

## §9 Exact stuck step (precise)

The missing implication, in the sharpest form:
> **"A Witten-complex small-eigenvalue datum of `S_W` on `A/G` ⟺ the 0⁺⁺ glueball mass `m₁`."**
> It never fires, for **four independent reasons**, each sufficient:
> 1. **(index stratum)** the topological datum is `g`-free ⇒ rigid ⇒ free value `0` (`INVENT-NEW-INDEX`);
> 2. **(splitting stratum)** the small eigenvalues are **Eyring–Kramers instanton/θ tunneling** `≈e^{−2·8π²/g²}` — wrong channel (⊥ neutral `tr F²`), **upper** bound not lower, and the formula needs index-1 Morse–Smale saddles that YM lacks (Savvidy is extensive-index);
> 3. **(gap-above stratum)** it is the **harmonic Hessian at the critical point = `−Δ` at the vacuum = free**, `λ_min(k→0)=0`, identical on free and interacting ⇒ vacuous;
> 4. **(structural)** no finite-dimensional **Morse–Smale complex** exists on 4D `A/G` (Uhlenbeck bubbling; non-compact scale-flat Morse–Bott); Floer/Donaldson are `g`-free TQFTs with no propagating mass.
>
> **The 0⁺⁺ mass lives in the "beyond-Morse" measure-variance `G″=E[H]−Var` (=H1), which the Morse/Witten complex discards by construction.** Pinning it means bounding the measure's quadratic form — which the Morse structure (critical points + tunneling + harmonic Hessian) cannot supply.

This is the **spectral-flow / index-rigidity wall** (`FORGE-SPECTRAL-FLOW`, `INVENT-NEW-INDEX`) re-derived at the **level of the low spectrum rather than the invariant**, plus the one genuinely-new observation: **Stratum II escapes rigidity but is the wrong channel and the wrong direction, with Savvidy in its prefactor.** So the spectral-flow no-go missed nothing that helps — the only thing it was silent on lands on the other two named walls.

---

## §10 Honest boundary & flip condition

- **This is NOT a proof of Clay and NOT a global no-go.** It is a **[NG] for the named class** (Morse/Witten low-spectrum + Atiyah–Bott equivariant cohomology + index/spectral-flow), plus an **[O]** exact missing step, plus **[J]**/"wrong-object" classification. `∀class ≠ ∀possible` (`KEY-EXISTS §5`).
- **Not under-selling Stratum II.** The Eyring–Kramers splitting *is* a continuous, non-rigid, non-perturbative, acid-test-clean datum — genuinely more than an index. But it is a *tunneling rate in the θ/0⁻⁺ channel*, bounded *above* in the wrong coordinate; that is a precise negative, not a dismissal.
- **Flip condition.** This flips to **genuine leverage** iff someone exhibits a coupling-dependent Morse/Witten datum that is (i) a **lower** bound on the **0⁺⁺** support edge (not an upper bound on near-zero eigenvalues, not the θ/condensate channel), (ii) **not** the free harmonic Hessian at the vacuum, and (iii) uniform in `(L,a)`. By §2 that datum would have to encode the measure variance `E[H]−Var` — i.e. it would *be* H1, not a Morse invariant. Equivalently: a finite-dim Morse–Smale complex whose small-eigenvalue splitting **equals** the glueball mass — which would need the twist/SUSY (§6) in the **neutral** channel (none exists) or the bubble-free 4D compactness (none exists, §5).
- **What it does buy (honest).** A clean **coordinate-system dictionary**: "index-rigidity" (Stratum I), "instantons≠glueballs / condensate wall" (Stratum II), and "Lichnerowicz = Savvidy / free Hessian" (Stratum III) are the **same wall** at three heights of the Witten spectrum. It localizes the missing step as *the measure-variance remainder the Morse complex cannot see*, and it certifies that the two live "proof" families in this area (energy-gap / gauge-fixing-Morse) target **Stratum III / gauge copies**, not the gap.

---

## §11 References

**External — verified/standard this session:**
- E. Witten, "Supersymmetry and Morse theory," *J. Differential Geom.* **17** (1982) 661–692. **[R]**
- B. Helffer, J. Sjöstrand, "Puits multiples en mécanique semi-classique IV. Étude du complexe de Witten," *Comm. Partial Differential Equations* **10** (1985) 245–340; and the small-eigenvalue/tunneling asymptotics (Helffer–Klein–Nier; Helffer–Nier; Eyring–Kramers `λ≈(h/π)|λ̂₁|√(det Hess(min)/det Hess(saddle))·e^{−2Δf/h}`). **[R/V-method]** (statement confirmed this session).
- M. F. Atiyah, R. Bott, "The Yang–Mills equations over Riemann surfaces," *Phil. Trans. R. Soc. Lond. A* **308** (1983) 523–615. **[R]**
- M. F. Atiyah, N. J. Hitchin, I. M. Singer, "Self-duality in four-dimensional Riemannian geometry," *Proc. R. Soc. Lond. A* **362** (1978) 425–461. **[R]**
- K. Uhlenbeck, "Removable singularities in Yang–Mills fields," *CMP* **83** (1982) 31; "Connections with Lᵖ bounds on curvature," *CMP* **83** (1982) 11. **[R]**
- A. Floer, "An instanton-invariant for 3-manifolds," *CMP* **118** (1988) 215; S. K. Donaldson, "An application of gauge theory to four-dimensional topology," *JDG* **18** (1983) 279. **[R]**
- N. Berline, E. Getzler, M. Vergne, *Heat Kernels and Dirac Operators*, Springer (1992); V. Mathai, D. Quillen, "Superconnections, Thom classes, and equivariant differential forms," *Topology* **25** (1986) 85. **[R]**
- G. 't Hooft (monopole/instanton); A. Belavin, A. Polyakov, A. Schwarz, Y. Tyupkin, *Phys. Lett. B* **59** (1975) 85; G. Savvidy, *Phys. Lett. B* **71** (1977) 133. **[R]**
- P. Feehan, "Morse theory for the Yang–Mills energy function near flat connections," arXiv:1906.03954. **[R]** (energy gap ≠ mass gap, §7).
- Kapustin–Witten / Haydys–Witten 5d–4d Morse–Smale–Witten (instanton) framework; SIGMA 2026 survey. **[R]** (g-free TQFT, §5).

**External — flagged unrefereed (audit targets §7, do not treat as established):**
- "On the Completion of Agawa's Proof…" (Cambridge *coe* working paper, 2025) — Morse on gauge-fixing functional / Gribov uniqueness. **[R, unrefereed]**
- Zenodo `records/18800920` ("Master4.pdf") and similar 2025–26 preprints applying Witten–Helffer–Sjöstrand / vorticity bounds to a YM "mass gap." **[R, unrefereed]**

**Internal (absolute paths):**
- `D:\Dev\math\CLAUDE.md`; `D:\Dev\math\docs\BRIEF-FOR-MODEL.md`.
- `D:\Dev\math\research\unconditional\PROOF-MORSE-FLOER.md` — Morse–Bott / topology-not-gap (the invariant level).
- `D:\Dev\math\research\unconditional\H1-T7-localization.md` — equivariant/BV localization = Q-cohomology.
- `D:\Dev\math\research\unconditional\INVENT-NEW-INDEX.md` — index/spectral-flow **rigidity** no-go; APS counts wrong channel.
- `D:\Dev\math\research\unconditional\FORGE-SPECTRAL-FLOW.md` — Hessian negative index frozen (0⁺⁺ PSD / UV Savvidy).
- `D:\Dev\math\research\unconditional\BPS-SUSY-PROTECTION.md` — SUSY/BPS = charge sector (⊥ 0⁺⁺).
- `D:\Dev\math\research\unconditional\PROOF-WITTEN-SPINOR.md` — Witten spinor consumes positivity; Lichnerowicz **is** Savvidy.
- `D:\Dev\math\research\unconditional\H1-T5-NCG.md`, `PROOF-THIMBLE.md`, `PROOF-STOCHASTIC-LOCALIZATION.md` — re-encodings / wrong-resource / flow=H1.
- `D:\Dev\math\research\unconditional\PROOF-SAVVIDY-EXACT-SIGN.md` (`H=−D²+2adF`, `N_neg`), `W6-BE-at-2.5` (`−1.12β`, `90–97%` orthogonality), `INVENT-VARIANCE` (`G″=S−D`), `PROOF-COARSE-POINCARE` (H1=`m_eff²`), `GAP-IN11`, `PROOF-NONINVERTIBLE-SYM`, `E2-REGRESS-CIRCULARITY.md`, `HISTORICAL-TOOL-SWEEP.md`.

---

*Status: **[NG] for the named class** (Morse/Witten low spectrum + Atiyah–Bott equivariant cohomology + index/spectral-flow; wrong-object for the low-spectrum splitting = instanton/θ/condensate, upper-bound; vacuous for the harmonic gap-above = free Hessian at the vacuum) + **[O] exact missing step** ("Witten small-eigenvalue datum ⟺ 0⁺⁺ glueball mass," blocked by the measure-variance remainder `G″=E[H]−Var` that the Morse complex discards) + **[J]** classification. **1/5** on the hard conditions. NOT a proof of Clay; NOT a global impossibility theorem. The spectral-flow no-go missed nothing that helps: the only datum it was silent on (the continuous tunneling splitting) lands on the instanton/condensate and free-Hessian walls. **Clay stays open.***
