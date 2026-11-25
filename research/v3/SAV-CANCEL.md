# SAV-CANCEL — the Savvidy negative-mode cancellation algebra, tested in the **physical 0⁺⁺ channel**: the rank-one/intertwiner algebra is exact and `D_Sav` is carried by a single Hessian matrix element; but for the **gauge-invariant 0⁺⁺ insertion** (tr F²) the negative mode is selected out (exactly at the symmetric saddle, suppressed generically) and enters only as a bounded **pair** — a strict bound that is a **decoupling, not a floor**, so the negative-mode program terminates in a no-go for the physical mass

> 2026-11-24 · Attack **SAV-CANCEL**. Task: test the Savvidy negative-mode cancellation algebra — couple
> the negative mode to positive modes, the **intertwiner signs**, the **physical 0⁺⁺ insertions**, and the
> **RG scale** — and look for a nonperturbative **Ward-type relation** making the negative-mode contribution
> to the physical `k=0` spectral state strictly bounded.
> **Reads:** `Q3-SCHUR-DEFECT.md` (inertia lemma, one-plaquette dichotomy, `G″=S−D`, `E[Σ]`),
> `H6-HEAD-SPECTRAL.md` (typical point vs essential infimum, `λ₁≥ess inf G″`), `SM2-RANK-ONE.md` (overlap
> lemma, margin identity `G″=m_soft+|D₊|`, ceiling `m_soft≤G″`), `T-WAVE-SYNTHESIS.md` (the four predicates
> `A,W,I,G`; the gap route `W∧H1`). **Also:** `SM1-SOFT-EXACT.md`, `SM4-CONSTANT-MODE.md`,
> `SM5-TCI-HEAD.md` (`TCI(strong)` false, `TCI(avg)=I(ρ)=H1`), `H2-DEFECT-BOUND.md` (`D=D₊+D_Sav`),
> `G4-GENERAL-BOUND.md` (`G″=A−Bβ²`), `CONSTRUCT-SAVVIDY-CELL.md` (scale separation + bounded softening),
> `SAVVIDY-WI-SYNTHESIS.md` (§2 WI routes all = A3 = H1), `research/unconditional/K0-WARD-LOWER.md` (the
> plain `k=0` Ward/Stein tower is acid-vacuous and `=` local `(◆)=H1`), `research/v3/WG-INVARIANT-ESTIMATOR.md`
> and `WG-FSS-FIT.md` (the **gauge-invariant** `χ_O` grows `L^{0.37}`, no clean plateau).
> **Witness:** `research/v3/w7/sav_cancel.py` (`py -3.12 -I`, numpy only, deterministic; output quoted §7).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (no QED; Clay 0%)

The algebra closes — and closing it is the finding. Five results decide the task:

1. **The cancellation algebra is one matrix element. [T]** With the rank-one Savvidy Hessian
   `K = K₊ − |λ₋|P₋` (`P₋=φ₋φ₋ᵀ`, `K₊⪰0`, `φ₋∈ker K₊`), the exact split
   ```
   ⟨K⁻¹c,c⟩ = ⟨K₊⁺c,c⟩ − |⟨φ₋,c⟩|²/|λ₋| ,   D = D₊ + D_Sav ,   D₊≤0 (BL),  D_Sav = E|⟨φ₋,c⟩|²/|λ₋| ≥ 0
   ```
   has the **intertwiner sign fixed**: the positive modes *stiffen* (`+`), the negative mode *softens*
   (`−`, the minus in front of the rank-one projector). On the exact one-plaquette the **entire**
   coarse↔negative coupling is the single Hessian element `c_r = g′(r)/r − g(r)/r² = sin(θ/2)·(stiffness)`:
   `c_r≡0` at `θ=0` (axial symmetry ⇒ `D_Sav=0`, `G″=S>0`) and maximal at `θ=π` (collective ⇒ `G″<0` past
   `β*=4.30661`). Verified exact to all digits (§7·1–2).

2. **The physical 0⁺⁺ insertion is a different carrier, and the negative mode is selected out of it.
   [T]+[V]** `D_Sav` is defined for a **linear** coarse mode `v` (through `c(v)=Π₀ Hess V v`). The physical
   `0⁺⁺` operator is the gauge-invariant **scalar composite** `O=tr F²` (the plaquette energy). On the
   Savvidy cell (4D `L=1` periodic `SU(2)`, constant field `U₁=e^{icσ₃/2}`, `U₂=e^{icσ₁/2}`) the negative
   mode is **exactly orthogonal** to `∇O` at the symmetric saddle (`φ₋·∇O = 0` to machine zero, both
   negative modes, two field strengths), while it has O(1) overlap (`0.354`) with the gauge-**dependent**
   Cartesian/uniform coordinate; `∇O` lies exactly on the *stiffest positive* eigenvector; and the negative
   mode reaches `O` only at **second order** (the pair term `φ₋ᵀHess(S/β)φ₋ = −0.15…−0.32`, strictly
   bounded). At generic configs the linear overlap is *suppressed* (`cos²`: `0.014` vs `0.085`, `≈6×`)
   rather than exactly zero. This is the colour + point-group selection: a single octet vector (TM) mode
   cannot interpolate a gauge-invariant scalar (`J=0`) at tree level. (§7·3–4.)

3. **The strict bound the task sought exists — and is a decoupling, not a floor. [NG]** Colour-singletness
   + angular momentum give `⟨0⁺⁺|O|1-gluon⟩=0` exactly in the continuum; the negative mode therefore
   contributes to the physical `0⁺⁺` spectral density only at the `φ₋φ₋` **pair threshold**
   `s ≥ 4m_stab² = 8g|B|` (`m_stab²=V″(φ₀)=2g|B|>0`), with a bounded pair coupling. So the negative-mode
   contribution to the physical `k=0` state **is strictly bounded** (zero linear, bounded quadratic). But
   this bound says **the negative mode cannot BE the atom**; it says nothing about whether the atom from the
   *positive* sector has `m₁>0`. Bounding the negative mode cannot lower-bound the physical mass.

4. **The RG-scale / trace-anomaly Ward identity confirms the split, not the gap. [T]-identity, [NG]-tool.**
   `O=tr F²` *is* the trace-anomaly density, `T^μ_μ=(β/2g)tr F²`. The selection of (2) says the negative
   mode does not feed the anomaly density at linear order; the physical `0⁺⁺` normalization runs by
   dimensional transmutation (an anomalous dimension), a fluctuation at scale `√(gB)` being decoupled from
   that running. This is *why* the selection exists (the anomaly is a scale-derivative of the positive-
   sector condensate), and it carries **no positivity** for `m₁`.

5. **Acid test passes, non-vacuously and without requiring the negative mode to vanish. [T]+[V]** Free
   massless: `F=0` ⇒ no Savvidy background ⇒ **no negative mode**, and no atom (`f₁²=0`); the selection
   returns "zero negative-mode contribution" and manufactures **no** gap. The negative mode is **not**
   required to vanish: it exists (`λ_min<0`) and is simply decoupled from the scalar (§7·3–5).

**One line.** *The Savvidy cancellation algebra is exact and is one Hessian matrix element with a fixed
intertwiner sign; in the physical `0⁺⁺` channel the negative mode is selected out (exact at the symmetric
saddle, `≈6×`-suppressed generically) and enters only as a bounded pair. That strict bound is a
**decoupling**: it proves the negative mode was never the physical atom's obstruction and relocates the
open step to the positive-sector atom — it is not the split `D=D₊+D_Sav`, not `H1`, not the weight `W`, and
**not the gap**. **No QED. Clay 0%.***

---

## 1 · The cancellation algebra and the intertwiner signs `[T]`

**Setup.** `H = Hess V`, coarse unit `v`, `a=⟨v,X⟩`, `K=Π₀HΠ₀`, `c=Π₀Hv`, `Π₀=1−vvᵀ`. The Savvidy Hessian
is rank-one negative: `K = K₊ − |λ₋|P₋`, `P₋=φ₋φ₋ᵀ`, `K₊⪰0`, `K₊φ₋=0`, `|λ₋|=g|B|`.

> **Lemma SC.1 (exact cancellation; rank-one sign). [T]** For any `c=c_∥+c_⊥`,
> ```
> K⁻¹ = K₊⁺ − P₋/|λ₋| ,     ⟨K⁻¹c,c⟩ = ⟨K₊⁺c,c⟩ − |⟨φ₋,c⟩|²/|λ₋| ,
> ```
> so with the exact marginal-curvature identity `G″=S−D`, `S=E[α]`, `D=Var(f)` (`Q3` Lemma Q3.2(a)):
> ```
> D = D₊ + D_Sav ,   D₊ := Var(f) − E⟨K₊⁺c,c⟩ ≤ 0  (Brascamp–Lieb on K₊⪰0),   D_Sav := E|⟨φ₋,c⟩|²/|λ₋| ≥ 0 .
> ```
> The **intertwiner sign** is the minus in front of `P₋`: the positive fibre *stiffens* the coarse
> curvature (`D₊≤0` on `log`-concave `K₊`), the negative mode *softens* it (`D_Sav≥0`). `[T]`; witness
> split error `0` (`SM5`, `H2`).

This is exactly `SM2`/`H2`/`SAVVIDY-DTERM`'s split — reproduced here as the **substrate**, not the result.
The point of the task is *which* `v` is physical.

**The one coupling, exactly, on the one-plaquette `[T]+[V]`.** In the flat-chart `SU(2)` one-plaquette
(`S2`/`G4` model: `V=h(ρ)+g(ρ)u`, `g=β sin(ρ/2)sin(θ/2)`) the coarse↔negative coupling in the
`(radial, e₃)` block is
```
c_r(r) = g′(r)/r − g(r)/r² = sin(θ/2) · (a β-stiffness in r) .
```
`c_r` is the **whole** off-diagonal Hessian entry `H_{re₃}`; `D_Sav ∝ c_r²/|λ₋|`. Hence:
- `θ=0`: `sin(θ/2)=0` ⇒ `c_r≡0` **exactly** ⇒ `D_Sav=0` ⇒ `G″=S=α>0` (fibre-confined; the negative mode,
  which *is* present in the radial fibre `K_rad<0` past the crossover, is decoupled from `e₃`).
- `θ=π`: `c_r` maximal ⇒ the negative mode is a **collective** `(radial,e₃)` mix ⇒ `G″=A−Bβ²<0` for
  `β>4.30661`.

Witness (§7·1–2): anchors `A=0.20511530`, `B=0.01105929`, `β*=4.30661` reproduced to all digits; `c_r`
exactly `0.00000` at `θ=0` and `−0.354`/`−0.533` (`r=1.5,2.5`) at `θ=π`. **The sign and magnitude of the
negative-mode contribution are a single intertwiner matrix element.**

---

## 2 · The physical 0⁺⁺ insertion and the selection `[T]` (continuum) + `[V]` (Savvidy cell)

**The carrier mismatch (the crux). `[T]`-structural.** `D_Sav` is a property of a **linear** coarse mode
`v` (an element of the tangent space, hence gauge-*covariant*, not invariant). The physical `0⁺⁺` operator is
the **gauge-invariant scalar composite**
```
O(x) = tr F²(x)   (lattice: O = Σ_{μ<ν}(1 − ½ tr P_{μν}) = S_W/β) ,
```
the standard lowest-dimension `0⁺⁺` glueball interpolator. Its coupling to the negative mode is governed by
a **selection rule**, not by a linear-mode overlap.

> **Lemma SC.2 (continuum selection; Ward-type global symmetry). [T]** The physical `0⁺⁺` state is a colour
> singlet (`R=0`) with `J^P=0⁺`. A single negative mode is a colour-**octet**, `J=1` (transverse-magnetic,
> helicity `±1`) fluctuation. Colour invariance and angular-momentum conservation give
> ```
> ⟨0⁺⁺ | O | 1-gluon⟩ = 0     ⇒     the negative mode contributes to the physical 0⁺⁺ spectral density
>                                         only at the φ₋φ₋ pair threshold s ≥ 4 m_stab² = 8 g|B| .
> ```
> Equivalently `⟨φ₋, c(v_phys)⟩=0` for the physical scalar, so `D_Sav(v_phys)=0`. [T] (global symmetries).

**Lattice realization on the Savvidy cell `[V]`.** 4D `L=1` periodic `SU(2)` (one site, four links, six
plaquettes `P_μν=U_μU_νU_μ†U_ν†`), `β=2.5`, constant field `U₁=e^{icσ₃/2}`, `U₂=e^{icσ₁/2}` (the minimal
carrier of `−D²+2F`; `SM1`). The Savvidy tachyon is the 2-fold transverse mode (`λ_min=−0.21` at `c=0.8`,
`−0.62` at `c=1.2`). Compute `φ₋·∇O` and the eigenbasis of `Hess(S+Φ)`:

| `c₁=c₂` | `λ_min` | `#neg` | `φ₋·∇O` (linear) | `φ₋ᵀHess(S/β)φ₋` (pair) | `\|φ₋·û_uniform\|` | `∇O` eigen-alignment |
|---|---|---|---|---|---|---|
| 0.8 | −0.2102 | 2 | `1.2×10⁻⁹` | −0.151 | 0.354 | `1.000·(eig +1.08)` |
| 1.2 | −0.6219 | 2 | `−1.6×10⁻⁹` | −0.317 | 0.354 | `1.000·(eig +1.55)` |

So at the symmetric saddle: **`∇O` lies exactly on the stiffest positive eigenvector**, the negative mode
is **orthogonal to `∇O` to machine zero**, the negative mode reaches `O` **only quadratically** (pair term,
strictly bounded `|·|<0.32`), and the negative mode **does** overlap the gauge-dependent uniform coordinate
(`0.354`). This is the lattice face of Lemma SC.2 (point-group/colour selection).

**Honest robustness caveat `[V]`.** The exact zero is a **high-symmetry identity** (`c₁=c₂`). Under generic
perturbations the orthogonality degrades to suppression: `max|cos(φ₋,∇O)|` `0.000` (saddle) → `0.010–0.018`
(pert `0.05`) → `0.073–0.090` (pert `0.15`), while `|cos(φ₋,û_uniform)|` stays `0.34–0.57`; over `N=1464`
generic configurations the mean squared overlaps are `cos²(φ₋,∇O)=0.0137` vs `cos²(φ₋,û)=0.0846`, a
**`≈6×`** suppression. The `L=1` cell (no continuous rotation group) can only realise the selection as
suppression, not exact orthogonality — consistent with Lemma SC.2, which is a continuum statement.

---

## 3 · The RG-scale Ward identity (trace anomaly): what it is, and what it is not `[T]`/`[NG]`

The physical `0⁺⁺` density `O=tr F²` is **exactly** the trace-anomaly operator:
```
∂_μ j_D^μ = T^μ_μ = (β(g)/2g) tr F²  (+ ...) ,   ⟨T^μ_μ⟩ = (β/2g)⟨tr F²⟩ = −(β²/2g²)·(gluon condensate) ,
```
with `β(g)<0` (asymptotic freedom). The dilatation Ward identity is thus a **global-symmetry identity
relating the physical `0⁺⁺` insertion to the generator of RG scale transformations**.

**What it gives `[T]`.** (i) The physical `0⁺⁺` normalization `M₀` and residue `f₁²` run with the anomalous
dimension set by this same `β` and the operator's `γ_O` (Callan–Symanzik) — a fixed, RG-consistent running.
(ii) It **explains** the selection of §2: the trace anomaly is a scale-derivative of the *positive-sector*
condensate, so a fluctuation at scale `√(gB)` decoupled from the condensate is decoupled from `T^μ_μ`, hence
from `O` at linear order. The negative mode is at the wrong (short) scale to feed the anomaly density.

**What it does not give `[NG]`.** A Callan–Symanzik identity fixes the **running** of the residue; it does
**not lower-bound the pole mass** `m₁`. The atom/residue `f₁²` remains the unprotected object (`K0-WARD-LOWER`
§5: `f₁²≥r₀` is a theorem only at `β<1/48`). So the RG-scale Ward identity closes the *normalization* but
carries no positivity for `m₁` — it is the third face of the same fact: bounding the negative-mode/anomaly
channel is not bounding the atom.

---

## 4 · The strict bound is a decoupling, not a floor — the no-go `[NG]`

Collecting §1–3: the negative-mode contribution to the **physical** `0⁺⁺` `k=0` spectral state is **strictly
bounded** — zero at linear order (Lemma SC.2; exact at the symmetric saddle), a bounded pair term
`|φ₋ᵀHess(S/β)φ₋|≲0.32` above the pair threshold `8g|B|`. This is the nonperturbative Ward-type relation
the task asked for, and it is genuine. **But it cannot make the physical mass strictly positive:**

> **No-go SC.3. [NG]** *A strict bound on the negative-mode contribution to the physical `0⁺⁺` state does
> **not** lower-bound the `0⁺⁺` mass.* Reason — the bound has the form "negative mode ⊥ (or suppressed
> in) the physical channel"; that is a statement that a *particular contribution is absent*, not that some
> *other* contribution is present and gapped. The mass `m₁` is supplied by the positive-sector residue
> `f₁²=ν_O({m₁²})` and its isolation (`T-WAVE`'s `A∧I`); the selection is blind to it.

Two independent cross-checks that this is the real state of affairs (not a tool artefact):

- **The corpus's `D_Sav` saturation is a gauge-dependent-coordinate effect. [T]/[V]** The measured
  `D/S≈0.97` (`H2`, `SAVVIDY-LATTICE-CELL`) is on the *softest/uniform linear mode* — where the overlap
  here is `0.354`, not zero. It does **not** transfer to the gauge-invariant insertion. Consistently, the
  gauge-**invariant** `χ_O` (`WG-INVARIANT-ESTIMATOR`, `WG-FSS-FIT`) does **not** plateau cleanly — it
  grows (`L^{0.37}`, `χ_O: 0.33→0.59` over `L=4…16`; a flat plateau is excluded at `χ²/dof=7.3`). The
  physical object behaves *differently* from the coarse mode, exactly as the carrier mismatch predicts.
- **`K0-WARD-LOWER` already showed the plain `k=0` Ward tower is one identity `=p″(0)<0`=`H1`, acid-vacuous.**
  The selection here is a *different* Ward identity (global colour/angular-momentum/parity, plus the trace
  anomaly) — it constrains *which* contributions enter, not the sign of `p″(0)`. Consistent: neither the
  old tower nor this selection proves `m₁>0`.

**Where the wall actually is.** Not the negative mode. The negative-mode `D_Sav`/Schur/TCI program
(`SM1–SM5`, `H2`, `CONSTRUCT-SAVVIDY-CELL`) all bound `D_Sav` for *linear* (gauge-covariant) modes — the
`H1` carrier. The physical `0⁺⁺` atom is a *different*, positive-sector object. So the negative-mode line
**terminates independently** (`[NG]` for "bound `D_Sav` ⇒ physical gap"), and the open step is unchanged:
the positive-sector atom mass floor (`T-WAVE`'s `A∧I`, equivalently `W∧H1`), `[O]`.

---

## 5 · Kill-condition audit — what this is and is not

The task's kill conditions, checked explicitly:

| kill if it reduces to… | verdict | reason |
|---|---|---|
| **defect split `D=D₊+D_Sav`** | **no** | The split is the *substrate* (§1), but the **result** is that for the *physical* insertion `D_Sav(v_phys)=0`/suppressed — the *opposite* regime to the split's `D_Sav=O(S)` for the gauge-dependent mode. It is a statement about *which carrier* the split applies to, not a re-derivation of the split. |
| **`H1`** (uniform coarse-mode curvature floor `G″_{k→0}≥λ`) | **no** | `H1` is a *floor* on the gauge-covariant linear mode; this is a *decoupling* of the negative mode from a *different* (gauge-invariant) operator, and it explicitly demonstrates `H1`'s carrier ≠ the physical atom's carrier. It does not prove `H1` and does not use it. |
| **source alignment `W`** (`w₁=f₁²/M₀`) | **no** | `W` is a *value* (the strong-coupling weight floor and its continuation); this is a *selection/decoupling* (zero linear coupling, bounded pair), fixing no weight. `W` concerns the positive-sector residue; the selection is silent on its magnitude. |
| **the gap** | **no** | Split-neutral, `H1`-neutral, `W`-neutral, and gap-neutral: it is a `[NG]` for the negative-mode route. |

**Positive content that is genuinely new here:** (i) the exact identification of the whole negative-mode
coupling with the single intertwiner element `c_r∝sin(θ/2)` on the exact model (§1); (ii) the **physical-
vs-linear carrier separation** with an exact lattice selection identity (`φ₋·∇O=0` at the symmetric saddle,
`∇O`=stiffest positive eigenvector) and its `≈6×` generic suppression (§2); (iii) the explicit demonstration
that this strict bound is a decoupling (`[NG]`, §4) and *why* (colour + angular momentum + trace anomaly).
None of the three is in the corpus.

---

## 6 · Acid tests (charter) `[T]`+`[V]`

1. **Free massless — must fail (no manufactured floor). [T]** `F=0` ⇒ no Savvidy background ⇒ **no negative
   mode** (`λ₋→0⁻`), and `f₁²=0` (no bound state). The selection returns "zero negative-mode contribution"
   and produces **no gap**; the witness's free row gives `G″=k²→0`. Passes.
2. **Cancellation must not require the negative mode to vanish. [T]** At the Savvidy cell the negative mode
   **exists** (`λ_min=−0.21, −0.62 <0`) and the cancellation holds anyway — it is decoupled from the scalar
   (`φ₋·∇O=0` at the saddle), not removed.
3. **Discriminating (the compact-measure trap). [T]** The one-plaquette `β=0` Haar floor `G″=0.2051` and
   the physical-insertion marginal curvature are *compact-measure* curvatures (bounded support), **not**
   free-field or nonperturbative floors — they do not certify an atom (`H1_does_not_imply_atom`, `P6`). The
   selection therefore makes no positivity claim that could be mistaken for a gap.
4. **`β<1/48` contrast. [T]** In the strong-coupling regime the full LSI holds (`SZZ`, `K_S=N(½−24|β|)>0`)
   and the marginal is log-concave — there the negative mode is absent/controlled and the selection is
   vacuous-but-true; the head `β≈2.5` is the disjoint regime where the selection is non-trivial and the
   atom is open. Consistent with `SM5` §2.6.

---

## 7 · Witness — `research/v3/w7/sav_cancel.py` (`py -3.12 -I`, deterministic)

Quoted output (abridged; full in `sav_cancel.out`):

```
(1) ANCHORS: A=0.20511530 B=0.01105929 beta*=4.30661
    beta=2.5 th=pi: S=0.20512 D=0.06912 G''=0.13599 E[Sigma]=0.05549   (A-B b^2 = 0.13599)
    beta=6.0 th=pi: S=0.20512 D=0.39813 G''=-0.19302 E[Sigma]=-0.65673 (A-B b^2 = -0.19302)
(2) INTERTWINER c_r = sin(th/2)*(...):   th=0 -> 0.00000/0.00000 ;  th=pi -> -0.35433/-0.53264 (r=1.5,2.5,beta=6)
(3) SAVVIDY CELL (4D L=1):
    c1=c2=0.8: lam_min=-0.2102 #neg=2 |gradO|=0.1538 ; gradO_phys = 1.000*(eig +1.080)
        negmode: phi.gradO=1.2e-9 [LINEAR=0] ; phi^T Hess(S/b) phi=-0.15148 [PAIR] ; |phi.uniform_s3|=0.3536
    c1=c2=1.2: lam_min=-0.6219 #neg=2 |gradO|=0.4202 ; gradO_phys = 1.000*(eig +1.554)
        negmode: phi.gradO=-1.6e-9 [LINEAR=0] ; phi^T Hess(S/b) phi=-0.31708 [PAIR] ; |phi.uniform_s3|=0.3536
(4) ROBUSTNESS: max|cos(phi-,gradO)| 0.000(saddle)->0.018(pert.05)->0.090(pert.15) ; uniform stays 0.34-0.57
    generic configs: mean cos^2(phi-,gradO)=0.0137 vs cos^2(phi-,uniform)=0.0846  (suppression ~6x)
(5) ACID: free massless G''=k^2->0, no negative mode, no atom; the negative mode need not vanish.
```

Evidence quality: (1) exact Gauss–Legendre quadrature (`n=2400–4000`), `A,B,β*` to all printed digits;
(2) exact closed-form `c_r`; (3) exact `12×12` finite-difference Hessian + eigen-decomposition, orthogonality
to `10⁻⁹`; (4) deterministic perturbation/generic ensembles (`N=1464`); (5) exact. **A witness, not a proof**;
the theorem content is Lemmas SC.1–SC.2 and No-go SC.3; the open step is the positive-sector atom (`[O]`).

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| exact cancellation + rank-one sign `K⁻¹=K₊⁺−P₋/|λ₋|`, `D=D₊+D_Sav`, `D₊≤0`, `D_Sav≥0` | substrate (Lemma SC.1) | **[T]** |
| one-plaquette: whole coupling = single element `c_r=sin(θ/2)·(…)`; `c_r≡0`@`θ=0`, max@`θ=π` | exact (§1, §7·2) | **[T]**+**[V]** |
| `G″=A−Bβ²`, `A=0.20511530`, `B=0.01105929`, `β*=4.30661` | exact quadrature | **[T]**+**[V]** |
| physical `0⁺⁺` = gauge-invariant scalar composite `tr F²`; `D_Sav` is a **linear-mode** object | carrier mismatch | **[T]**-structural |
| continuum selection `⟨0⁺⁺|O|1-gluon⟩=0` ⇒ negative mode only at pair threshold `s≥8gB` | colour + angular momentum | **[T]** |
| Savvidy cell: `φ₋·∇O=0` exactly at symmetric saddle; `∇O`=stiffest positive eigvec | exact §7·3 | **[T]**-in-model+**[V]** |
| generic configs: linear overlap suppressed `≈6×` in `cos²` vs gauge-dependent coordinate | §7·4 | **[V]** |
| negative mode couples to `O` only at 2nd order, pair term `|φ₋ᵀHess(S/β)φ₋|<0.32` | §7·3 | **[V]** |
| trace-anomaly Ward identity `T^μ_μ=(β/2g)tr F²`; fixes running, not the pole mass | §3 | **[T]**-identity / **[NG]**-tool |
| **strict bound ⇒ physical `m₁>0`** | **false**: bound = "negative mode absent", not "atom gapped" | **[NG]** (No-go SC.3) |
| corpus `D/S≈0.97` is the softest **linear** mode; gauge-invariant `χ_O` grows `L^{0.37}` (no plateau) | carrier confirmed by `WG-FSS-FIT` | **[T]**+**[V]** |
| atom mass floor in the positive sector (`A∧I` / `W∧H1`) | the open step, unchanged | **[O]** |
| free massless: no negative mode, no atom, no manufactured gap | acid §6 | **[T]**+**[V]** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 9 · References (internal)

- `Q3-SCHUR-DEFECT.md` §1–§2 (inertia lemma, one-plaquette dichotomy, `G″=S−D`, `E[Σ]`), §4.
- `SM2-RANK-ONE.md` §1 (overlap lemma `⟨φ₋|c⟩=⟨H_cfφ₋|v⟩`, helicity projection), §0/§3 (ceiling `m_soft≤G″`).
- `H6-HEAD-SPECTRAL.md` §3 (typical point vs essential infimum); `SM5-TCI-HEAD.md` §0–§3 (`TCI(strong)` false,
  `TCI(avg)=I(ρ)=H1`); `H2-DEFECT-BOUND.md` §1 (`D=D₊+D_Sav`), §2.4 (IR-marginal, Landau pole).
- `G4-GENERAL-BOUND.md` §2–§3 (`G″=A−Bβ²`, `e₃` barrier); `SM1-SOFT-EXACT.md` §2 (`L=1` Savvidy tachyon,
  transverse 2-fold); `SM4-CONSTANT-MODE.md`; `CONSTRUCT-SAVVIDY-CELL.md` §2–§4 (scale separation + bounded
  softening); `SAVVIDY-WI-SYNTHESIS.md` §2 (WI routes = A3 = H1).
- `research/unconditional/K0-WARD-LOWER.md` §2–§5 (the plain `k=0` Ward/Stein tower is one identity
  `=p″(0)<0`=`H1`, acid-vacuous; `f₁²≥r₀` unprotected).
- `T-WAVE-SYNTHESIS.md` §1–§5 (predicates `A,W,I,G`; `I⇒A`, `W⇒A`, gap via `W∧H1` or `I∧(ρ²≤C²)`).
- `research/v3/WG-INVARIANT-ESTIMATOR.md` §0–§6, `research/v3/WG-FSS-FIT.md` §0–§3 (gauge-invariant `χ_O`,
  `L^{0.37}` growth, ceiling-not-edge).
- Standard (textbook, cited from corpus, not re-fetched): Savvidy, Phys. Lett. B **71** (1977) 133;
  Nielsen–Olesen, Nucl. Phys. B **136** (1978) 45; Brascamp–Lieb (1976); Prékopa (1971); Haynsworth/Schur
  inertia; Callan–Symanzik / trace anomaly (Adler–Bardeen–Bell–Jackiw; Collins); Shen–Zhu–Zhu CMP **400**
  (2023) 805, arXiv:2204.12737.
- Witness `research/v3/w7/sav_cancel.py` → `sav_cancel.out` (`py -3.12 -I`, numpy only, deterministic;
  one process). No third-party sources downloaded.

*No QED. Clay 0%. The wall for the physical mass is unmoved; the negative-mode line is shown to bound a
different carrier, and closes on a strict decoupling rather than a floor.*
