# Q-WAVE-SYNTHESIS — attacking H1: seven routes, **four equivalent missing lemmas**, and the **carrier hazard** that demotes H1 from finish line to stepping stone

> **VERDICT (no QED, Clay 0%).** Seven routes attacked H1, the last open input of the ultraviolet
> face: `χ = Var_μ(ξ_{k→0}) ≤ 1/λ` uniformly in `(L,a)` at `β=2.5`. **None proves it.** The wave's
> positive content is structural: the standard method classes fail for **directional** reasons, and
> their four escape hatches are the **same inequality in four faces** (a variance bound, a connected
> four-point bound, a fibre gap, and a single-well/gate-density bound). The wave's negative content is
> sharper: Q6 shows H1 is **strictly weaker than the gap**, and that on the gauge-invariant composite
> `tr F²` the H1 ceiling is **vacuous** (`Var ∝ 1/V`), while on the gluon carrier it is
> gauge-dependent. So H1 is a necessary stepping stone, not the finish line; the irreducible target
> is the **atom** on the `0⁺⁺` Källén–Lehmann measure.

---

## §0 · The seven verdicts

| Route | Verdict | Result |
|---|---|---|
| **Q1** non-Abelian infrared bound | `[NG]` | The class {FSS Gaussian domination, chessboard/RP, compact transfer, heat-kernel/LSI, Dobrushin} cannot reach `β=2.5`. FSS is **vacuous at `k=0`** (the free massless reference has `Ĉ(0)=∞`); the free noncompact photon is RP with `χ=∞`; compactness gives only `Var≤π²N_links=Θ(L⁴)`; the interaction reweighting is `e^{2βN_p}=e^{Θ(L⁴)}`; Dobrushin `α(2.5)=180`. Missing lemma: `Var_μ(a) ≤ Var_Haar(a)+C₀`, `C₀` independent of `(L,a)`. |
| **Q2** loop bootstrap | `[NG]` | The finite-N Makeenko–Migdal equation closes in `W` only at `N=∞`, leaving an `O(1/N²)` connected two-loop term. Writing `χ=χ_G+χ_conn`, equivalently `G″=S−D`, the identity checks to 6 digits. Reflection positivity gives `D≥0`, which only **lowers** `G″`—the wrong direction. UV asymptotics fixes only large-`s` moments. Missing term: `D ≤ S−λ` (the Savvidy rank-one part `D_Sav`). |
| **Q3** integrate out the negative mode | `[T]`/`[NG]` | Exact inertia lemma (Haynsworth): for `H=[[α,cᵀ],[c,K]]`, `Σ=α−cᵀK⁻¹c`, `n₋(H)=n₋(K)+n₋(Σ)`, so `Σ≻0` iff every negative mode is fibre-confined. One-plaquette `θ=0`: `c_r=0`, `K` indefinite but `Σ=α>0` — excision works. `θ=π` barrier: `K≻0` yet `Σ<0`; the negative mode is **collective** `(radial,e₃)`, aligned with the coarse direction — no fibre mode to excise, positivity not restored. The rank-one lemma `χ≤(1+η)/ρ_rest` reduces to `χ≤1/ρ_*`; the `(1+η)` is cosmetic and small `D/S` does not imply `ρ_rest>0`. Missing 4D estimate: a uniform fibre gap `K⪰γI`, i.e. `E⟨K⁻¹c,c⟩<S`. |
| **Q4** Bakry–Émery with defect | `[T]`/`[C]`/`[NG]` | `Ric_{SU(N)}=(N/4)g>0` is exact (witness 1e-15), but `Ric+Hess S_W≈−2.3<0`: the Savvidy Hessian dominates. The `k=0` projection is lossless, so the marginal BE constant **is** `ess inf G″`; group Ricci is silent on χ. Weighted rescue via exact 1D Cheeger `λ₁≥h²/4`, so `χ≤4/h²`; finite uniform χ iff the `k=0` marginal is single-well with gate density bounded below. 1D SU(2): `λ₁=0.2645`, `χ≤3.78`; `β=2.5` floor `0.222` gives `χ≤4.5`. Defect measure is the **wrong control**: `ν(G″<0)` shrinks 41× while `λ₁` collapses 170× (Kramers). Missing estimate: a uniform lower bound on the Savvidy gate density, i.e. single-wellness of the 4D `k=0` marginal. |
| **Q5** 4D numerical plateau | `[V]`/`[O]` | Pilot at `β=2.5`: `χ=Var(a)=4.17±0.18` (`L=4`) and `3.61±0.24` (`L=6`), ratio `0.87`; a free massless field would rise by `(6/4)²≈2.25`. Growth is excluded; a plateau at the Haar value `≈3.7` is consistent. Polyakov correlator decays within `ξ≲1–2` sites. Caveats: the coordinate is gauge-invariant only infinitesimally (≤4% drift under local transforms); the exact force-moment `G″(k=0)=+0.03±0.08` is unresolved (broad marginal, `sd=2.15`); composites are vacuous (`Var∝1/V`). Decisive run specified: Landau+FP+Gribov, `L=8,12(16)`, `N≥2000`, one job at a time. |
| **Q6** logical status of H1 | `[T]`/`[NG]` | H1 is **strictly weaker than the gap**: `gap⇒H1` is a theorem; the converse is false and machine-checked (`H1_does_not_imply_atom`; thin atom `ν=wδ_ε+(1−w)δ_{s₂}`, `w=ε/2λ`, has `χ=1/λ` yet `inf supp ν=ε→0`). The three claimed equivalences (clustering / transfer gap / `G″≥λ`) are false as biconditionals and hold only in the Gaussian-marginal limit. H1 is a theorem in seven regimes (strong coupling `β<1/48`, `β=0` Haar `χ=3.7198`, 1D chain, finite volume non-uniformly, 2D, large-N strong coupling, UV face `[C]` on a different object) and open in the eighth (4D `SU(N)`, `β=2.5`). The free massless `U(1)` field violates H1. **Carrier hazard:** on the gauge-invariant composite `tr F²`, `Var(q̄)=χ/V→0`, so H1 is vacuously true; on the gluon `a=⟨v,log U⟩` it is gauge-dependent. |
| **Q7** sixth-species forge | `[T]`/`[NG]` | Seven candidate families; five collapse onto banked objects. The two survivors are computed no-gos. **Direction theorem:** analytic-duality/Fisher–Rao functionals give `E_ρ[G″]=I(ρ)≥1/Var`, i.e. `Var≥1/I`—a **variance floor**, the logical opposite of H1; Bernstein gives `G″≤τ²` (a ceiling), and a band-limited `ρ_ε=1+εcos(τa)` has `G″=ετ²/(1+ε)→0` at fixed bandwidth. **Frozen-outlier:** the Savvidy soft mode is an exact deterministic eigenvalue `−g|B|`, not a BBP fluctuation; on `G″=S−D` with `D/S≈0.97`, the gap is a 3% residue, below every available RMT relative error (`≥10%`). |

---

## §1 · The four faces of the missing lemma

The wave's main positive result is that the four escape hatches are the **same inequality**. Each is a
different coordinate on H1:

```
   (i)   Var_μ(a) ≤ Var_Haar(a) + C₀              (Q1: finite-mass Gaussian domination)
   (ii)  D ≤ S − λ                                 (Q2: connected 4-point / softened H1)
   (iii) K ⪰ γ I   (equivalently E⟨K⁻¹c,c⟩ < S)    (Q3: uniform fibre gap)
   (iv)  ρ(gate) ≥ c > 0  (single-well marginal)   (Q4: gate density / Cheeger capacity)
```

Under the stated normalizations these are equivalent to `χ≤1/λ`; the equivalence is `[T]` where the
marginal is Gaussian and `[C]` otherwise. This is a **robust characterization**, not a proof: it says
that a proof of any one face would be a proof of H1, and that the campaign now knows four independent
ways to state the same missing input. It also explains why the routes fail together: they are all
attacking the same object from directions that the Savvidy tail blocks.

---

## §2 · The direction theorems: why the standard methods fail

Four failures are **directional**, not quantitative. This is the wave's sharpest structural lesson.

1. **Reflection positivity is the wrong direction at `k=0`.** FSS controls non-constant modes; the
   free massless reference has `Ĉ(0)=∞`; the free photon is itself RP with `χ=∞`. Any RP bound that
   also bounds the free photon is vacuous (acid test). On the loop side, RP gives `D≥0`, which only
   lowers `G″`.
2. **Fisher information is a variance floor.** `E[G″]=I(ρ)≥1/Var` is equivalent to `Var≥1/I`: a lower
   bound on the variance, the opposite of H1's upper bound. Bernstein's inequality caps `G″` above;
   the band-limited family `ρ_ε=1+εcos(τa)` drives `G″→0` at fixed bandwidth. The exact anchors
   `0.2051`, `0.222` are therefore **ceiling-side** values, not floor-certificates.
3. **Compact transfer and Dobrushin are volume-blind at `β=2.5`.** Compactness gives only
   `Θ(L⁴)`; the interaction reweighting is `e^{Θ(L⁴)}`; Dobrushin `α=180`. None is uniform.
4. **The negative mode is collective where it matters.** The inertia lemma shows that a fibre
   negative mode can be excised, but at the `θ=π` barrier the negative mode is aligned with the coarse
   direction and `K≻0` while `Σ<0`; there is no fibre mode to remove. The rank-one defect rescue is
   cosmetic: `(1+η)` does not control `ρ_rest`, and small `D/S` does not imply `ρ_rest>0`.

Two further warnings:

- **Defect measure is the wrong control.** `ν(G″<0)` can shrink while `λ₁` collapses (Kramers);
  `D/S≈2%` is an averaged variance ratio, not a gate-density bound.
- **RMT is blind at the required precision.** The Savvidy outlier is deterministic, and the 3% residue
  `G″/D` is below the `≥10%` relative error of every available RMT/free-probability result.

---

## §3 · The 4D numerical status, and why `G″≥λ` is not the right face

Q5 gives the first controlled 4D pilot: `χ=4.17±0.18` (`L=4`) and `3.61±0.24` (`L=6`), ratio `0.87`,
against a free-field growth of `2.25`. Growth is excluded; the plateau is consistent with the Haar
value `≈3.7`. But the exact force-moment estimator gives `G″(k=0)=+0.03±0.08`, while `1/χ≈0.24`:
`R=G″·χ≈0.13`, far from the Gaussian value `1`. **In 4D the marginal is strongly non-Gaussian**, so
the equivalence `H1 ⇔ G″≥λ` fails quantitatively. H1 must be stated as the **variance bound**
`χ≤1/λ`, not as a pointwise curvature bound, unless Gaussianity is proved. This is a correction to
the campaign's shorthand.

The pilot also inherits the carrier caveat: the flat-chart coordinate is gauge-invariant only
infinitesimally (≤4% drift under moderate local transforms), so it is evidence, not a gauge-invariant
measurement. The decisive run is specified (Landau+FP+Gribov, `L=8,12(16)`, `N≥2000`, one job).

---

## §4 · The carrier hazard (Q6): H1 alone is not the finish line

Q6 identifies the deepest problem. H1 can be read on three carriers:

| carrier | H1 status | problem |
|---|---|---|
| coarse field `ξ_{k→0}` (gluon log) | non-vacuous, but gauge-dependent | not a physical observable; acid test ambiguous |
| gauge-invariant composite `tr F²` | **vacuously true** (`Var∝1/V→0`) | carries no gap content |
| gauge-invariant non-self-averaging carrier | would be the right object | **its linear existence is the structural theorem that does not hold** |

So "H1 ⇔ gap" presupposes a carrier that the campaign has not constructed. On the composite, H1 is
empty; on the gluon, it is gauge-dependent; the third carrier is exactly the missing structural
object. This is why H1 is strictly weaker than the gap, and why the irreducible target is the
**atom**:

> **Target `H*`:** `(f₁² ≥ r₀ with isolation I₂) ∧ (χ ≤ 1/λ)` on a total-normalized gauge-invariant
> carrier. Given the pair, the gap follows `[T]`; the single irreducible object is the `0⁺⁺`
> Källén–Lehmann atom. H1 supplies only the ceiling; the atom is the edge.

The falsification battery `F1–F5` is the guard: F1 free field (must fail), F2 thin atom (χ alone must
not certify), F3 deep double well (no barrier-uniform λ), F4 carrier/vacuity (must be non-vacuous), F5
dimension (`G″_a=a²m_phys²`, a ratio).

---

## §5 · What the wave won

1. **An exact inertia lemma** (Q3): `n₋(H)=n₋(K)+n₋(Σ)`, which decides when integrating out a negative
   mode can restore positivity. It is a reusable structural tool, and it explains the one-plaquette
   `θ=0` success versus the `θ=π` failure.
2. **A direction theorem** (Q7): Fisher–Rao/duality functionals are variance floors, not ceilings;
   Bernstein caps curvature; band-limited densities drive `G″→0`. This corrects the corpus's use of
   the Fisher identity as a positive H1 tool.
3. **A defect-measure warning** (Q4): measure of the negative set is not a proxy for the spectral gap;
   the Kramers collapse is a counterexample to the natural intuition.
4. **A 4D pilot plateau** (Q5): growth is excluded at `β=2.5` for `L=4,6`; the plateau is consistent
   with the Haar value. This is evidence, not proof, and it is carrier-limited.
5. **A four-face equivalence** (§1): the missing lemma has four independent statements, all H1.
6. **A reusable falsification battery** `F1–F5` (Q6).

---

## §6 · What the wave did not win

- **No proof of H1** in 4D at `β=2.5`. Every route ends at one of the four faces of §1.
- **No mass gap.** H1 is strictly weaker (Q6); the atom is untouched.
- **No carrier.** The gauge-invariant non-self-averaging carrier that would make H1 non-vacuous is not
  constructed; on the composite H1 is empty, on the gluon it is gauge-dependent.
- **No sixth species.** Q7's census collapses seven candidate families; five reduce to banked
  objects, two to sharp no-gos. The sixth species is not found.
- **Clay 0%. No QED.**

---

## §7 · Strategic consequence

H1 remains a legitimate **stepping stone**: it is necessary, non-circular, and a proof at `β=2.5`
would be a genuine theorem. But the Q wave shows it is not the headline target. The campaign should:

1. **State H1 on a definite carrier** and prove the carrier is non-vacuous (F4); otherwise the
   statement is empty.
2. **Run the decisive 4D measurement** (Q5 protocol) to settle the plateau on a gauge-invariant
   observable.
3. **Shift the headline to the atom**: `f₁²≥r₀` with isolation `I₂` on the `0⁺⁺` Källén–Lehmann
   measure, with `χ≤1/λ` as the necessary ceiling. The atom is the irreducible object; H1 is its
   ceiling, not its substitute.
4. **Use the four faces of §1 as the checklist** for any proposed proof of H1, and `F1–F5` as the
   falsification battery for any proposed proof of the gap.

**One line.** *Seven routes to H1 all end at the same four-faced lemma; the standard methods fail for
directional reasons; the 4D pilot plateaus but is carrier-limited; and H1 is strictly weaker than the
gap, with a carrier hazard that makes the atom, not H1, the irreducible target. No proof, no QED,
Clay 0%.*

*Refs:* `Q1-INFRARED-BOUND.md`, `Q2-LOOP-BOOTSTRAP.md`, `Q3-SCHUR-DEFECT.md`, `Q4-BAKRY-EMERY.md`,
`Q5-4D-CHI.md`, `Q6-H1-STATUS.md`, `Q7-TOOL-FORGE.md`; witnesses `w7/q1_*.py` … `w7/q7_*.py`;
`P-WAVE-SYNTHESIS.md`, `STAR2-CLASSIFICATION.md`, `H6-HEAD-SPECTRAL.md`, `D6-DECISION.md`.
