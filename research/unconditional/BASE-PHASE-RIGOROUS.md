# BASE-PHASE-RIGOROUS — can rigorous statistical mechanics certify a phase/gap **at** the crossover β≈2.5?

> **VERDICT (no QED).** **No.** Every rigorous mass-gap / phase-certificate machinery for lattice gauge theory is an **absolute (small-β / high-temperature) convergence** method, and the crossover lies **outside the convergence radius of all of them**. The peer-reviewed mass-gap radius ends at **β<1/48** (Shen–Zhu–Zhu: λ>96); the crossover is λ≈7.3 (β_W≈2.47) — **≈13× in λ (β_W 0.19 vs 2.47), an analyticity wall, not a constant to optimise**. The one route that crosses to large β (loop-equation / area-law, Cao–Nissim–Sheffield) stops **one hard step short**: it proves **area law, not the transfer-matrix mass gap**. Pirogov–Sinai is inapplicable (no discrete zero-T order parameter); Dobrushin/disagreement-percolation thresholds sit at β = O(0.1–0.2), an order of magnitude below 2.5; a Peierls/contour surface-tension bound at β=2.5 would be **equivalent to the target**; Balaban controls the **UV** (β→∞) but not the IR gap. **The crossover is genuinely beyond all these machineries** — `[NG]` for the class. The *base* is therefore: **two solid ends (β<1/48 and β→∞) with a notch in the IR middle exactly where the creature must stand.** `Clay 0%`.

---

## §1 · The honest distance (numbers, not adjectives)

| datum | value | source |
|---|---|---|
| proven mass-gap radius | **β<1/48** ⟺ `λ≡g²N>96` (all d, all N; d=4: 1/(16·3)) | Shen–Zhu–Zhu, arXiv:2204.12737 / *CMP* **400** (2023) 805 |
| crossover (Hooft) | `λ≈7.3` (SU(3)), `β_W≈2.47`; in SZZ units `β≈0.27` | measured scaling onset, `ATOM-strong-coupling-extend` §1 |
| distance to cross | **≈13× in λ** (`β_W` 0.19 → 2.47) | `ATOM-strong-coupling-extend` §1 |
| character of the barrier | **analyticity radius** of absolute methods; **not** a thermodynamic transition | `ATOM-strong-coupling-extend` §1, §5 |
| best 2026 reach (unrefereed) | β<195/7168≈0.027 | gauge-orbit-curvature preprint, **unrefereed** |

Even a rigorous proof at β=2.5 would not finish Clay: the lattice `β` is the *bare* coupling, and the continuum theory needs the **β→∞** limit along the scaling curve (the second, separate gap = `I2`). This file addresses the **first** gap (1/48 → 2.5).

---

## §2 · Machinery by machinery: what it gives, where **exactly** it fails at β≈2.5

**(a) Cluster / log-Sobolev (Shen–Zhu–Zhu, Osterwalder–Seiler).** `[T]` for β<1/48: LSI + Poincaré + exponential clustering ⇒ volume-/L-uniform mass gap, together with the area law. **Fails** at β≈2.5 because the cluster/expanding-map radius of convergence is finite and ~two orders of magnitude below the crossover; SZZ's own constant is `1/(16(d−1))`. **No known variant** enlarges it — the radius is structural for the absolute-convergence class.

**(b) Reflection positivity (Osterwalder–Schrader).** `[T]`: gives OS axioms ⇒ Hilbert space/transfer matrix (G1 side). **Not a gap certificate.** The massless free scalar is reflection-positive with `Δ=0`; RP cannot exclude a gapless atom. `[NG]` for gap certification.

**(c) Dobrushin uniqueness / disagreement percolation (Dobrushin 1968; van den Berg–Maes; Dobrushin–Shlosman).** `[T]` for weak coupling: a unique Gibbs state and exponential clustering, hence a gap, **when the Dobrushin coefficient `α_D<1`.** Estimate for 4D SU(N): each link touches `2(d−1)=6` plaquettes; `α_D ≲ 6β·c` with `c=O(1/2)` ⇒ threshold `β = O(0.1–0.2)`. **At β=2.5, `α_D≈O(7)≫1`** — fails by roughly an order of magnitude, and worsens as N grows (Hooft scaling). `[NG]` at the crossover.

**(d) Pirogov–Sinai / contour (PS 1975–76).** `[T]` for models with a **discrete order parameter and ≥2 coexisting ground states** (Potts, lattice Higgs, vortex phases). **Inapplicable at zero T**: pure SU(N) gauge has no discrete order parameter for the `0⁺⁺` sector; the centre `Z_N` breaks only **at finite temperature** (deconfinement), not at the zero-T mass gap. PS could at best speak to a finite-T transition, not `G2`. `[NG]` (wrong object).

**(e) Area-law / loop-equation route (Cao–Nissim–Sheffield).** `[T]`: Wilson **area law** for β<c_d (absolute, N-independent; arXiv:2505.16585) and into the Hooft regime β~N (2509.04688). **This is the only route that reaches large β** — but it proves **confinement (area law), not the transfer-matrix gap**. `area-law ⇒ mass gap` is **unproven**, and the missing bridge `m_{0⁺⁺}=c√σ` is a measured transmutation, in a **different (0-form vs 1-form) channel**. Stops **one hard step short** = the Clay step. `[O]` at exactly that bridge.

**(f) Balaban renormalisation group (CMP 109/116/122).** `[T]`: UV stability, renormalisation and a controlled effective action as **β→∞** (continuum). **Controls the UV, not the IR gap.** It does not produce a positive `m₁`; it guarantees the continuum construction is well-defined (the G1 side). `[NG]` for the IR gap.

**(g) Adhikari–Cao (Ann. Prob. 53 (2025) 140).** `[T]`: weak-coupling results — but only for **finite gauge groups**; **no transfer to continuous SU(N)**. `[NG]`.

**(h) Peierls / contour surface tension at β=2.5.** The conceivable positive variant: a contour/Peierls estimate needs a **rigorous lower bound `σ(β)>0` on the surface tension at β=2.5**. But `σ(β=2.5)>0` (perimeter-law lower bound on a Peierls contour) is essentially **confinement at β=2.5**, i.e. **the target itself**. `[O]` — circular.

---

## §3 · Is the crossover rigorously "in a phase"?

Two distinct questions, often conflated:

- **"Is there a thermodynamic phase transition at β≈2.5?"** Believed **no** (zero-T pure SU(N) is a single confined phase for all β>0; the "crossover" is a smooth scaling-onset, not a singularity). So "the crossover is in a phase" is believed **yes** — but this is *not* what Clay needs.
- **"Is the `0⁺⁺` mass gap positive at β≈2.5 (and uniform as β→∞)?"** **Open.** No rigorous machinery certifies it there.

**The base shape.** The rigorous frontier is solid at **both ends** — deep strong coupling (`β<1/48`, SZZ) and the UV continuum (`β→∞`, Balaban) — but **the IR middle (the scaling region's gap) is a notch exactly where the creature must stand.** That notch is not a missing constant; it is the absence of a small parameter at `α_s~O(1)`.

---

## §4 · Could any variant work?

The only two conceivable positives, and why each is not a construction:
1. **A rigorous surface-tension lower bound at β=2.5** (contour/Peierls) — **equivalent to the target** (§2h). Circular.
2. **An exact solvable large-N / 't Hooft (or Gross–Witten) corner** — the large-N gap is rigorous only at fixed **small** β (strong coupling, already gapped); `N=∞→3` is asymptotic with an unbounded `c/9` remainder, and `β=2.5, g≈1.5` needs `N<0.07` for the strong-coupling angle — the regions are **disjoint** (`PROOF-LARGE-N-RIGOROUS`, `ABC2-DUALITY`). Not a construction.

**No variant of the absolute-convergence or contour class is known to cross the analyticity wall at finite N.** This is a structural no-go, not an optimisation gap.

---

## §5 · Verdict

- **`[NG]`** for the **entire class** of rigorous phase/mass-gap certificates (cluster/LSI, RP, Dobrushin, disagreement percolation, PS, Peierls/contour, Balaban) applied **directly at β≈2.5**. Each fails at a named step (§2).
- **`[O]` (exact missing step)**: the IR gap in the scaling region at `α_s~O(1)`, whose only "condition" is the absence of a small parameter — the same wall as `PROFILE`/`SLOT-SYNTHESIS`, reached from the rigorous-stat-mech side.
- **What the BASE actually is:** a plinth solid at `β<1/48` and at `β→∞`, with a notch in the IR middle. The notch is *not* where a constant is missing — it is where **every absolute method's analyticity radius ends and no non-absolute method reaches the gap** (the loop-equation route reaches area law but not the gap).
- **Honest scope:** not a theorem that *no conceivable* rigorous method works — it is that **every named, standard machinery fails at a precisely located step**, and the one non-absolute escape stops one hard step short. `Clay 0%; no QED.`

*Labels:* `[T]` proven machinery (β<1/48 gap; RP axioms; area law β<c_d; Balaban UV); `[NG]` the class at β≈2.5; `[O]` exact missing step (IR gap at scaling; area-law→gap bridge; surface tension at 2.5).
*Refs.* Shen–Zhu–Zhu, arXiv:2204.12737 / *CMP* **400** (2023) 805; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Cao–Nissim–Sheffield, arXiv:2505.16585 & 2509.04688; Adhikari–Cao, *Ann. Prob.* **53** (2025) 140; Balaban, *CMP* **109/116/122** (1987–89); Dobrushin (1968); van den Berg–Maes (1994); Dobrushin–Shlosman; Pirogov–Sinai (1975–76); Osterwalder–Schrader (1973/75); Brydges–Fröhlich–Seiler, *CMP* **83** (1982) 123; Fröhlich–Spencer. Internal: `ATOM-strong-coupling-extend`, `PROOF-P3-UNCONDITIONAL`, `PROOF-LARGE-N-RIGOROUS`, `ABC2-DUALITY`, `PROFILE`, `SLOT-SYNTHESIS`.
