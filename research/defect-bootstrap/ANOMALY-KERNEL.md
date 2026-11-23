# How the `Z_N^{(1)}`–θ anomaly enters the crossing kernel

> CPB-attack probe (sub-task 7). Verdict: the anomaly is an **irremovable 2-cocycle phase**;
> it is a **selection rule**, not a gap proof.

1. **Anomaly.** Electric `Z_N^{(1)}`; mixed with θ: a `Z_N^{(1)}` rotation shifts
   `θ → θ + 2π/N` (Witten effect), encoded by a 5d BF-type symmetry TQFT (GKSW 2015;
   Cordova–Tachikawa).

2. **On the line.** A line of charge `q∈Z_N`; its worldline theory carries a 0-form `Z_N`
   symmetry (defect operators grade by `q_i`). Inflow gives a **1d 't Hooft anomaly** for it,
   classified by `H²(Z_N,U(1)) = Z_N` with level `k` (θ enters as the level/background).

3. **`Z_N` grading.** `H_def = ⊕_{q_i} H^{(q_i)}`; the line of charge `q` fuses `q_i → q_i+q`.

4. **Modified crossing.** The anomaly is an **irremovable 2-cocycle phase** in defect
   fusion/crossing: `K_{ij} → ω(q_i,q_j) K_{ij}`, `ω(q_i,q_j) = e^{2πi k q_i q_j/N}`
   (projective `Z_N` rep), not removable by defect phase redefinitions. The bulk–defect
   crossing for a neutral bulk `φ`,

   ```
   ⟨φφ⟩_def(q) = Σ_O C^{φφ}_O 𝒢_O(ξ) = Σ_{q_i} a_i² 𝔽_i(ξ),
   ```

   keeps its form but the defect sum runs over the `Z_N`-graded sectors with **relative
   phases fixed by `ω`**.

## Honest boundary

The exact phase (`e^{iθ/N}` vs `e^{2πi k q_iq_j/N}`, and the θ-normalisation) is
BF-convention-dependent and must be fixed against the bootstrap convention. This is a
**selection rule**, not a gap proof.
