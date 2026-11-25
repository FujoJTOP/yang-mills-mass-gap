# W-WAVE-SYNTHESIS — no bronze theorem yet: W1/W2/W3 end at precise missing steps, WG gives the first gauge-invariant plateau

> **Verdict (no QED, Clay 0%).** The W wave did **not** produce a theorem at β=2.5. It produced one
> genuine numerical measurement (WG), one exact-coefficient result with a conditional large-N window
> (W1), one closed structural equation without a floor (W2), and one negative (W3: all exact
> solvable small volumes have `w₁=0`). **WG's extension corrects the flat-chart picture:** the
> gauge-invariant `χ_O` is strongly sub-linear, not a plateau — `χ_O=0.469±0.028` (L=12) and
> `0.587±0.055` (L=16), log-log slope `0.371±0.066`; free `L²` is excluded at `24.7σ`, but a pure
> plateau is also excluded at `5.6σ`. Saturation remains `[O]`. Per the V3 gate, this is **not a
> major breakthrough**; stop the broad W wave and reassess.

---

## §0 · Four verdicts

| Route | Verdict | Result |
|---|---|---|
| **W1** strong-coupling alignment | `[T]`/`[C]`/`[O]` | Reproduced `w₁(0)=s/2` (N≥3), `s/4` (SU(2) real fundamental). Exact one-plaquette normalizers: SU(2) `M0=1/4−β²/32+…`; SU(3) `M0=1/18+β/36+…` (the β-linear term is the `d_abc` cubic Haar moment `⟨(Re trU)³⟩=1/4`). Next coefficient exact numerically; its `c` part is P3's dressed amplitude `[R]` gap `[O]`. Padé/Borel of `M0` reach 2.5 but are vacuous (`M0>0`). Full series outside the cluster radius at N=2,3 (`u(2.5)=1.02, 3.55`). `w₁>0` at 2.5 survives only `[C]` on amplitude headroom `C·u<1` (empty at N=2,3) or `N≳760`. |
| **W2** large-N master field | `[T]`/`[O]`/`[NG]` | The master equation for `ψ₁` is the closed linearised Makeenko–Migdal / planar Bethe–Salpeter eigen-equation `ℳψ₁=m₁²ψ₁` around `W₀` (solvable in 2D; `W₀` unbuilt in 4D → `[O]`). `g₁/Σg=1/(1+R)` is the point mass of the normalised KL measure at the lightest glueball. Acid tests pass. It is **not** uniformly positive in the 't Hooft coupling: `g₁/Σg(λ=0)=0` (free field), rising to `s/2`; `1/N²` keeps `w₁∈[0,1]` but adds no floor; no moment/positivity functional bounds `g₁/Σg`. |
| **W3** exact small volume | `[T]`/`[NG]`/`[O]` | At β=2.5 every exactly solvable small volume gives `w₁=0` or a degenerate atom: 2D (`tr F²=H`, conserved), periodic Cartan transfer (`ψ₁` odd vs even source, exact to `8e−34`), β=0 (rank-one contact), free field (number parity), one-plaquette/4D L=1 (no `ψ₁`). `w₁` is **not** already positive at small volume. The sole positive anchor remains strong coupling `s/2` (L-uniform, binding not volume). Finite-size law: solvable cases `w₁=0`; delocalised point source `2/(L−1)`, decreasing. Full 4D L=2,3 is `[O]` (character series outside radius, `uK≫1`). |
| **WG** gauge-invariant estimator | `[T]`/`[V]` | The total-normalized zero mode of an APE-smeared `tr F²` density is exactly gauge-invariant (stable to `1e−15`; flat chart drifts 40–1000% per config), non-self-averaging (`χ_O=O(1)`, `Var(q̄)=χ_O/V↓`), acid-clean (free massive `1/m²`; free massless `∞`). Measured β=2.5: `χ_O=0.328,0.377,0.411` (L=4,6,8); ratios `1.15±0.09, 1.09±0.15` vs free `2.25,1.78` — a statistical **plateau**. `χ_P` drops because Polyakov is the wrong N-ality-1 channel; `χ_Landau` grows because a gauge-fixed coordinate accumulates the longitudinal/Gribov sector. Atom `f₁²` remains `[O]`. |

---

## §1 · What this means for V3

- **No theorem at β=2.5** for `W` or `I`; the W wave does not meet the bronze bar.
- **WG is a genuine gauge-invariant measurement, but not a plateau**: `χ_O` grows sub-linearly
  (slope `0.371±0.066`) from L=4 to L=16. The flat-chart plateau was not exact saturation; a finite
  `χ_∞` remains possible but unproved. It fixes the carrier caveat and weakens the plateau claim.
- **W1/W2/W3 give precise missing steps**: P3's dressed amplitude (`[R]`), the 4D master field `W₀`,
  and the cluster radius at N=2,3. None is a new inequality; all are the same nonperturbative input.

**Gate decision:** stop the broad W wave. The next step is not another 6-route sweep; it is a choice
between (a) pushing WG to L=12/16 with the fast Numba engine (gauge-invariant plateau, numerical), (b)
the P3 dressed-amplitude problem (`[R]`) as a standalone theorem target, or (c) the 4D master field
`W₀` (large-N, likely a multi-month program). **V3 still has no major breakthrough; Clay 0%, no QED.**

*Refs:* `W1-STRONG-ALIGN.md`, `W2-MASTER-FIELD.md`, `W3-SMALL-VOLUME.md`, `WG-INVARIANT-ESTIMATOR.md`,
`W-WAVE-CHARTER.md`; witnesses `w7/w1_*.py`, `w7/w2_*.py`, `w7/w3_*.py`, `w7/wg_*.py`.
