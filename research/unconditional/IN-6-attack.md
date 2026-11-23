# IN-6 — IR limit with a unique vacuum on ℝ⁴

**Verdict.** IN-6 bundles two claims of opposite status. **Vacuum uniqueness is standard** — a consequence of OS4/ergodicity. **Uniqueness of the limit** (all boundary conditions and subsequences reach one Gibbs state) is **open**, and is *not* the gap. The exact open content is per-sector boundary screening `W1_0 ∧ W1_1`.

**(1) Precise statement.** `{μ_L}` on `Λ_L↑ℤ⁴` (fixed-volume `a↘0` first, then `Λ→∞`) has a unique weak limit `μ` satisfying OS0–OS4. Split:
- **VAC** (each limit point has one vacuum): translation-invariant subspace of `L²(μ)` = `ℂ·1`.
- **UNI-lim** (one limit): every bc/subsequence converges to the same `μ`.

**(2) Status / evidence.** **VAC: standard.** OS4(=MIX) ⇒ ergodic (mean ergodic theorem) ⇒ trivial invariant subspace (REPORT L4). No gap needed: critical 2D Ising is mixing with a unique vacuum, power-law decay, no gap (Example 2). So vacuum uniqueness is **independent of, and strictly weaker than, the gap**. **UNI-lim: open.** MG on one state does not bound it — low-T Ising pure states cluster yet coexist (Example 1). It needs `W1_0` (neutral) ∧ `W1_1` (string) screening; only `W1_0` is served by the gap.

**(3) Advance / missing step.** *Proved piece:* `VAC ⇐ OS4` (L4: Choquet simplex, `μ` extreme ⇔ ergodic ⇔ trivial invariant subspace). An implication, not new input; OS4 is weaker than the gap (Ex. 2). *Missing step, pinned:* `UNI-lim ⟺ W1_0 ∧ W1_1` (REPORT L2: DLR + boundary screening ⇒ unique Gibbs state). `TIGHT ∧ W1_0 ∧ W1_1 ⇒ UNI` is proven (Thm A); `W1` is unproven for 4D YM. **Transfer-matrix Perron–Frobenius cannot bridge it:** on `T⁴_L`, Krein–Rutman yields a unique positive ground eigenspace *only if* the top eigenvalue is isolated — that isolation **is** the gap (circular). And finite-volume PF is not bc-uniform, so distinct bc carry distinct ground states into the limit (Ising). PF gives VAC-from-gap, not UNI-lim.

**(4) Flip condition.** IN-6 (UNI-lim) **fails** iff two distinct limit points exist — two bc/subsequences landing in different Gibbs states: e.g. distinct 't Hooft sectors, or center-unbroken phases sharing `𝒜_0` but with different string tension (`W1_1` breaks).

*Refs:* `docs/conditional-theorem.md` §IN-6; `research/uniqueness/{REPORT.md §5,L4,L2, §5.5,R32.md}`.
