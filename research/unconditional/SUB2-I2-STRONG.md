# SUB2-I2-STRONG — can the cost grading close I2? No: it grades *amplitudes*, not *eigenvalues*; the failure is a `t`-uniform `O(u)`

> 2026-11-23 · Adjudicates the documented contradiction: PROOF-P1P2 §3 asserts the integer cost grading yields I2 (`λ₂/λ₁→0`, `δ=(p'−p)ln(1/u)/a`); PROOF-P3-UNCONDITIONAL (corrected to residue-only) and REDTEAM-FORGE retract it. **Verdict: the retraction is correct; PROOF-P1P2 §3 is the error.** The cost grading closes **I1** (residue `f₁²≥r₀>0`) — the single tube is the unique *leading weight*, its `f⊗f̄→1` recoupling is `+`, KP is `V,L`-uniform. It **cannot** close **I2** (isolation, `spec(H)∩(m₁,(1+δ)m₁)=∅`) because it grades *polymer amplitudes* (an integer plaquette cost `|P|`, weight `≤(Cu)^{|P|}`) while I2 is a datum about *transfer-matrix eigenvalues* (`λ_n=e^{−am_n}`, continuous rates, NOT sector costs); the bridge is the second-level gap, and the KP/OS remainder is a `t`-uniform `O(u)`, never `O(e^{−δt})`. Regime `β<1/48`, fixed `a`, uniform in `V,L`.

## What the grading does deliver — I1 (+ classical gap)
For `β<β₁=min(1/48,c/N)` the KP-graded connected `tr F²` expansion converges uniformly in `V,L`. Every polymer linking slices `0,t` has cost `|P|≥t/a`; equality forces a single fundamental tube, weight `a₀u^{4t/a}`, `a₀>0`. Hence `e^{m₁t}C_L(t)=f₁²+O(u)`, `m₁=4 ln(1/u)/a` 〔factor-4 corr. 2026-11-24，见 P3-AUDIT-CORRECTED.md〕, `f₁²≥r₀(β)>0`, uniform in `L,t`: I1 and the classical OS strong-coupling gap `m_gap=m₁` 〔结论不变：因子 4 为公共正因子〕.

## The theorem that *would* close the target (if I2 held)
Theorem A (`mass-spectrum/STATUS`): `I1∧I2 ⟺` sharp decay. With I2 at `δ=(p'−p)ln(1/u)/a→∞`, uniformly in `L`: `ν_O=f₁²δ_{m₁²}+ν_rest`, `supp ν_rest⊂[(m₁+δ)²,∞)`, `m₁=Δ` — the regime-limited **gap+residue** theorem. **It does not follow**, and not for a weak constant.

## The exact failing estimate
I2 needs a two-exponential bound
`|e^{m₁t}C_L(t)−f₁²| ≤ B e^{−δt}`, some `δ>0`, uniform in `L,t`.
The grading supplies only (PROOF-P3 §4) `|e^{m₁t}C_L(t)−f₁²| ≤ Cu` (t-uniform, `C=O(1)`). Two faces of one failure:
- **(precision)** `Cu` is a *constant*; `Be^{−δt}` decays. The `O(u)` is the *un-resummed self-energy/residue decoration* of the leading tube (a `t`-independent multiplicative `1+O(u)`), sitting at the **same order in `u`** as the excited-sector amplitude. A `t`-uniform `O(u)` bound is satisfied by *any* second level with `λ₂/λ₁∈[1−Cu,1)` — `δ` arbitrarily small, even a cut opening just above `m₁`. It yields **no positive `δ`**.
- **(placement)** integer grading fixes the excited sector only up to an *integer* multiple of `c₁`. The physical `0⁺⁺*` at `≈1.5m₁`, strictly **below** the `2m₁` two-glueball continuum (operative `Λ`, PROOF-ISOLATION §4.4), is *sub-leading*; the leading integer cost (a boxy `2×1` tube) sits at `2m₁` — on/above the continuum, with **no isolation gap**. The `1.5m₁` binding is exactly the sub-leading content the grading cannot resolve.

## Why `λ₂=c₂u^{p'}` (P1P2 §3) is the error
Writing `λ₁=c₁u^p+O(u^{p+1})`, `λ₂=c₂u^{p'}` treats an *asymptotic* series as if it had a uniform remainder. But the `O(u^{p+1})` self-energy correction to `λ₁` is the *same order* as `λ₂` when `p'=p+1`, so `λ₂` is never resolved as a distinct eigenvalue — the "two-gap" statement is the one-gap (`m₁`) statement renamed. This is PROOF-EDGE-E2's bottom-functional blindness: every controlled datum sees `λ₁`; none pins `λ₂`. Lemma SC (coarse Poincaré, SZZ `β<1/48`) does not help — a bulk-bottom inequality fixing `m₁`, silent on `λ₂`.

## 5-point verdict (target = I2)
1. **Nonperturbative — PASS for I1, FAIL for I2.** The `λ_n=c_nu^{p_n}` identification is a partial strong-coupling series with no uniform remainder.
2. **Sign-controlling — FAIL.** Leading recoupling is `+` (gives I1); the `O(u)` decoration/excited remainder is sign-indefinite (SU(N) character alternation + self-energy).
3. **Non-circular — FAIL at closing.** Asserting `λ₂` at `u^{p'}` presupposes the eigenvalue structure that *is* I2.
4. **Non-SUSY + fails free field — PASS.** Massive free scalar binds no tube, `f₁²=0`; strong coupling binds, `a₀>0`.
5. **Uniform `m_eff²(k→0)≥λ` — gap passes, isolation fails.** `m_gap²=c₁²ln²(1/u)/a²` is uniform (classical OS); the *isolation* `δ>0` is **not** — the missing estimate.

## Sharpest line
Cost grading, I1 and I2 all "look like" whole powers of `u`; only I1 is. I1 is a statement about a *weight* (a leading amplitude, positive at the recoupling), which the integer grading makes unique and positive, `V,L`-uniform. I2 is a statement about a *rate* (`λ₂`), and a weight-grading whose remainder is a `t`-uniform `O(u)` cannot distinguish a second atom at `(p'−p)ln(1/u)/a` from one at `ε→0`, nor place the `0⁺⁺*` below the `2m₁` continuum. The "costs" coincide only in the formal series; rigorously `p'−p` is a **weight** gap, not yet a **spectral** gap. Closing I2 means converting one into the other — 4D non-abelian bound-state isolation (no free-field anchor, no termwise sign), which no cluster/LSI/Poincaré method supplies. **I2 stays open; the exact failing estimate is `Cu` vs `Be^{−δt}`.**

*Refs.* PROOF-P1P2 §3 (retracted claim); PROOF-P3-UNCONDITIONAL §4 (the `O(u)` remainder, its own correction); PROOF-ISOLATION §3–4.4 (`1.5m₁`=operative `Λ`); PROOF-EDGE-E2 (bottom-functional blindness, exact two-mass estimate); PROOF-IN11 (`a_{0⁺⁺}=a_{2⁺⁺}=1`); PROOF-COARSE-POINCARE (Lemma SC silent on `λ₂`); REDTEAM-FORGE; mass-spectrum/STATUS (Theorem A, I1/I2); PROOF-F1SQ (gapped free scalar, `f₁²=0`); Osterwalder–Seiler 1978; Shen–Zhu–Zhu CMP 400 (2023); Glimm–Jaffe–Spencer 1974; Kotecký–Preiss 1986.
