# G2-CORRELATION-GEOMETRY — geometric near-orthogonality is the correct *diagnosis* of why `X≈0`, not a *tool*: it is free-field-blind, and every route to a usable "`|X|` small" constant is gauge-only, Ward-forbidden, or circular ⇒ `|X|<W_geo` is still `VarG`=H1

> **判词 (verdict).** The mechanism is fixed: `X=−2Cov(G,Φ′|0)`; measured `corr(a,∇S)≈0` (force nearly independent of the coarse amplitude) but `corr(a,Φ′)≈+0.85` (amplitude tracks the Jacobian), so `X≈−0.04` and `G″≈W_geo`. *Proving* near-orthogonality cannot become the bound. **(a)** Elitzur kills only the **longitudinal** force (`∇_{Dω}S=0` exactly); the sign-carrier is the **transverse** Savvidy force `∝|F|`, untouched — and once `a,G` are both gauge-**in**variant, Elitzur is silent. **(b)** "many plaquettes vs one mode" decoupling is **Ward-forbidden**: `∇_VS` is the conjugate force, `G″·Var(a)=1`, so decoupling it makes `Var(a)=∞` — the free field. **(c)/(d)** transfer-matrix/Schmidt constants *are* the spectral gap. And the acid test decides: `corr(a,∇S)=0` is the **massless free-field** value with `G″=0`; YM is the inverse regime (`Var∇S≈2` large, `corr≈0`), and exploiting it needs `Var∇S` bounded = `VarG` = H1. **Near-orthogonality does not give the bound; `|X|<W_geo` is still `VarG`=H1.**

## §1 The target
`G″=Δ_dyn+W_geo+X`, `Δ_dyn=E[H]−VarG`, `W_geo≥N/12` proven. A **lower** bound needs `X≥−ε` with `ε<W_geo−Δ_dyn` (a constant) **or** `Δ_dyn` controlled below. Near-orthogonality addresses only `|X|`.

## §2 (a) Elitzur — real but the wrong half
Gauge invariance ⇒ zero longitudinal force (a true theorem, part of `corr≈0`). But the negative, sign-carrying mode is **transverse magnetic** `∝|F|` (Savvidy, W4), carries the action, and gauge invariance never zeroes it. The physical covariance is between two gauge-**invariant** objects → Elitzur is silent.

## §3 (b) many-vs-single — contradicted by Ward
`a_{k→0}` is the **soft** mode (`Var(a)~O(1)`, volume-stable), not CLT-suppressed `1/√V`. `G″·Var(a)=1` makes `∇_VS` conjugate to `a`; decoupling them ⇒ `Var(a)=∞` ⇒ no gap — so (b) predicts the free field and proves too much. Bounding the true response = H1.

## §4 (c)/(d) transfer matrix & Schmidt — circular
`Cov(a,G)=Σ_{n≥1}⟨0|a|n⟩⟨n|G|0⟩` needs `Δ_n` (the gap); PF gives the vacuum but `a` lives in `Ω^⊥`, the RP-cone-degenerate face (PROOF-GAUSS-INVARIANT-PF). Schmidt of a local Gibbs measure: the entanglement spectrum's low edge is the physical spectrum (Li–Haldane). Both routes *are* the gap.

## §5 Acid test: `corr=0` is free-field-satisfied
Massless free field at `k→0`: `m=0` ⇒ `∇S` has no soft part ⇒ `Var∇S=0` ⇒ `corr=0` trivially, `G″=0`, `f₁²=0`. YM at β=5.7: `Var∇S≈2.0` **large** yet `corr≈0` — a **nonzero-variance cancellation**, the opposite limit. Using the YM orthogonality *requires* showing `Var∇S` large-but-orthogonal ⇒ controlling `Var∇S`=`VarG`=H1. A tool built on `corr≈0` is sign-blind and free-field-vacuous — it does **not** fail the acid test.

## Five-point verdict
1. **Nonperturbative — PASS** (gauge/Ward/Prékopa exact).
2. **Sign-controlling LOWER bound — FAIL**: `|X|` small is magnitude; `Δ_dyn` uncontrolled; per-mode `X∈[−0.09,+0.16]`, no uniform sign.
3. **Non-circular — FAIL**: quantifying "small" = `VarG` (C-S=FH) or T/Schmidt spectrum (gap) = H1.
4. **Non-SUSY + fails free field — FAIL (decisive)**: `corr=0` is the massless-free value with `G″=0`.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL**: no `λ`.

## Exact failing estimate
```
Ward: G″·Var(a)=1.  G″ = Δ_dyn + W_geo + X,  Δ_dyn=E[H]−VarG, X=−2Cov(G,Φ′).
corr≈0 ⇒ |X| small, but G″ ≥ W_geo + Δ_dyn − |X|, Δ_dyn uncontrolled below.
β=2.5: VarG≈2, VarΦ′≈0.02, W_geo≈0.32; C-S |X|≤0.40>W_geo ⇒ G″≥−0.08.
Need |X|≤W_geo−λ (⇔ VarG≤W_geo²/(4VarΦ′)≈0.8–1.3, vs ≈2) OR sign X≥0 (absent).
Free: ∇S=0 ⇒ Var∇S=0 ⇒ corr=0, G″=0 — acid-test-satisfied, sign-blind.
```

**Sharpest line.** The near-orthogonality is a diagnosis, not a certificate: the amplitude is slaved to the flat-Haar Jacobian, not the Wilson force. Elitzur kills only the action-free longitudinal force; decoupling the force from `a` is Ward-forbidden (would give `Var(a)=∞`, the free field); T/Schmidt constants are the gap. Decisively, `corr(a,∇S)=0` is the massless free-field value *with* `G″=0` — the YM case is the inverse regime (`Var∇S≈2` large yet `corr≈0`), and exploiting it needs `Var∇S` bounded = `VarG` = H1. **Near-orthogonality does not give the bound: `|X|<W_geo` is still `VarG`=H1.**

*Refs.* C4-COV-LASTMILE; B4-VARIANCE-MAG (`VarG≈2` L-stable); C5-NEG-ASSOC / C6-DEPENDENCE-CLASS (tight `|X|`=H1); B3-TYPICAL-COERCIVITY (`ρ≈−0.55`, C-S short 1.8×); COMPUTE-GPP-SU3 (`E[H]≈Var∇S≈2.0`, `W_geo≈0.32`, β=5.7 `Cov` sign-varying); PROOF-JACOBIAN-SAVES (identity, `W_geo≥N/12`, Ward); PROOF-GAUSS-INVARIANT-PF (`Ω^⊥` edge); PROOF-SAVVIDY-WGEO-DOMINATION (transverse negative mode); SUB1-COV-BOUND (Stein=`Var(a)`=H1); Ward/Slavnov–Taylor; Elitzur 1975; Prékopa 1971; Fréchet–Hoeffding; Li–Haldane 2008; Savvidy PLB 71 (1977) 133.
