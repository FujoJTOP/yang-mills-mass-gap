# D-BBS-INTEGRATE-NEG — explicitly integrating out the Savvidy negative modes to leave a positive-complement effective theory: exact ⇒ vacuous, one-loop ⇒ complex; the surviving *fast* Hessian is provably positive (a real lemma), but it never reaches the *slow* `G″`

> **判词 (verdict).** The mechanism has a genuine, bankable piece and a fatal one, and they sit on opposite sides of one identity. **(Lemma NS — new, true):** eliminating a *negative-definite* block from a quadratic form strictly *raises* the surviving block's curvature, `Σ=A−BC⁻¹Bᵀ⪰A` (verified), so on the exact positive eigenspace the reduced fast Hessian is `⪰λ_+(F)=+g|B|>0`. This is the positive complement, and it re-establishes BBS's *fast* one-step contraction after the negatives are removed — a refinement over the corpus line "`H_ff` indefinite ⇒ Schur undefined." **(Fatal):** the reduced theory is only *clean* under the **Gaussian/one-loop** truncation, and there the elimination of the negative block yields `i^{N_neg}(2π)^{N_neg/2}|detH_-|^{−1/2}` — the **Savvidy determinant**, whose imaginary part is the physical Nielsen–Olesen instability. The exact (non-Gaussian, Fubini) integration *changes no physics* but confers no log-concavity (Prékopa needs joint log-concavity, which fails). Either way the Schur gain lives entirely in `E[PᵀHP]`; the H1 wall is the **variance** term `−Var(Pᵀ∇S)`, which the block elimination never touches. The reduced theory admits the trap iff `Σ_fast−Var≥λ ⟺ H1`. Does not close. Exact residual = H1.

## (a) The negative-mode subspace, per shell
Second variation on transverse 1-forms is `H_j=−D²+2adF_{B_j}`; its lowest Landau rung (`n=0`, spin anti-aligned) is `λ_min=−g_j|B_j|` with degeneracy `g_j|B_j|A_{12}/2π`. Per shell the negative subspace `V_-^{(j)}` is the lowest-Landau-level, transverse-**magnetic** set, rank `N_neg^{(j)}=c_N(g_j|B_j|)²V_j/8π²` (`c_{SU2}=1`, `c_{SU3}=3`), depth `−g_j|B_j|`. These are **orbit-transverse** (PROOF-PHYSICAL-SUBSPACE-CONVEXITY): physical, in `𝒜/𝒢`, not pure-gauge — so no BRST choice deletes them. `N_neg` *shrinks* in the UV (`g_j→0`) but is *large* at β=2.5 (`g≈1.5`).

## (b) Integrate out — exact vs one-loop, and the Schur lemma
Split `H=H_+⊕H_-`, `H_-≺0`.
- **Exact (Fubini):** `∫_{V_-}e^{−S}dφ_-` over a finite-dim subspace. Physics-preserving and well-defined, but `S_eff(φ_+)=−log∫e^{−S}dφ_-` is **not** log-concave — Prékopa requires the *joint* `S_W` log-concave, which the Savvidy Hessian violates (PROOF-COARSE-POINCARE (a)(iv)). Exact elimination buys nothing.
- **One-loop (Gaussian):** each `λ_k<0` gives `∫e^{−½λ_kx²}dx=√(2π/|λ_k|)·i`, so `∫_{V_-}e^{−½φH_-φ}=i^{N_neg}(2π)^{N_neg/2}|detH_-|^{−1/2}` — Re real, Im `=πN_neg/2`: **the Savvidy determinant, complex.** Its modulus is computable; its phase is the instability.
- **Lemma NS (bankable):** in the quadratic form `[[A,B],[Bᵀ,C]]`, `C≺0 ⇒ Σ=A−BC⁻¹Bᵀ=A+B|C|⁻¹Bᵀ⪰A`. On the exact `+`-eigenspace (`B=0`) `Σ=A⪰λ_+`; with residual coupling it is strictly larger (numerically `Σ−A⪰0`, min `+0.16`). **Eliminating a negative block only helps the survivor.**

## (c) BBS trap on the reduced theory; is physics changed?
- **Fast half: yes.** `Σ⪰λ_+(F)=+g|B|>0`, so the reduced fast Hessian is positive and the one-step contraction `e^{−cλ_+}` survives after removing the negatives. Genuine, not in the corpus.
- **Slow half: no.** `G″=E[PᵀHP]−Var(Pᵀ∇S)+W_geo+X`. Lemma NS improves only `E[PᵀHP]`; the binding term `−Var` (H1) is entropic, orthogonal to the block-Schur. So trap ⟺ `Σ_fast−Var≥λ ⟺ G″≥λ ⟺ H1`.
- **Physics:** exact Fubini leaves `Z` invariant (no change). One-loop reduction *does* change it: real-indefinite → **complex-positive** (metastable). The Im is physical; "covariant elimination" cannot delete it, because the negatives carry no gauge nullspace to quotient. The reduced theory is a *quasimode*, not a genuine probability measure.

## (d) Does the determinant substitute for φ⁴ convexity? No.
φ⁴'s large field is *more* convex (`+12λφ²`); BBS's LF weight *decreases*. Here (i) `|detH_-|^{−1/2}` **grows** with `N_neg∝(gB)²V` — enhancement, not damping; (ii) `log|detH_-|=½TrlogH` is exactly the Savvidy `(11/48π)(gB)²log` term, whose background-Hessian **is** the indefinite object being removed — a re-expression, not a convexification; (iii) the phase `i^{N_neg}` kills positivity, so log-concavity/PF/Brascamp–Lieb have no carrier globally (only inside the metastable well).

## Exact residual
```
Lemma NS: Σ_fast ⪰ λ_+(F)=+g|B| > 0     (true, fast mode)
missing : Var(Pᵀ∇S|0) ≤ Σ_fast − λ  ⟺  G″ ≥ λ  ⟺  H1   (slow mode, unchanged)
one-loop : Z_eff complex, Im = (π/2)N_neg  (Savvidy, physical)
```

**Sharpest line.** Integrating out the negatives is *exact and vacuous* (Fubini, no convexity) or *one-loop and complex* (Savvidy determinant) — never a *real, log-concave* reduced theory. The genuine gain, Lemma NS, is real and new, but it certifies only the **fast** half of BBS; the **slow** wall `−Var` is entropic and block-orthogonal, so the trap lands back on H1, and the determinant itself is a complex, field-growing prefactor that neither damps large fields as φ⁴'s convexity does nor carries a positive measure. Since the negatives are orbit-transverse, no gauge-covariant choice removes them; their imaginary part is the physical instability, and removing it is the gap.

*Refs.* PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, `N_neg=c_N(gB)²V/8π²`); PROOF-BORN-OPPENHEIMER §2 (one-loop/Schur complex, `H_ff` indefinite); PROOF-COARSE-POINCARE (a)(iii)–(iv), H1; PROOF-PHYSICAL-SUBSPACE-CONVEXITY (negatives orbit-transverse, 7.1% after Q⊥); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`); C-BBS-MARGINALLOG §3–5 (poly LF convexity breaks; same lock); FORGE5-SL-UV-CONE-EMPTY (β=2.5 cone high-momentum, `osc=Θ(N_□)`); Savvidy, PLB **71** (1977) 133; Nielsen–Olesen, NPB **144** (1978) 376; Leutwyler, NPB **179** (1981) 129; Bauerschmidt–Brydges–Slade, LNM **2242** (2019).
