# REDUCTION-AUDIT — adversarial audit of the central "⟺" chain: two definitions, one conjecture, one false biconditional

> 2026-11-23 · Adversarial audit. Claim under test:
> **Clay gap ⟺ `f₁²=ν({m₁²})>0` is a KL atom ⟺ `m_eff²(k→0)≥λ` ⟺ `Var_μ(a_{k→0})≤1/λ` (H1).**
> **判词 (verdict, up front):** the chain is **not** a chain of equivalences. Read against the corpus's *own* definitions it is **(iv) FALSE → (iii) conjecture → (ii) definition**, plus one genuine **(i)** one-directional theorem binding them. The load-bearing equivalence — gap ⟺ H1 — holds in **one direction only**. The chain manufactures necessity as sufficiency and renames the **inverse susceptibility** as "mass." It is information-losing, not wrong in direction, and the corpus's own files already concede every gap below.

## Link-by-link

**(A) Clay gap ⟺ `f₁²>0` atom — verdict (iv) FALSE as ⟺, (i) one way.**
`Clay gap` (the Lean statement, `lean/YMStatement/Basic.lean`: uniform exponential clustering of the generating functional at rate `m>0`, i.e. `spec(H)∩(0,m)=∅`) does **not** imply any particular operator has an atom. The corpus proves this itself: `G1-UNIQUENESS-VACUUM.md` — "the generalized free field satisfies all OS axioms with `Δ>0` yet **no atom** (Theorem B′)"; `PROOF-F1SQ.md` — the massive free scalar `m>0` with `O=:φ²:` has `f₁²=0` (its KL measure is a cut on `[4m²,∞)`). So **gapped theory with `f₁²=0` exists** ⇒ the forward `⟹` is false. The reverse (isolated atom at the bottom of a channel ⇒ that channel clusters at rate `m₁>0`) is a theorem, but becomes the *mass* gap only if that channel carries the lightest state overall — an unproved ordering input. Correct statement: **atom (isolated, at `inf supp ν`) ⇒ gap in that channel**; `gap ⇏ atom`. The written `⟺` is false.

**(B) atom ⟺ `m_eff²(k→0)≥λ` — verdict (iii) CONJECTURE, not theorem.**
The atom is a weight of the composite `tr F²` KL measure; `m_eff²(k→0)=G″=V_eff″` is the curvature of the *coarse zero-mode* effective potential of the fundamental connection. `G2-ASSEMBLE-100.md §7` says it plainly: the unification `H*=H1` is "**a structural identification** … One closure would likely close both faces — but that is **a conjecture, not a theorem**." The rigorous bridge (`GAP-G2-edge` residue–level-gap lemma) only trades `f₁²` against a *second-level threshold* `Λ` given moments; it never produces the curvature `G″`. The `⟺` is a gloss over an identification.

**(C) `m_eff²(k→0)≥λ` ⟺ `Var_μ(a)≤1/λ` — verdict (ii) DEFINITION (hence tautology), or (iv) FALSE if `m_eff` is read as the lightest mass.**
The corpus defines `m_eff²:=Γ_T(0)=1/χ=1/Var(a)` (`D1-BLOCK-MARTINGALE.md`, `D3-SUSCEPTIBILITY-VARIATIONAL.md`). With that definition both sides are literally the same inequality — a relabeling, zero content. In the other reading (`m_eff²=G″=V_eff″`), only `G″≥λ ⟹ Var(a)≤1/λ` is a theorem (Brascamp–Lieb 1-D, one direction); the converse and the identity `G″·Var(a)=1` hold **only for a Gaussian marginal** (`ABC-NEW-INEQ.md`; measured "2–5%," `D3`). So the `⟺` is either a tautology or false.

## The three targeted questions

**(a) gap ⟺ atom?** No. As above: gapped generalized free field and gapped free `:φ²:` both have `f₁²=0`. Necessity of the atom fails; only sufficiency (atom+isolation ⇒ gap) is a theorem.

**(b) mass or susceptibility?** **Susceptibility.** `a` is the zero-momentum mode, so `Var(a)=χ=∫s^{-1}dν` (zero-momentum connected propagator). The lightest mass is `inf supp ν`; the two are different functionals. Counterexample on `[0,∞)`: `ν=wδ_ε+(1-w)δ_S` with `w/ε≤1/λ` gives `χ≤1/λ` while `inf supp ν=ε→0`. So `Var(a)≤1/λ` is a **screening/susceptibility** bound — a residue-weighted *harmonic* mean of the masses, dominated by the lightest residue, **not** the lightest mass. It equals the pole mass only in the single-pole (Gaussian) case. (Also `G″` vs `1/χ` is only 2–5% even numerically, `FORGE5-SL-IR-SINGLE.md`.)

**(c) tautology that loses content?** Not a pure tautology — it is non-vacuous (`Var(a)=∞` on the free massless field, so it fails the acid test), and `gap ⇒ Var(a)≤1/λ`(×normalization `M₀`) is a real theorem via `χ=∫s^{-1}dν≤M₀/λ`. But it **loses the content**: it keeps a susceptibility bound and discards the hard lightest-mass/edge datum. Getting back to a mass gap needs the residue concentration `I1` (`f₁²` a controlled fraction of `M₀`) *and* isolation `I2` — precisely `H*`, which the corpus itself says is **strictly stronger** than bare `m>0` (`RESIDUE-CONDITIONAL.md`, `G2-RESIDUE-ATOM.md`: "atom strictly stronger than gap").

## What actually survives (the honest skeleton)

- **(i) Theorem:** isolated atom at `inf supp ν` of a channel ⇒ exponential clustering in that channel (`Theorem A`).
- **(i) Theorem:** `gap ⇒ χ≤M₀/λ` (i.e. `Var(a)` bounded, up to normalization).
- **(i) Theorem:** `G″≥λ ⇒ Var(a)≤1/λ` (Brascamp–Lieb, one way).
- **(iii) Conjecture:** `H*=H1` "one object."
- **(iv) False:** `gap ⟺ atom`; `H1 ⟺ gap`.

**Bottom line.** Written as `⟺`, the chain is unsound. Its true logical shape is a **diamond of implications converging on the necessary condition `Var(a)≤1/λ` (H1 = susceptibility)**, with the *sufficient* condition being the strictly stronger `H* = I1∧I2` (residue + isolation). Proving H1 alone does **not** yield the Clay gap (tiny-residue counterexample); proving H* does, but H* is not H1 — the unification is the conjecture, not the result. The reduction to `Var(a)≤1/λ` preserves a genuine inequality but **relocates, rather than removes, the residue content**, and silently equates inverse susceptibility with lightest mass. No circularity at the level of the *identities* (Ward, BL, residue–level-gap are real); the circularity is logical: the chain's final `⟺` presumes the very residue/ordering input it claims to have reduced away.

*Refs.* `PROOF-SKELETON.md` §1; `G2-ASSEMBLE-100.md` §1–2,§7; `GAP-G2-edge.md`; `PROOF-EDGE-E2.md`; `PROOF-ISOLATION.md` §1,§4; `CONTINUUM-CP.md`; `RESIDUE-CONDITIONAL.md`; `G2-RESIDUE-ATOM.md` §9; `PROOF-F1SQ.md`; `G1-UNIQUENESS-VACUUM.md`; `D1-BLOCK-MARTINGALE.md`; `D3-SUSCEPTIBILITY-VARIATIONAL.md`; `FORGE5-SL-IR-SINGLE.md`; `ABC-NEW-INEQ.md` §1; `lean/YMStatement/Basic.lean` (`MassGapStatement`).
