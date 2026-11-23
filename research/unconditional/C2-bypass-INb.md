# C2 — Bypassing IN-11 (IN-b): can the conditional theorem drop the channel selection?

> 2026-11-23 · Companion to `GAP-IN11.md`, `IN-11-attack.md`, `W4-cert-vs-gap.md`.
> **Verdict: IN-11 drops from the theorem's axiom list (12→11), but not for free — it survives as a supplier note of the defect-bootstrap (DB) route, where it is provably necessary.**

## 0. Decisive fact: IN-11 is not a hypothesis of the formal theorem

In `Handoff.lean` the certificate is a **scalar** function `cert : ℝ → ℝ` and the bridge is
`hcert : ∀μ, cert μ ≤ Δ/Λ` (`Handoff.lean:77`) — **channel-blind**. `handoff_gap` / `mass_gap_of_certificate` carry no channel index and never consume IN-11. IN-11 enters only as the **supplier obligation of DB**: DB bounds two channels, so `cert = min(cert_{0⁺⁺}, cert_{2⁺⁺})`, and the identity `min(Δ_{0⁺⁺},Δ_{2⁺⁺}) = Δ` *is* IN-11. The skeleton already omits it; IN-11 is needed only to instantiate `cert` from DB.

## 1. Three reformulations

**(a) Channel-blind certificate.** Needs `cert_blind ≤ inf_c Δ_c = Δ` with `cert_blind>0` — a *lower* bound on `inf spec(H|_{1^⊥})` with no channel selection, i.e. **gap-strength with zero slack**. Constructible? The known channel-blind mechanisms (Osterwalder–Seiler strong coupling, Balaban RG, SLC⇒Poincaré) are exactly the IN-3 wall: rigorous only below the crossover. Not constructible today; no easier than IN-3, and *strictly harder* than DB's two-channel certificate (which at least has a machine).

**(b) Total density of states.** `ν_tot := ⊕_c ν_c` is the KL measure of a channel-summed local operator — distinct `J^{PC}` sectors are orthogonal, so cross terms vanish and `ν_tot ≥ 0`, with `inf supp ν_tot = Δ` (Reeh–Schlieder, complete channel set). But moments `∫s^k dν_tot` give *upper* bounds on the edge, `Δ ≤ (M_k/M_0)^{1/k}` — **wrong direction** (moment lemma). Only edge control `ν_tot([0,δ))=0` gives `Δ ≥ δ` = the gap verbatim. As a bootstrap it needs `0⁻⁺`/odd-C/high-spin blocks, which `GAP-IN11` shows do not exist (`q=∂·K` unprotected). So (b) is **gap-equivalent at best, strictly harder** than the two-channel bootstrap.

**(c) Count fold.** `#{states with 0<m<δ} ≤ K` is insufficient for any `K≥1`: place one unseen state at `ε→0`; `Δ=ε` is arbitrarily small while the count stays ≤ K. Only `K=0` transfers the bound, and `K=0` *is* `(b′)` = IN-11. Separation model kills every `K≥1`; no quantitative `δ/K` gap exists (the lightest of K states is unbounded below).

## 2. Proposition (relative necessity — proven)

Fix the supplier to a certificate valid only in `C ⊊ {channels}`. Then any hypothesis set concluding `Δ ≥ min(δ_W,δ_S)` with `cert>0` must imply `inf_{c∉C} Δ_c ≥ min_{c∈C} Δ_c ≥ cert`.

*Proof.* Otherwise a model has some `c₀∉C` with `Δ_{c₀}=ε < min_{c∈C} Δ_c`; then `Δ=ε` while every `C`-certificate hypothesis holds; pick `ε ≤ min(δ_W,δ_S)/2` to falsify the conclusion. The weakest monotone repair is `Δ_{c₀} ≥ cert`, i.e. IN-11 weak form `(b′)`. ∎

## 3. Verdict

- **The theorem:** a channel-blind reformulation **exists and is no harder** — state `cert` channel-blind (the Lean already does); **IN-11 drops** from the axiom list.
- **The cost:** the only channel-blind substitutes are gap-strength (a)/(b) or insufficient (c). Dropping IN-11 from the *statement* moves its content into "supply a channel-blind `cert`," which is ≥ gap.
- **The DB route:** by §2 it **strictly needs IN-11**; it cannot be removed or count-weakened there.
- **Honest resolution:** demote IN-11 from an *input* to a *supplier note* of the DB route. Do **not** claim IN-11 is necessary for *all* gap proofs — false; a channel-blind route evades it by *being* the gap.

*No fabricated proof: §2 is the only new theorem; (a)–(c) are separation-model obstructions, not claims that the gap is impossible.*
