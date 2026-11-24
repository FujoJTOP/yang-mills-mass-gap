# Transmutation: does σ>0 ⇒ m_{0⁺⁺} ≥ c√σ, non-circularly?

*Independent first-principles analysis. No other project files consulted.*

## The two objects are in different sectors

- **String tension σ**: defined by the asymptotic decay of the fundamental Wilson loop, ⟨W(C_{R,T})⟩ ∼ exp(−σRT) as R,T→∞. It is a property of a **1-form** loop operator, living in the **charged** (static-fundamental-source) sector, and stated **asymptotically** in large loop size.
- **Glueball mass m_{0⁺⁺}**: the gap above the vacuum in the **neutral, gauge-invariant, local-operator (0-form)** sector: the least E_n with ⟨0|O_{0⁺⁺}(0)|n⟩ ≠ 0 for a local scalar operator.

Different sector (charged vs neutral), different operator (1-form line/loop vs 0-form local), different scale (asymptotic large-R vs short distance). Any derivation must bridge all three.

## The direction-of-inequality problem (the core)

The area law is an **asymptotic upper bound** on the loop correlator: ⟨W⟩ ≤ C e^{−σA}. By the transfer-matrix representation ⟨W(C_{R,T})⟩ = Σ_n w_n e^{−E_n T}, a bound on this decay gives a **lower bound on the lightest energy coupled to W**, namely E₁ ≥ σR + o(R). But that E₁ is the **flux-tube/string state**, in the charged sector, with energy ∝ R.

Now m_{0⁺⁺} ≥ c√σ is a **lower bound on a neutral-sector energy**. To derive a lower bound on a mass one needs an **upper bound on the 0⁺⁺ vacuum correlator** at some distance (decay cannot outrun e^{−mx}) — equivalently a **lower bound on the 0⁺⁺ spectral weight**. The area law supplies neither: it controls a different correlator, in a different sector, and only asymptotically. So the implication is not merely unproven; **the available inequality points the wrong way**. This — wrong-direction inequality plus sector-and-scale mismatch — is the precise obstruction. It is *not* circularity *per se*.

## (a) Flux-tube modes → a glueball of mass c√σ?

Heuristic: an open string of length R has transverse modes k = πn/R, E ≈ σR + π(n−…)/R; a *closed* string of size R ∼ 1/√σ gives m ∼ √σ, its ground state being scalar → 0⁺⁺. This produces an **estimate** m ≃ c√σ (lattice c ≈ 3–4), **not an inequality**. The closed-string ground state ↔ local 0⁺⁺ operator identification is assumed, not proved; and at R ∼ 1/√σ the effective-string description is outside its controlled regime (R comparable to the string thickness) — exactly where control is needed. The lightest closed string is the stable ground state and cannot decay, so its existence yields no bound; only the *ratio* is numerically observed.

## (b) Sector bridge: does one exist?

Honestly: **no rigorous operator-to-state bridge exists here.** The Wilson loop naturally creates charged/string states; the 0⁺⁺ glueball is a neutral closed-string excitation. Identifying the closed-string ground state with a specific local 0⁺⁺ operator — with controlled overlap, hence a spectral-weight inequality — *is* the missing step. In lattice language it would be a reflection-positive bound of the form "area-law decay in the charged sector ⇒ a lower bound on |⟨0|O_{0⁺⁺}|lightest⟩|²," i.e. a **no-hiding theorem**. None is known.

## (c) String breaking

In pure SU(N) with no dynamical fundamentals, the fundamental string between static sources **does not break**; the area law is genuinely asymptotic. Closed strings can decay into glueball pairs only above threshold E > 2m_g, and the lightest 0⁺⁺ sits by definition below the two-particle continuum in its channel. So string breaking does not bound the lightest state; it constrains excited/large-R states — i.e. it lives in the asymptotic large-R regime where the area law is stated, while the glueball is short-distance. String breaking is a **symptom of the same scale mismatch**, not an independent obstruction.

## (d) Circularity

The string-picture derivation assumes √σ is the mass-per-length of a flux tube *and* that its closed version is a 0⁺⁺ glueball of mass c√σ. That already presupposes the confinement→glueball identification one wants to prove. Using the breaking threshold 2m_g to reason about the string would be outright circular. Absent (b), the non-circular content is empty.

## (e) Transport from strong coupling

At strong coupling σ > 0 (Osterwalder–Seiler area law) and the full mass gap (Osterwalder–Seiler / cluster expansion) are **sibling theorems from the same expansion**. One cannot read σ ⇒ m off σ alone: there the gap comes from the cluster expansion, not from σ. Transporting to the continuum (β→∞) requires uniform control (Balaban-type RG) of both, with constants surviving the limit — open.

## Verdict

- **Proved (given):** σ > 0 at strong coupling. Full mass gap at strong coupling (separate theorem).
- **Conditional:** m_{0⁺⁺} ≥ c√σ holds at strong coupling **only if** the sector/spectral-weight bridge is added; c ≈ 3–4 numerically — consistent, not a theorem.
- **Open:** a uniform, non-circular σ ⇒ m_{0⁺⁺} inequality toward the continuum.

**The implication "σ>0 ⇒ m_{0⁺⁺} ≥ c>0" is NOT theorem-able from σ alone as stated.**

**One precise missing step:** a reflection-positive **sector/scale bridge** — an upper bound on the neutral 0⁺⁺ correlator at x ∼ 1/√σ (equivalently a lower bound on its spectral weight) derived from the charged-sector area-law decay. String breaking, circularity, and transport all reduce to this. The obstruction is a sector mismatch plus a wrong-direction inequality plus asymptotic-vs-bound — not circularity.
