# W4 — Certificate vs. gap: is the “≡” real?

**Verdict.** The certificate is a *sufficient* route, not a necessary one. “cert ≡ gap” is false as written. `conditional-theorem.md` §4 is correct; `CROSSOVER.md` §2 is wrong.

## (1) Is gap ⇒ cert argued anywhere?

No. Nothing in the docs or literature argues it (the “crossover certificate” is project-internal). Worse, for the project’s own `cert` the implication is **false**, not merely unproven.

- The only candidate, `CROSSOVER.md` §2(⟸): “if Δ>0, take the analysis that proves Δ>0 as cert.” This quantifies over *all* lower bounds — cert₍∃₎ — and merely restates Δ>0.
- But `cert μ` (`IN-10`, `Handoff.lean:77`) is the specific scale-analysis functional: a **one-sided lower bound**, `cert ≤ Δ/Λ`. A one-sided bound need not be tight: a gapped theory (Δ/Λ>0) is consistent with `cert ≤ 0` — the analysis just fails to certify. With `IN-10`, **gap ⇏ cert>0**. `IN-11` sharpens it: a gap can sit in a channel `0⁺⁺/2⁺⁺` does not see, so `cert ≤ 0` while Δ>0.
- So gap⇒cert holds only for cert₍∃₎ (vacuously), and fails for cert. The “≡” proof equivocates the two. The program already believes the theory gapped yet the crossover uncertifiable — presupposing this failure.

## (2) What the gap needs, independent of the certificate

- **Universal:** a rigorous nonperturbative IR statement — mass generation by transmutation, Δ∼Λ, exponentially small, invisible to finite-order β-expansions (T4). No rigorous 4D continuum pure-YM exists (Jaffe–Witten 2000; 2D: Osterwalder–Seiler 1978).
- **Route-specific:** the “crossover wall” — Pringsheim radius 1/48 vs β≈2.47, Savvidy negative modes, expansion window λ∈[~1,7.3] — belongs to the **certificate/expansion/RG/domain-decomposition family**; all nine mapped methods are this type.
- **Bypasses exist; they do not cross this crossover:** reflection-positivity + transfer-matrix uniform gap (Osterwalder–Schrader 1973/75; Nachtergaele–Sims 2006); measure-geometric CD(m,∞) contraction (Bakry–Émery 1985; Otto 2001; project W3-D); exact/deformation (SUSY: Seiberg–Witten; holography) — non-rigorous (T2). Each avoids β≈2.47 yet replaces it with an equally-or-more-open problem (uniform-in-(a,L) gap; Savvidy concavity → curvature→0); none is rigorous or excluded.
- Hence “**every** gap proof must cross the crossover” is **not established**: it conflates “all *known* methods fail here” with “all *possible* methods must.”

## (3) Honest statement

Replace `CROSSOVER.md`’s “≡” with: *the crossover certificate is one sufficient route (cert>0 ∧ IN-10 ∧ IN-11 ⇒ gap); necessity is false, so “stuck at cert ⇒ stuck at gap” does not follow, and certificate-free routes are not excluded.* Keep §4; demote “≡” everywhere it appears.
