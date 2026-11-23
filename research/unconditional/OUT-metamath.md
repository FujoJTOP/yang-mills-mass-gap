# OUT — Metamathematics of the 4D YM mass gap

> **Verdict: a normal open analytic statement.** The single rigorous metamathematical handle is
> *negative*: a coding/absoluteness argument makes the statement Σ¹₂, hence invariant under forcing.
> That rules out CH-style independence and makes `choice` harmless — but it builds no measure. There
> is no large-cardinal, new-axiom, or relative-consistency lever: the difficulty is analytic *width*,
> not consistency-strength *height*.

## What the statement actually is (correcting "Π₂")

`BREAK-IN-8` writes Clay as `∀G∃Δ` ("Π₂"). That is quantifier-skeleton shorthand, **not** the
arithmetic hierarchy. The matrix quantifies over a generating functional `Z` on `𝒮(ℝ⁴)`
(`YangMillsInfra/Basic.lean`: OS0–4 predicates on `Z`, exponential-clustering bound `∃C ∀a`), i.e. over
**measures on a Polish space**. Coding such a `μ` (and `m, C`) by reals, with `∀G` ranging over the
classified, computably-enumerable list of compact simple groups, **Clay-YM is Σ¹₂** (third-order
arithmetic). Not arithmetic Π₂: the matrix is analytic, not decidable.

## The one real lever: Shoenfield absoluteness

Σ¹₂ statements are absolute between V, L, and every forcing extension (Shoenfield 1961). Hence:

- **(a) Independent-by-forcing: impossible.** No CH/Whitehead-style forcing proof can separate ZFC from
  ¬YM; the truth value is fixed across the generic multiverse. Non-forcing independence is not excluded,
  but no mechanism is known, and YM is *not* Gödelian — it encodes no Con-sentence.
- **`choice`/forcing is harmless.** The Lean base (`propext/choice/Quot.sound`) cannot be an independence
  source. AC-heavy constructions may be used freely without weakening the result.

## (b) New axiom / large cardinal: no handle

Large cardinals buy *consistency strength* and (via determinacy) decide projective statements — but Σ¹₂
is already L-absolute and forcing-invariant, so they change nothing here, and they *construct no measure*.
A "YM axiom" carries no consistency strength over ZFC: `Con(ZFC) → Con(ZFC + ∃ massive YM measure)` is
trivial if true and otherwise *is* the construction. Relative-consistency tricks collapse to the analytic
problem.

## (c) Equivalence to a known hard statement: none, and it would not help

The project's reductions (G1 ↔ constructive QFT, G4 ↔ 4D CFT uniqueness) are to **open analytic**
statements in the same subject, not to known-independent ones (`T2`, `W3-adjacent`). No coding to a
logical statement exists; by Σ¹₂-absoluteness such an equivalence would carry zero independence content.
There is no independence result for any QFT or statistical-mechanics statement — **zero precedent**.

## Honest boundary

The Σ¹₂ claim assumes standard coding (measures on Polish spaces as reals; countable dense test families).
A non-standard coding demanding a full AQFT functor could reach Σ¹₃, where Shoenfield does not apply — but
no independence mechanism is known there either. Net: metamathematics yields one theorem
(forcing-invariance, `choice`-harmlessness) and **no construction lever**. G1/G2 stay analytic walls; the
logical framing removes worries, not difficulty.
