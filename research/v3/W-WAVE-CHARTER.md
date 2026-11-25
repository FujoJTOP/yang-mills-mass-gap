# W-WAVE-CHARTER — the next wave: the weight/source-alignment problem, stated precisely enough to attack

> **Status.** This is a charter, not a result. It fixes the problem statement for the next wave after
> U1–U6, so that the next round does not drift back to variance/positivity tools. **No QED. Clay 0%.**
> **Labels:** [T]/[C]/[NG]/[O]/[V].

---

## §0 · The two targets (independent; U6)

```
   W  (weight)      w1 = f1^2 / M0 >= w0 > 0          uniformly at beta=2.5
   I  (isolation)   rho^2 = mu2/mu1 <= C^2,  Lambda* > m1^2   (window theorem)
```

with the master identity (U1)
```
   w1 = |<psi1 | o-1>|^2 / Var(o) ,
```
where `o = (1/4) tr F^2 / <O>` is the local source and `psi1` the first excited transfer eigenvector.
Neither target implies the other; the gap needs `W ∧ H1` or `I ∧ (rho^2 <= C^2)`.

**The W-wave attacks `W` and `I` separately, and never conflates them.**

---

## §1 · The six routes (each with a deliverable and an acid test)

### W1 · Strong-coupling exact alignment and continuation
Compute `w1` exactly in the strong-coupling character expansion (U3 gives `w1=s/2`), including the
first correction. Determine whether the series has a positive radius and whether a Padé/Borel
continuation reaches β=2.5 with `w1>0`. **Deliverable:** the exact coefficients, the continuation, and
the first sign change or the missing estimate. **Acid test:** free field `w1=0`.

### W2 · Large-N master field
Construct the large-N master field / first excited eigenvector and compute
`w1∞=g1/Σg` analytically. Determine whether `g1/Σg` is bounded below uniformly in the 't Hooft
coupling, and whether the `1/N^2` correction preserves positivity. **Deliverable:** the master-field
equation for `psi1`, the ratio `g1/Σg`, and the N-correction. **Acid test:** 2D `tr F^2` `w1=0`.

### W3 · Exact small-volume alignment
Compute `w1` exactly on L=2,3 at β=2.5 by transfer matrix / character expansion (no Monte Carlo).
Determine the L-scaling and whether `w1` is already positive at small volume. **Deliverable:** the
exact table and the finite-size law. **Acid test:** free field and one-plaquette.

### W4 · Source optimization
`w1` depends on the source `o`. For a variational family of smeared sources (APE/gradient flow at
fixed physical radius), compute the optimal `w1` and determine whether a **provable** lower bound
follows from a Rayleigh/min-max argument. **Deliverable:** the variational principle and the best
bound. **Acid test:** free field; a parity-misaligned source must give `w1=0`.

### W5 · Operator mixing and the renormalized source
`tr F^2` mixes with other scalar operators under RG. Determine whether the physical atom couples to a
combination that improves the alignment, and whether the mixing matrix can be controlled at β=2.5.
**Deliverable:** the mixing matrix and its effect on `w1`. **Acid test:** `gamma_w1=0` (U4) must be
recovered in the no-mixing limit.

### W6 · Adversarial: separate W from I and G
Construct models with `W∧¬I`, `I∧¬W`, `W∧I∧¬G`, and `G∧¬W` (the `:φ²:` gapped atomless case).
Determine the strongest no-go and the weakest hypothesis under which `W` implies `I` (or vice versa).
**Deliverable:** the separation table and the falsification battery. **Acid test:** free field.

---

## §2 · The acid-test battery (every W route must pass all)

1. **Free massless field:** `w1=0` exactly (no manufactured atom).
2. **2D pure Yang–Mills:** `tr F^2=H` conserved ⇒ `w1=0` at all N (the charged/neutral blindness).
3. **One plaquette:** no time direction; the atom is degenerate; the route must say so, not invent one.
4. **Strong coupling:** `w1=s/2`, N-uniform (U3); the route must reproduce it.
5. **Large N:** `w1=O(N^0)`, `rho^2=O(N^0)` (U3); a route that gives `N^{-1}` or `N^{-2}` is wrong.
6. **Continuum:** `gamma_w1=0` (U4); a route that predicts a running ratio is wrong.

---

## §3 · What a successful W-wave would and would not prove

- `W` alone proves the atom exists with a uniform weight; it does **not** prove isolation or the gap
  (U6).
- `W ∧ H1` proves `m1^2 ≥ f1^2/chi ≥ w0 λ` at the lattice coupling; it is a **lattice mass gap**, not
  the continuum Clay statement.
- `I ∧ (rho^2 ≤ C^2)` proves the gap from the support side; it is independent of `W`.
- A Clay paper needs **both** routes controlled uniformly in `(L,a)`, plus the continuum/OS
  reconstruction. The W-wave is one input, not the whole proof.

---

## §4 · What not to do

- Do not attack `W` with a variance, positivity, spectral, Perron, RG, RMT, duality, center-vortex,
  Wilson-loop, entanglement, or Harnack tool: U1/U2/U5/T1/T5/Q1/Q2/Q7 already show each is either
  diagonal, operator-blind, direction-reversed, or blind at the 3–40% level.
- Do not identify `W` with isolation (U6); the two are incomparable.
- Do not claim the gap from `W` alone; the double well has `W` and `I` with no gap.
- Do not launch a route without the six acid tests in §2.

*Refs:* `U-WAVE-SYNTHESIS.md`, `U1-HARNACK-BOUND.md`, `U3-LARGE-N-BINDING.md`,
`U4-CONTINUUM-ATOM.md`, `U6-BINDING-LOGIC.md`, `T3-TP2-ISOLATION.md`, `T6-BINDING-ADVERSARIAL.md`.
