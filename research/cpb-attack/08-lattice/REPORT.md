# Method 8 — finite-lattice rigorous gap + rigorous FSS

> CPB-attack probe. Verdict: **dead for the crossover certificate**.

## (a) Fixed-lattice certification: feasible only for tiny L, and certifies the wrong quantity

At fixed `β≈2.5` the lattice gap is `O(1)` in lattice units (`ξ/a ≲ 0.3 ⟹ ma ≳ 3`) —
a UV/lattice-scale state, **not** the physical glueball. Rigorous certification
(interval arithmetic / transfer-matrix spectral bounds on the gauge-invariant Hilbert
space) is technically possible at `L=4` with an electric cutoff, but the state space
grows like `exp(cL³)` and is infeasible at `L=6–8`. A gap at fixed `a` is **not** the
certificate: the certificate is `m/Λ`, reached only as `a→0`.

## (b) The FSS obstruction is the mass gap

The certificate needs the gap `≥c>0` uniformly in `L` and in `a→0`. Volume-uniformity
at `β=2.5` has no known proof: cluster expansions are volume-uniform only at small
`β`, and SZZ/IJGMP-type results do not reach `β=2.5`. The continuum limit `ma→0` must
be controlled against `a` — i.e., the physical mass — which *is* the mass gap. So (b)
is not a technical extension of (a); it is the target itself.

## Verdict

**Dead for the crossover certificate.** Needs-X: a volume- and spacing-uniform gap
(the RG / single-scale certificate).

## Honest boundary

A finite-L bound would be a decision input certifying a UV state, not the physical mass.
