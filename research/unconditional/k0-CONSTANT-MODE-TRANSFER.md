# k0-CONSTANT-MODE-TRANSFER — diagonalize the transfer operator on the zero-momentum mode: no invariant subspace, a quartic/flat valley, a gap that is the target respelled

> **Verdict (判词).** "Restrict `T` to the constant (zero-momentum) mode, diagonalize it as a 1-D Schrödinger operator, read off the ground/second gap" is **not a well-posed reduction**, and its gap is **not** an independently computable physical gap. **(A) No invariant subspace:** the electric term preserves momentum, but the plaquette (magnetic) term is a **momentum convolution** — at cubic order it maps `k=0` into every `k`, so `P₀TP₀` is no sub-operator; diagonalizing "on the constant mode" means *eliminating* `Θ(N_□)` modes, i.e. the full theory. **(B) Neither 1-D nor harmonically gapped:** the mode is `4(N²−1)`-dimensional; the tree potential `V₀(a)=Σ_□(1/2N)tr[a_μ,a_ν]²` is **quartic** transverse and **exactly zero** on the `4(N−1)`-dim Cartan valley (Wilson-line moduli), so `V₀″(0)=0` everywhere. **(C) Any curvature is the target:** elimination gives `V_eff″(0)=G″(0)=m_eff²(k→0)`, whose one-loop value is Savvidy-indefinite (`H_ff` indefinite). Born–Oppenheimer in mode space: H1 relocated by zero distance. **Fails #1,#2,#3,#5; #4 only formally.**

## (A) No constant-mode sub-operator
`H=(g²/2)ΣE²+(1/2g²)Σ_□(2−Re tr U_□)`: the electric term is mode-diagonal, the plaquette a **momentum convolution**, diagonal only at quadratic order (`k²|A_k|²`). So **only the free theory** decouples the mode (as the exact zero mode); at cubic order — the mass-generating nonabelian part — it creates `k≠0` quanta from it. Hence `[T,P₀]≠0`; a self-adjoint operator on `a₀` alone needs the other modes *eliminated* (Schur/Feshbach) — nonlocal, and H1.

## (B) Non-1-D; tree curvature zero
`a_μ∈g` is `4(N²−1)` real parameters. `V₀(a)=(V/4g²)Σ_{μ<ν}tr[a_μ,a_ν]²` is **quartic** near `a=0` (so `V₀″(0)=0`) and vanishes iff the `a_μ` commute — a `4(N−1)`-dim flat valley = the flat-connection/holonomy-torus moduli on `T⁴`, a centre order-parameter space. An exactly flat valley is not harmonically gapped; the compact-group measure/Vandermonde adds only a `β`-independent curvature `≥N/12` at `β=0` (k0-CONSTRAINT-EFFPOT) — entropic compactness, not dynamical `m_eff²`.

## (C) The curvature read off is G″(0)
`V_eff″(0)=E_μ[H_vv|0]−Var_μ(∂_nS|0)=G″(0)=m_eff²(k→0)` (PROOF-BORN-OPPENHEIMER; k0-CONSTRAINT-EFFPOT) — "read the gap off the operator" *is* "compute `G″(0)`." The one closed-form route, one-loop/Schur, is **undefined** (Savvidy-negative direction on the `k→0` slice, `H_ff` indefinite; PROOF-COARSE-POINCARÉ (a)(iii)). Exactly solvable models give only the Gaussian slice `1/S`, with `Var_red≤Var_true ⇒ gap_red≥gap_true` — the **wrong direction** (k0-EXACT-ONEMODE).

## Five-point verdict
1. **Nonperturbative — FAIL.** The only explicit operator is tree (quartic, curvature 0) or one-loop (Savvidy-indefinite); the nonperturbative object is the eliminated free energy = H1.
2. **Sign-controlling — FAIL.** Tree curvature `=0` on a flat valley; one-loop indefinite. The sign of `m_eff²(0)` is the open binary (`G″(0)≈−0.165` lean, VAR-FINAL).
3. **Non-circular — FAIL.** `V_eff″(0)=G″(0)` identically = the target; zero slack, a rename.
4. **Non-SUSY + fails free field — PASS, formally only.** Noncompact massless field decouples the mode → exact zero mode, `V₀≡0`, gap 0, `f₁²=0`. It passes because it *is* `m_eff²`, not via a discriminating mechanism.
5. **Uniform lower bound — FAIL.** (i) flat valley ⇒ tree curvature 0; (ii) the "second eigenenergy" is **not the glueball** — gauge-invariantly the constant mode is the holonomy/centre modulus, orthogonal to the flux sector carrying the gap; (iii) it is *one* mode: `Var(a₀)≤1/λ` fixes a single covariance entry, not the spectral infimum over all soft modes — flooring the physical gap is H1 on *every* soft mode, not this one alone.

## Exact failure
```
[H_B, P_{k=0}] ≠ 0  (plaquette = momentum convolution, cubic order) ⇒ P₀TP₀ not a sub-operator;
V₀″(0)=0 everywhere  (quartic transverse, zero on the flat Cartan valley);
V_eff″(0)=G″(0)=m_eff²(k→0),  one-loop value Savvidy-indefinite.
```
The premise dies: no invariant constant-mode sector exists to diagonalize, and the eliminated operator's curvature is exactly `G″(0)` — whose only computable value is the Savvidy-indefinite wall.

**Sharpest line.** Diagonalizing `T` on the zero-momentum mode is neither 1-D nor decoupled: the interaction that would *create* the operator is exactly what destroys the reduction; the tree potential is quartic with a flat holonomy valley; and the curvature a genuine elimination yields is `G″(0)` — the target respelled, undefined at one loop. Fails #1,#2,#3,#5; #4 a formal pass by identity.

*Refs.* PROOF-BORN-OPPENHEIMER (`V_eff″=G″`, fast≠slow); k0-CONSTRAINT-EFFPOT (`β=0` `N/12`); k0-EXACT-ONEMODE (wrong-direction, spurious free-point gap); PROOF-COARSE-POINCARÉ (H1, `H_ff` indefinite); PROOF-GAUSS-INVARIANT-PF (constant mode = holonomy/centre); HARDFORGE-SL-FULL (Savvidy concavity); PROOF-P3 (`β<1/48`); W6-BE-at-2.5; VAR-FINAL (`−0.165`); Gross–Pisarski–Yaffe, RMP **53** (1981) 43; Lüscher, CMP **54** (1977) 283.
