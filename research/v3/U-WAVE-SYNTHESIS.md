# U-WAVE-SYNTHESIS — binding is **source alignment**, not positivity; **weight and isolation are independent**; the continuum running is free; large-N binding survives

> **VERDICT (no QED, Clay 0%).** Six routes attacked the nonperturbative binding theorem and the
> continuum limit left by the T wave. The wave produces one master identity, one correction, two
> positive structural results, and three no-gos:
>
> 1. **Master identity (U1).** `w₁ = |⟨ψ₁|o−1⟩|² / Var(o)`: the atom fraction is the squared
>    direction cosine of the local source `o−1` onto the first excited eigenvector. Binding is
>    **source alignment**, not a property of the transfer kernel. PF, Harnack, the node theorem, and
>    TP₂ control the kernel, not the alignment; none gives `w₁≥w₀` unconditionally.
> 2. **Weight vs support (U6).** The uniform residue floor `W: w₁≥w₀` is a **weight** datum; it gives
>    `W⇒A` and nothing else. It is **incomparable** with isolation `I`: a threshold atom has weight
>    but no isolation; a thin isolated atom has isolation and a gap but `w₁→0`. The T wave's
>    ``w₁≥w₀ = isolation = H*'' is a weight-vs-support conflation, false both ways.
> 3. **The gap has two independent sufficient routes:** `W ∧ H1` (squeeze,
>    `m₁²≥f₁²/χ≥w₀λ`) and `I ∧ (ρ²≤C²)` (the window theorem, `μ₁≥δ/(C²−1)`). Neither route alone is
>    the gap; `W` and `I` are independent inputs.
> 4. **The continuum running is free (U4).** `γ_{w₁}=0` exactly: `f₁²` and `M₀` carry the same
>    `Z_O²`, so the ratio is RG-invariant. `w₁` has a finite `a→0` limit at fixed physical volume;
>    positivity is `[O]`, not an RG statement.
> 5. **Large-N binding survives (U3).** `w₁=O(N⁰)`, `ρ²=O(N⁰)` bounded; the `1/N²` cancels in the
>    ratio. `w₁∞=g₁/Σg>0` iff the lightest `0⁺⁺` is an isolated single-trace glueball (generic at
>    `N=∞`); the master field is unconstructed, so this is evidence, not proof.
> 6. **The diagonal tools fail (U2/U5).** The Wilson loop sees the diagonal energy density `σ`, not
>    the off-diagonal residue; `σ>0` gives only an upper bound. Entanglement is operator-blind and
>    anti-correlated with the atom: the exact single-mode link is symmetric under `w₁↔1−w₁`, and the
>    bound is vacuous exactly in the atom regime.

---

## §0 · The six verdicts

| Route | Verdict | Result |
|---|---|---|
| **U1** Harnack/PF lower bound | `[T]`/`[NG]`/`[O]` | Master identity `w₁=|⟨ψ₁|o−1⟩|²/Var(o)`. The vacuum Doob transform is the identity (ground state constant); the `ψ₁`-transform is not a Markov semigroup (needs the absent second cone). Harnack controls `|ψ₁|` shape on a nodal domain, never alignment; its only bound is conditional on source-sign alignment. Sturm/TP₂/node order is absent on `𝒜/𝒢`; even a strictly TP₂ Gaussian kernel gives `w₁=0` for a parity-misaligned composite. Exact no-gos: source direction (`cos²θ→0`), free parity (`w₁=0`), delocalization (`2/(L−1)→0`). Missing: `⟨ψ₁|o−1⟩²≥c₁Var(o)` uniformly (E-U1). |
| **U2** Wilson-loop OPE | `[T]`/`[NG]` | The area/shape derivative of `W(C)` equals the `F`-insertion and its contraction the `tr F²`-insertion on the worldsheet — exact (mod contact). Normalized by the loop, this is the **diagonal** on-loop energy density `σ`; the two-`tr F²` OPE gives the continuum `M₀` (a lower bound on the denominator). The overlap `w₁=f₁²/M₀` is the **off-diagonal** vacuum-to-glueball residue, which neither sees. `σ>0` gives at most an upper bound. 2D YM: `σ_R∝C₂(R)>0` but `tr F²` diagonal ⇒ `f₁²=0`. |
| **U3** large-N binding | `[T]`/`[C]` | `w₁=O(N⁰)`, not `N^{-1}` or `N^{-2}`: numerator and denominator are built from the same single-trace `O`, so the `1/N²` source normalization cancels. Large-N factorization makes `M₀` single-trace-dominated; multi-trace/continuum leakage is `1/N⁴`; `w₁→g₁/Σg_n+O(1/N²)`. `ρ²=μ₂/μ₁=O(N⁰)` bounded. Binding survives `[C]` iff the lightest `0⁺⁺` is an isolated single-trace glueball (generic at `N=∞`). 2D `tr F²` carrier: `w₁=0` at all `N` (acid pass); holonomy `ρ²=2N²/(N²−1)→2`. Strong coupling: `M₀=1/(2N²)`, `w₁=a₀/M₀=s/2`, `N`-uniform. |
| **U4** continuum limit | `[T]`/`[O]` | `γ_{w₁}=0` exactly: `f₁²` and `M₀` are both `Z_O²` times one insertion of `tr F²`, so `d ln w₁/d ln μ=2γ_O−2γ_O=0`. After renormalization at a fixed physical scale (gradient flow at fixed `t`, or point-splitting at fixed `ε`), `w₁` is finite at fixed physical volume; the bare ratio's apparent `→0` is the `M₀~a^{-8}` contact artifact. `ρ` is also `a`-independent. Positivity is `[O]`, not an RG statement: `γ=0` fixes the value, not its sign. Free field and 2D give `w₁=0`; 2D is the counterexample (`ρ²=8/3>1` yet `w₁=0`). |
| **U5** entanglement | `[T]`/`[NG]` | Single-mode reduction gives `w₁=Z`, `S=h(Z)`, so `w₁≥h^{-1}(S)`; but `h(Z)=h(1−Z)` is symmetric, so `S=0` fits both `w₁=0` and `w₁=1`. The bound is vacuous exactly in the atom regime and strengthens only as `w₁→½` — anti-correlated with the target. No universal floor: free massless chain has `S∝log l` with `w₁=0`. MI is never negative (SSA); an entanglement deficit must be relative entropy against a reference tube = the atom = circular. `w₁` is not a function of `ρ_A` (same `S`, two operators: `w₁=0.20` vs `0.82`). 2D: `σ>0` but `w₁=0`. |
| **U6** binding logic | `[T]`/`[NG]` | `W: w₁≥w₀` is a weight datum; it gives `B⟹A` (strict) and nothing else. `W` and `I` are incomparable: a threshold atom has weight but `Λ★=m₁²` (`W∧¬I`); a thin isolated atom has `I+G` and bounded `ρ²=2` but `w₁→0` (`I∧¬W`). The T wave's ``w₁≥w₀ ⟺ isolation'' is a weight-vs-support conflation. The window's lower bound `ρ²≥1+c` is isolation, not binding. Double well: `w₁≈0.97` uniform, `Δ` large, but `λ₁→0`, `χ→∞`, `ρ²→∞` (`W∧I∧¬G`). Gap only via `W∧H1` (squeeze) or `I∧(ρ²≤C²)` (window). Strongest no-go U6-NG1: weight data independent of `I,G`; cannot force `ρ²≤C²`. |

---

## §1 · The master identity: binding is source alignment

For a local source `o(x)=¼tr F²(x)/⟨O⟩` and the transfer matrix `T` with ground state `ψ₀` and first
excited state `ψ₁`,

```
   w₁ = |⟨ψ₁|o−1⟩|² / Var(o) .                                      [T]
```

This is a **direction cosine**: `w₁` measures how much of the source `o−1` lies in the first excited
channel. It is not a property of the kernel `T` alone. Three consequences:

1. **PF/Harnack/node/TP₂ cannot give a lower bound.** They control the *shape* of `ψ₁` (or the
   ordering of eigenvalues), not its *alignment* with the source. The Harnack bound is conditional on
   a source-sign hypothesis; the node theorem needs an order that is absent on `𝒜/𝒢`; even a
   strictly TP₂ Gaussian kernel gives `w₁=0` for a parity-misaligned composite.
2. **The exact no-gos are source-side.** The direction cosine can vanish at fixed kernel (source
   direction), on the free field (parity), or in the delocalized limit (`2/(L−1)→0`).
3. **The missing estimate is a source-alignment bound:** `⟨ψ₁|o−1⟩² ≥ c₁ Var(o)` uniformly at
   β=2.5 (E-U1). This is not closable by positivity; it needs binding *into the source's channel*.

---

## §2 · Weight and support are independent

The T wave identified the missing object with ``uniform `w₁≥w₀` = isolation = `H*`''. U6 shows this
is a conflation of two different data:

| datum | definition | implication |
|---|---|---|
| **A** (visibility) | `f₁²>0` | — |
| **W** (weight floor) | `w₁=f₁²/M₀≥w₀>0` | `W⇒A` |
| **I** (isolation) | `Λ★=inf supp(ν−f₁²δ)>m₁²` | `I⇒B_spec⇒A` (T6) |
| **G** (gap) | `m₁²>0` | — |

`W` and `I` are **incomparable**:

- threshold atom: `W` holds, `Λ★=m₁²`, so `¬I`;
- thin isolated atom: `I` and `G` hold, `ρ²=2` bounded, but `w₁→0`, so `¬W`.

The gap has **two independent sufficient routes**:

```
   W ∧ H1        ⇒  G      (squeeze: m₁² ≥ f₁²/χ ≥ w₀λ)
   I ∧ (ρ²≤C²)   ⇒  G      (window theorem: μ₁ ≥ δ/(C²−1))
```

Neither `W` nor `I` alone gives `G`; the double well satisfies `W∧I` yet has `λ₁→0`,
`χ→∞`, `ρ²→∞`. So the campaign's target is the **pair `(W, I)` with the two side conditions**,
not a single binding inequality.

---

## §3 · The diagonal/off-diagonal split

The Wilson-loop identity is exact but diagonal: the shape derivative sees the on-loop energy density
`σ`, and the two-`tr F²` OPE sees `M₀`. The atom fraction `w₁` is the off-diagonal
vacuum-to-glueball residue. Hence `σ>0` cannot lower-bound `w₁`; the 2D model (`σ>0`, `f₁²=0`)
is the counterexample. This is the same charged-vs-neutral blindness as the center-vortex route
(T5), now at the level of the loop OPE.

---

## §4 · Entanglement is operator-blind

The single-mode identity `w₁=Z`, `S=h(Z)` looks like a bound, but `h(Z)=h(1−Z)` makes it symmetric:
`S=0` is compatible with both `w₁=0` and `w₁=1`. The bound is vacuous in the atom regime and
anti-correlated with it. Entropy is a basis-independent aggregate; the atom is a basis-dependent
coherence. The free massless chain (`S∝log l`, `w₁=0`) kills any universal floor. An entanglement
deficit would have to be measured against the atom itself, which is circular.

---

## §5 · The continuum limit is free; positivity is not

`γ_{w₁}=0` exactly, because both `f₁²` and `M₀` carry the same `Z_O²`. After renormalization at a
fixed physical scale, `w₁` is finite at fixed physical volume and `ρ` is `a`-independent. The RG
half of the continuum problem is therefore solved: the ratio does not run. What remains is
positivity, `inf_{L,a,t} w₁≥w₀>0`, which is the same `[O]` as before.

---

## §6 · Large-N: binding survives, with a flat N-scaling

`w₁=O(N⁰)` and `ρ²=O(N⁰)`; the `1/N²` cancels in the ratio. Large-N factorization makes `M₀`
single-trace-dominated, and `w₁∞=g₁/Σg>0` iff the lightest `0⁺⁺` is an isolated single-trace
glueball. This is evidence for the N-structure of a uniform binding theorem, not a proof: the master
field is unconstructed. The strong-coupling test gives `w₁=s/2`, N-uniform.

---

## §7 · What the wave won, and the next object

**Won:**
- The master identity: binding is source alignment.
- The weight-vs-support correction: `W` and `I` are independent; the gap has two sufficient routes.
- The diagonal/off-diagonal split: Wilson loops and entanglement are blind to `w₁`.
- The continuum running is free (`γ_{w₁}=0`).
- Large-N binding survives with a flat N-scaling.

**Not won:**
- No proof of `W` or `I` at β=2.5.
- No source-alignment bound `⟨ψ₁|o−1⟩²≥c₁Var(o)`.
- No reduction of Clay: the gap still needs `W∧H1` or `I∧(ρ²≤C²)`.

**The next object** is a **source-alignment + isolation pair**: (i) a uniform lower bound on the
direction cosine of the local `tr F²` source onto the first excited transfer eigenvector, and
(ii) the TP₂/window isolation bound. They are independent; neither is a variance, positivity,
spectral, Perron, RG, entanglement, or Wilson-loop statement.

**One line.** *Binding is the alignment of the `tr F²` source with the first excited state; weight and
isolation are independent; the continuum ratio does not run; large-N binding survives; Wilson loops
and entanglement are diagonal or operator-blind. No proof, no QED, Clay 0%.*

*Refs:* `U1-HARNACK-BOUND.md`, `U2-WILSON-OPE.md`, `U3-LARGE-N-BINDING.md`, `U4-CONTINUUM-ATOM.md`,
`U5-ENTANGLEMENT-BINDING.md`, `U6-BINDING-LOGIC.md`; witnesses `w7/u1_*.py` … `w7/u6_*.py`;
`T-WAVE-SYNTHESIS.md`, `R-WAVE-SYNTHESIS.md`.
