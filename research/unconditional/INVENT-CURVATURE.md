# INVENT-CURVATURE — thermodynamic-geometry (Fisher–Rao / Ruppeiner) curvature is **not** a new strut: it is *flat* (1 parameter), *constant −1* on every Gaussian (acid-vacuous), *homogeneous in χ* (2 parameters), and *embedding-dependent* (not a functional of μ alone)

> **VERDICT (no QED).** The bare-handed round tried to invent an information-geometric species: the scalar
> curvature `R(β)` of the Fisher–Rao / Ruppeiner metric on the family `{μ_β}`. **It is [NG] as the strut.**
> Four exact obstructions, any one fatal: **(A)** a 1-parameter family is flat (`R≡0`), and pure YM at fixed
> volume has one coupling ⇒ the naive curvature is vacuous; **(B)** every Gaussian family has scalar
> curvature `R≡−1` (Gauss `K=−1/2`), identical for massless and massive free fields ⇒ acid-**vacuous**;
> **(C)** the metric component in the source direction *is* the susceptibility, `g_hh = ∂²_h log Z = Var(tr F²)
> = χ` exactly, so any 2-parameter curvature is a function of `χ` and its slopes ⇒ **homogeneous** (same wall
> as Theorem BE); **(D)** `R` depends on the *chosen* 2-parameter family (two Gaussian embeddings give
> `R=−1` vs `R≈−6.9`), so it is **not a functional of `μ` alone** (`PROFILE` P1 fails) and its drift is
> READ non-Gaussianity, not a sign-definite force. **Exact failing step: `R` inherits `χ` (§C) and needs a
> variable choice (§D); bounding `R` ⇔ bounding `χ` = `H1`. Clay 0%.**

---

## §0 · Labels, iron laws, acid

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the class · **[V-num]** new numeric.
**IL-i** non-positivity (must fail on free massless). **IL-ii** nonabelian. **IL-iii** no small parameter.
**Acid**: `≤0`/`=0` on free massless *and discriminatingly* on the gapped-atomless `:φ²:`.

---

## §1 · The candidate species

Thermodynamic geometry attaches a Riemannian metric to a family of measures. For the gauge family
`{μ_β}` (grand-canonical at coupling `β`, plus optionally a source `h` for `tr F²`):

- **Fisher–Rao:** `g_{ij} = E[∂_i log p · ∂_j log p] = ∂_i∂_j (−log Z)` (the Hessian of the free energy),
- **Ruppeiner:** the same metric up to the sign convention, `g_{ij} = −∂_i∂_j (free energy)`,
- and `R(β)` = its scalar curvature — a *coordinate-free* functional, which is why it looked promising as a
  "phase-sensitive" strut.

Witness: `w7/invent_curvature.py` (numpy, `py -3.12 -I`, clean).

---

## §2 · The four tests

### (A) One parameter ⇒ flat `[T]`
A 1D Riemannian manifold has no curvature tensor: `R ≡ 0`. Pure YM at fixed volume has **one** coupling
`β`; the naive curvature of `{μ_β}` is **identically zero**. To have any curvature at all one must add a
second parameter (a source `h`, or the volume). So the *obvious* object is vacuous.

### (B) Gaussian ⇒ constant `[T]/[V-num]`
With the second parameter taken as `(μ, σ)` (mean, variance of the coarse mode), the Fisher metric is the
Poincaré half-plane form `g = diag(1/σ², 2/σ²)`. Bonnet's formula gives **Gauss curvature `K=−1/2`**,
i.e. **scalar curvature `R=−1`**, *independent of `σ`*. Witness (`σ=1, 0.1, 10`):
```
Gaussian sigma=1 : scalar R = -1.000
Gaussian sigma=0.1: scalar R = -1.003
Gaussian sigma=10 : scalar R = -1.000   (approach massless: chi->inf)
```
⇒ **every Gaussian (massless AND massive) gives `R=−1`.** The curvature does **not** distinguish gapless
from gapped free field. **Acid-vacuous** in the discriminating sense: it does not fire differently on the
two adversary classes.

### (C) `g_hh = χ` exactly ⇒ homogeneous `[T]/[V-num]`
For the family `μ_h ∝ e^{−S + h·O}` with `O = tr F²`, the second derivative of the log-partition is the
variance of `O`:
```
g_hh = ∂²_h log Z = Var(O) = χ     (the susceptibility / inverse G″(0))
```
Witness (Gaussian `a=1,b=0` and quartic `b=0.08`):
```
g_hh = 2.000000 = Var(O)=chi = 2.000000   (diff 2e-8)
g_hh = 1.250702 = Var(O)=chi = 1.250702   (diff 1e-8)
```
⇒ the metric **literally contains the susceptibility `χ`** as its source-direction component. Any
2-parameter curvature `R = R(g_{ββ}, g_{βh}, χ; ∂g)` is therefore a functional of `χ` and its slopes;
bounding `|R| ≥ c` requires controlling `χ`. **Homogeneous** — precisely the failure mode of Theorem BE.
The wall, in geometric clothing.

### (D) Embedding-dependent, not a `μ`-functional `[T]/[V-num]`
`R` depends on *which* 2-parameter family you embed the measure in. Two embeddings of a **Gaussian** give
different curvature:
```
(mu, sigma) Gaussian family        : scalar R = -1
(a, b) quartic family at its Gaussian point b=0 : scalar R ≈ -6.9
```
So `R` is **not a functional of `μ` alone** — it needs a choice of thermodynamic variables. This violates
`PROFILE` P1 (must be a functional of the field measure). And the drift of `R` with the quartic coupling `b`
(`b=0→0.5`: `R` from `−6.9` to `−2.0`) tracks **non-Gaussianity** — the `η`/3rd–4th-cumulant defect — which is
a **READ** quantity (it reports the measure), not a **FORCE**, and whose sign is set by the interaction
(`Ruppeiner`: `R<0` vs `R>0` for different interaction signs), not by any general principle.

---

## §3 · Why it cannot be the strut (and the tie to B2)

The strut must be a **κ-independent, sign-definite source** on the field measure (`PROFILE` P4), acid-clean,
inhomogeneous. Thermodynamic curvature is none of these:

1. **Flat** (1 parameter) or **homogeneous in `χ`** (2 parameters, §C) — never a κ-independent *additive* source.
2. **Constant `−1` on all Gaussians** (§B) — cannot tell gapless from gapped free field.
3. **Embedding-dependent** (§D) — not a functional of `μ` alone (fails P1).
4. **A transition detector, not a gap detector.** Ruppeiner's heuristic is `|R| ~ correlation volume`,
   and `R` *diverges at a critical point* `β_c`. But **(B2, roof round)** zero-T 4D pure YM is a **single
   phase** — only the strong-weak **crossover**, no thermodynamic critical point. So `R` has **nothing to
   signal** at `β≈2.5`: it is smooth through the crossover, like every other smooth thermodynamic datum.
   The one regime where curvature *would* fire (a critical point) does not exist.

The only content `R` adds is a geometric repackaging of **non-Gaussianity = `η`**, already known as a
READ-not-FORCE object (`Ward/Stein`). It is a **relabelling, not a new species**.

---

## §4 · Exact failing step and honest bottom line

- **Failing step `[O]`:** §C (`g_hh = χ`) + §D (embedding dependence). `R` inherits `χ` and needs a variable
  choice; hence **bounding `R` is equivalent to bounding `χ` = `H1`**, and it is not even a well-defined
  functional of `μ`.
- **Verdict `[NG]`:** thermodynamic-geometry curvature is **not** a κ-independent sign-definite source; it is
  *flat* ∨ *homogeneous in `χ`* ∨ *embedding-dependent* ∨ *constant on Gaussians*. The information-geometric
  route does not escape Theorem BE — it is Theorem BE in Riemannian dress.
- **What would still be needed (unchanged):** a **non-order-parameter** (B2), **non-Markov** (BE),
  **non-topological** (rigidity), **euclidean-nonlocal**, **sign-definite**, **`μ`-alone** functional — i.e.
  exactly the "new species" the roof round named, of which curvature is one more **non**-instance.

**Acid:** free massless ⇒ `R=−1` (the Gaussian constant); `:φ²:` ⇒ also a Gaussian (`R=−1`) — no
discrimination. The witness is regression-safe: `g_hh=χ` verified to `1e−8`, Gaussian `R=−1` to `1e−3`.

**Clay 0%. No QED. Not committed (parent handles).**

*Files:* `research/unconditional/w7/invent_curvature.py` (witness, runs clean).
*Refs:* `ROOF-SYNTHESIS.md` (B2, Theorem BE), `BASE-FIELDMEASURE.md` (Theorem BE), `PROFILE.md` (P1, P4).
