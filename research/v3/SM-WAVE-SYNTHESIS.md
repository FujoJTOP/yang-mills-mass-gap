# SM-WAVE-SYNTHESIS — the soft-mode margin is positive at small volume (min +0.32), decays ≈1/L (not the box gap), with an L-independent Haar floor; decompactification is the remaining question

> **VERDICT (no QED, Clay 0%).** Six angles attacked the Savvidy soft mode, the one mode on
> which the UV descent breaks. The result is **positive at small volume and sharply conditional at
> large volume**:
>
> 1. **SM1 `[T]`/`[V]`:** the conditioned soft-mode curvature is **positive in every small-volume
>    SU(2) model** (1×1, 1×2, 1×3, 1×4, 2×2 cross, 4D single site), minimum `+0.32 > N/12=0.167`.
>    Its scaling is `≈1/L` (`1.01→0.65→0.42`), **not** the `4π²/L²` box-gap collapse; it increases
>    with `β` (`0.49→0.83→1.60` for `β=1→4`); the `β=0` Haar floor is `0.23–0.25` and
>    **`L`-independent**. The full Hessian is indefinite already at small volume
>    (`frac λ_min<0`: `0.08→1.00` with `L`), so `γ_j<0`, `η_j≥1`, `D=+0.19>0`: the positivity is
>    **measure-carried, not a Hessian gap**.
> 2. **SM4 `[T]`/`[V]`:** the constant mode has an exact closed-form marginal; at `β=2.5` its floor
>    is `+0.36…+0.69` (Cartan and off-diagonal/Savvidy classes), twist-robust, **no `1/L²`
>    softening** on `L=1,2,3`.
> 3. **SM3 `[T]`:** anisotropic blocking makes `η=0.28` at `β=2.5` (vs `22` isotropic) and gives a
>    volume-uniform β-dependent curvature floor `≈0.40`; but the soft mode is a blocking fixed
>    point, so no blocking lifts `G″_soft` above its bare value.
> 4. **SM2 `[T]`:** exact overlap lemma and margin identity `G″=(S−D_Sav)+|D₊|=m_soft+|D₊|`; the
>    refined bound `D_Sav≤g²Var(gB,λ₄)` corrects the corpus (`2/λ₄` is an asymptote, peak
>    `2.21/λ₄`); the ceiling `m_soft≤G″≈3.5%·S` shows positivity is **BL-carried**.
> 5. **SM5 `[NG]`/`[T]`:** TCI(strong) is false at the Savvidy barrier (it would imply the gap);
>    TCI(typical) is true but pointwise; TCI(avg) is the Fisher identity, uniform margin `=H1`.
> 6. **SM6 `[V]`:** the existing lattice data are gauge-dependent plane waves; the sampled soft
>    mode softens with volume (`+0.200→+0.164→+0.142`; `L=8` exact `+0.080±0.189`, zero-consistent),
>    so positivity is **not certified**. The minimal decisive run is specified.

---

## §1 · What is now established

- A **positive soft-mode margin at small volume**, with an explicit floor (`min +0.32`) and a slow
  `≈1/L` decay. This is the first positive, quantitative soft-mode statement in the campaign.
- The positivity is **measure-carried** (Haar/BL), not a Hessian gap; the negative mode only thins
  it (ceiling `≈3.5%·S`).
- The **constant mode** (the actual `k=0` IR mode) has an exact floor with **no `1/L²` softening**;
  it is the most `L`-stable object found.
- Anisotropic RG controls the descent (`η<1`) but cannot lift the fixed-point soft mode.

## §2 · The decisive open question

The soft-mode floor decays `≈1/L` at small volume, while the `β=0` Haar floor is `L`-independent.
The remaining question is **decompactification**:

> Does the `≈1/L` decay tend to zero (then the uniform floor is `H1`), or does it saturate at the
> `L`-independent Haar floor (`0.23–0.25`)?

The two scenarios are distinguishable by a constant-mode calculation at `L=5–8` and by separating
the constant mode from the lowest nonzero mode. The lattice run specified in SM6 (Landau gauge +
Gribov copies, block-Hessian diagonalization at `β=2.5`, `L=8,12`) would test the same question
statistically.

## §3 · What this changes

- **Not** a proof of the gap: the uniform margin is still `H1`, and the strong forms are false.
- **A real positive signal:** the IR soft mode is not a Hessian gap; it is a measure-carried floor
  that is positive at small volume and slowly decaying. This is exactly where a proof could live if
  the Haar floor survives decompactification.
- **A sharper target:** decompactification of the constant mode, not another broad no-go wave.

**Clay = 0%. No QED.**

*Refs:* `SM1-SOFT-EXACT.md`, `SM2-RANK-ONE.md`, `SM3-ANISOTROPIC-RG.md`, `SM4-CONSTANT-MODE.md`,
`SM5-TCI-HEAD.md`, `SM6-SOFT-LATTICE.md`; witnesses in `w7/`; `STAR2-CLASSIFICATION.md`.
