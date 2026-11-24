# ROOF-RESURGENCE — resurgence pins the *obstruction* (the u=2 ambiguity **is** a resurgence invariant, computable from perturbation theory), **not** the condensate (the transseries parameter, free)

> **VERDICT (no QED).** Resurgence does **not** pin `⟨tr F²⟩` scheme-independently and does **not** cancel the `u=2` IR-renormalon ambiguity. What it *does* do — and this is the sharpest possible reformulation of the wall — is split the object into two, one **computable** and one **not**: the **ambiguous part is a resurgence invariant** (the Stokes constant = the Borel residue at `u=2`, readable from the large-order coefficients `a_n ~ n!`), while the **condensate itself is the transseries parameter `σ`**, a free integration constant that perturbation theory (with or without resurgence) never fixes. `[T]` for the dictionary; `[NG]` for resurgence-as-pinning. Equivalently: **perturbation theory + resurgence knows the ambiguity exactly, and the condensate not at all.** **Clay 0%.**

---

## §0 · Labels & iron laws

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric · **[R]** inherited. Acid: free massless ⇒ `=0`, `:φ²:` ⇒ `=0`. IL-i non-positivity, IL-ii nonabelian, IL-iii no small parameter (all carried; not re-derived here — this report is a *reformulation* of the obstruction, not a construction).

---

## §1 · The object, in resurgence language

The gluon condensate `⟨tr F²⟩` is defined by the OPE as the coefficient of the unit operator; its perturbative estimate is an asymptotic series in `α_s` whose **Borel transform has an IR renormalon pole at `u = d/2 = 2`** (`d=4`). In transseries language the physical quantity is

```
F(s) = [ perturbative Borel sum ]  +  σ · e^{−2s} · [ series ],
```

where `s ~ 1/α_s` (so `e^{−2s}` is the `u=2` nonperturbative scale), and

- **`σ` = the transseries parameter = the condensate magnitude** (the nonperturbative input);
- the **Borel residue at `u=2`** controls the **large-order behaviour** of the perturbative coefficients and the **Stokes constant** across the `u=2` singularity.

This is the standard `u=2` structure (Mariño; Aniceto–Schiappa–Vonk; Costin). The question the roof needs: **does the resurgent structure fix `σ`?**

---

## §2 · What resurgence gives, exactly `[T]`

Three exact statements, all verified on the canonical toy (Euler series, Borel `1/(1−t)`, pole at `t=1`, nonperturbative scale `e^{−1/x}=e^{−2s}`):

1. **The ambiguity is a resurgence invariant.** The lateral Borel sums `S_±` differ by
   `S_+ − S_− = 2πi · Res_{u=2} = 2πi · (−e^{−2s})`, i.e. `|S_+ − S_−| = 2π e^{−2s}` — **exactly** the residue (residue theorem). The residue = the **Stokes constant**, **computable from the large-order asymptotics** `a_n ~ n!` (witness C: `a_n/n! ≡ 1`).
2. **The median (PV) sum is real and canonical.** `(S_+ + S_−)/2 = PV` is real and scheme-canonical (witness A: PV real at every `x`).
3. **The alien calculus closes on the perturbative series.** The Stokes automorphism `Δ_{θ}` relating the two sectors is fixed by the residue — so the *ambiguity* is fully determined by perturbation theory.

⇒ **Resurgence gives a complete, exact dictionary for the ambiguity.** The `u=2` obstruction `ROOF-CONDENSATE-FLOOR` found by hand (Borel pole, ambiguity `~Λ⁴` = condensate's own size) is here identified as **the Stokes constant**, a computable resurgence invariant. That is a genuine sharpening: the obstruction is not "unknown," it is **an exact invariant of the perturbative series.**

---

## §3 · What resurgence does NOT give `[NG]`

1. **`σ` is free.** `F = PV + σ e^{−2s}` is real for **any** real `σ` (witness B: `σ ∈ {0, −0.5, 1.7, 10}` all give real `F` from the **same** PV and the **same** perturbative data). Reality does not pin `σ`. There is no equation for it inside perturbation theory.
2. **The median sum is canonical but not physical.** The median corresponds to `σ = 0` (or to a specific symmetric prescription); the *physical* condensate has `σ ≠ 0`, and the median sum's value is **not** the physical value unless `σ` is supplied from outside.
3. **`σ` does not appear in `{a_n}`.** The large-order coefficients encode the *residue*, never the *transseries parameter* (witness C). So no amount of higher-order perturbation theory, and no resurgence resummation of it, can read `σ`.

⇒ **Resurgence is the sharpest statement of the wall, not a way through it:** the condensate is *by construction* the one datum perturbation theory cannot know. This is `ROOF-CONDENSATE-FLOOR`'s conclusion, now exact: the ambiguity **is** the Stokes constant; the condensate **is** the transseries parameter; they are structurally different objects.

---

## §4 · The only loopholes (and why they close for pure YM)

Resurgence *can* pin `σ` in special cases, always via a **second, independent nonperturbative handle**:

- **(L1) A known "instanton" amplitude.** If the nonperturbative sector is a genuine saddle with a computable fluctuation determinant, the saddle fixes `σ`. For the `u=2` **renormalon** there is **no saddle** — the "instanton" is the **scale anomaly / the condensate itself** (the `u=2` singularity is an IR renormalon, not an instanton), so there is no independent amplitude to fix `σ`. The identifying-the-ambiguity-with-the-condensate is exactly `ROOF-CONDENSATE-FLOOR`'s transseries statement, confirmed here.
- **(L2) A reality / positivity constraint on the full transseries.** If the physical `F` must be the *unique* real resummation, one might hope to pin `σ`. But `PV + σ e^{−2s}` is real for *every* real `σ` (witness B), so reality is **one equation short** of pinning `σ`. Positivity of `⟨tr F²⟩` (if it held) is an **inequality** on `σ`, not an equality — it gives a *range* `σ ≥ σ_min`, not the value. (And `⟨tr F²⟩ > 0` is itself not proven.)
- **(L3) A resurgence triangle closing `σ` to a Stokes constant of a *different* sector.** This requires a *second* nonperturbative scale whose transseries parameter is independently known and which couples to the same Stokes constant. Pure YM at `θ=0` has no such second handle at the `u=2` scale (the `u=1` and `u=2` renormalons are separate, and neither is a physical saddle).

⇒ **All three loopholes need a nonperturbative input that pure YM does not supply.** Resurgence alone cannot close `σ`.

---

## §5 · Concrete numbers from the witness (`[V-num]`)

Toy Euler series, `x = 1/(2s)` so the nonperturbative scale is `e^{−2s}`:

| `x` | `s` | PV (median) | `Res = −e^{−1/x}` | `2π|Res|` (discontinuity) |
|---|---|---|---|---|
| 0.5 | 1.0 | 6.699e−01 | −1.353e−01 | 8.503e−01 |
| 0.1 | 5.0 | 1.131e−01 | −4.540e−05 | 2.853e−04 |
| 0.02 | 25.0 | 2.042e−02 | −1.929e−22 | 1.212e−21 |

- Discontinuity `= 2π|Res|` exactly (residue theorem) — the **invariant**.
- `σ ∈ {0, −0.5, 1.7, 10}` at `x=0.1` → `F ∈ {0.1131461, 0.1131234, 0.1132233, 0.1136001}`, all real — **`σ` free**.
- `a_n/n! ≡ 1` — the pole is read off perturbation theory; `σ` is not in `a_n`.

---

## §6 · Verdict and relation to the roof

- **`[T]`** Resurgence gives the exact, computable **dictionary**: the `u=2` renormalon ambiguity **is** the Stokes constant (Borel residue), readable from large-order perturbation theory.
- **`[NG]`** Resurgence-as-pinning: it does **not** fix the condensate `σ`, cancel the ambiguity, or yield a canonical *physical* value. The median sum is canonical but not physical; reality is one equation short; positivity (unproven) gives only a range.
- **Relative to `ROOF-CONDENSATE-FLOOR`:** the "not well-posed / `u=2` ambiguity = condensate's size" result is **confirmed and sharpened** — in resurgence language the obstruction is *an invariant* (`[T]`), and the condensate is *a free transseries parameter* (`[NG]`). The roof is still missing the same strut.
- **What would still be needed:** a **nonperturbative, scheme-independent characterization of the `0⁺⁺` anomaly** that fixes `σ` without passing through perturbation theory — i.e. exactly the "new species" (non-order-parameter, non-Markov, renormalon-free functional on `μ`) that `ROOF-SYNTHESIS` §4 named. Resurgence **relocates** the wall to its exact invariant form; it does not cross it.

**Clay 0%; no QED.** Labels: `[T]` dictionary (Stokes = residue, computable), `[NG]` pinning, `[V-num]` toy numbers, `[O]` = fix `σ` nonperturbatively (the same `H1`).

*Witness:* `research/unconditional/w7/roof_resurgence_witness.py` (runs clean, `py -3.12 -I`, numpy/scipy).
*Refs:* `ROOF-CONDENSATE-FLOOR.md`, `BASE-FIELDMEASURE.md`, `ROOF-SYNTHESIS.md`; standard resurgence (`Mariño`, `Aniceto–Schiappa–Vonk`, `Costin`).
