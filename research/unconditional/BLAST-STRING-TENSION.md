# BLAST-STRING-TENSION — can `σ` (the 1-form string tension) DERIVE `H1` (`χ ≤ 1/λ`)? A genuine `[C]` sufficient condition `λ ≤ 16σ/M` exists; it does NOT close at β≈2.5

> 2026-11-24 · final blast on `H1` from the string tension. **Verdict: `[C]` a real sufficient condition — the string tension derives `H1` **provided** the coarse Poincaré constant satisfies `λ ≤ 16σ/M` (equivalently `λ ≤ m₁²/M`), via the exact spectral squeeze + the strong-coupling bridge `m₁ ≥ 4√σ`. But this does NOT close at the crossover: (i) the bridge is 1-form→0-form and only rigorous at `β<1/48`; (ii) `M` (total 0⁺⁺ weight) is UV-sensitive there; (iii) `λ ≤ 16σ/M` is unproven and is a *separate* coarse-RG input, not given by `σ`. The bound is **tight at the atomless (gapless) configuration**, so it proves `H1` (the inverse susceptibility) — the *weaker* statement — and can never prove the atom. `[NG]` for a proof of `H1` at β=2.5. Clay 0%; no QED.**

Distinct from `TRANSMUTE-TENSION-GAP` (which killed `σ ⇒ m_{0⁺⁺}`, 2/5). Here the target is the **susceptibility** `H1`, not the mass.

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The exact squeeze (the only rigorous `σ → χ` link)

**Lemma G `[T]`** (PNT-CHI-CERTIFICATE; the spectral squeeze of the paper). For the 0⁺⁺ KL measure `ν` with atom `f₁²` at `m₁² = inf supp ν`, total weight `M = ∫dν`, susceptibility `χ = ∫s⁻¹dν`:
```
f₁² / m₁²  ≤  χ  ≤  M / m₁² .
```
The upper bound `χ ≤ M/m₁²` is the whole content for `H1`: it converts a mass scale into a susceptibility **upper** bound.

## §2 · The string-tension bridge

**Bridge `[C]`, β<1/48** (ROOF-CONFINEMENT-GAP): area law and the 0⁺⁺ gap share the polymer weight `u`, giving `m₁ ≥ 4√σ`. (At β>1/48 no rigorous bridge; TRANSMUTE-TENSION-GAP shows `σ ⇒ m_{0⁺⁺}` itself is open there. We use only the strong-coupling bridge.)

## §3 · The derived `[C]` sufficient condition for `H1`

Chain: `χ ≤ M/m₁²` (Lemma G) and `m₁² ≥ 16σ` (bridge) ⇒
```
χ  ≤  M / (16σ) .
```
`H1` is `χ ≤ 1/λ`. Since both are upper bounds on `χ`, `H1` follows **iff** the string-tension bound is at least as strong:
```
[Lemma G] ∧ [bridge m₁ ≥ 4√σ] ∧ ( λ ≤ 16σ/M )   ⟹   χ ≤ M/(16σ) ≤ 1/λ   ⟺   H1 .
```
The exact sharp form (no bridge needed) is `λ ≤ m₁²/M`. **This is a genuine `[C]` theorem: strong confinement with `λ` below the string-tension scale `16σ/M` implies `H1`.** It is non-vacuous (both sides are same-dimension `O(1)` comparisons; `m₁²/M` and `λ` both scale as one-over-mass²).

## §4 · Why it does NOT close at β≈2.5

- **(channel, 1-form vs 0-form)** `σ` is a 1-form datum (area law); `χ` is 0-form (trace-anomaly channel). `σ` reaches `χ` **only through `m₁`**, via the strong-coupling bridge. No direct channel link (TRANSMUTE §1: a point carries no 1-form charge; no Ward identity transports `σ` into `ν_s`).
- **(bridge regime)** `m₁ ≥ 4√σ` is proven only `β<1/48`. At β=2.5 the bridge is a measured identity (`c≈3.9`), not a theorem — `σ ⇒ m_{0⁺⁺}` is itself open (TRANSMUTE 2/5; `inf_L E_closed = 0`).
- **(M is UV-sensitive)** `M = ∫dν` is the full 0⁺⁺ spectral weight; in the continuum it is the (renormalized) coincidence limit — not uniformly bounded at the crossover.
- **(λ is independent)** even granting `m₁`, the criterion `λ ≤ 16σ/M` involves `λ`, the **coarse Poincaré constant** — an RG object `σ` does not touch. `σ ⇒ H1` holds **iff** `λ ≤ 16σ/M`; this is the "size of `σ` vs `1/χ`" question made precise: `H1` is obtained when `16σ/M ≥ λ`, i.e. when the string-tension scale **exceeds** the coarse gap. Large `σ` (strong confinement) helps; small `σ` relative to `λM` gives no `H1`.
- **(tight at the wall)** the atomless (bump-at-floor) measure **saturates** `χ → M/m₁²` (§2 witness). So the bound `χ ≤ M/m₁²` is *tight exactly at the gapless configuration*. `σ` therefore cannot distinguish gapped from gapless: it proves `H1` (inverse susceptibility, the weaker half) and **never the atom `f₁²>0`** — consistent with `H1 ⊊ gap`.

## §5 · Verdict

- `[T]` Lemma G (`χ ≤ M/m₁²`).
- `[C]` **`H1 ⟸ (β<1/48 bridge) ∧ (λ ≤ 16σ/M)`** — a genuine new sufficient condition: *strong confinement + coarse gap below the string-tension scale ⇒ `H1`.*
- `[NG]` this does **not** derive `H1` at β≈2.5 (bridge regime; `M`; `λ ≤ 16σ/M` unproven).
- `[O]` exact missing inputs: the bridge at 2.5 (1-form→0-form dynamical bridge) **and** `λ ≤ 16σ/M` (coarse Poincaré small vs string-tension scale).
- The route is **susceptibility-only** — `σ ⇒ H1`, never `σ ⇒ atom`.

**Clay 0%.** `σ` derives `H1` only under `λ ≤ 16σ/M`, which the crossover does not provide.

*Refs.* Lemma G: PNT-CHI-CERTIFICATE. Bridge: ROOF-CONFINEMENT-GAP (`m≥4√σ`, β<1/48). `σ⇒m` open: TRANSMUTE-TENSION-GAP. Channel: PROOF-1FORM-ANOMALY. H1 def: lean/YMStatement/Reduction.lean, PNT-CHI-CERTIFICATE.
