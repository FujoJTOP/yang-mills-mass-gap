# FALSIFY-WARD-NONLINEAR — nonlinear combinations of the exact YM identities do NOT force the atom; the only-exit claim (YM-structure or H1) STANDS

> **VERDICT: NOT falsified (this route).** I tried to break the "only-exit = YM-structure or `H1`" claim with a
> **nonlinear combination** of the exact identities — Ward `κχ=1−η`, trace anomaly `A=∫s⁻²dν`, positivity, AF
> tail — that might exclude the atomless measure. **It does not.** The structural reason is exact: every such
> identity is a **moment / integral** datum (Ward is a tautology — it *defines* `η:=1−κχ`, SLOT-WARD-STEIN), and
> an **atomless flat bump** `ρ_h = f₁²/h` on `[E₁,E₁+h]` matches **every moment** (all `k`, incl. `M₀, χ, A`) to
> `O(h)` while having `f₁² = 0`. Any nonlinear combination is a **continuous function of these matched moments**,
> so it takes the **same value** for the atom and the atomless bump. Witness: `R_CS=χ²/(M₀A)`, `detH=M₀A−χ²`,
> `R_inv`, `R_sq=M₁M₋₁/M₀²` all converge to the atom's values as `h→0` (max `|ΔM_k| ~ h`). **No forcing.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The atomless adversary (acid-tested, valid KL measure)

Take the atom-bearing measure `ν_atom = f₁² δ_{E₁} + g·1_{[E₂,E₂+L]}` and replace the floor atom by a flat
**bump** `ρ_h = f₁²/h` on `[E₁, E₁+h]` (same total mass, same floor `E₁`, continuum unchanged). Then:
- **`M₀ = ∫dν`, `χ = ∫s⁻¹dν`, `A = ∫s⁻²dν`, and every `M_k = ∫s^k dν`** match the atom's to `O(h)` — because
  `(f₁²/h)∫_{E₁}^{E₁+h} s^k ds → f₁² E₁^k =` the atom's moment (the `k=−1` log case likewise).
- **positivity** `ρ_h ≥ 0` ✓; **floor** `inf supp = E₁` ✓; **`f₁² = 0`** (no atom at the floor) ✓;
- **Ward** `κχ=1−η` holds by defining `η:=1−κχ` (the identity constrains nothing — SLOT-WARD-STEIN); **AF tail**
  untouched (bump at low `s`).
⇒ The bump is a **valid atomless KL measure** satisfying **all** the exact identities, with `f₁²=0`.

---

## §2 · Nonlinear combinations fail (witness `w7/falsify_ward_nonlinear.py`, clean)

The strongest nonlinear candidates (all functions of the moment data, all from positivity):
- Cauchy–Schwarz `R_CS = χ²/(M₀A) ≤ 1` (i.e. `M₁² ≤ M₀M₂`), `detH = M₀A − χ² ≥ 0` (2×2 Hankel PSD),
  `R_inv = A M₀/χ² ≥ 1`, `R_sq = M₁M₋₁/M₀² ≥ 1`, and higher 3×3 Hankel PSD.

Witness (atom `f₁²=0.5` at `E₁=1` + uniform continuum): as `h→0`, `R_CS→0.730382`, `detH→0.446667`,
`R_inv→1.369146`, `R_sq→1.278053` — **identical** for atom and atomless bump (max `|ΔM_k|` falls as `1.0×h`),
while `f₁²_bump=0`. **Every such combination is blind** by construction: it is a continuous function of matched
moments.

---

## §3 · Structural reason (and the exact escape)

A nonlinear combination of the exact identities is a **continuous function of finitely many moments**. The
atomless bump matches **all** moments to `O(h)` with `f₁²=0`; hence **any** continuous function of them is
matched ⇒ cannot force `f₁²>0`. To distinguish, one needs a **non-continuous** functional — a **singular-kernel**
integral (e.g. the `L²` route `N(ε)=∫|G|²`, TOOL-HARDY) — but that is **not** a combination of the exact data,
and even it yields *presence*, not the *weight* (weight = A3). Equivalently, excluding the bump requires the
**YM-specific admissibility** (`ν` is the KL measure of `tr F²` in *interacting* YM, not an arbitrary positive
measure) — which **is** the claimed escape.

**Conclusion.** The "only-exit = YM-structure or `H1`" claim is **confirmed, not falsified**, on the
nonlinear-exact-identity route. This does **not** prove no YM-structure tool works (that is exactly the open
escape) — it proves the exact **integral** identities, however combined, are continuum-compatible.

*Labels:* [T] moments are integrals / Ward tautological; [V-num] the convergence table; [NG] nonlinear
moment-combinations as atom-forcers; [O] escape = YM-specific admissibility or a singular-kernel functional =
`H1`. **Clay 0%; no QED.**
