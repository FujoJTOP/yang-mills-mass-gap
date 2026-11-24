# V3-WAVE-SYNTHESIS — 3 no-gos, 1 data correction, 2 structural advances: the monotone-rail and N_f shortcuts are dead; only the (★) conditioned-variance bound survives

> **VERDICT (no QED, Clay 0%).** Six agents attacked the three literature-identified surfaces,
> two angles each. **A and B are dead as proof routes** (three clean no-gos). **C survives, but
> only as a reduction**: the Balaban step is no longer blocked by volume growth, and SZZ's
> constants are polynomial in `N`; what remains is one named inequality `(★)`. The wave also
> forces a correction to the campaign's own "one species" picture: the UV and IR faces have the
> **same form but are not the same object** (TCI is IR-native and does not close L3). Clay is still
> 0%.

---

## §0 · Labels

`[T]` theorem · `[C]` conditional · `[NG]` no-go · `[O]` exact missing step · `[V]` numerical.
All six witnesses run under `py -3.12` (numpy 2.4.6), exit 0.

---

## §1 · The six verdicts

| Route | Verdict | The load-bearing statement |
|---|---|---|
| **A1** `a`-anomaly → residue | `[NG]` | `a` is a **3-point** `⟨TTT⟩` coefficient; `f₁²` is a **2-point** `⟨ΘΘ⟩` atom. No kernel connects the channels. `Δa ≥ c·f₁²` only reads as `f₁² ≤ Δa/c`, an **upper** bound — the wrong direction. `a_IR ≥ F(f₁²,m₁)` is false (gapped IR has `a_IR=0`). Acid-vacuous on free field and `:φ²:`. |
| **A2** condensate → residue | `[NG]` | The only `u=2`-ambiguity-cancelling object is the full correlator `Π=∫dν/(Q²+s)`, which is **linear in ν** ⇒ atom-blind (FALSIFY-OPE). `G2 = c·M₋₁+subtr.`, atom term `f₁²/m₁²`. Numerically `m₁⁴/G2 ≈ 696`: a natural residue overshoots by ~700×, so a match needs `~1.4×10⁻³` of natural — renormalon near-cancellation, not atom dominance. Residual `[O]` = Stieltjes inversion = H1. |
| **B1** `N_f`-monotonicity | `[NG]` | Proved `[T]`: `β₀,β₁` strictly decreasing, sign change at `306/38`, `a` strictly decreasing, `β(g)` pointwise increasing. But **any order parameter vanishes on the open conformal window**, so it is constant on an interval and cannot be strictly monotone; non-strict monotonicity gives only `O(0)≥0` (vacuous). A rigorous BZ IRFP is itself only `[C]`. Counterexample: LSD finds `⟨ψ̄ψ⟩/F³` enhanced at `N_f=6`. The missing lemma is **impossible**. |
| **B2** lattice `N_f` data | `[V]` | 8 values have usable data. Monotonicity is **coarse only** (gapped 0–8, conformal 16, edge bracketed). `γ_m` is not monotone (0.45–1.0 at `N_f=12`). Sharpest: "`N_f=10` conformal (2023)" and "`N_f=12` χSB (Fodor)" are the **unique non-monotone assignment**; if the edge is unique/monotone, they cannot both hold. Window sentence corrected (see §4). |
| **C1** SZZ `N`-scaling | `[T]`(+`[V]`) | SZZ's Bakry–Émery constant is **exactly linear in `N`**: in `d=4`, `K_S = N(1/2 − 24|β|)` (verified against the paper's source: `(N+2)/2 − 1 − 8(d−1)N|β|`). The strong-coupling radius `|β|<1/(16(d−1))=1/48` is **`N`-independent**. The mass-gap rate is polynomial (`c_N ≍ N^{-3}` from the proof's constants; witness ratio `2^{-3}`), not exponential. **The large-`N` route is alive**; the wall is the `N`-independent `β`-transport/`a→0`, not `N`-degradation. |
| **C2** `(a,Λ)`-uniform Balaban | `[T]`+`[O]` | Balaban's theorem is `a`-uniform at fixed volume, at the level of the **action**, not the **force Hessian** L3 needs. **Lemma B `[T]`:** the large-field entropy channel is volume-free per unit volume; convergence is the constant inequality `c > c* = d ln L/(2β₀)` (witness `c*≈20–30`). "Intrinsic volume growth" is refuted; the real `Λ`-dependence is the IR tower bottom = the gap. Reflection positivity is **not** the wall (exact RG preserves the OS cone; GL4 transports RP from lattice measures). **TCI/two-block is IR-native and does not close L3 — same form, different object.** Residual = `(★)∧(★dec)`, exactly L3. |

---

## §2 · What this does to the figure at the centre of the paper

The v2 figure (and Chapter 5) presents one inequality with two faces (H1 and L3) and a dashed
"same shape". The wave sharpens the status of that dashed arrow:

1. **A1 and A2 kill the hope of a monotone rail.** The only exact monotone available (`a`) is in the
   wrong channel; the condensate is the right channel but linear in `ν` and numerically 700× too
   small to be dominated by the atom. So the transport machine still has **no rail**.
2. **B1 kills the `N_f` shortcut.** No strictly monotone universal order parameter can exist,
   because any order parameter vanishes identically on the open conformal window.
3. **C2 shows "same form" ≠ "same object".** TCI is IR-native; it does not close L3. The
   identification `L3 ⇔ H1` is therefore **more conjectural than v2's figure suggests**, not less.
4. **C1 and C2 leave exactly one live target:** a `Λ`-tempered, per-shell, one-sided upper bound on
   the conditioned force variance,
   \[
   (\star)\quad
   \operatorname{Var}_{\mu(\cdot|U')}\bigl(\langle P^{\mathsf T}\nabla S_j,v\rangle\bigr)
   \le E\bigl[\langle P^{\mathsf T}H_jP\,v,v\rangle\bigr]+W_j(v)-\delta_W ,
   \]
   uniform in `j,U',a` (with `(\star dec)` for the decimated version). This is the same species as
   H1, but the wave established that proving the IR instance does **not** automatically give the UV
   instance. The "one species" claim must be rewritten as **"one form, two instances; the
   dictionary is the open problem."**

---

## §3 · The updated attack surface (only one left)

**C — `(★)`.** The wave reduced the UV side to:
- **Lemma B `[T]`** (done): the large-field count is volume-free per unit volume; no intrinsic
  exponential volume growth.
- **`(★)`** (open): the conditioned force-variance bound, per shell, `Λ`-tempered, `N`-uniform.
- **`(★dec)`** (open): the decimated version needed for the RG map.
- **RP** (not the wall): exact RG preserves the OS cone; GL4 transports RP.

The honest next move is to attack `(★)` directly, with Lemma B as the base. A and B are no-gos and
should be recorded as such, not retried.

**A2 residual.** The full-correlator route is not a new proof route: it is exactly Stieltjes
inversion of the full spectral density, i.e. H1. Record it as an equivalence, not a tool.

**Infrastructure.** Formalize the acid test (free field has no gap) and the corrected figure claim
("same form, identity open") in Lean. Useful, not progress on the gap.

---

## §4 · Corrections to the campaign's own claims

1. **"One species" wording.** v2 says L3 and H1 are "one species". Correct to: *same form; the
   identification is a conjecture; TCI is IR-native and does not close L3* (C2). The figure's
   dashed arrow should be labelled `conjecture`, and H1 should carry the visual weight.
2. **Residue-floor direction.** Never use `Δa ≥ c·f₁²` as a floor; it is an upper bound
   `f₁² ≤ Δa/c` (A1). The useful reverse inequality `Δa ≤ c′f₁²` is not available from the KS
   integral (UV-sensitive).
3. **Condensate as residue probe.** The condensate is not an O(1) probe of `f₁²`; the atom can be
   the lightest IR feature but cannot dominate `G2` (A2: `m₁⁴/G2≈696`).
4. **Conformal-window sentence.** Replace "9≤N_f≤16 is gapless-conformal" with: *rigorously only
   the two-loop statement holds (BZ zero for `N_f>306/38`, AF for `N_f<16.5`); lattice supports an
   IRFP at `N_f=16` and χSB for `N_f≤8`; `N_f=12` is contested; the boundary is bracketed, not
   located* (B2). Note the `N_f=10` vs `N_f=12` inconsistency.
5. **SZZ `N`-scaling.** `K_S` is linear in `N` and the radius is `N`-independent `[T]`; the gap
   *rate* is polynomial (`N^{-3}` from the proof's constants, non-optimal) — do not state an
   `N`-uniform gap constant, but do state that there is no exponential `N`-degradation (C1).
6. **L3 vs Balaban.** Balaban is `a`-uniform at fixed volume for the action, not the force Hessian;
   the volume-growth objection is removed by Lemma B, but `(★)` remains (C2).
7. **Sussman–Alpöge 2026** is `[R]`: unrefereed, conditional, gap excluded (C2). Do not cite as a
   theorem.

---

## §5 · What is now dead (do not retry)

- **`a`-anomaly → `f₁²` bridge** (A1, `[NG]`): wrong channel, wrong direction, acid-vacuous.
- **Condensate → `f₁²` finite sum rule** (A2, `[NG]`): the only renormalon-clean object is linear
  in `ν`; numerically the atom cannot dominate.
- **`N_f`-monotonicity theorem** (B1, `[NG]`): no strictly monotone universal order parameter can
  exist; non-strict is vacuous.
- **"Volume growth is intrinsic to Balaban"** (C2, refuted): Lemma B shows the entropy channel is
  volume-free per unit volume.
- **"TCI closes L3"** (C2, refuted): same form, different object.
- **"`N`-degradation blocks the large-`N` route"** (C1, refuted): `K_S` grows linearly in `N`.

---

## §6 · What survives

1. **`(★)`** — the conditioned force-variance bound, per shell, `Λ`-tempered, `N`-uniform. This is
   the only live route among the three, and it is the transport machine's UV half.
2. **Lemma B** — a genuine `[T]` base for `(★)`; the large-field count is volume-free.
3. **SZZ `N`-scaling** — the large-`N` route is not blocked; `K_S` is linear in `N`.
4. **B2's data table** — the corrected conformal-window statement and the `N_f=10`/`N_f=12`
   inconsistency, as a concrete lattice target.

---

## §7 · Honest bottom line

The wave did what the literature map predicted: the two "shortcut" surfaces (a monotone, an
`N_f` inequality) are **provably dead**, and the one structural surface (the UV conditioned-variance
bound) is **alive but reduced** to a single named inequality. The campaign's own "one species"
picture is corrected: same form, not the same object, until `(★)` is proved at both scales.

**Clay = 0%. No QED.** The wall did not move; the map around it got sharper, and one of its three
doors is now the only one left open.

*Refs:* `A1-SUMRULE-A.md`, `A2-SUMRULE-CONDENSATE.md`, `B1-NF-MONOTONE.md`, `B2-NF-LATTICE.md`,
`C1-SZZ-NSCALING.md`, `C2-BALABAN-UNIFORM.md`; witnesses in `w7/`; `../lit-v3/V3-LIT-SYNTHESIS.md`.
