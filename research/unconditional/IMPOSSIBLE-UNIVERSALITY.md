# IMPOSSIBLE-UNIVERSALITY — the strongest universality match (4D SU(2) finite-T deconfinement ⟺ 3D Ising) is real but channel-mismatched and finite-T; zero-T YM has no critical point, so universality cannot force the 0⁺⁺ gap

**Impossible base.** Argue 4D YM at the crossover sits in a *gapped universality class*, matched to a known gapped system by an *exact universal quantity*, forcing the gap as a phase property. **Result: the strongest match is real but off-target; universality does not force the 0⁺⁺ gap.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric.

---

## §1 · The strongest genuine match `[T]`

4D SU(2) pure YM at **finite temperature** has a second-order **deconfinement transition** whose order parameter is the **Polyakov loop** (a `Z₂` magnetization) in the **3D Ising universality class** — `β=0.326`, `ν=0.630`. This is a genuine, exact universality match of a 4D gauge theory to a 3D spin model, and the strongest available. **But it is the finite-`T` transition, in the center/Polyakov (1-form) channel, and it disappears as `T→0`.**

---

## §2 · Why it does not force the 0⁺⁺ gap — the exact failure

1. **No critical point at zero `T` `[T]`.** 4D YM at `T=0` is a **smooth crossover** (no phase transition, `B2`) ⇒ there are **no critical exponents** ⇒ critical-exponent universality is **vacuous**. The Ising match exists only at `T_c`.
2. **Channel mismatch `[T]`.** The match is the Polyakov/center (1-form) order parameter; the 0⁺⁺ gap is the 0-form `tr F²` scalar. Gauge-invariant Hilbert-space decomposition makes them **Schur-orthogonal** (`BASE-ORDER-PARAMETER`, `ROOF-CONFINEMENT-GAP`): the universality fixes the *confinement transition*, not the 0⁺⁺ mass.
3. **The universal quantities that exist are gap-blind `[T]`.** The 4D **`a`-anomaly** (a-theorem) is universal and RG-monotone, but it is a **UV** datum: a free massless (gapless) and a free massive (gapped) scalar have the **same** `a=1/90`; a gapless CFT has `a≠0`. So `a` does not force a gap.
4. **Edge statistics presuppose the gap `[T]`.** Spectral-edge universality (RMT, Tracy–Widom) **requires an edge to exist** — circular; and the free density-of-states edge exponent `(D−2)/2` is **phase-space**, `D`-dependent (witness: `0, 1/2, 1` for `D=2,3,4`), **not** a gapped-class classifier.
5. **Large-N ratio `[C]`.** `m_{0⁺⁺}/√σ→3.30` is a universal *leading-order* ratio, but it **assumes** the gap — it does not force it.

---

## §3 · Verdict

**`[NG]`.** Universality does **not** force the 0⁺⁺ gap. The zero-`T` crossover is **genuinely non-universal** (smooth, no critical exponents); the one real universality match (SU(2) finite-`T` ⟺ 3D Ising) is a **confinement-transition** statement, **1-form**, **`T≠0`**, Schur-orthogonal to the 0⁺⁺ scalar; and the universal quantities that exist (`a`-anomaly, edge exponents, large-N ratio) are either UV/gap-blind or presuppose the gap. **The impossible base fails at exactly the criticality/channel step.** Genuine yield: the finite-`T` Ising match is a real theorem but off-target. **Clay 0%; no QED.**

*Witness:* `w7/impossible_universality.py`. *Refs:* `BASE-ORDER-PARAMETER.md`, `ROOF-CONFINEMENT-GAP.md`, `KAPPA-LARGEN.md`.
