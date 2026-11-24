# CONSTRUCT-DUAL-STABILIZED — the trapped vortex ring IS a genuine KL atom for `B>B_c`, but the **minimal** dual GL has `B<0` (Lüscher), so the atom still requires non-minimal input `= H1`

> **VERDICT (no QED).** Building the `+Bξ²/R`-stabilized dual ring settles the atom question cleanly.
> **(1) For any `B>0` the ring energy is strictly convex with a unique interior minimum `R*`, a *stable
> discrete ring* ⇒ a genuine δ in `ρ_{trF²}` ⇒ a real KL atom `f₁²>0`** (a theorem *within the model*).
> **(2) The atom does not survive `B→0` as a physical state**: `R*→(ξ/8)e^{3/4}≈0.265ξ<ξ`, below the
> thin-core validity; a physical atom (`R*>ξ`) needs `B>B_c≈1.33A`, an order-one trapping. **(3) The
> minimal dual GL / YM has NO positive `B`**: its only power-`1/R` term is the transverse-fluctuation
> (Lüscher) term `B_⊥=−(D−2)π/6<0`, which keeps the ring **monotone ⇒ no atom** (a neutral closed flux
> tube *shrinks*). Positive `B` is non-minimal input. **(4) Even in model B the atom weight `f₁²` is set
> by the same non-minimal `tr F²↔|φ|²` coupling**, so its magnitude is model-dependent. **The atom is
> real in the modified model but is `H1` in the dual coordinate** — the same pattern as the rigid-string
> master field (`m/√σ=2π√(2κ)`, `κ=H1`).

---

## §0 · Labels

**[T-model]** theorem within the explicit model · **[V-num]** new numeric · **[C]** conditional (a model,
not YM) · **[O]** exact missing step · **[NG]** no-go for the minimal class.
Witness: `w7/construct_dual_stabilized.py` (runs clean, numpy only, `py -3.12 -I`).

---

## §1 · The model and why `B>0` gives a genuine atom `[T-model]`

Units `ξ=1`, `A=1`. Trapped vortex-ring energy

```
E(R) = A·R·( ln(8R/ξ) − 7/4 ) + B·ξ²/R ,
E′(R) = A·( ln(8R/ξ) − 3/4 ) − B·ξ²/R² ,
E″(R) = A/R + 2B·ξ²/R³ .
```

For `B>0`: `E″>0` on `R>0` ⇒ `E` strictly convex ⇒ **at most one** critical point; and
`E′(R→0⁺)→−∞`, `E′(R→∞)→+∞` ⇒ **exactly one** root `R*`, a stable minimum (`E″(R*)>0`). A stable
discrete ring of definite radius `R*` is an **isolated level** `m*=E(R*)` ⇒ a **δ** in the spectral
density ⇒ **`f₁²=ν({m*²})>0`**. (Witness §1.)

| `B` | `R*` | `R*>ξ`? | `m*=E(R*)` | `m*/√σ` | atom? |
|---|---|---|---|---|---|
| 0.5 | 0.711 | no | 0.695 | 0.695 | no (R*<ξ) |
| 1.0 | 0.903 | no | 1.313 | 1.313 | no |
| 1.329 | 1.000 | — | 1.659 | 1.659 | threshold |
| 2.0 | 1.163 | yes | 2.278 | 2.278 | **YES** |
| 3.0 | 1.355 | yes | 3.072 | 3.072 | **YES** |
| 5.0 | 1.652 | yes | 4.400 | 4.400 | **YES** |
| 10.0 | 2.178 | yes | 7.004 | 7.004 | **YES** |

So **the atom is a theorem of the model for `B>0`** — there is a real, isolated, stable ring. Reaching the
physical ratio `m*/√σ≈3.9` needs `B≈4.3` (interpolating `B=3→3.07`, `B=5→4.40`).

---

## §2 · It does NOT survive `B→0` as a physical atom `[V-num]`

As `B→0⁺`, `R*` retreats: `R*(B=0.5)=0.711ξ`, `R*(0.2)=0.534ξ`, `R*(0.05)=0.377ξ`,
`R*(0.01)=0.297ξ`, limit `R*→(ξ/8)e^{3/4}=0.2646ξ`. **The minimum moves *below* the core radius `ξ`**,
where the thin-core ring energy is invalid. A *physical* atom (`R*>ξ`) requires

```
B > B_c = A·(ln8 − 3/4) = 1.3294·A  ≈ 1.33 A      (order-one trapping).
```

⇒ The trapped-ring atom is **not a small-`B` limit of a physical state**; it needs `B` of order `A`
(i.e. of order the line tension `σ`). It is a genuine atom only in a **non-minimal** model with an
`O(σ)` trapping term.

---

## §3 · The minimal dual GL / YM has `B<0` (no atom) `[NG]`

The minimal (London-limit, massless-transverse) dual-string ring energy is

```
E_ring(R) = 2πR·σ + (curvature/log) + c/R + …
```

The **only** power-`1/R` term of the minimal theory is the free transverse-fluctuation (Lüscher) term

```
B_⊥ = −(D−2)π/6 = −(4−2)π/6 = −1.047   (D=4)   STEP [T],
```

which is **negative**. A negative `B` makes `E′(R)=A(ln(8R/ξ)−3/4)+|B|ξ²/R²` **even more positive** on
`R≥ξ` ⇒ still strictly monotone ⇒ **no interior minimum ⇒ no atom**. Physically: **a neutral closed flux
tube shrinks** — there is no stable neutral ring in minimal YM. (A positive contribution could come only
from *heavy* transverse modes, `∝ e^{−2m_V R}/R`, exponentially small for `R>ξ=1/m_H` in the relevant
`κ_GL≫1` regime, and not a clean stabilizing `+1/R`.) ⇒ **`B_minimal ≤ 0`: minimal YM has no atom.**
(Witness §3.)

---

## §4 · Even in model B, the atom *weight* `f₁²` is model-dependent `[O]`

In model B the atom is a genuine δ at `m*=E(R*)` with weight `f₁²=w(R*)²`, where `w(R)=⟨0|tr F²|ring(R)⟩`
is the `tr F²` overlap of the ring. `w` is set by the **same non-minimal `tr F²↔|φ|²` coupling** that
would supply a positive `B`; the minimal AHM fixes *masses*, not this overlap (exactly the
`PNT-DUAL-CONDENSATE` obstruction). So **the atom's mass is definite but its weight is not fixed by the
minimal theory** — `f₁²` magnitude is model-dependent. (Witness §4.)

---

## §5 · Verdict

- **Atom real?** **Yes, in the trapped model (`B>B_c≈1.33A`)** — a genuine isolated stable ring ⇒ δ in
  `ρ_{trF²}` ⇒ `f₁²>0` `[T-model]`.
- **Survives `B→0`?** **No** — `R*` drops below `ξ` (0.265ξ), outside thin-core validity; a physical atom
  needs order-one `B` `[V-num]`.
- **Does minimal dual GL have `B`?** **No positive `B`** — only the Lüscher term `B_⊥=−1.047<0` ⇒ no atom;
  the neutral tube shrinks `[NG]`.
- **Does it advance the wall?** It **converts** "no atom in the dual picture" into an **explicit
  one-parameter family** `m*/√σ=f(B)` with a threshold `B_c`: the atom exists iff `B>B_c`. This is the
  dual-coordinate analogue of the rigid-string `m/√σ=2π√(2κ)` — **the missing order-one input
  (`B`, resp. `κ`) is `H1`**. The minimal theory gives the *wrong sign* (`B_⊥<0`); supplying a positive
  `B` and the overlap `w` is the exact missing step.

**Clay 0%; no QED.** Labels [T-model]/[V-num]/[C]/[O]/[NG].

*Refs:* `CONSTRUCT-DUAL-LATTICE.md` (continuum baseline), `PNT-DUAL-CONDENSATE.md`,
`ROOF-DUALSUPERCONDUCTOR.md`, `CONSTRUCT-MASTERFIELD.md` (rigid string, the same one-parameter pattern),
`PNT-CHI-CERTIFICATE.md` (H1 as the certificate). Witness: `w7/construct_dual_stabilized.py`.
