# TOOL-NEVANLINNA — Nevanlinna–Pick / Carathéodory interpolation does NOT give a non-integral forcing lower bound; it splits into the two known blind camps

> **VERDICT `[NG]` for the missing tool.** Complex-function-theory interpolation is a genuine, rigorous
> handle — but it does **not** fill the blind spot. On point **values** it is a finite **integral**
> functional (OPE-RHO applies → **upper** bounds only, acid-dirty); on **zeros** it is the **A3**
> read-not-force class. The blind spot (a non-integral **forcing** lower-bound atom detector) is
> *robustly* closed within this family too. **Clay 0%; no QED.**

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

## §1 · Construction & witness (`w7/tool_nevanlinna.py`, runs `py -3.12 -I`)

`Π(z)=∫dν/(s−z)` is a Stieltjes / Pick function (`Im Π/Im z ≥ 0`); the atom is `f₁²=−Res_{s=1}Π`.
Nevanlinna–Pick on point values `Π(z_j)` gives the feasible set of measures matching them.

**(1) An ATOMLESS bump mimics the atom.** Truth = atom `0.5` at `s=1` + continuum on `[2,10]`. Replace
the atom by a narrow *atomless* bump near `s=1` (total mass 0.5, width `h`). Its point evals match:

| `h` | `max |ΔΠ|` over 4 complex `z_j` |
|---|---|
| 1e−2 | 1.33e−4 |
| 1e−3 | 1.33e−6 |
| 1e−5 | 1.33e−10 |

⇒ as `h→0` an **atomless** measure reproduces the atom's `Π(z_j)` exactly ⇒ **`f₁²=0` is feasible**.

**(2) LP max atom weight** from the same evals: `f₁²_max = 0.5000` (= the true atom weight). So the
NP-feasible range is `f₁² ∈ [0, 0.5000]` — **includes 0** ⇒ **NOT forcing**; NP gives only an **upper** bound.

## §2 · Why it fails (structural, not technical)

Point evals `Π(z_j)=∫dν/(s−z_j)` are finite measurable **integral** functionals (bounded kernels
`1/(s−z_j)`). **OPE-RHO's theorem applies verbatim**: an atomless bump at the atom location matches every
`Π(z_j)`, and more generally a floor-bump of weight `CΔ²` matches with `Δ→0`. Hence no finite point-eval
bound can **lower**-bound `f₁²`: NP yields `Δ_max` (upper), never `Δ_min` (lower). **Acid:** it returns
`f₁²∈[0,f_max]` for an atomless truth too ⇒ non-discriminating ⇒ **acid-dirty**.

## §3 · The non-integral variant is A3

Interpolating **zeros** of `Π` (or Padé-pole positions) — not values — *is* non-integral (zero location),
but it is exactly the **WI/𝒜 class = A3 (read-not-force)**; already shown to fail for atom+continuum
(`WI-EMBEDDED`) unless reweighted into the Stieltjes residue (=A3 again).

## §4 · Classification (the `[T]` yield)

Complex-analysis interpolation splits cleanly into the two known camps, **neither** filling the blind spot:
- **VALUES → integral class** (OPE-RHO; upper bounds; acid-dirty);
- **ZEROS → A3** (read-not-force).

**The blind spot — a non-integral, forcing (lower-bound), acid-clean atom detector — is not filled by
complex-function-theory inequalities; it is robustly closed here as well.** This *narrows* the search: an
escape cannot be a standard analytic-function inequality on `Π`'s values.

*Refs:* `PNT-OPE-RHO.md` (integral no-go), `WI-ESCAPE.md`/`WI-EMBEDDED.md` (zeros=A3), `A3`/`PROOF-F1SQ.md`.
