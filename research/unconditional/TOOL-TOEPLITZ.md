# TOOL-TOEPLITZ — the Christoffel/Szegő (Toeplitz–Hankel) construction gives an UPPER bound only, no forcing lower bound exists (Δ_min=0), and the limit is A3: it does NOT escape the OPE-RHO blind spot

> **VERDICT (no QED).** Tested whether a non-integral Szegő/Christoffel lower bound on `f₁²` can
> escape the OPE-RHO moment-blindness. **It cannot.** Three exact facts:
> **(i)** the Christoffel function is an **UPPER** bound (`f₁² ≤ λ_n(m²)`, `λ_n ↓ f₁²`), not a lower one;
> **(ii)** both `λ_n` and the Gauss node-weight are **functions of finitely many moments** (`M_0..M_{2n-2}`)
> — i.e. **integral data**; **(iii)** **no lower bound on `f₁²` from finitely many moments can exist** — over
> the truncated-moment set the atom weight attains 0 (`SLOT-BOOTSTRAP` `Δ_min=0`). The Szegő "non-integral
> rewording" is the moment recurrence in disguise (Verblunsky = moments). The n→∞ limit reads the atom =
> **A3 (read-not-force)**; attaining it at controlled `n` needs the continuum = `H1`. **Same wall, Toeplitz
> coordinates.**

---

## §0 · Labels

**[T]** 真定理 · **[V-num]** new numeric · **[O]** 精确缺失步 · **[NG]** no-go. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The three facts (witness `w7/tool_toeplitz_witness.py`)

Measure A: atom `f₁²=0.5` at `s=1` + continuum `0.5` on `[4,6]` (gap 3).

### (i) Christoffel is an UPPER bound `[T]`
`λ_n(m²) = 1/K_n(m²,m²) = 1/(v^T H^{-1} v)`, `H_ij=M_{i+j}`, `v=(1,m²,…,m²{}^{n-1})`. Standard bound
`f₁² ≤ λ_n`, and `λ_n ↓ f₁²`. Witness: `λ_n(1) = 1.000, 0.510, 0.500, 0.500, …` — **decreasing to 0.5**.
So the Christoffel function **upper**-bounds the atom (it can *exclude*: if `λ_n` small, no atom; it cannot
*force*: `λ_n ≥ c` does not imply `f₁² ≥ c`).

### (ii) The Gauss node-weight also converges from ABOVE `[V-num]`
The n-node Gauss rule (`exact deg ≤ 2n−1`) has a leftmost node `x_1 → m²` and weight `w_1 → f₁²`. Witness:
`w_1 = 0.529, 0.5009, 0.50002, 0.500000, …` — **decreasing to 0.5 from above**, not increasing from below.
So no from-below ("forcing") sequence materialises in the natural Gauss construction either.

### (iii) No finite-moment lower bound can exist `[T]`
`λ_n` (and the Gauss rule) are determined by `M_0..M_{2n-2}` (`M_0..M_{2n-1}`) — **integral data**. Over the
truncated-moment set `{ν : M_k fixed, k≤2n}`, the atom weight `ν↦ν({m²})` is a linear functional ranging in
an interval whose **minimum is 0** (`SLOT-BOOTSTRAP`/`PNT-OPE-RHO` `Δ_min=0`: a continuum with the same
moments and no atom exists). Hence **no functional of finitely many moments can lower-bound `f₁²` by a
positive constant** — the atom can always be "hidden" into matched moments. This is the OPE-RHO obstruction,
restated for the atom weight.

---

## §2 · Acid-dirty at finite n `[V-num]`

Atomless continuum `[0.3,6]`: `λ_n(1) = 1.00, 0.371, 0.302, 0.300, …, 0.135` (`n=8`) — decreasing **slowly**,
not to 0. A genuine atom of `f₁²=0.134` gives `λ_n(1)→0.134`. So at `n=8` the **atomless continuum
(0.135) is indistinguishable from an atom of 0.135** — the finite-n Christoffel is acid-dirty (a gapless
continuum mimics a small atom), exactly the moment-blindness.

---

## §3 · Why the Szegő rewording does not escape

The "non-integral" route (Szegő recursion / Verblunsky coefficients on the unit circle / generating function
of `Π`) is the **moment recurrence in disguise**: the generating function's Taylor coefficients are the
moments; the Verblunsky coefficients are moment-derived. The analytic structure of `Π` (Stieltjes transform)
has the atom as a **pole residue** `f₁²=−Res Π` = **A3 (read-not-force)** — the same object, not a new one.
So the Toeplitz/Szegő construction lives entirely in the integral/A3 class.

---

## §4 · Honest bottom line

- **No forcing lower bound.** The Christoffel function is an upper bound; no finite-moment lower bound on
  `f₁²` exists (`Δ_min=0`); the Gauss weights come from above. The Toeplitz/Szegő tool does **not** manufacture
  a `f₁² ≥ c` certificate.
- **No escape from OPE-RHO.** The construction is moment-determined (integral class); the acid is dirty at
  finite `n` (atomless continuum mimics a small atom); the limit reads the atom = A3.
- **What it does give:** a clean, convergent **upper-bound** sequence `λ_n(m²) ↓ f₁²` and a Gauss rule
  `x_1→m²`, `w_1→f₁²` — useful as a *read* of an atom whose existence is already assumed (an A3-style
  diagnostic), not a forcing tool.
- **The blind spot remains.** Reading `f₁²` from finite data is blocked by the truncated-moment
  underdetermination (`Δ_min=0`); the only reader is the limit = A3 = `H1`. **Clay 0%.**

*Refs:* `w7/tool_toeplitz_witness.py`; `PNT-OPE-RHO.md` (the integral no-go), `SLOT-BOOTSTRAP.md`
(`Δ_min=0`), `WI-MODIFIED-FUNCTIONAL.md` (Re-slope = A3), `KEY-DIRECT-DERIVATION.md`.
