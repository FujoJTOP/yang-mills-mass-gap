# PNT-IR-ANALYTIC — the IR-only analytic functional EXISTS, and it is the **third renormalon-free reformulation**: the nearest Lee–Yang/Fisher zero of the **0⁺⁺-projected** `Z(τ)` in complex time gives `π/m₁` (acid-clean, atom-positive, renormalon-free) — but its evaluation is `H1`

> **VERDICT (no QED).** LEEYANG's failure was using the complex **coupling** `β=1/g²`, whose nearest
> singularities are the **kinematic `β=0`** branch points (acid-dirty: fires on the free field). The fix is
> to use the complex **time extent `τ`** at fixed real coupling, **projected to the 0⁺⁺ sector**:
> `Z_{0⁺⁺}(τ)=Tr_{0⁺⁺}e^{−τH}`. Its nearest Lee–Yang/Fisher zero in `τ` satisfies **`Im τ_zero = π/m₁`**
> iff an **isolated 0⁺⁺ atom** sits below the continuum. So the functional
> **`𝒜[μ] := π / Im(τ_zero)`** is (i) **acid-clean** (free massless and `:φ²:` have no atom zero), (ii)
> **atom-positive** (fires iff `f₁²>0`), (iii) **renormalon-free** (`π/m₁` is a physical pole — the same gain
> as the large-N mass and the dual Higgs). **BUT it is a re-encoding of the gap** (`𝒜` *is* `m₁`), and
> computing `Z_{0⁺⁺}(τ)` needs the **non-local 0⁺⁺ projector** and the **`a→0, L→∞` atom-isolation** — the
> transport wall (`H1`). **The type exists; the content is `H1`. Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** witness. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The construction

The 0⁺⁺ sector's partition function is a Boltzmann sum over the 0⁺⁺ spectrum `{E_n}` with weights `{g_n}`:

```
Z_{0⁺⁺}(τ) = Tr_{0⁺⁺} e^{−τH} = g₀ + Σ_n g_n e^{−τ E_n}      (vacuum E₀=0, gap E₁=m₁)
```

In the **complex time-extent plane** `τ ∈ ℂ0`, this has zeros. For a single isolated atom below the continuum,
keeping the two lightest levels `g₀ + g₁ e^{−τm₁}=0`:

```
τ_zero = (1/m₁)[ ln(g₁/g₀) − iπ(2k+1) ]        ⇒   Im τ_zero = π(2k+1)/m₁,  k∈ℤ
```

so the **nearest** zero (k=0) has `Im τ_zero = π/m₁`. **Define**

```
𝒜[μ] := π / Im(τ_zero)      (= m₁ = the 0⁺⁺ gap, when an atom is the lightest isolated state)
```

**Why the 0⁺⁺ projection is essential.** The *full* `Z(τ)` over all states has its nearest zeros set by the
**lightest excitations across all channels** — on a lattice these are **kinematic free/box modes** (the
LEEYANG acid problem, in the `τ` plane). The **0⁺⁺-projected** `Z_{0⁺⁺}` removes them; the projector
`P_{0⁺⁺}` is non-local, and that non-locality is exactly what buys acid-cleanliness.

---

## §2 · The three required properties (with witness)

Witness `w7/pnt_ir_analytic.py` (exact analytic + constrained scan, `py -3.12 -I`, clean):

```
(0) analytic pure atom g₀+g₁e^{−τm}=0 (g₁/g₀=2,m=1):  τ = (0.693, 3.142), (0.693, 9.425), …
(A) atom m₁=1.0 + continuum s₀=2.0:   Re τ=0.688  Im τ=3.140  π/Im=1.000  = m₁   ✓
(B) atom m₁=0.6 + continuum s₀=2.0:   Re τ=1.164  Im τ=5.234  π/Im=0.600  = m₁   ✓
(C) :φ²: (no atom, continuum s₀=4):   nearest |Z|≈1.00 — NO true zero; π/Im=5.48 has no atom meaning
(D) free massless (no atom, cont s₀=0): nearest |Z|≈0.99 — NO true zero
(E) atom weight g₁: 2.0→0.05 keeps Im τ≈π/1=3.14; g₁=0 (no atom) destroys the zero (|Z|≈1)
```

- **(i) Acid-clean `[T]/[V-num]`:** free massless (0⁺⁺ continuum from threshold 0) and `:φ²:` (continuum from
  `4m²`) have **only a cut, no isolated level below it** ⇒ **no zero at `π/m_atom`** ⇒ `𝒜` does not fire.
  Witness (C),(D): the nearest `|Z|` stays `≈1` (no zero). *Contrast LEEYANG: there `D=2.5>0` fired on the free field.*
- **(ii) Atom-positive `[T]`:** the isolated zero exists **iff** an atom sits below the continuum; as `g₁→0`
  it merges into the cut and disappears (witness E). So `𝒜` fires iff `f₁²>0` with isolation — exactly the
  KL atom, not a mere gap.
- **(iii) Renormalon-free `[T]`:** `m₁` is a **physical spectral pole**, so `π/m₁` is scheme-independent — **no
  `u=2` renormalon**. This is the **genuine gain over the gluon condensate** (`ROOF-CONDENSATE-FLOOR`: ambiguity
  = its own size), and it **matches** the logic of the two other renormalon-free routes (`ROOF-LARGEN-MASTER`,
  `ROOF-DUALSUPERCONDUCTOR`): *target a physical mass, not a condensate.*

---

## §3 · Why this is an IR-only fix of LEEYANG (and what LEEYANG got right)

| | LEEYANG | PNT-IR-ANALYTIC |
|---|---|---|
| analytic plane | complex **coupling** `β=1/g²` | complex **time** `τ` (fixed real `g²`) |
| nearest singularity | **kinematic `β=0`** branch point (every theory) | 0⁺⁺-projected Lee–Yang zero |
| acid test | **FAIL** (`D>0` on free field) | **PASS** (no zero on free/`:φ²:`) |
| what it detects | transfer-matrix positivity (universal) | the **0⁺⁺ atom** `f₁²>0` |
| renormalon-free | (n/a) | **YES** (`π/m₁` physical) |

LEEYANG's *type* (analytic-structure functional) was right; the plane and the projection were wrong. This
constructs the right plane (`τ`) and the right projection (0⁺⁺).

---

## §4 · Why it still does not close — the exact missing step `[O]`

`𝒜[μ] = π/Im(τ_zero)` is **a re-encoding of the gap** (`𝒜 ≡ m₁`); it is **not a new computable handle**.
Evaluating it requires both:

1. **The 0⁺⁺ projector `P_{0⁺⁺}`** — a non-local operator. The *full* `Z(τ)` zeros are kinematic (free/box
   modes) = acid-dirty; only the projected `Z_{0⁺⁺}` is clean, and the projection is non-local.
2. **The thermodynamic / continuum limit with the atom isolated below the continuum** — `L→∞`, `a→0`, and the
   isolation `inf supp ν_below < threshold`. In finite volume the "continuum" is *discrete*, so "atom vs
   continuum" is a **thermodynamic-limit** distinction; the clean zero `Im τ=π/m₁` with `m₁` the *infinite-volume
   atom mass* exists only in that limit.

Both are **exactly `H1`** (the transport / continuum wall), in the analytic-structure coordinate. So the
functional satisfies all three type-spec properties *as a definition*, but its **evaluation is the same wall**.

---

## §5 · Honest bottom line

- **YES, such an IR-only analytic functional exists** and meets all three requirements at the level of
  definition: `𝒜[μ] = π/Im(τ_zero)` of the 0⁺⁺-projected `Z(τ)` is **acid-clean, atom-positive, and
  renormalon-free** — a genuine, clean object that LEEYANG did not have.
- **It is the third renormalon-free reformulation** (with large-N master field and dual superconductor): *all
  three work by targeting a physical mass, avoiding the `u=2` condensate renormalon.* This convergence is the
  round's structural content.
- **But it is a re-encoding, not a new handle.** `𝒜 ≡ m₁`; extracting it needs the non-local 0⁺⁺ projector +
  `a→0, L→∞` atom-isolation = `H1`. It joins `A3` (Fredholm determinant) and the flow-response as a
  **renormalon-free but `H1`-cost** object.
- **What it adds:** a *clean, acid-passing, renormalon-free* statement of the strut — the sharpest form the
  analytic route can take — and a concrete numerical target (0⁺⁺-projected Lee–Yang zeros) if a lattice
  transfer-matrix computation in the 0⁺⁺ sector can be pushed.
- **Clay = 0%.** No QED.

*Files:* `PNT-IR-ANALYTIC.md`; witness `w7/pnt_ir_analytic.py`.
*Refs:* `INVENT-LEEYANG.md`, `ROOF-CONDENSATE-FLOOR.md`, `ROOF-LARGEN-MASTER.md`, `ROOF-DUALSUPERCONDUCTOR.md`, `INVENT-FLOWRESP.md`, `PROOF-F1SQ.md` (A3).
