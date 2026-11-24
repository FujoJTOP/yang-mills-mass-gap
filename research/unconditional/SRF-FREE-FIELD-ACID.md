# SRF-FREE-FIELD-ACID — the decisive free-Gaussian acid test of SRF 2.0 (Soft-mode Retaining Feshbach Flow): **PASS on the faithful global-V reduction — inf spec coarse is EXACTLY 0 (massless) and m² (massive)** — and the test's teeth: the *local single-block* criterion spuriously gaps by `2−2cos(π/b)>0` (vacuous), the `"52"` fiber inequality is satisfied at m=0 and hence free-field-blind, and no L-uniform free gap exists

> 2026-11-23 · Acid test of SRF 2.0 against the iron rule of the charter: **a gap tool MUST fail on the free massless field; if it produces a strictly-positive soft gap there it has secretly deleted the soft mode and is vacuous** (`docs/BRIEF-FOR-MODEL.md` §6; `CLAUDE.md` §1). SRF 2.0 = block-average coarse `V` (soft modes RETAINED) + fiber `η∈ker C` + exact Feshbach/Schur `S(z)=A − z − Σ(z)`, `Σ(z)=W(D−z)⁻¹Wᵀ ⪰ 0`, target = a "soft-spectrum lower bound" on `V` (the "SRF-CUBE local criterion"). **判词: PASS — but ADMISSIBLE ≠ USEFUL.** The faithful reduction is *exactly* gapless on the free field (inf spec coarse `= 0`, machine-zero `≤1e-16`; massive `= m²`, exactly) because the k=0 mode is block-constant (∈V), is annihilated by the fiber coupling `Wᵀ` and hence by `Σ(z)` for all `z` and **all sign conventions** — the self-energy cannot touch the massless soft mode at all. The −Σ term (correct sign) only *lowers* the non-soft coarse energies (measured `{2/3,2/3,4/3} → {0.316,0.316,1.09}`); it does **not** manufacture a gap. The test's teeth bite two *variant* implementations: (i) a **local single-block** SRF-CUBE (freezing the inter-block k=0 mode) spuriously certifies `2−2cos(π/b)>0` on the free massless field → that incarnation is **VACUOUS**; (ii) the **"52" fiber-Hessian inequality** (hexagonmath `2610.00341`, **unreviewed**) is *satisfied by the free massless field* (`γ_fib>0` already at `m=0`), so fiber-control is blind to the gap — sharp free constant is `1/γ_fib = O(1)`, not 52. Passing is **necessary, not sufficient**, and earns *zero* positive evidence: exactness forces it. The load-bearing unknown — an L-uniform, non-loop, non-circular bound crossing the Savvidy crossover — is untouched and remains **H1**.

---

## §0. What the acid test decides — and the two readings that split the verdict

The charter's free-field test has one job: a proposed gap tool must **produce no gap** on the free massless field, and must **not** claim a positive L-uniform floor that is false there (`BRIEF-FOR-MODEL.md` §2 item 4, §6). Applied to SRF 2.0 there are exactly two readings of "the coarse `V` is what SRF-CUBE bounds," and they get **opposite** verdicts:

- **Reading G (global-V, faithful).** The coarse space `V` is the **full** block-average space (all block-constant functions, *including the inter-block k=0 constant mode*), and SRF-CUBE bounds `inf spec S(z)` over that whole space. → the k=0 soft mode is **in** `V`. Verdict below: **PASS** (exactly gapless).
- **Reading L (local single-block).** SRF-CUBE is a *per-block* / single-cell Hessian positivity criterion (the literal "local criterion"), i.e. the coarse energy is read off **one isolated block** with its interior `k→0` motion frozen. → the k=0 mode is **not** in the single block and is **deleted**. Verdict below: **FAIL / vacuous** (spurious gap `2−2cos(π/b)>0`).

The whole point of running the test is to *separate* these. The framework is admissible **iff** it commits to Reading G. Below: G in §§1–5, L in §6.

Throughout, **[真定理]** = exact free-Gaussian fact (per-block- and dimension-independent unless noted); **[no-go]** = a variant provably fails the acid test; **[缺失步]** = the untouched load-bearing step.

---

## §1. Setup: free Gaussian, block-average map, exact orthogonal split

Free scalar on a periodic lattice (spacing `a=1`), `d`-dimensional, sites `x`, field `φ∈ℝ^N`, `N=L^d`:

```
[M φ]_x = (−Δ φ)_x + m² φ_x = Σ_μ (2φ_x − φ_{x−e_μ} − φ_{x+e_μ}) + m² φ_x ,   φ periodic.
```
`M` is the Hessian; its spectrum *is* the free physical spectrum. Fourier: `λ(k)=Σ_μ 2(1−cos 2πk_μ/L)+m²`, so `inf spec M = m²` (constant mode `e₀ = N^{−1/2}(1,…,1)`) and the first nonzero rung is `2(1−cos 2π/L)=4 sin²(π/L)` (one direction excited).

**Block-average coarse map `C`** (`d=1`, block size `b`; `d>1` blockwise): `(Cφ)_I = b^{−1}Σ_{x∈block_I} φ_x`. Its row space is the **block-constant** subspace; its kernel is the **fiber**. These are orthogonal complements, giving the SRF split

```
ℝ^N = V ⊕ F ,   V = row(C) = {block-constant},   F = ker(C) = {block-sum-zero (fiber η)} ,
dim V = n = N/bᵈ ,   dim F = N − n ,   P_V = Cᵀ(CCᵀ)⁻¹C ,   P_F = 1 − P_V .
```
In any orthonormal basis `U_V (N×n)`, `U_F (N×(N−n))` of the two spaces, the Hessian blocks are

```
M = [[A, W],[Wᵀ, D]] ,   A = U_Vᵀ M U_V (n×n) ,  D = U_Fᵀ M U_F (f×f) ,  W = U_Vᵀ M U_F (n×f) ,
S(z) = A − z − Σ(z) ,    Σ(z) = W (D − z)⁻¹ Wᵀ        (Feshbach/Schur, exact) .
```
**Exactness [真定理].** `det(M−z) = det(D−z)·det S(z)` (Schur; verified numerically to `8.9e-16`, §table). So the physical eigenvalues are exactly the fiber eigenvalues (`det(D−z)=0`) together with the roots of `det S(z)=0`. The map is **lossless** — this is the entire reason a faithful SRF passes the free test trivially.

> **Note on the cited inequality.** SRF 2.0 cites a block-average fiber Hessian inequality `‖η‖² ≤ 52‖dη‖²` from hexagonmath `2610.00341`. The project has already graded that preprint **"AI-heavy, venue suspect, only a lead, not a basis"** and flags a co-cited `2602.0054` as **retracted** (`research/physics_handoff.md` §风险 4). We **do not rely on it** and recompute the free fiber gap in §5.

---

## §2. FREE MASSLESS (`m=0`) — the decisive case: inf spec coarse `= 0` EXACTLY, the soft mode survives

**[真定理 — the decisive result, all `d`, all block sizes.]** At `m=0` the constant map `e₀` is block-constant so `e₀∈V`; and `M e₀ = 0`. Therefore

```
Wᵀ e₀ = P_F M e₀ = P_F 0 = 0            ⇒   e₀ ∈ ker(Wᵀ)  ⇒   Σ(z) e₀ = W(D−z)⁻¹Wᵀ e₀ = 0  for EVERY z,
A e₀ = P_V M e₀ = 0                       ⇒   S(0) e₀ = (A − 0 − Σ(0)) e₀ = 0 .
⇒  inf spec S(z)=0-problem  =  0   EXACTLY  (e₀ is a genuine eigenvector, eigenvalue 0).
```
So the Feshbach-reduced **coarse operator is gapless, exactly, on the free massless field**. The k=0 soft mode lives in `V` and survives the fiber elimination **identically**: it has *zero* overlap with the fiber (`Wᵀ e₀=0`), so no self-energy — of **either sign** — can move it off 0.

**Exact numbers** (1D `L=12`, `b=3` ⇒ `N=12`, `n=4`, `f=8`; single process, dense, `numpy`):

| quantity | value | meaning |
|---|---|---|
| `inf spec M` (full) | `6.2e-17` | `=0` (truth), the k=0 mode |
| `inf spec A` (coarse block) | `1.2e-17` | `=0` already on `A` |
| `inf spec (A − Σ(0))` (coarse op) | `9.4e-18` | `=0` — Feshbach coarse operator exactly gapless |
| `‖e₀ − P_V e₀‖` | `1.9e-16` | `=0` — the soft mode is **in V** |
| `‖Wᵀ e₀‖ = ‖P_F M e₀‖` | `0.0` | soft mode **exactly decoupled** from fiber |
| `‖Σ(z) e₀‖` at `z=0` | `3.9e-17` | self-energy annihilates the soft mode |
| `min eig Σ(0)` | `5.8e-17` | `Σ ⪰ 0` (PSD) — sign confirmed |
| Schur residual `max|logdet(M−z) − logdet(D−z) − logdet S(z)|` | `8.9e-16` | Feshbach is **exact** |
| first nonzero `M` | `0.2679 = 4 sin²(π/12)` | matches dispersion |
| fiber gap `γ_fib` (§5) | `1.0` | `>0` — UV sector gapped even at `m=0` |

2D check (`4×4`, block `2×2` ⇒ `N=16`, `n=4`): `inf spec M = −5e-16 (=0)`, `inf spec(A−Σ(0)) = −6.7e-16 (=0)`, `‖P_F M e₀‖ = 0.0`, `γ_fib = 2.0`. Same verdict in `d=2`.

> **Decisive check = PASS on Reading G.** inf spec(coarse) is `0` to machine precision (≤ `1e-16`), **not** small-positive. The massless soft mode survives the Feshbach map with zero self-energy.

---

## §3. FREE MASSIVE (`m>0`) — coarse spectrum `= m²` exactly, no shift

**[真定理.]** Now `M e₀ = m² e₀`. Since `e₀∈V`, `Wᵀ e₀ = P_F M e₀ = m² P_F e₀ = 0` (because `P_F e₀=0`), and `A e₀ = P_V M e₀ = m² e₀`. Hence `Σ(z) e₀ = 0` for all `z`, and

```
S(m²) e₀ = (A − m² − Σ(m²)) e₀ = (m² − m² − 0) e₀ = 0 .
```
So `z = m²` is an **exact** coarse eigenvalue, and `inf spec` of the coarse problem `= m²`. The Feshbach map reproduces the **physical mass with no shift**.

**Exact numbers** (1D `L=12,b=3`, `m²=0.5`): `inf spec M = 0.4999999999999996`; `inf spec A = 0.4999999999999998`; `min eig S(z=m²) = −0.0` (i.e. `0`: `z=m²` is a coarse eigenvalue); `‖Wᵀ e₀‖ = 4.8e-17`; `γ_fib = 1.5 = m² + γ_fib(0)`; Schur residual `6.2e-15`. 2D (`4×4, 2×2`, `m²=0.5`): `inf spec M = 0.5`, `γ_fib = 2.5`. **[真定理: massive coarse spectrum = `m²`, Feshbach accurate, no shift.]**

---

## §4. Self-energy sign check — `Σ⪰0`, the `−Σ` *lowers*, it does not gap; soft-mode capacity `= 0`

**[真定理.]** For `z < γ_fib` (in particular `z=0` and `z=m²<γ_fib`), `D − z ⪰ D − 0 ≻ 0` (since `D ⪰ γ_fib·1` when `m=0`, §5), so `(D−z)⁻¹ ⪰ 0`, hence `Σ(z) = W(D−z)⁻¹Wᵀ ⪰ 0` (PSD). The coarse operator carries `A − z − Σ(z)`: the **minus** sign in front of `Σ` makes it **negative-semidefinite** in the coarse operator, i.e. it **lowers** the coarse energies. It **cannot** create a gap.

Measured (massless 1D `L=12,b=3`), coarse-block eigenvalues before/after the self-energy:

```
A                eigenvalues : {  0.000000 , 0.666667 , 0.666667 , 1.333333 }
(A − Σ(0))       eigenvalues : {  0.000000 , 0.315789 , 0.315789 , 1.090909 }
per-mode change  (A − (A−Σ)) : {  0.000000 , 0.350877 , 0.350877 , 0.242424 }   ≥ 0
Σ(0)             eigenvalues : {  0.000000 , 0.242424 , 0.350877 , 0.350877 }   ⪰ 0
```
Two exact facts, both load-bearing:

1. **`−Σ` lowers every non-soft coarse mode** (Cauchy interlacing / Courant–Fischer min-max: `λ_k(A−Σ) ≤ λ_k(A)` since `Σ⪰0`; here `{2/3,2/3,4/3} → {0.316,0.316,1.09}`). It pushes the coarse spectrum *toward* zero — the **opposite** of creating a gap. **[真定理]**
2. **On the soft mode `Σ e₀ = 0` exactly** (the `0` entries above), so the soft eigenvalue is *unmoved* at `0` (massless) / `m²` (massive). The **relative spectral capacity of the soft mode for the fiber is exactly 0** — the physical soft mode is orthogonal to the fiber coupling, so no fiber elimination (no `z`, no sign) can move it. **[真定理 — the optional "relative capacity for one Gaussian block": capacity(soft) `= 0`; capacity(non-soft) `= 0.35/0.67 ≈ 0.53`, `0.24/1.33 ≈ 0.18`.]**

Consequence: a naive sign error (`+Σ` instead of `−Σ`) would *inflate* the non-soft coarse spectrum and **manufacture a spurious gap**; the correct sign removes that failure mode. And even the wrong sign cannot move the soft mode (capacity 0), so the massless PASS is sign-robust.

---

## §5. Fiber gap at `m=0`, and why the `"52"` inequality is free-field-blind

**[真定理 — the fiber is UV and gapped for free.]** The fiber `F = ker C` consists of within-block zero-sum (`k≠0`-within-block) modes. Even at `m=0`, `D = P_F(−Δ)P_F` is **positive definite**: `γ_fib := inf spec D > 0`. Exact values (1D `L=12`, `m=0`):

| `b` | `γ_fib` | sharp Poincaré `1/γ_fib` | local single-block spurious `2−2cos(π/b)` (§6) |
|---|---|---|---|
| 2 | 2.000000 | 0.500 | 2.000 |
| 3 | 1.000000 | 1.000 | 1.000 |
| 4 | 0.840533 | 1.190 | 0.586 |

Fiber spectrum for `b=3`: `{1.0, 1.78, 1.78, 3.0, 3.0, 3.55, 3.55, 3.67}` — all `>0`. The exact Poincaré inequality on the fiber is `‖η‖² ≤ (1/γ_fib)‖dη‖²` with the **sharp** constant `1/γ_fib = O(1)` (`0.5–1.2` for `b=2,3,4`).

**The acid-test kill [no-go for the fiber-only route].** `γ_fib > 0` **already at `m=0`** means the free *massless* field satisfies `‖η‖² ≤ 52‖dη‖²` identically (and even the sharp `1/γ_fib` version). So:

- the `"52"` inequality, and *any* bound resting on the fiber Hessian gap, **controls the UV/fiber sector — which is already gapped on the free massless field** — and is **silent on the soft sector `V` where the mass lives**. By the charter's rule it is **free-field-blind: vacuous as a gap tool**.
- the magnitude `52` is not reproduced for the free scalar (sharp constant is `O(1)`); per `research/physics_handoff.md` §风险 4 the source `2610.00341` is unreviewed and must not be relied on. Whether `52` is a loose 4D-block instance or a different normalization is **not** the point: a fiber-control inequality **satisfied at `m=0`** cannot pass a discriminating acid test.

**[真定理.** Fiber control is necessary hygiene (it keeps `η` bounded so `Σ(z)` exists) but contributes **nothing** to the gap: the gap, if any, is entirely in `V` (§2), and on the free field `V` is exactly gapless.**]**

---

## §6. The decisive fork: local single-block SRF-CUBE spuriously gaps; no L-uniform free gap exists

### §6a. Reading L (local single-block) — **FAIL / vacuous**

**[no-go.]** If SRF-CUBE reads the coarse spectral floor off a **single block** `b` sites wide, the block's own constant is a zero mode — but that local constant is **not** the global soft mode (which is the *coherent* pattern `V_I = const` across **all** blocks; a single block cannot host it). Whatever boundary the block is given, its first **relative/gradient** mode is strictly positive at `m=0`:

```
isolated free block  (open path Laplacian):  first relative mode = 2 − 2cos(π/b)
      b=2: 2.000 ,  b=3: 1.000 ,  b=4: 0.586                 (all > 0)
block with neighbours pinned to its mean (Dirichlet):
      lowest mode = 2 − 2cos(π/(b+1))  ;  b=2: 1.000 ,  b=3: 0.586 ,  b=4: 0.382  (all > 0)
```
Both modelings are **strictly positive on the free massless field** — a **spurious gap**. It arises precisely because Reading L **deletes the inter-block k=0 mode**: the global constant is the *coherent* pattern across blocks, which no per-block criterion can see. By the decisive check this incarnation "has secretly deleted the soft mode → the whole framework is VACUOUS" (`BRIEF-FOR-MODEL.md` §2(4)). The verdict is boundary-model-independent: **Reading L fails the acid test** whichever way the single block is closed. (Verified numerically: open `→ {2.000, 1.000, 0.586}`, Dirichlet `→ {1.000, 0.586, 0.382}` for `b=2,3,4`.)

### §6b. No L-uniform free gap exists — a hard ceiling on what any correct SRF-CUBE may claim

**[真定理.]** The first nonzero eigenvalue of the free massless lattice Laplacian is `4 sin²(π/L)` (one direction excited; same in every `d`), which **→ 0 as `L→∞`**:

| `L` (1D) | first nonzero `= 4 sin²(π/L)` | | 2D `L×L` | first nonzero |
|---|---|---|---|---|
| 12 | 0.267949 | | 4 | 2.000 |
| 24 | 0.068148 | | 6 | 1.000 |
| 48 | 0.017110 | | 8 | 0.586 |
| 96 | 0.004282 | | →∞ | → 0 |

So the free massless field has **no L-uniform spectral gap** — inf spec `= 0` exactly, first excited `→0`. Hence **any** framework that claims an **L-uniform positive** soft floor `m_eff²(k→0) ≥ λ > 0` is **already false on the free field** — a non-theorem. A correct SRF-CUBE must therefore evaluate to **exactly `0`** (or `→0`) on free. Reading G does (§2). This is the acid test's second set of teeth: it forbids the tempting "local positivity ⇒ uniform gap" upgrade.

---

## §7. Scaling to 4D — every decisive statement is dimension-independent

The free-field verdict does not change from `d=1` to `d=4`, because each decisive fact is a *structural* property of the constant mode vs. the block map, or a per-direction dispersion:

- **PASS (Reading G), `d=4` [真定理].** `e₀ = const` is block-constant (∈V) in **any** `d`; `M e₀ = m²e₀`; `Wᵀ e₀ = m²P_F e₀ = 0`; `Σ(z)e₀ = 0`; coarse inf spec `= 0` (massless) / `m²` (massive), block- and dimension-independent. The 4D k=0 soft mode survives identically.
- **Fiber gap `d=4` [真定理].** The 4D fiber (blockwise zero-sum) is UV; `γ_fib = γ_fib(0)+m² > 0` at `m=0` (verified `d=1,2`; structural in `d=4`). Fiber control stays free-field-blind.
- **Local spurious gap `d=4` [no-go].** The lowest mode of the open `dD` block Laplacian excites one direction: `2−2cos(π/b) > 0`, independent of `d` (others flat). Reading L is vacuous in 4D too.
- **No L-uniform free gap `d=4` [真定理].** first nonzero `= 2(1−cos 2π/L) = 4 sin²(π/L) → 0`; same ceiling.

So the entire acid-test verdict, and both trap-kills, transfer verbatim to the physical `d=4`, `SU(N)` setting insofar as the **free Gaussian** is concerned. (What does **not** transfer is the nonabelian interaction — that is §8.)

---

## §8. Scorecard (charter 5-point + acid-test discipline) and verdict

| # | Constraint | Verdict | Exact fact |
|---|---|---|---|
| 1 | Nonperturbative | **PASS (on free)** | Feshbach/Schur exact, no loop; but the test only exercises the Gaussian. |
| 2 | Sign-controlling lower bound `m_eff²≥λ` | **VACUOUS-PASS (free)** | On free it returns `0` (massless)/`m²` (massive) — correct, *not* a positive `λ`; the SU(N) `λ>0` is not produced. |
| 3 | Non-circular | **UNTOUCHED** | Free test says nothing about the SU(N) circularity (cf. `PROOF-COARSE-POINCARE`, H1). |
| 4 | **Fails free field** | **PASS (Reading G) / FAIL (Reading L)** | G: correctly gapless, inf spec `=0`. L: spurious gap `2−2cos(π/b)>0` → vacuous. |
| 5 | L-uniform `≥λ`, reach β≈2.5 | **FAIL (free ceiling)** | first nonzero `4 sin²(π/L)→0`; no L-uniform free gap. |
| — | Self-energy sign `Σ⪰0`, `−Σ` lowers | **PASS [真定理]** | `Σ⪰0`; per-mode `−Σ` *lowers* `{2/3,2/3,4/3}→{0.316,0.316,1.09}`; soft-mode capacity `=0`. |

**VERDICT — PASS, ADMISSIBLE, NECESSARY-NOT-SUFFICIENT.**
- **Free-massless: PASS in the faithful (Reading G) reduction.** inf spec(coarse) `= 0` **exactly** (machine zero `≤1e-16`), **not** a spurious gap. The k=0 soft mode is retained in `V` and survives the Feshbach map with zero self-energy. Free massive: inf spec `= m²` exactly, no shift.
- **The pass is forced, not earned.** Exactness of the Schur map plus the fact that the block-average `V` contains `k=0` make *every* faithful exact reduction pass. So the free test yields **no positive evidence** toward the gap; it only certifies the soft mode is **not deleted** and rules out sign errors.
- **The test kills the two live traps.** (i) Reading L (local single-block) **FAILS**: spurious gap `2−2cos(π/b)>0` at `m=0` → vacuous. (ii) The `"52"` fiber inequality (hexagonmath `2610.00341`, unreviewed) is **satisfied at `m=0`** → free-field-blind, cannot be the gap mechanism.
- **Admissibility condition.** SRF 2.0 may proceed to `SU(2)`/`SU(N)` **only** with SRF-CUBE defined as the **global-V** Feshbach criterion (Reading G), **and** with the explicit acknowledgement that on free it returns `0` — its SU(N) content must supply a genuine positive `λ` that the free field lacks.

---

## §9. Sharpest line, exact stuck step, refs

**Sharpest line.** On the free field the massless soft mode `e₀` sits in the **kernel of the fiber coupling itself** (`Wᵀe₀ = P_F M e₀ = 0`), so the Feshbach self-energy cannot touch it — for any `z`, any block size, any sign convention: with the correct sign `−Σ` merely **lowers** the non-soft coarse energies (`{2/3,2/3,4/3}→{0.316,0.316,1.09}`) and the coarse operator is gapless **exactly** (`inf spec = 0`, machine-zero). SRF 2.0 therefore **passes the free-massless acid test** — but the pass is a consequence of *exactness*, not of any mechanism, so it buys nothing: it is the **local single-block** criterion (which freezes the inter-block `k=0` mode and spuriously gaps by `2−2cos(π/b)>0`) that fails, and it is the **`"52"` fiber inequality** that is *satisfied at `m=0`* and thus blind to the gap. SRF is **admissible iff it is the global-V criterion**; its whole scientific burden is then exactly the piece the free field cannot touch.

**Exact stuck step [缺失步].** The free acid test reduces to two exact trivialities — "the Schur map of a Gaussian is exact" and "the block average keeps `k=0`" — so it certifies nothing about `SU(N)`. The load-bearing step is **unchanged**: an **L-uniform**, **non-loop** (one-loop is Savvidy-indefinite, `PROOF-COARSE-POINCARE` (a)(iii)), **non-circular** lower bound on the `SU(N)` coarse-mode spectrum that reaches `β≈2.5` — i.e. **H1** (`G″(k→0)=1/χ`, `Var(a_{k→0})≤1/λ`). On the free field that bound is `0`; SRF-CUBE must *supply* the positivity the free field lacks, and that supply is H1, not a consequence of the Feshbach/block structure.

---

*Refs.* Internal (ground, not invented): `docs/BRIEF-FOR-MODEL.md` §2(4),§6 (free-field acid test; a tool must fail on free massless) and §7 (the wall); `CLAUDE.md` §1 (acid-test iron rule), §2 (labeling); `research/physics_handoff.md` §风险 4 (**`2610.00341` AI-heavy/suspect, "only a lead, not a basis"; `2602.0054` retracted**) — source of the `"52"` inequality, *not relied on*; `k0-CONSTANT-MODE-TRANSFER` (constant k=0 mode = soft mode; free theory decouples it as the exact zero mode; any "constant-mode" reduction that keeps it is exactly gapless); `D1-BLOCK-MARTINGALE` (block decomposition; `Var(a)=χ`, `m_eff²(k→0)≥λ` is H1); `PROOF-COARSE-POINCARE` (H1; one-loop `H_ff` Savvidy-indefinite (a)(iii)); `FORGE-FESHBACH-ANGLE` (Feshbach/Schur sign, Haynsworth inertia; the cross-block sign on the Savvidy complement); `docs/yang-mills-prior-work.md` §1.4 (clustering ⇏ gap; hexagonmath counterexample). Code: `research/unconditional/srf_free_acid.py` (exact dense, single process, `N≤256`). External: Schur complement / Feshbach projection (Feshbach, *Ann. Phys.* **5** (1958) 357; *ibid.* **19** (1962) 287); Haynsworth, *Lin. Algebra Appl.* **1** (1968) 73 (inertia additivity); Cauchy interlacing / Courant–Fischer min-max, from which `λ_k(A−Σ) ≤ λ_k(A)` for `Σ⪰0` (eigenvalue monotonicity under a PSD subtraction); standard discrete-Laplacian spectrum `2Σ_μ(1−cos 2πk_μ/L)`.
