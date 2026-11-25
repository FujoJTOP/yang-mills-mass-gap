# FOSP-A — the representation-fusion transfer kernel `K_{αβ}=⟨α|T|β⟩`: TP₂ holds for SU(2) *only because* its fusion graph is a 1-D path (Jacobi/Mehler = T3 in disguise), and **fails robustly for SU(3)**; the confining Casimir's positive cross term makes the SU(3) kernel *submodular*, so no representation order (1-D or 2-D/MTP₂) supplies the missing second cone — **kill**

> **2026-11-24 · FOSP-A — the first elimination test of the "Fusion-Ordered Spin-network Positivity"
> species.** Question: on gauge-invariant spin-network/fusion states, ordered by highest-weight
> dominance, is the transfer kernel `K_{αβ}=⟨α|T|β⟩` TP₂ (all `2×2` minors `≥0`) with a *stable*
> positive structure? Reads `T3-TP2-ISOLATION.md`, `R2-ISOLATION.md`, `X2-FACTORIZATION.md`,
> `T4-CARRIER.md`, `docs/paper-v2/ch04-gap.tex`. Witness `research/v3/w7/fosp_a.py`
> (`py -3.12 -I`; numpy+scipy+mpmath; **35 checks, 0 FAIL**; output `w7/fosp_a.out`).
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

FOSP proposes to beat T3's "no order on `𝒜/𝒢`" no-go by using a **gauge-invariant** order — the
representation labels with their highest-weight (dominance) order — instead of the orderless
position basis. This is a genuine new attempt; the answer is **negative**, and sharply so:

1. **The structural criterion is Gantmacher–Krein, and it forces 1-D.** `K=e^{-H}`, where
   `H=(g²/2)C₂(R) − (1/g²)M` (electric Casimir + magnetic fusion adjacency, the gauge-invariant
   single-plaquette Kogut–Susskind truncation). **If the generator `−H` is a Jacobi (tridiagonal)
   matrix with nonnegative off-diagonals, then `K` is totally nonnegative and strictly TP₂**
   (Gantmacher–Krein; discrete Sturm–Liouville). Tridiagonality `⇔` the fusion graph is a **path**
   `⇔` the order is **1-D**. **[T-standard]**

2. **SU(2) passes — but it *is* T3, not something new.** The SU(2) fusion adjacency `j→j±½` is
   tridiagonal, so `K` is TP₂ (verified to **40 digits**: `0` negative minors of `1296` at
   `g²=0.5,1.0,2.4`; `0/6084` at `j_{max}=6`). This is exactly the 1-D Mehler/Sturm–Liouville total
   positivity T3 already banked — **the transfer matrix in disguise**. **[T]+[V]**

3. **SU(3) fails — robustly, at every coupling, truncation, and order.** The SU(3) fusion graph is
   the **2-D weight lattice** (`bandwidth(M)=4–5`, not tridiagonal), so Gantmacher–Krein gives
   nothing. Direct computation: **≈26–29 % of all `2×2` minors are negative**, to 40 digits, in
   *every* linear extension of the dominance order (`−0.19` at `g²=1`; `−0.019` at physical
   `g²=2.4`). The **2-D rescue** — MTP₂ (log-supermodularity on the dominance poset) — **also
   fails** (`−0.044` at `g²=1`). **Analytic cause:** the SU(3) Casimir has a **positive cross term**
   `∂_p∂_q C₂ = 1/3 > 0` in Dynkin coordinates, making `log K` **submodular** — the exact 2-D
   analogue of T3's anti-correlated Gaussian (§2.3). **[V]+[NG]**

4. **The kernel is order-1 positive but order-2 blind.** `K` is entrywise `>0` and positive-definite
   (Perron positivity holds), yet fails TP₂. FOSP reproduces, *in representation labels*, the
   R2/T3 **tier wall**: order-1 (Perron) is free, order-2 (TP₂) is where the physical group dies.
   This is a genuine refinement of the T3 no-go — the obstruction is not "no order" but "no **1-D**
   order" — but it is still a no-go. **[NG]**

**Kill verdict.** Two of the three kill clauses fire: **TP₂ fails** (SU(3), the physical group) and
the object **is the transfer matrix in disguise** (`K=e^{-H}`, SU(2) TP₂ = Jacobi/TN). The
"sign-indefinite" clause does **not** fire (`K` is order-1 positive); the failure is strictly at the
**second order** — which is precisely the tier the campaign cannot cross. **No QED. Clay 0%.**

---

## 1 · The FOSP object, and why it is a fair test

T3's no-go (`No-go T3-C`) says the coarse transfer kernel on the orderless quotient `𝒜/𝒢` carries
**no global order**, so the second invariant cone (TP₂) cannot even be *stated* there. FOSP's
proposal: the **representation labels are gauge-invariant**, and they carry a natural **dominance
order**; use that as the order and test TP₂ on `K_{αβ}=⟨α|T|β⟩` in the fusion basis. This is a
legitimate, non-strawman attempt — it looks for the missing cone in *group-theoretic* data rather
than geometric position data. We test it exactly.

**The gauge-invariant truncation.** A single spatial plaquette with dynamical holonomy `U_p`; the
gauge-invariant Hilbert space is spanned by the character states `|R⟩` (plaquette in irrep `R`) — a
spin-network loop. The Kogut–Susskind generator in the `|R⟩` basis is
```
     H(g) = (g²/2)·C₂(R)  −  (1/g²)·M ,     K = e^{−H} ,
     M_{RR'} = ⟨R| (χ_fund + χ_fundbar)/2 |R'⟩  = ½ · (fund × antifund fusion adjacency) .
```
`M` comes from the fusion rules
```
   SU(2):  ½ ⊗ j        = (j−½) ⊕ (j+½)
   SU(3):  [1,0]⊗[p,q] = [p+1,q] ⊕ [p,q−1] ⊕ [p−1,q+1]      (drop negative)
           [0,1]⊗[p,q] = [p,q+1] ⊕ [p−1,q] ⊕ [p+1,q−1]      (drop negative)
```
(dimension-counted in the witness). `H` has the two structural ingredients of every lattice gauge
transfer matrix — an **electric** diagonal Casimir term and a **magnetic** off-diagonal fusion term —
so it is the minimal honest testbed. **Caveat, stated plainly:** this is a `0+1`-D truncation (the
single-plaquette QM), *not* the 4-D YM transfer matrix; it tests the **structural** question FOSP
poses (does the representation order make the kernel TP₂-stable?), which is exactly FOSP-A's
remit. The full 4-D generator is *further* from tridiagonal (§5), so the kill extends.

**Ordering.** "Highest-weight dominance" is a **partial** order. `C₂` is strictly monotone under it
(`λ ⪰ μ ⇒ C₂(λ) ≥ C₂(μ)`), so we use Casimir-ascending as the primary linear extension and test
alternatives (level `p+2q`, `p`, `q`, conjugate flips). This matters: **TP₂ is only defined for a
total order**, so its very statement presumes a linear extension — an immediate foundational
objection to FOSP on a partial order.

---

## 2 · The structural criterion: TP₂ `⇔` a tridiagonal Jacobi generator (so `⇔` 1-D)

The decisive mechanism is classical and sharp.

> **Theorem FOSP-A-1 (Gantmacher–Krein; `[T-standard]`).** Let `A` be a symmetric **Jacobi**
> (tridiagonal) matrix with **nonnegative** off-diagonal entries. Then `e^{tA}` is **totally
> nonnegative** for all `t ≥ 0`, and **strictly totally positive** if the off-diagonals are `>0`.
> Equivalently (Karlin, via T3 Theorem T3-A): a strictly TP₂ kernel is the heat kernel of a regular
> **Sturm–Liouville** operator — a **1-D** object. TP₂ of `e^{−H}` is therefore *equivalent to the
> generator `−H` being a 1-D Jacobi operator*.

Here `−H = −(g²/2)·diag(C₂) + (1/g²)·M`. Its off-diagonal is `+(1/g²)M ≥ 0`, its diagonal is
`−(g²/2)C₂ < 0`. So:
- **`M` tridiagonal ⇒ `−H` Jacobi ⇒ `K` TP₂** (the theorem applies).
- **`M` not tridiagonal ⇒ no TP₂ guarantee** (the theorem is silent; a specific kernel may still or
  may not be TP₂ — the actual object decides).

The witness confirms both the positive direction (a tridiagonal generator gives a TP₂ kernel; a
banded generator can break it) **and** the concrete conclusion for the real kernels (SU(2) passes,
SU(3) fails). The key identity: **TP₂ `⇔` 1-D order.** That is the whole story.

---

## 3 · SU(2): TP₂ holds — and it is T3's 1-D Mehler kernel, not a new cone `[T]+[V]`

For SU(2) the fusion adjacency is `j ↔ j±½`, i.e. the **path** (bandwidth 1). Hence `−H` is a
tridiagonal Jacobi matrix, and `K=e^{−H}` is TP₂ by Theorem FOSP-A-1. Because `C₂(j)=j(j+1)` is
**quadratic** in `2j`, `H` is a **discrete harmonic oscillator**; `K` is a **Mehler-type** heat
kernel — the exact operator T3 §2.3 already identified as the canonical TP₂ kernel. Witness
(40-digit):

| `g²` | min `2×2` minor | #negative / 1296 |
|---|---|---|
| 0.5 | `+4.08×10⁻⁶` | 0 |
| 1.0 | `+6.47×10⁻¹²` | 0 |
| 2.4 | `+1.97×10⁻²⁰` | 0 |

and at `j_max=6` (13 states), `0/6084` negative at every `g²`. **SU(2) is TP₂** — but the minors
become *vanishingly* positive in the UV corner (the kernel is near low-rank there, `K_{jj'} ~
e^{−c(j+j')}`), i.e. TP₂ is **marginal, not uniformly strict** — exactly T3's caveat that "TP₂
gives strictness, not the constant `c`." **This is not a new second cone; it is the known 1-D
Jacobi/Sturm–Liouville positivity of the transfer matrix.** The transfer matrix, in disguise.

---

## 4 · SU(3): TP₂ fails, and so does the 2-D rescue `[V]+[NG]`

For SU(3) the fusion adjacency is the **2-D weight lattice** in the dominant Weyl chamber
(bandwidth 4–5, not tridiagonal), so Theorem FOSP-A-1 does **not** apply. The actual kernel decides,
and it decides against TP₂.

**1-D test (linear extensions of dominance).** Direct, 40-digit computation of all `2×2` minors:

| order | `g²=1`: min minor | #negative / total |
|---|---|---|
| Casimir ↑ | `−1.913×10⁻¹` | 38224 / 142884 |
| Casimir ↓ | `−1.913×10⁻¹` | 38224 / 142884 |
| level `p+2q` | `−1.913×10⁻¹` | 34958 / 142884 |
| `p` | `−1.913×10⁻¹` | 36246 / 142884 |
| `q` | `−1.913×10⁻¹` | 36246 / 142884 |

At small truncation and 40-digit precision: `p+q≤3` → `536/2025` negative (`26 %`) at `g²=1`,
`580/2025` (`29 %`) at `g²=2.4`; `p+q≤4` → `2950/11025` (`27 %`), `3184/11025` (`29 %`). The
failure is **not** a numerical artifact (40-digit) and **not** an ordering accident (all five
extensions fail).

**2-D rescue (MTP₂ / log-supermodularity on the dominance poset).** The natural order-invariant
notion for a *partial* order is multivariate TP₂ / log-supermodularity: for comparable `x≼x'`,
`y≼y'`, require `K(x,y)K(x',y') ≥ K(x,y')K(x',y)`. It **also fails**:

| `g²` | MTP₂ min | #negative / 3025 |
|---|---|---|
| 1.0 | `−4.38×10⁻²` | 534 |
| 2.4 | `−2.89×10⁻⁴` | 630 |

So FOSP fails under *both* the 1-D and the 2-D readings. There is no representation order that
rescues SU(3).

**The analytic cause — the confining Casimir is submodular.** In Dynkin coordinates the quadratic
Casimir is `C₂(p,q) = (p²+pq+q²)/3 + (p+q)`, so its mixed second difference is
```
     ∂_p ∂_q C₂ = 1/3 > 0        (discrete: C₂(2,2)−C₂(2,1)−C₂(1,2)+C₂(1,1) = 1/3) .
```
Because `log K ⊃ −(g²/2)C₂`, the cross-derivative is `∂_p∂_q log K = −(g²/2)(1/3) = −g²/6 < 0`:
`log K` is **submodular**, so MTP₂/TP₂ fails. **This is the exact 2-D analogue of T3's
anti-correlated Gaussian** (`∂_x∂_y log K = −b < 0`, a `2×2` minor `−0.54`, §2.3). The confining
potential of SU(3) — the positive `pq` cross term — is what breaks the order-2 positivity. A
second, structural obstruction compounds it: **the dominance order only relates same-triality
(same centre/N-ality) irreps**, while the magnetic term connects different-triality sectors, so `K`
necessarily mixes sectors that the order cannot order.

---

## 5 · The three elimination tests

**(1) Free massless must not give a positive gap — PASS.** In the KS generator the only parameter
that can produce a gap is the electric Casimir spacing `(g²/2)·ΔC₂`; the free-massless limit is
`g²→0`, and the gap `→0` (`(g²/2)(¾) = 0.375 → 0.00375` as `g²: 1 → 0.01`). In the
magnetic-dominated regime the `−M` **band edge** is approached continuously, so the bulk gap `→0`
with truncation (`2.84×10⁻¹ → 6.79×10⁻³` as `j_max: 4 → 32`) — a continuum, not a gap. FOSP
certifies **no** positive gap on the free/massless limit. **[T]+[V]**

**(2) 2D YM `ρ²=8/3` must not be misread as the `trF²` atom (`w₁=0`) — PASS.** The 2D YM transfer
matrix in the character basis is `e^{−t C₂(R)}` (diagonal); the holonomy ratio is
`ρ² = C₂(adj)/C₂(fund) = 8/3` exactly. But `trF²` in 2D **equals the Hamiltonian** and is conserved,
hence diagonal in the *same* basis, with eigenvalue `C₂(R)`; so the atom
`w₁ = |⟨0|trF²|fund⟩|² = 0` exactly. The `8/3` is the **holonomy tower**, not the `trF²` isolation;
FOSP must not advertise it as 4-D `trF²` evidence. It does not. **[T]+[V]** (matches T3 §2.1 / R5).

**(3) SU(3) minors must be stable at physical coupling, not only strong — FAILS (this is the kill).**
At `β=2.5` (`g²=6/β≈2.4` for SU(3)), the minors are **stably negative**, not stably positive:
- **Truncation:** min minor `−1.910×10⁻²` for *every* cutoff `p+q≤3,4,5,6,7,8` (10→45 states) —
  identical, because the minimal minor is attained among the lowest irreps present from the start.
- **Coupling:** min minor negative across the whole physical window — `−6.61` (`g²=0.5`), `−1.91×10⁻¹`
  (`1.0`), `−1.91×10⁻²` (`2.4`), `−3.53×10⁻⁴` (`6.0`) — no fine-tuned sign change, robustly negative.

So the SU(3) `2×2` minima are *stable* — but stably **negative**. The criterion demanded by FOSP-A
(stable *positive* structure at physical coupling) is **not met**. **[V]+[NG]**

---

## 6 · "Is it the transfer matrix in disguise?" — yes `[NG]`

By construction `K = e^{−H}` is the transfer matrix of the single-plaquette QM. The witness confirms
`eig(K) = exp(−eig(H))` to `≤2×10⁻¹⁶`. Its positivity is therefore a statement about the **transfer
generator**, not a new object:

- **SU(2):** `TP₂ ⇔ −H` is Jacobi `⇔` the SU(2) fusion graph is a path `⇔` **1-D**. This is T3's
  Sturm–Liouville/Mehler theorem, re-derived. **Renamed, not new.**
- **SU(3):** the generator is not Jacobi, and `K` fails TP₂ — the transfer matrix is *not* TP₂.

And the **tier point** crystallises the whole picture: SU(3) `K` is **order-1 positive** (entrywise
`>0`, positive-definite — Perron holds) but **order-2 blind** (TP₂ fails). FOSP thus reproduces, in
representation labels, exactly the R2/T3 tier wall: **order-1 (Perron) is free; order-2 (TP₂) is
where the physical group fails.** The "second cone" FOSP hoped the representation order would supply
is *the same* missing cone, now shown to be absent for SU(3) — the orderless-quotient obstruction
T3 identified is in fact a **no-1-D-order** obstruction, which representation labels do not remove.

---

## 7 · Exact missing step `[O]`, and what would flip the verdict

The obstruction is now exact and two-headed, both **analytic**:

1. **Structural (`[T]`):** TP₂ of `e^{−H}` requires a **tridiagonal** (1-D Jacobi) generator
   (Gantmacher–Krein). SU(3)'s fusion graph is 2-D — no 1-D order exists. `[NG]` for the class
   "TP₂ via a linear order on SU(3) irreps."
2. **Dynamical (`[T]`):** the confining SU(3) Casimir has `∂_p∂_q C₂ = 1/3 > 0`, forcing `log K`
   **submodular** — no MTP₂ either. `[NG]` for the class "MTP₂ on the dominance poset."

**Missing step `[O]`.** A **genuinely 2-D order-2 positivity** that is (i) *non-submodular* in the
confining direction and (ii) coarse-graining-stable — i.e. an order-2 cone on a 2-D weight lattice
that survives the positive Casimir cross term. No such structure is known; TP₂/MTP₂ are blocked
analytically, and the only surviving 1-D case (SU(2)) is the transfer matrix. **Flip condition:**
exhibit a 2-D, non-submodular, coarse-graining-stable order-2 kernel on SU(3) fusion states whose
minors are stably positive at `β=2.5` and that reduces to the YM transfer matrix — equivalently
neutralise the `+1/3` cross term without changing the physics. That is a new species of positivity,
not currently in hand.

**Reach at `β=2.5`:** zero. The representation order — the *one* gauge-invariant order FOSP could
appeal to — is 2-D for SU(3) and the minors are stably negative. The toy is `0+1`-D, but the full
4-D generator's magnetic part is *higher*-dimensional still (the spin-network configuration graph),
so the Jacobi route is even more foreclosed, and the positive-Casimir-cross-term obstruction is
coupling- and dimension-independent.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| FOSP object | `K=e^{−H}`, `H=(g²/2)C₂−(1/g²)M` (electric Casimir + magnetic fusion) on single-plaquette spin-network states | **[T]** (definition) |
| structural criterion | TP₂ of `e^{−H}` `⇔` `−H` a tridiagonal Jacobi matrix (Gantmacher–Krein) `⇔` 1-D fusion graph | **[T-standard]** |
| SU(2) adjacency | path (`j→j±½`), tridiagonal ⇒ `K` TP₂ (Mehler/discrete-Sturm–Liouville) | **[T]+[V]** |
| SU(2) verdict | TP₂ holds to 40 digits (`0/1296`, `0/6084` negative), but = T3's 1-D cone; **transfer matrix in disguise** | **[T]+[V]** |
| SU(3) adjacency | 2-D weight lattice, **not** tridiagonal (`bandwidth 4–5`) | **[T]+[V]** |
| SU(3) 1-D TP₂ | **fails** in every linear extension; ≈26–29 % of minors negative, 40-digit robust | **[V]+[NG]** |
| SU(3) 2-D MTP₂ | **fails** (`−0.044` at `g²=1`); `log K` **submodular** since `∂_p∂_q C₂=1/3>0` | **[T]+[V]+[NG]** |
| triality | dominance orders only same-triality irreps; magnetic term mixes sectors ⇒ harms TP₂ | **[T]** |
| elimination (1) | free massless: gap `→0` (electric `→0`, magnetic band gap `→0`); no positive gap | **[T]+[V]** (PASS) |
| elimination (2) | 2D YM `ρ²=8/3` (holonomy) coexists with `trF²` atom `w₁=0` | **[T]+[V]** (PASS) |
| elimination (3) | SU(3) minors **stably negative** at physical `g²≈2.4` (all cutoffs/couplings) | **[V]+[NG]** (FAIL) |
| tier | SU(3) `K` order-1 positive (entrywise, PD) but order-2 blind (TP₂ fails) | **[NG]** |
| kill | TP₂ fails (SU(3)) **and** transfer matrix in disguise; "sign-indefinite" clause does *not* fire | **[NG]** |
| missing step | a non-submodular, coarse-graining-stable 2-D order-2 cone on SU(3) fusion states at `β=2.5` | **[O]** |
| Clay mass gap | unbuilt | **[O] · 0 %** |

---

## 9 · Sharpest line

**FOSP's idea — order gauge-invariant spin-network states by highest-weight dominance and demand a
stable TP₂ kernel — is testable and fails, for a sharp reason.** TP₂ of `K=e^{−H}` is equivalent
(Gantmacher–Krein) to the generator being a **1-D Jacobi** matrix. SU(2)'s fusion graph is a path, so
its kernel is TP₂ — but that is precisely T3's 1-D Mehler/Sturm–Liouville total positivity, the
**transfer matrix in disguise**, not a new second cone. SU(3)'s fusion graph is the **2-D weight
lattice**, and its kernel fails TP₂ **and** MTP₂ robustly (≈26 % of `2×2` minors negative, 40-digit,
every order, every truncation, every coupling through `β=2.5`), with the analytic cause being the
confining Casimir's **positive cross term** `∂_p∂_q C₂ = 1/3 > 0`, which makes `log K` submodular —
the exact 2-D analogue of T3's anti-correlated Gaussian. The kernel is order-1 positive (Perron)
but **order-2 blind**: FOSP reproduces the R2/T3 tier wall in representation labels. The three
acid tests behave correctly (free massless gapless, 2D `ρ²=8/3` ≠ `trF²` atom, SU(3) stably
non-positive), so the test is not vacuous — it is a genuine **kill**: TP₂ fails for the physical
group, and for SU(2) it is the transfer matrix renamed. **No QED. Clay 0%.**

---

*Refs.* **Order-2 / TP₂:** `T3-TP2-ISOLATION.md` (TP₂ ⟺ Sturm–Liouville; Mehler; anti-correlated
Gaussian; window theorem); `R2-ISOLATION.md` (order-1 vs order-2 tier; `ρ²=8/3`; `I1↛I2`). **Carrier
honesty:** `T3` §2.1, `R5-EXACT-ANCHORS.md`, `T4-CARRIER.md` (2D `trF²=H` conserved ⇒ atom `0`).
**Factorization:** `X2-FACTORIZATION.md` (LC = transfer matrix in disguise — the same kill shape).
**Standard tools as such:** Gantmacher–Krein, *Oscillation Matrices and Kernels* (tridiagonal ⇒ TN
`e^{tA}`); Karlin, *Total Positivity* Vol. I (TP₂ ⟺ Sturm–Liouville heat kernel); MTP₂ / multivariate
log-supermodularity (Karlin–Rinott). **Gauge model:** Kogut–Susskind Hamiltonian; Migdal 2D YM
(character/heat-kernel transfer matrix); SU(3) quadratic Casimir `C₂(p,q)=(p²+pq+q²)/3+(p+q)`
(standard). **Witness:** `research/v3/w7/fosp_a.py` (numpy/scipy/mpmath, `py -3.12 -I`,
deterministic, 35 checks, 0 FAIL; output `w7/fosp_a.out`). **No third-party downloads; no `[C]`
upgraded to `[T]`; no QED; Clay 0%.**
