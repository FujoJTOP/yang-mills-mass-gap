# BASE-INFINITE-VOL — can the strong-coupling residue be carried toward the crossover by a β-uniform tool? **No: β\* is the KP radius, unmoved**

> **VERDICT (no QED).** The P3 residue floor `f₁²_L ≥ r₀(β)>0` (L-uniform, I1) is a **strong-coupling** theorem whose radius `β*` is set by the polymer/Kotecký–Preiss convergence criterion. We test the four β-uniform structural tools named in the directive — **reflection positivity across `L`, Gaussian domination (Fröhlich–Simon–Spencer), monotone coupling, correlation inequality** — for whether any extends the residue floor past `β*` toward the crossover `β≈2.5` *without crossing it*. **None does.** The residue is a *magnitude* (a KL atom weight); the four tools are *structural* (positivity, monotonicity, domination of expectations). The only object that carries the magnitude is the convergent polymer expansion, whose radius is
> `β* = 3.4×10⁻³` (crude KP, `a=1`), `≤ 8.9×10⁻³` (grading `a` optimised), benchmark `1/48 = 2.08×10⁻²` (Shen–Zhu–Zhu) — all **≥ 2 orders of magnitude below** `β≈2.5`. **The argument stops exactly at the KP boundary `λ_aΦ(a)=½`; the only exit is a rigorous 4D-gauge RG, which does not exist.** The residue is a base that cannot be grown toward the crossover by any β-uniform inequality.

---

## §0 · Labels, iron laws, acid test (carried from P3)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-here]** new this round.
**IL-i** non-positivity (fail on free massless: `f₁²=0`). **IL-ii** nonabelian. **IL-iii** no small parameter (relevant only once one tries `β≈2.5`).
**Acid** (`P3-DRESSED-MASS`): free field `a₀=0 ⇒ f₁²=0`; our `r₀(β)→0` as `u→0` reproduces it (witness §5).

---

## §1 · The base piece: tube–vacuum factorization of the residue `[V-here]`

The P3 residue is produced by the **four-wall tube** (`P3-FOURWALL-DERIVATION`: `|P|≥4t/a`, `m₁=4ln(1/u)/a`). Its contribution to `C_L(t)` factorizes **exactly** into a *tube weight* times a *vacuum reweighting*:

```
  C_L(t) = W_tube(β) · R_vac(β) · e^{−m₁ t} (1 + O(u)),        f₁²_L = lim_t e^{m₁ t} C_L(t) = W_tube·R_vac(1+O(u)),
  W_tube(β)  = a₀ u^{4/a},     a₀>0  explicit for ALL β,
  R_vac(β)   = Z_{tube}/Z_0 ∈ (0,1],   a positive ratio of partition functions.
```

**Why this is the right base object.** `W_tube` is *β-explicit and monotone* for all `β` — it is **not** the obstruction. **Every bit of β-dependence that decides whether the residue survives is inside the single scalar `R_vac`**, a positive partition-function ratio. So the question "extend the residue past `β*`" becomes, sharply:

> **Is there a lower bound `R_vac(β) ≥ c(β)>0` that is uniform in `L` and valid up to `β≈2.5`?**

The strong-coupling expansion answers it only for `β<β*`: `R_vac ≥ 1 − O(u)` (the tube perturbs the vacuum by `O(u)`), which is positive iff `λ_aΦ(a)<½`, i.e. `β<β*`. **Past `β*` this lower bound goes negative (vacuous), and no β-uniform tool replaces it.** That is the whole content of the wall, isolated to one positive scalar.

`[V-here]` is the *factorization*, not an evaluation of `R_vac` — it relocates the open step from "the residue" to "a vacuum partition-function ratio," which is the object the four tools are designed to bound. §2 shows they cannot.

---

## §2 · Tool-by-tool: none extends `β*`

**(A) Reflection positivity across `L` `[T]` — vacuous for the floor.**
RP (spatial reflection on `T³_L×ℝ`) is a *structural* property of the Wilson action, valid for **all `β≥0`**, independent of the coupling. It is exactly the imported line of `P3-CYLINDER-RESIDUE` (gives OS/KL, `dμ_O≥0`, self-adjoint `H≥0`). Applied to `R_vac` it yields only `0 ≤ R_vac ≤ 1` — the *trivial* interval, with **no `β`-dependent lower bound**. RP is β-uniform but *magnitude-blind*: it certifies the spectral representation, never the residue weight. **Does not move `β*`.**

**(B) Gaussian domination (Fröhlich–Simon–Spencer) `[T]` — upper bound, wrong direction.**
FSS Gaussian domination gives `⟨e^{f}⟩ ≤ exp(½⟨f,f⟩_{Gauss})` — a bound on the *expectation of an exponential*, i.e. an **upper** bound on `R_vac`, and it degenerates (the Gaussian covariance grows with the tube volume) as `β` grows. For a *lower* floor on `R_vac` it gives nothing. (Gaussian domination is what makes the strong-coupling *mass gap* rigorous — but *inside* the small-`β` regime, combined with the polymer expansion, not replacing it.) **Does not move `β*`.**

**(C) Monotone coupling `[NG]` — unproven and physically wrong for a KL atom.**
A clean extension would follow if `R_vac(β)` (hence `f₁²(β)`) were monotone in `β`: then `f₁²(β) ≥ f₁²(β*_0)>0` for all `β≥β*_0`. But (i) **no monotonicity of a KL atom weight in `β` is proved** for nonabelian lattice gauge; (ii) the mechanism fails: `f₁²` is a spectral *weight*, not an expectation, and the corpus's own measurement (`β=1.0` residue "dissolves into the cutoff") shows it is **not** monotone even numerically (a lattice-resolution effect, but enough to kill a naive monotone-coupling proof). **Dead.**

**(D) Correlation inequality (FKG / Griffiths / Ginibre) `[NG]` — wrong object.**
FKG-type positivity requires the action to be a sum of *ferromagnetic* pair interactions; the Wilson action expanded over **all** irreps (`P3-FULLREP-KP`) has sign-indefinite higher-character contributions — it is **not FKG-positive**. And correlation inequalities bound *expectations of local observables* (`⟨W(C)⟩`), not a Källén–Lehmann **atom weight** `f₁²=ν({E₁})`. No correlation inequality reaches a spectral weight directly. **Dead.**

**Summary `[V-here]`.** The four tools are β-uniform in *form* but none carries the residue *magnitude*. The lower bound `f₁² ≥ r₀(β)` remains available **only** from the convergent polymer expansion, hence only for `β<β*`.

---

## §3 · Exact `β*` and where the argument stops `[T]` (mod constants)

KP criterion (`P3-FULLREP-KP` §2.3, a=1): `Σ_{P∋x}|w(P)|e^{|P|} ≤ λ_aΦ(1)/(1−λ_aΦ(1)) ≤ 1` **iff `λ_aΦ(1) ≤ ½`**, with
```
  Φ(a) = e^{β}(e^{βe^{a}}−1),      λ_a = eΔ = 20e ≈ 54.4   (crude, unoptimised; Δ=20).
```
Solving `λ_aΦ(1)=½` (witness §1):
```
  β*(crude, a=1) = 3.36×10⁻³,           log10(2.5/β*) = 2.87  orders below the crossover;
  β*(optimised a) = 8.9×10⁻³,           log10(2.5/β*) = 2.4;
  β*(benchmark, SZZ) = 1/48 = 2.08×10⁻², log10(2.5/β*) = 2.08.
  at β=2.5:  λ_aΦ(1) = 5.9×10⁵  — the KP series is violated by six orders.
```

**Where it stops.** The argument stops *at the KP boundary `λ_aΦ(a)=½`*: beyond it the activity sum `Σ|w(P)|e^{a|P|}` diverges (animal + rep-label entropy overtakes the factorial-suppressed activity), so the expansion gives **no finite bound on `R_vac`**, and §2 shows no β-uniform inequality supplies one. **The gap `β* → 2.5` (≥2 orders in `β`, ≥5–6 orders in `λ_aΦ`) is not closable by constant-optimisation** (witness §3: optimising `a` and `λ_a` buys a factor `≲3`, i.e. `<0.5` orders). It is a **method** wall, not a physical transition — the residue is expected to persist to `β≈2.5`, but no known β-uniform *certificate* does.

---

## §4 · What the base can and cannot be

- **Can be:** the factorization `f₁² = W_tube·R_vac` with `W_tube` explicit for all `β` and `R_vac∈(0,1]` positive — a clean reduction of the extension problem to **one positive scalar** with no spectral/representation content. This is a reusable base piece: any future tool need only lower-bound `R_vac`.
- **Cannot be:** a β-uniform extension of the residue floor. The four named tools are provably insufficient (§2). The residue radius `β*` is the Kotecký–Preiss radius (§3).
- **The only exit** (`[O]`): a **rigorous renormalization group for 4D lattice gauge** reaching the crossover (Brydges–Slade-level machinery adapted to gauge theory). Such an RG **does not exist** for 4D nonabelian gauge. Building it *is* the Clay-scale tool, not a β-uniform inequality. So the base stops here: **the residue is a strong-coupling object; carrying it to the crossover is the wall itself.**

---

## §5 · Honest bottom line

- **β\* computed:** `3.4×10⁻³` (crude KP), `≤8.9×10⁻³` (optimised grading), `1/48` (benchmark) — **≥2 orders below the crossover `2.5`.**
- **None of the four β-uniform tools** (RP, Gaussian domination, monotone coupling, correlation inequality) gives an L-uniform residue floor past `β*`: RP is magnitude-blind (`≥0`), Gaussian domination is an upper bound, monotone coupling is unproven (and fails for a KL atom), correlation inequalities bind expectations not spectral weights.
- **The base piece** is the tube–vacuum factorization `f₁²=W_tube·R_vac`: the whole `β`-dependence lives in one positive scalar `R_vac`, which no known tool bounds below past `β*`.
- **Acid test passes** (`β→0 ⇒ u→0 ⇒ f₁²→0`, free field). **Labels:** the factorization `[V-here]`; the four non-extensions `[T]`/`[NG]` as marked; the radius `[T]` mod crude constants. **No QED; Clay 0%.**

*Witness:* `research/unconditional/w7/base_infinite_vol.py` (numpy/scipy, <1 s, prints β*, the KP table, the optimisation floor, and the four-tool verdict).
*Refs:* `PROOF-P3-UNCONDITIONAL.md` §1–§5; `P3-FULLREP-KP.md` §2.3, §4; `P3-FOURWALL-DERIVATION.md`; `P3-CYLINDER-RESIDUE.md`; `P3-DRESSED-MASS.md`; Kotecký–Preiss CMP **103** (1986) 491–498; Osterwalder–Seiler Ann. Phys. **110** (1978) 440; Fröhlich–Simon–Spencer (Gaussian domination / RP, 1976); Shen–Zhu–Zhu CMP **400** (2023) 805 (radius `β<1/48`); Brydges–Slade (rigorous RG, absent for 4D gauge).
*Status: **β-uniform extension of the residue — [NG]; base piece = factorization + exact β\* — [T]/[V-here].** I2 / `a→0` / `β→2.5` remain open. No QED.*
