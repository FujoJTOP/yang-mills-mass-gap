# KEY-DIRECT-DERIVATION — a genuine direct assault on `f₁² ≥ c > 0` ( ⟺ `G″(k→0) ≥ λ`) at the crossover, and the exact step each attack breaks

> 2026-11-24 · **Direct constructive attempt**, per charter §1 ("造工具 + 硬算", do not stop at "hard"). Five concrete attacks on `f₁² = ν({m₁²}) ≥ c > 0` (residue face) and/or `m_eff²(k→0)=G″(k→0) ≥ λ > 0` (curvature face) for the lightest `0⁺⁺` KL atom of `O=¼tr F²`, **uniform in `(L,a)`, at `β≈2.5` (`α_s~O(1)`), failing on the free massless field**. Each attack is carried through its arithmetic (witness `research/unconditional/w7/key_direct_derivation.py`, `[V-here]`), acid-tested, and pushed to the **precise** iron law / estimate it breaks.
>
> Labels (charter §2): **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a *named class* · **[O]** exact missing step · **[V-here]** computed this round · **[R]** external as-reported. **This does NOT prove Clay. It is an honest map of the wall from a fresh angle. No QED on anything unproven.**

---

## VERDICT (up front)

> ### **No attack gives a valid non-positivity, nonabelian, crossover lower bound. The best attack (I, Savvidy) passes iron laws (i)+(ii) and the acid test but dies on the *minorant inequality* `Φ ≤ m₁` — it is the scale of an unstable saddle, not a spectral floor. Every other attack dies on a precise, named line: (II) the variance split is an identity that cancels → H1; (III) variational machinery needs a small parameter or `E₂`; (IV) the `0⁺⁺` spectral density is sign-definite → upper bounds only; (V) the one positive nonabelian invariant (condensate) gives the wrong direction / is the gap (SGL).**
>
> **What genuinely survives (partials, §6):** (a) `[T]` the Savvidy negative mode is **exactly nonabelian-only** (U(1): 0, SU(2): ≥1) — a clean discriminator, but an instability (wrong side); (b) `[T]` an exact **U(1) curvature theorem** (`D≤0 ⇒ G″≥S≥0` at all couplings) that *coexists with `m₁=0`* — the sharpest demonstration that the curvature/positivity face cannot gap; (c) `[T]` the exact **cancellation** `S−D ≡ E[h]−Var(f)`, so the discrepancy decomposition carries no inequality — the gap is literally "stiffness beats force-variance"; (d) `[O]` the single named missing inequality (§7).

---

## §0 Targets, the `H*` vs `H1` subtlety, and the iron laws made operational

**Two faces of the target, and why "equivalently" is a Gaussian idealization.**
- *Residue face* (the true Clay object): `f₁² = ν({m₁²}) = |⟨0|O|G⟩|² ≥ c > 0`, `m₁=inf supp ν`, with isolation `I2: E₂−m₁²≥δ`. This is `H* = I1∧I2`. **[T, definition]**
- *Curvature face*: `G″(k→0) = 1/χ`, `χ=∫s⁻¹dν` (inverse zero-momentum susceptibility). **[T]** The exact squeeze is `f₁²/χ ≤ m₁² ≤ M₀/χ` (H1-IS-IT-THE-GAP). Hence `G″>0` and `f₁²>0` are **equivalent only for a single-atom (Gaussian) `ν`**; for the gauge-invariant composite they are *different functionals* — `G″` is a residue-weighted inverse harmonic moment that only *upper*-bounds `m₁²` (harmonic mean ≥ min). So the task's "equivalently" is the Gaussian-LM limit; in general a valid derivation needs **both** faces. We attempt both and flag where they diverge. **[T]** (this is the corpus's `H1 ≠ H*`; §6b)

**The three iron laws, operationalized.** Any derivation of a lower bound must be:
- **(IL-i) non-positivity:** not reducible to an inequality that the free massless field also satisfies (`f₁²=0`, `G″=0` there). Positivity/convexity/LSI/transport/BL are *Schur-concave*: they cap totals, hence only *upper*-bound the edge atom — vacuous. **[T]** (SCHUR-CONVEX)
- **(IL-ii) uses the nonabelian interaction:** else the Abel no-go — 4D pure U(1) is massless. **[T]**
- **(IL-iii) no small parameter:** must hold at `α_s~O(1)` (`β≈2.5`); no weak-coupling, no `β→0`. **[T]**

**The free-field acid test.** A candidate `Φ` fails the acid test if `Φ(free massless) > 0` (vacuous), and is *vacuously acid-clean* if it also gives `≤0` for the free **massive composite** `:φ²:` (gapped, `f₁²=0`) — the discriminating case. We test against both.

---

## §1 Attack I — a **non-positivity minorant from the Savvidy negative mode** (the best attack)

**The idea (charter's hint #1).** The one genuine *signed / indefinite* structure of nonabelian YM — not available to U(1) — is the spin–field coupling `2 ad F_B` in the gluon fluctuation Hessian `−D² + 2 ad F_B`. In a background field `B` it produces a **negative mode** whose magnitude `|ω₋²| ∝ g²·C₂(adj)·B` is fixed in sign by a *positive* nonabelian invariant (`C₂(adj)=N`, `b₀=11N/48π²`). Take `Φ_Sav(B) := √(2gB)` (the Savvidy/Nielsen–Olesen "magnetic mass") as a candidate minorant.

**Arithmetic (exact, `[V-here]` §A).** Charged massless spin-1 in a constant `B`, lowest Landau level: `E² = 2|q|B(n+½) + 2qB\,s_z`, `s_z∈{−1,0,+1}`.
- **U(1) pure gauge:** the photon is *neutral* (adjoint = trivial rep, `q=0`). No charged vector exists → the only gluon is the massless photon `−∂²` → spectrum `{0,…}`, **0 negative modes.** **[T]**
- **SU(2):** background along `T₃`; the two off-diagonal gluons carry charges `±g` (verified: `ad(T₃)` eigenvalues `{0,±i}`). The `s_z=−1` mode has `E² = gB − 2gB = −gB < 0` → **2 negative modes** at `gB=1`. **[T]** **[V-here]**

So the negative mode is an **exact, nonabelian-only** object (IL-ii satisfied) and it is **signed, not a positivity** (IL-i satisfied *in form*). Acid test: `Φ_Sav(0)=0`; `Φ_Sav→0` as the background → 0. The **free field / U(1) / free-massless all give `Φ_Sav=0`** — acid-clean. **`[V-here]` §B; [T]**

**Why it is NOT a minorant — the exact failing step.** The candidate must satisfy `0 < c ≤ Φ ≤ m₁`. It fails `≤ m₁`:
1. As a **function** of the background, `Φ_Sav(B)=√(2gB)` is monotone increasing and `→0` at `B→0`. The physical (confining) vacuum has **no uniform background** (`⟨B⟩=0`) yet `m₁>0`; so `Φ_Sav(0)=0 ≪ m₁`. The function is not a lower bound at the physical point. **`[V-here]` §B**
2. Evaluated at a putative one-loop minimum, `m_Sav=√(2gB_min) ~ Λ` — but this is the **same dimensional-transmutation scale** as `m₁` (an *identity of scales*, not an inequality), and it is the **curvature of an unstable saddle** (it has the negative mode), i.e. a *width / instability* scale, not a *spectral floor* of the ground state. There is **no operator inequality** `m₁ ≥ √(2gB)`. **[T]**

**Which iron law / estimate breaks.** IL-(i) and IL-(ii) **pass**; IL-(iii) is not the issue. The break is the **minorant estimate `Φ ≤ m₁` itself**: the Savvidy scale belongs to the *wrong vacuum*. To convert "background-field scale" into "true-vacuum spectral floor" requires the nonperturbative **Savvidy → physical transport**, which is *exactly* the missing crossover machine (H1/(CP)/confinement; BRIEF §4 "transport ≡ confinement"). So Attack I is the best attack and it dies on the one step the whole campaign dies on — but now stated in the cleanest possible form: **the single signed nonabelian structure we have is an instability, and an instability scale decreases precisely where the gap is largest.** **[O]**

---

## §2 Attack II — the exact variance decomposition `G″ = S − D`, with a direct bound on the defect `D` (charter's hint #2)

**The idea.** `G″ = S − D`, `S=E[h_n−⟨K⁻¹c,c⟩|0]` (Gaussian–Schur curvature, `K=Π₀HessV Π₀`), `D=Var(f|0)−E[⟨K⁻¹c,c⟩|0]` (fiber BL defect; `D>0 ⇔ K indefinite ⇔ Savvidy`, and `K` is indefinite **only** for nonabelian — `[T]` from §1: `2adF_B=0` for U(1)). The hope: bound `D<S` using the *nonabelian* indefinite structure that U(1) lacks.

**The arithmetic that kills it — the Schur term cancels identically.** Add and subtract `⟨K⁻¹c,c⟩` and expand:
```
S − D = [E[h] − E⟨K⁻¹c,c⟩] − [Var(f) − E⟨K⁻¹c,c⟩] ≡ E[h] − Var(f) = G″.
```
The `K⁻¹` term appears with the **same sign** in `S` and in `D` and **cancels exactly**. Therefore:
- `D < S ⟺ Var(f)<E[h]` — **tautologically** `G″>0`, i.e. **H1 itself**. The split is a *rearrangement*, not an inequality. **[T]**
- The **nonabelian indefinite `K` shifts `S` and `D` by the same `−E[⟨K⁻¹c,c⟩]`**; it cannot separate them. The Savvidy defect is *adversarial*, not helpful: `D>0` (indefinite fiber) **reduces** `G″ = S−D`. So the nonabelian interaction is simultaneously **necessary** (U(1) has no indefinite fiber, and is massless) and **adversarial** (in the curvature face it pushes `G″` down). **[T]**

**Toy corroboration (`[V-here]` §D).** For `V=½m²x²+½ky²+gxy²+¼λy⁴` the exact marginal gives `G″(0)=m² − g²Var(y²|0)` (verified to 5 digits against quadrature), the identity `S−D=E[h]−Var(f)`, and:
| `m²` | `g` | `λ` | `G″` (exact) | reading |
|---|---|---|---|---|
| 1.0 | 0 | 0 | **1.000** | free Gaussian; `→0` as `m²→0` (**acid passes**) |
| 1.0 | 0.6 | 0.5 | **0.818** | nonabelian coupling reduces `G″` |
| 0.0 | 0.6 | 0.5 | **−0.182** | free-massless marginal mode driven **negative** (Savvidy marginal instability) |

So `G″` is literally `stiffness(m²) − (nonabelian) force-variance`; proving `G″>0` is `m² > g²Var`. That is a **variance upper bound** `Var(f|0) < E[h|0]`.

**Which iron law / estimate breaks.** It breaks on **IL-(i)**: any bound on `Var(f|0)` is a Schur-concave positivity/variance inequality; the free field **saturates** it (`D=S=0`, `G″=0`), so it cannot be *strict*, and proving it for YM is exactly `H1`. It also fails to supply any *new* lever because the nonabelian `K` cancels (§ above). **Exact failing step: `Var_μ(∂_nV|a=0) < E_μ[h_n|a=0]`, uniform in `(L,a)` — this is `H1` verbatim.** **[O]**

---

## §3 Attack III — a **variational construct** whose value is `≤ m₁²` and `>0` only for the interacting theory (charter's hint #3)

**The idea.** Build `Φ = ` a variational lower bound on `m₁` via `Φ≤m₁`, `>0` only for YM. The corpus's lower-bound machinery is exactly five methods (E2-FLIP-MINORANT): enclosure (Gershgorin/Temple/Weinstein), determinant/trace, cluster/KP, comparison/domination, transfer-positivity. All were tested there. We re-derive the two that are *dynamical* (not positivity) and could in principle use the interaction:

**Arithmetic (`[V-here]` §C).**
- **Rayleigh / min–max:** `min` over 2×10⁵ random trials of `R=⟨ψ,Hψ⟩` on `(vac)^⊥` stays `≥E₁`; it *equals* `E₁` only at the exact eigenvector. Every trial value is an **upper** bound on `E₁`; a "minorant `Φ=R≤E₁`" is **false**. **[T]** (witness C1)
- **Temple:** `b(ρ)=R−Var/(ρ−R)` is sound **iff `ρ≤E₂`**. With *consistent* 2-level data, `ρ=E₂` reproduces `b=E₁` exactly; `ρ=1.3E₂` and `ρ=3` give `b=1.0095, 1.0150 > E₁` — **overclaims** (false too-high "lower" bounds). A certified `ρ≤E₂` **is the second level (isolation `I2`)**. **[T]** (witness C2; `⇒` Theorem ENCL)
- **Bogoliubov:** `m₁ ≥ gap(H₀) − ‖V‖ = 0 − 0.7 = −0.7 ≤ 0`: the free Yang–Mills `H₀` is **gapless** (`gap(H₀)=0`) and the nonabelian interaction is `O(1)`. Vacuous. **[T]** (witness C3)

**Which iron law / estimate breaks.** Two distinct precise failures:
- Rayleigh/Weinstein is an **upper** bound (direction wrong — not a minorant at all).
- Temple needs `ρ≤E₂`: this is **isolation `I2 = H*`**, a strict strengthening of the gap — circular by one rung.
- Bogoliubov needs `gap(H₀)>‖V‖`, i.e. **`α_s` small** — breaks **IL-(iii)** at `α_s~O(1)` (and the free `H₀` is gapless so no small-parameter expansion has a positive base). **Exact failing step: there is no variational lower bound that does not consume `E₂` (isolation) or a small coupling; the only dynamical lever available (`‖V‖ < gap(H₀)`) is false at the crossover.** **[O]/[NG]** for the class.

---

## §4 Attack IV — a **superconvergence / Weinberg-type sum rule** for the residue (a fresh attempt)

**The idea (not in the death map in this exact form).** In QCD the first Weinberg sum rule `∫[ρ_V−ρ_A]ds = f_π²>0` *lower*-bounds a residue, because the integrand is a **sign-indefinite difference** whose integral is fixed positive. Can an analog force `f₁² ≥ c` in the `0⁺⁺` channel?

**Arithmetic.** Two exact anchors, both sign-*definite*:
1. **Trace anomaly:** `T^μ_μ=(β/2g)trF²`, `β<0` (AF, `b₀=11N/48π²>0`). The dilation low-energy theorem gives `⟨0|¼trF²|G⟩ = f_G m_G² g/(2β)` — an *exact relation* but it contains `f_G` and `m_G` (the unknowns). Homogeneous in `f₁²`; cancels. **[T]** (INEQ-T4)
2. **Spectral moments:** `O≥0` pointwise and `ν≥0` give `∫s^k ρ_0(s)ds = ⟨O(−∂²)^kO⟩_c = ` moments; and the anomaly correlator `⟨T^μ_μ T^ν_ν⟩` is a **positive** spectral density. Hence **every exact sum rule in the `0⁺⁺` channel is a positive functional of `ν`**, i.e. a **moment**, and by Schur-concavity **upper**-bounds the edge atom: `f₁² m₁^{4} ≤ ∫s²ρ_0`, etc. **[T]**

**The missing carrier.** To get a *lower* bound one needs a **sign-indefinite** superconvergence, i.e. a difference of two channels with a common OPE. In pure glue the natural partner of `0⁺⁺` (`trF²`) is `0⁻⁺` (`trF\tilde F`), but `trF\tilde F = ∂_μ K^μ` is a **total derivative (topological)** — its correlator is **contact-only** at leading order, so it carries no matching continuum and no superconvergence lever. There is **no `V−A`-type signed channel** in the scalar glue sector (no axial current in `0⁺⁺`).

**Which iron law / estimate breaks.** **IL-(i)**: the `0⁺⁺` spectral density is **sign-definite** (`O=¼trF²≥0`, anomaly correlation positive), so the exact sum rules are moments → **upper** bounds; the needed sign-indefinite difference has **no carrier** in this channel. This is the **signed-Fourier death** (SIGNED-FOURIER-CORE, score 1/5) re-expressed in sum-rule language: *Positivity of `ρ` gives upper bounds; a residue lower bound needs a signed cancellation, and the only candidate channel (`0⁻⁺`) is a decoupled total derivative.* **Exact failing step: no sign-indefinite superconvergence relation exists for `O=¼trF²`.** **[NG]** for the sum-rule class. **[T]**

---

## §5 Attack V — the **gluon condensate as a positive nonabelian invariant** (charter's hint #1, second form)

**The idea.** `⟨¼trF²⟩>0` is a positive nonabelian invariant whose **sign is fixed** (by `b₀>0`, `β<0`: `⟨T^μ_μ⟩<0 ⇒ ⟨trF²⟩>0`), and it is **absent for U(1) pure gauge** (`b₀=0`). Use it to anchor `m₁` or `f₁²`.

**Arithmetic.** The only dimensionally consistent relation is `m₁⁴ ≳ κ⟨¼trF²⟩`, `κ>0` — the **SGL inequality** (MS-BLIND §4). It is **not known**; positivity gives only the **reverse** (`m₁² ≤ moment ratio`, an upper bound). For the residue: `f₁² ≤ Var(O) = Σ_{n≥1}|⟨0|O|n⟩|²` — again an **upper** bound (the atom is one term of a sum of nonnegative terms). No route from `⟨O⟩` to a *lower* bound on `f₁²`: the condensate is `⟨0|O|0⟩`, the residue is `|⟨0|O|G⟩|²`, different matrix elements, and there is no general inequality `|⟨0|O|G⟩|² ≥ ϑ(⟨0|O|0⟩)` without knowing `|G⟩`.

**Which iron law / estimate breaks.** The needed step **is** the SGL inequality (`m₁⁴ ≥ κ⟨¼trF²⟩`), which the corpus proves *is the gap* (MS-BLIND); and the residue direction is a moment **upper** bound. So it breaks on the **direction (IL-i)** and reduces to the **SGL = gap** identification. **[O]**

---

## §6 What genuinely survives (partials — stated as loudly as the failures)

**(a) `[T]` The Savvidy negativity is exactly nonabelian-only.** U(1) pure gauge has **0** negative modes (neutral photon); SU(2) has **≥1** (`E²=−gB`, verified exactly). This is the cleanest possible realization of IL-(ii): the one signed nonabelian structure is *provably absent* for U(1). **But it is an instability** (right object, wrong side — §1). Partial worth banking as a sharp **discriminator**, not a bound. `[V-here] §A`

**(b) `[T]` Exact U(1) curvature theorem that coexists with masslessness.** For pure U(1) in its (Savvidy-free, `B=0`) vacuum, the fiber Hessian is `K=−∂²≥0` (no spin–field coupling, §1), so by Brascamp–Lieb `D≤0` and `G″≥S` hold; and `S≥0` is the Gaussian–Schur curvature. Yet 4D pure U(1) is **massless** (`m₁=0`, `G″→0` as the mode goes soft). This is the sharpest possible demonstration that the entire `G″≥S` curvature machinery — even as an **exact theorem** — is compatible with `m₁=0`: *curvature cannot gap* (H1-IS-IT-THE-GAP, made exact). It also shows the nonabelian indefinite fiber (`D>0`) is **necessary for a nonzero curvature to be "fighting" something**, while entering with the **adversarial** sign. `[V-here] §A`

**(c) `[T]` The variance split carries no inequality.** `S−D ≡ E[h]−Var(f)`: the discrepancy decomposition is an identity that cancels its `K⁻¹` term. Its *value* is diagnostic (measured `D/S≈2%`, Gaussian-exactness), not a lever. The gap is, exactly and only, "**stiffness > force-variance**" (`m_eff² > Var(∂_nV)`), i.e. `H1`.

**(d) `[O]` The single named missing inequality, in its sharpest curvature-face form.** Combining (b)+(c): the exact missing bound is
> **`W_geo > D − Δ_dyn − X`**, i.e. the **proved-positive nonabelian Haar curvature `W_geo≥0` (Prékopa, `[T]`) must outrun the nonabelian fiber defect `D`** (margin measured `D/S≈2%`).
Here `W_geo` is β-independent and provably `≥0` on every mode incl. `k→0` (`[T]`, PROOF-JACOBIAN-SAVES), and `D>0` is the Savvidy fiber. Proving this **uniformly is `H1`**. **Caveat that makes it two-faced:** `W_geo` is a *curvature* quantity ⟹ per H1-IS-IT-THE-GAP it certifies `G″≥λ` but only **upper**-bounds `m₁`; so even granting it one still needs `f₁²≥c` **separately** (the squeeze). The residue face is *not* free.

---

## §7 The exact step, and the honest verdict

**The single exact missing step** (all five attacks funnel here):
> A **non-positivity, nonabelian, `(L,a)`-uniform lower bound** `f₁² ≥ c > 0` (residue face) **and** `G″(k→0) ≥ λ > 0` (curvature face) at `β≈2.5`, that (i) is not a Schur-concave positivity (else the free massless field — which it must *fail* — is also certified), (ii) uses the nonabelian interaction essentially (the negative mode / `2adF_B` / `b₀>0`), and (iii) needs no small parameter. Equivalently: **`Var_μ(∂_nV|0) < E_μ[h_n|0]`** (curvature face, `=H1`) **and** **`ν({m₁²}) ≥ c`** (residue face, `=H*`), together, uniformly.

**Why each attack lands there, in one line:**
| Attack | IL-(i) | IL-(ii) | IL-(iii) | exact failing step |
|---|---|---|---|---|
| **I Savvidy** | ✓ | ✓ | — | **`Φ≤m₁` unprovable** — scale of an unstable saddle; `Φ→0` where the gap is largest ⇒ Savvidy-transport = H1/CP |
| **II `S−D`** | ✗ | ✓(cancels) | — | **`S−D` is an identity** ⇒ needs `Var(f)<E[h]` = H1; nonabelian `K` cancels |
| **III variational** | — | ✓ | ✗ | Rayleigh=upper; **Temple needs `E₂`(=I2)**; Bogoliubov needs `α_s` small |
| **IV sum rule** | ✗ | ✓ | ✓ | **`0⁺⁺` `ρ` sign-definite ⇒ moments ⇒ upper**; no signed (`V−A`) channel |
| **V condensate** | ✗ | ✓ | ✓ | **SGL = gap**; residue direction is a moment upper bound |

**Honest verdict.** This is **the Clay wall**, seen from five fresh angles. The one genuinely nonabelian, non-positivity, acid-clean object we possess — the Savvidy negative mode — is an **instability** whose scale **vanishes precisely where the gap is largest**; the task's "sign fixed by a positive nonabelian invariant" is realized (§1, §5) but the fixed sign is an **instability sign / a condensate sign**, and both give the **wrong direction** for a *lower* bound on a *spectral floor*. The curvature face is fully closed modulo `H1` (=`Var(f)<E[h]`); the residue face needs, in addition, `H*` (`f₁²≥c`), which positivity cannot supply. **No partial closes the crossover.** The genuine partials (§6) are: the exact **nonabelian-only-ness** of the negativity, the exact **U(1)-massless-but-`G″≥S`** theorem, the **cancellation** identity, and the named missing inequality `W_geo>D−(Δ_dyn+X)` with measured margin `≈2%`. **Clay remains open; 0%.**

---

## §8 References

**Witness:** `research/unconditional/w7/key_direct_derivation.py` (`py -3.12`, single small job; no lattice MC — exact linear algebra + 1D quadrature). Parts A (Savvidy U(1)/SU(2) spectrum), B (Savvidy scale monotonicity), C (Rayleigh/Temple/Bogoliubov direction & inputs), D (exact `G″=m²−g²Var` identity).

**Internal (absolute paths):** `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall = Savvidy concavity, key = `G″(0)`/H1); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five conditions; §4 transport≡confinement; §7 wall); `D:\Dev\math\research\unconditional\H1-IS-IT-THE-GAP.md` (squeeze `f₁²/χ≤m₁²≤M₀/χ`; `G″=1/χ`; H1≠H*); `PROOF-JACOBIAN-SAVES.md` (`W_geo≥0`, Prékopa; exact split `G″=Δ_dyn+W_geo+X`); `INVENT-VARIANCE-INEQ.md` (`G″=S−D`, `D>0⇔K indefinite`, `D/S≈2%`); `SIGNED-FOURIER-CORE.md` (modulus-square kills sign; Bochner=RP=KP; score 1/5); `E2-FLIP-MINORANT.md` (five lower-bound methods; Temple needs `E₂`; `Φ_sc` on `F_sc`); `E2-FLIP-NOGO.md` (Theorem 0; four-horn `C0`; `∀class≠∀possible`); `GOLDS-POS-MASS.md` (positivity pins only `m=0`); `PROOF-P3-UNCONDITIONAL.md` + `P3-FOURWALL-DERIVATION.md` (strong-coupling residue `m₁=4ln(1/u)/a`, `β<1/48`); `TRANSPORT-FIVE-WAY-VERDICT.md` (transport≡confinement; all five routes land on H1); `WASSERSTEIN-OLLIVIER.md` (W₁ dialect, sharpness deficit); `SCHUR-CONVEX.md`, `MS-BLIND.md` (SGL = gap), `GAP-G2-edge.md`, `PROOF-F1SQ.md` (free `f₁²=0`).

**External `[R]`:** Savvidy, *Phys. Lett. B* **71** (1977) 133 (negative mode / magnetic mass); Nielsen–Olesen, *Nucl. Phys. B* **144** (1978) 376 (charged-vector instability); Temple, *Proc. London Math. Soc.* (1928); Weinstein (1932); Kato (1949/1966); Bogoliubov (1960s); Källén (1952)/Lehmann (1954); Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Shen–Zhu–Zhu, *CMP* **400** (2023) 805 (strong-coupling gap); Kotecký–Preiss, *CMP* **103** (1986) 19; Brascamp–Lieb, *Adv. Math.* **22** (1976) 366; Prékopa (1971). Shifman–Vainshtein–Zakharov (glueball sum rules, approximate — not exact inequalities).

*Status: **5 direct attacks attempted and carried to their exact failing step. Best = Attack I (Savvidy signed minorant): passes IL-(i)+(ii) and acid test, dies on `Φ≤m₁` (unstable-saddle scale, Savvidy-transport = H1/CP). Genuine partials: nonabelian-only negativity `[T]`; U(1)-massless-but-`G″≥S` `[T]`; `S−D` cancellation `[T]`; named missing inequality `W_geo>D−(Δ_dyn+X)` `[O]`. Clay NOT proved; 0%. No QED on anything unproven.***
