# ROOF-DEFECT-CFT — defect CFT / boundary bootstrap gives **renormalon-free, protected** data, but in the **wrong regime** (requires `m=0` = gaplessness) and the **wrong channel** (1-form/defect-local, not the bulk `0⁺⁺` gap): **not the roof**

> **VERDICT (no QED).** The `0⁺⁺` operator `tr F²` is `'t Hooft`-neutral (`0`-form, center-singlet), and the
> only extended magnetically-charged defect is the `'t Hooft` loop `B(C)`, a **1-form** object. Defect CFT
> is a *genuinely renormalon-free, protected, in-principle bootstrap-computable* framework — but it fails the
> `0⁺⁺` gap for **two independent structural reasons**. **(R1) Regime** `[T]/[NG]`: a conformal defect CFT
> requires a **conformal bulk**; 4D pure YM is asymptotically free in the UV and **gapped in the IR**, with
> **no IR fixed point** — so applying defect CFT in the IR **presupposes the gaplessness the gap is about**.
> **(R2) Channel** `[T]`: the `'t Hooft` loop carries **1-form center charge**; its defect data (displacement
> dimension `Δ_D=d−p+1`, defect operators, bremsstrahlung) live in the **1-form / defect-local** sector, and
> the **bulk `0⁺⁺` gap is not a defect operator dimension**. So defect CFT inherits *both* the channel
> mismatch (`BASE-ORDER-PARAMETER`) *and* a regime obstruction unique to it. **Clay 0%.**

---

## §0 · Labels, iron laws, acid test (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**IL-i** non-positivity · **IL-ii** nonabelian · **IL-iii** no small parameter. **Acid**: must be `=0`/vacuous on free massless *and* fail on the gapped-atomless `:φ²:`.
The gap object is the lightest `0⁺⁺` KL atom `f₁²=ν({m₁²})=|⟨0|tr F²|G⟩|²` of the **local 0-form** operator `tr F²`.

---

## §1 · The question, stated precisely (and the precise channel)

In 4D pure `SU(N)`, the classical global structure is the **1-form center symmetry** `Z_N`. Its electric
probe is the Wilson loop `W(C)`; its magnetic dual is the `'t Hooft` loop `B(C)` (the disorder/twist defect
of the 1-form symmetry). Both are **line operators** = codimension-3 defects. `tr F²` is a **local 0-form**
operator, **center-neutral** (adjoint trace, zero `N`-ality).

**Question.** Defect CFT assigns to `B(C)` (viewed as a conformal defect) a set of **protected, renormalon-free**
data: the **displacement operator** `D^i` (dimension `Δ_D = d−p+1 = codim+1`), the **defect local operators** and their
anomalous dimensions, defect anomaly coefficients, and a **bremsstrahlung function**. Does any of this give a
**renormalon-free, bootstrap-computable** statement about the **bulk `0⁺⁺` gap**? Be precise about the channel.

---

## §2 · (R1) The regime obstruction `[T]/[NG]` — defect CFT needs a conformal bulk; pure YM's IR is gapped

A conformal defect CFT is the statement that the **bulk is conformal** and the defect preserves a residual
`SO(p+1,1)×SO(d-p)` symmetry; the whole framework (defect OPE, protected `Δ_D`, bootstrap) is built on
**conformal symmetry**. 4D pure `SU(N)` YM has:

- **UV**: asymptotically free (`β(g)<0`), approximately free — not a nontrivial CFT (free + logs);
- **IR**: **gapped** (the very statement `Δ>0`) — no fixed point, no conformal symmetry;
- **no Banks–Zaks window** (pure YM has no zero of `β` other than the trivial one; no fermions to open a window).

So there is **no conformal regime in which defect CFT applies** to the gap. Worse, the logic is circular:
**a conformal defect CFT in the IR would imply scale invariance = no mass gap**. Applying defect CFT to the
`0⁺⁺` gap **presupposes its negation**. The witness (`W1`) makes this concrete: a defect one-point function
(`∝` the bulk propagator) obeys the **defect-CFT power law only at `m=0`**; at `m>0` it decays **exponentially
with rate `m`** (effective rate `-d log G/dr → m`). The "defect CFT" scaling *is* the gaplessness.

> **Contrast.** The 3D Z₂ / 3D Ising **defect CFT** is solvable — but only **at the critical point `β_c`**
> (where it is a CFT, i.e. **gapless**), not in the gapped ordered/disordered phases. The gap lives *away*
> from `β_c`, where the defect is **massive** (area law), not conformal. Same acid: defect CFT describes the
> critical point, never the gap.

---

## §3 · (R2) The channel obstruction `[T]` — defect data are 1-form/defect-local, not the bulk `0⁺⁺` gap

Even granting a conformal defect, the channel is wrong. The `'t Hooft` loop `B(C)` is **1-form charged**
(the disorder operator of the center `Z_N`). Its defect CFT data are, by construction, **1-form-sector /
defect-local**:

- **Displacement operator** `D^i`: the response of the defect to transverse deformation; its dimension
  `Δ_D = d−p+1 = codim+1` is **protected by the transverse conformal symmetry** and is **independent of the
  bulk mass / coupling** (`W2` table: line-in-4D `Δ_D=4`, surface-in-4D `Δ_D=3`, …). → **renormalon-free but
  gap-blind** (the `'t Hooft` loop, a `p=1` line in `D=4`, has `Δ_D=4`).
- **Defect local operators** `ô(τ)`: operators living *on* the defect; their dimensions are **defect**
  anomalous dimensions, not bulk operator dimensions.
- **Bremstrahlung** `B`: the strength of defect emission; a defect-local, protected quantity.

Now put `tr F²` next to `B(C)`. Because `tr F²` is **center-neutral** (`W3`: `lk(point,loop)=0` ⇒ `'t Hooft`
phase `=1`), it has a well-defined **one-point function in the defect background** `⟨tr F²⟩_{B(C)}` — a
genuine `0`-form observable. **But** (i) this is a *response/condensate*, not a **spectral edge**; it does not
floor the gap (exactly `BASE-ORDER-PARAMETER`'s "one-point-in-defect-background ≠ gap"); and (ii) the
**operators it mixes with** in the defect OPE are **defect operators**, whose dimensions are §3's defect
data — **not** the lightest **bulk** `0⁺⁺` state. The bulk `0⁺⁺` gap is the **pole of the bulk two-point
function** `⟨tr F²(x) tr F²(0)⟩`, a bulk spectral quantity that **no defect operator dimension equals**.

So defect CFT **inherits the channel mismatch**, and adds one of its own (the defect-local/bulk distinction).
The magnetic sector stays 1-form; the `0⁺⁺` gap stays 0-form bulk.

---

## §4 · What defect CFT *does* give (and why it is still not the roof)

Defect CFT is the **right language for the 1-form sector** and delivers genuinely renormalon-free,
bootstrap-computable data — a real asset, just in the wrong place:

- **displacement multiplet** (`Δ_D = codim+1`, protected) — renormalon-free;
- **defect anomalous dimensions** and **defect anomaly coefficients** — renormalon-free;
- (in SUSY / holographic cases) the **bremsstrahlung function** `B` — exactly computable.

None of these is the bulk `0⁺⁺` gap. They are the **1-form-sector / defect-local** analogue of the "protected
data" that `BASE-DUALITY` already showed to be **`β=2.5`-uncontrolled** and orthogonal to `0⁺⁺`. So defect
CFT **relocates** the same wall: it makes the 1-form sector crisp and renormalon-free, and leaves the `0⁺⁺`
gap untouched.

---

## §5 · Witness (`W1–W3`) and acid test

`w7/roof_defect_cft.py` (numpy + `scipy.special.k0`, `<1 s`, EXIT 0):

- **W1 (regime):** effective decay rate `-d log|G|/dr` of the 2D propagator → `0.527` (`m=0.5`), `1.028`
  (`m=1.0`), `~0.02` (`m=0`, log-law). Defect-CFT power law **only at `m=0`**.
- **W2 (channel, protected):** `Δ_D = d−p+1` for line/surface in `D=4,3` (e.g. `'t Hooft` line-in-4D `Δ_D=4`);
  **independent of bulk mass** ⇒ gap-blind.
- **W3 (0-form neutrality):** `exp(2πi·lk(point,loop)/N)=1` for `N=2,3`; the `0⁺⁺` operator is 1-form-neutral.

**Acid test.** Defect CFT as a *tool for the gap* **fails the gapped adversary**: at `m>0` (`W1`) the
defect-CFT scaling is absent (exponential, not power-law), so it is **vacuous exactly where the gap lives**
(and the `:φ²:` adversary, being gapped, kills it the same way). It is **not** acid-clean for the gap.

---

## §6 · Verdict and consequence for the roof

- **Defect CFT does NOT give a renormalon-free, bootstrap-computable statement about the `0⁺⁺` gap.** Two
  independent obstructions: **(R1) regime** — it requires a conformal (gapless) bulk that pure YM does not
  have in the IR; its application to the gap presupposes gaplessness. **(R2) channel** — its data are
  1-form/defect-local; the bulk `0⁺⁺` gap is not a defect dimension. `[T]/[NG]`.
- **It does give** genuine, protected, renormalon-free **1-form-sector** data (displacement `Δ_D=d−p+1`, defect
  anomalous dimensions, bremsstrahlung) — a real tool, but orthogonal to `0⁺⁺`, consistent with
  `BASE-ORDER-PARAMETER` / `ROOF-0FORM-DUAL`.
- **Consequence:** "continue along existing theory" with defect CFT confirms — by yet another independent
  route — that the roof is **neither an order parameter, nor a defect-CFT datum**: it is a **bulk, 0-form,
  renormalon-free measure functional** on `μ` (`ROOF-SYNTHESIS` §4). Defect CFT sharpens the 1-form sector and
  leaves the `0⁺⁺` gap exactly where it was: `H1`.

**Clay 0%; no QED.**

*Refs:* `ROOF-0FORM-DUAL.md` (0-form 't Hooft impossible; `tr F²` is the anomaly, not an algebra),
`BASE-ORDER-PARAMETER.md` (channel mismatch; defect one-point ≠ gap), `ROOF-SYNTHESIS.md` §4 (the surviving
type), `ROOF-TRANSFER.md` (B2: single phase ⇒ no order parameter). Conformal-defect displacement dimension
`Δ_D = d−p+1` (M. Billò, V. Gonçalves, E. Lauria, M. Meineri, "Defects in conformal field theory",
*JHEP* **04** (2016) 091, arXiv:1601.02883 — verified). Witness: `w7/roof_defect_cft.py`.
