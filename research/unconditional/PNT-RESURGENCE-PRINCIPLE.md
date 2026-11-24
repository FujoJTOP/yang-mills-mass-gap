# PNT-RESURGENCE-PRINCIPLE — no physical principle pins the transseries parameter `σ`; the freedom is genuine for a single IR renormalon (no saddle, no second scale)

> **VERDICT (no QED).** The `u=2` IR-renormalon ambiguity is the **computable Stokes constant** (fixed by
> the large-order coefficients `a_n ~ n!`); the transseries parameter `σ` (the gluon-condensate / OPE power-
> correction coefficient) is a **free integration constant**. I audited every candidate selection principle
> — positivity, analyticity, OPE-consistency, large-`N` matching, 't Hooft anomaly, trace anomaly, a
> resurgence web, a solvable deformation — and **none pins `σ` to a value**. The exact obstruction is
> structural: **an IR renormalon has no saddle point** (unlike a UV renormalon / instanton, where the median
> sum *is* physical), and pure YM at `θ=0`, `T=0`, `V=∞` has **no second scale** to form a resurgence web.
> Positivity gives only a **range** (`|σ| ≤` ambiguity width), never a point. **Freedom genuine.** And —
> scope note — even a pinned `σ` would **not** lay the roof: the condensate is **not** the KL atom
> (`BASE-FIELDMEASURE`, `ROOF-CONDENSATE-FLOOR`). **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric (toy).

---

## §1 · Where `σ` lives

For an asymptotically-free observable `F(α) = Σ a_n α^n` whose Borel transform has an IR pole at
`u = 2` (dimension-4), the **transseries** is

```
F(α) = PV[a_n] + σ · e^{−2s} · (1 + …),        s = 1/(2β₀α),
```

where:
- the **Stokes constant** `Res` (the Borel residue) is **fixed by perturbation theory**: `a_n ∼ Res·Γ(n+1+ν)c^n`;
- the lateral Borel sums differ by `|S₊ − S₋| = 2π·Res·e^{−2s}` (the **ambiguity**);
- `σ` is a **real integration constant** multiplying the power correction `e^{−2s}`.

`ROOF-RESURGENCE` established that `{a_n}` fixes `Res` but **never** `σ`: `σ` does not appear in any `a_n`.
The question here: **is there a physical principle (beyond perturbation theory) that fixes `σ`?**

---

## §2 · The principle audit

| principle | what it gives for `σ` | why |
|---|---|---|
| **Positivity** `ρ(s) ≥ 0` | a **range** `\|σ\| ≲ π·Res·e^{−2s}`, **not a point** | `ρ(s)` is renormalon-free, but the OPE *split* `PV + σ e^{−2s}` is prescription-dependent; two prescriptions give the same `ρ` with `σ` differing by `2π Res` — positivity shifts with the split |
| **Analyticity** (cut plane) | nothing | already used to define the dispersion relation; does not touch the subtraction constant |
| **OPE consistency** | nothing | the OPE is asymptotic; it *defines* the coefficients it would need to fix `σ` (circular) |
| **Large-`N` matching** | nothing *yet* (circular) | would fix `σ` iff the 4D master field existed; but constructing it **is** the open target (`ROOF-LARGEN-MASTER`) — one cannot use it to pin `σ` without already having it |
| **'t Hooft anomaly** | nothing at `θ=0` | the center 1-form symmetry is **unbroken** in the confined `θ=0` vacuum (the anomaly acts at `θ=π`, forcing deconfinement — not the `θ=0` condensate) |
| **Trace anomaly** `θ^μ_μ=(β/2g)G²` | nothing | the identity ties `⟨G²⟩` to the vacuum-energy density, which carries the **same `Λ⁴` renormalon** — it transports the ambiguity, does not remove it (`ROOF-OPE`) |
| **Resurgence web / cross-constraints** | nothing (no web) | a web pins `σ` only if **≥2 coupled singularities** constrain each other; pure YM at `θ=0`, `T=0` has a **single** IR renormalon and no second nonperturbative scale |
| **Solvable deformation** (large-`β₀`, large-`N_f`, finite `T`, `θ≠0`, finite `V`) | fixes `σ` **in the deformed theory only** | the deformation supplies the missing second scale; but `σ_{(deformed)}` does **not** equal `σ_{pure YM}` — decompactifying / `N_f→0` is itself the wall |

**Eight candidates; zero pin `σ`.** Every one either gives a range, is circular, or works only in a
deformed theory whose limit back to pure YM is the original wall.

---

## §3 · The exact obstruction (why the freedom is genuine)

The decisive structural distinction:

- **UV renormalon / instanton (`u<0`)**: tied to a **saddle point** of the action. The instanton action fixes
  the *scale*, and in many theories (e.g. the double-well, the Painlevé equations) the **median resummation
  is the physical answer** — the Stokes constant *is* the physical nonperturbative amplitude. A principle can
  pin `σ` because `σ` is a physical saddle contribution.
- **IR renormalon (`u=2`)**: **no saddle point** — the renormalon is *not* an instanton; it is the perturbative
  image of an OPE power correction, i.e. of a *scheme-dependent* coefficient. There is no configuration to
  fix `σ`; the "ambiguity" is the scheme-dependence of the subtraction.

⇒ For the **IR** renormalon — which is exactly the gluon-condensate case — **`σ` is not a physical number**;
it is a subtraction convention. The physical content is `ρ(s)` (renormalon-free), and `ρ` is **floor-blind**
(`ROOF-OPE` = `SLOT-BOOTSTRAP`: `Δ_min(K)=0 ∀K`). **The freedom is genuine, not a gap in our cleverness.**

---

## §4 · Scope: even a pinned `σ` would not lay the roof

An important honesty point, independent of the above. `σ` is the **OPE coefficient of the gluon condensate**.
`BASE-FIELDMEASURE` / `ROOF-CONDENSATE-FLOOR` established that the condensate is **not** the KL atom
`f₁² = |⟨0|tr F²|G⟩|²`: the condensate can be nonzero while `f₁² = 0` (the `:φ²:` adversary). So even if a
principle pinned `σ`, it would give a **condensate**, not the **atom** — and the gap needs the atom
(`H*`). **Pinning `σ` is neither necessary nor sufficient for the roof.** This bounds how much the resurgent
route could ever deliver: at most a clean renormalon-free *scale*, not the mass gap.

---

## §5 · Witness (`w7/pnt_resurgence_principle_witness.py`, runs clean `py -3.12 -I`)

Transparent toy with a `u=2` pole (normalized `c=1/2`, `a_n = Γ(n+1)c^n`):
- **(A)** `a_n / (Γ(n+1)c^n) ≡ 1` for all `n` — the coefficients fix the residue `Res=1`, and **no `σ`
  appears**; adding `σ e^{−1/(cα)}` changes no `a_n`.
- **(B)** lateral-sum jump `= 2π Res e^{−1/(cα)}` (e.g. `2.85×10⁻⁴` at `α=0.2`, `0.85` at `α=1`); `F = PV +
  σ e^{−1/(cα)}` is **real for every real `σ ∈ {0, −0.5, 1.7, 10}`** — reality selects nothing.
- **(C)** positivity band `|σ| ≤ π Res e^{−1/(cα)}` (e.g. `1.4×10⁻⁴` at `α=0.2`) — a **range, not a point**.

This is a toy (not YM); it isolates the logical structure: `{a_n}` → residue; positivity → range; reality →
free; **single IR renormalon → `σ` free.**

---

## §6 · Bottom line

- **No physical principle pins `σ`.** Eight candidates audited; all fail (range / circular / deformed-only).
- **Freedom is genuine** for a single IR renormalon: no saddle (unlike instantons), no second scale at
  `θ=T=0, V=∞`. `σ` is a subtraction convention, not a physical number.
- **What would pin it:** a **second nonperturbative scale** forming a resurgence web (finite `T`, `θ≠0`,
  finite `V`, or a solvable deformation) and a **controlled limit** back to pure YM — but that limit is
  itself the transport wall (`H1`). Or the **4D master field** (`ROOF-LARGEN-MASTER`), which is the open
  target, not an available pin.
- **Scope:** even a pinned `σ` gives the condensate, not the KL atom — neither necessary nor sufficient for
  the roof.
- **Net:** this closes the "use resurgence to fix the condensate" hope cleanly. The resurgent route is
  **dead as a roof-strut**; its surviving value is diagnostic (the ambiguity = computable Stokes constant).
  **Clay 0%; no QED.**
