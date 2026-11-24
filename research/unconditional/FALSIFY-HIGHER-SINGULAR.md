# FALSIFY-HIGHER-SINGULAR — falsification attempt FAILS: higher-order singular kernels do NOT force the atom weight from finite data; they recover `f²` only given the pole LOCATION + a limit (= A3)

> **VERDICT (no QED).** Attempted to falsify "only exit = YM structure or `J→∞` (`H1`)" via stronger
> singularities: `∫|G|^p` (p≠2), `∫|G'|²`, multi-point kernels. **All fail to produce a forced weight lower
> bound from finite data.** Every higher-order kernel behaves the same as TOOL-HARDY's `∫|G|²`: near a pole it
> **diverges as `|x−m₁|^{-q}`** (`q=p` for `|G|^p`, `q=4` for `|G'|²`); over a **fixed** region the integral is
> **`+∞` (presence only, no `f²` value)**; the **divergence RATE** does carry `f^{2p}`, but recovering it
> requires **(i) centering at the pole location `m₁`** and **(ii) the `ε→0` limit** — i.e. **the location (the
> gap) + a limit (= A3)**. Finite point data is matched by an atomless measure (OPE-RHO). **The claim is
> CONFIRMED, not falsified**, for this route.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Witness: `w7/falsify_higher_singular.py` (numpy/scipy, `py -3.12 -I`, clean).

Model: `ν = f²δ(s−m₁) + smooth continuum ρ` on `[m₁+D, E_max]`, `f²=0.30`, `m₁=1`.

---

## §1 · The four facts (`§A–§F` in the witness)

- **(A) `p=2`**: `N₂(ε)=∫_{ε<|x−m₁|<δ}|G|²dx`, `ε·N₂ → 2f² = 0.180` (measured `0.1833→0.1811→0.1805→0.1803`). Recovers `f²` — **given `m₁`**.
- **(B) `p=3`**: `ε²·N₃ → f^{2p}·2/(p−1) = f³ = 0.027` (measured `0.02722→0.02707→0.02706`). **Stronger singularity, same story** — recovers `f³` given `m₁`.
- **(C) fixed region, `+∞`**: `N_p(ε)` over fixed `δ=0.95` grows like `ε^{−(p−1)}` (p=2: `18.3, 181, 1800`; p=3: `272, 2.7×10⁴, 2.7×10⁶`). **The integral's value is `+∞` — presence only; only the rate is finite, and the rate needs centering + `ε→0`.**
- **(D) atomless**: `N₂(ε) → 0.955` finite, no divergence.
- **(E) finite data blind**: a **single complex value** `G(z₀)`, `z₀=1.5+0.3i`, is matched to `3.7×10⁻¹³` by an **atomless** continuum (`c=0.228`, edge `1.124`, no δ at `m₁`). One point value does not force the atom (OPE-RHO / `TOOL-NEVANLINNA`).
- **(F) log-offset route**: `log|G(m₁+u)| + log u → log f² = −1.204` (measured `−1.222→−1.206→−1.204`). Again **given `m₁`**.

---

## §2 · Why higher singularity cannot help

The atom weight is the **residue** — **local** data at the pole. Any functional that recovers it must **locate the pole** (find `m₁`), because:
- **integrating near the pole requires knowing where it is** (`∫_{ε<|x−m₁|<δ}`);
- the **singularity strength `p` only changes the power** `ε^{−(p−1)}`, not the need for the center or the limit;
- over a **fixed** region the value is `+∞` (presence), never a quantitative `f²`;
- **finite point values are integrals** matched by atomless measures (OPE-RHO).

So the higher singularity is **exactly TOOL-HARDY**: a **presence detector** whose **weight recovery = location (= the gap) + `ε→0` (= A3)**. It is an instance of the "`J→∞` (infinite local boundary data) + location = `H1`" exit — **not a counterexample**.

---

## §3 · Verdict

- **Falsification FAILED.** No higher-order singular kernel yields a **forced weight lower bound from finite data**.
- **Sharpening:** the singularity order `p` is **irrelevant** to the forcing question — every `∫|G|^p` (`p>1`) and `∫|G'|²` needs the **pole location** (the gap) for the weight. The residue is intrinsically local-at-the-pole.
- **Claim survives**: to FORCE (not merely detect presence of) the atom weight, one needs **a nonlinear constraint using the YM structure**, or **infinite data `J→∞` together with the pole location** = `H1`.
- **What this adds (real):** a clean **family statement** — the entire higher-order singular-kernel class is `presence-only`; the weight is always A3-given-location. This **closes** the "stronger singularity" escape hatch, narrowing the only-exit claim to exactly {YM structure, `J→∞`+location = `H1`}.

**Clay 0%; no QED.** Files: `research/unconditional/FALSIFY-HIGHER-SINGULAR.md`, `research/unconditional/w7/falsify_higher_singular.py`.
