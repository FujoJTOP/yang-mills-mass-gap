# Yang–Mills Mass Gap — an evidence-first research log

> **Status: OPEN.** The 4D pure `SU(N)` Yang–Mills **mass gap** (Clay Millennium Problem) is **not solved here**. `Clay = 0%`. This repository is a **research log / ledger**, not a proof. Every claim is labelled and acid-tested; read it that way.

*(中文摘要见文末。)*

## What this is

An **AI-assisted, human-directed** campaign attacking the Clay problem *"Yang–Mills existence and mass gap"* for 4D pure `SU(N)` gauge theory. It is deliberately a **logbook**: ~450 report files recording what was tried, what died, and precisely *where* it died. The value is in the negative coordinates as much as the positive theorems.

## What is actually established (the positive ledger)

- **Conditional reduction.** The Clay gap (G2) is reduced to a single inequality **H1** (coarse-mode Poincaré / inverse susceptibility `G″(k→0) ≥ λ`) plus a finite list of named inputs — see the accompanying paper and `docs/PROOF-SKELETON.md`.
- **G1 (existence)** is, separately, essentially solved: a **machine-checked Lean 4 construction** (OpenAI, 2026) builds the 4D Wightman QFT via OS reconstruction — **existence only; it explicitly excludes the gap.**
- **Strong-coupling gap**: classical (Osterwalder–Seiler 1978; Shen–Zhu–Zhu, *CMP* 400 (2023), `β < 1/48`).
- **P3** — a finite-`L` strong-coupling **residue** (I1) theorem, `[T]` modulo the standard OS/Källén–Lehmann representation (corrected 2026-11-24; see `research/unconditional/PROOF-P3-UNCONDITIONAL.md`).
- **Structural theorems**: no-dissolution (kinematic), `W_geo ≥ 0` (Prékopa), gauge-invariant Perron–Frobenius, Feshbach–Bogoliubov two-block lemma, exact variance decomposition `G″ = S − D`, Ward/Stein identities.
- **Tools forged**: Fredholm-determinant identity (A3), multiplicative residue transport, gauge-invariant IR mass, gauge `W`-entropy, and a documented *death map* of ~11 species / 90+ routes.

## The wall (stated honestly)

The gap `= a non-positivity, nonabelian, no-small-parameter **lower bound** on the lightest `0⁺⁺` Källén–Lehmann **atom** at the strong-coupling crossover (`α_s ~ O(1)`) = **confinement**`.

- Higher ways to say the same: **H1** (`G″(k→0) ≥ λ`) is only the *inverse susceptibility* — it **upper**-bounds `m₁²`, strictly weaker than the gap; the gap is the edge atom `f₁² > 0` **with isolation** (`H*`).
- The **atom is upper-semicontinuous**; positivity / lower-bound (lsc) tools **cannot print it** — only exact identities *read* it, and those reduce to "locate the zero in the continuum" = the same wall.
- At the crossover there is **no small parameter**; every scale-separation tool (`(ℓ/L)^d`, UV→IR decoupling) vanishes there.

**The missing object** (never constructed) is a *nonlinear, non-positivity, IR-support-sensitive dynamical functional that takes different values on gapped vs gapless* — equivalently a stable sign-definite gauge-invariant quadratic minorant of the Yang–Mills Dirichlet form at `k→0`. It does not exist yet.

## Repository map

```
docs/
  YANG-MILLS-MASTER.md            # the master ledger
  PROOF-SKELETON.md               # Clay = G1 + G2, node-by-node status
  new-tools-registry.md           # every tool / no-go / structural theorem
  YANG-MILLS-CONDITIONAL-REDUCTION.pdf   # the paper, v1 (see below)
  YANG-MILLS-CONDITIONAL-REDUCTION-V2.pdf  # the paper, second edition (90 pp)
  YANG-MILLS-UV-FACE-V3.pdf       # focused paper, the UV face (17 pp)
  paper/                          # LaTeX source of v1
  paper-v2/                       # LaTeX source of the second edition
  paper-v3/                       # LaTeX source of the focused UV-face paper
research/unconditional/           # ~450 reports (the corpus) + w7/ witness scripts
lean/                             # Lean 4 formalization of the conditional reduction
CLAUDE.md                         # the project charter (rules of engagement)
```

## The paper

`docs/YANG-MILLS-CONDITIONAL-REDUCTION.pdf` — a **conditional-reduction** paper (deposited on Zenodo, DOI **10.5281/zenodo.24486316**; submitted). It is *not* a proof of Clay: it reduces the mass gap to H1 modulo explicit hypotheses and corrects the common conflation of the susceptibility `H1` with the spectral edge `H*`.

`docs/YANG-MILLS-CONDITIONAL-REDUCTION-V2.pdf` — the **second edition** (91 pages; LaTeX source in `docs/paper-v2/`). Same verdict (Clay = 0%), expanded with Part III: the wall drawn in **two independent coordinates** — measure/weak-* (the atom is a weak-* discontinuity; finite data are continuous) and β/RG/conformal (pure YM outside the conformal window; the window is real) — plus the machine-checked H1 sharpening, the Savvidy cleanup, the tool wave, the lattice diagnostics, a verification map, an errata record (including a corrected central-charge normalization, `a/c = 31/18`), and a new ultraviolet-face chapter classifying `(★)`, the soft mode, decompactification, the periodic constant-mode anchor, and the seven-route H1 attack with its carrier hazard. **Still not a proof.**

`docs/YANG-MILLS-UV-FACE-V3.pdf` — the **focused paper** (17 pages; LaTeX source in `docs/paper-v3/`). It decides the campaign's last ultraviolet inequality `(★)`: the strong form is **false** (exact SU(2) one-plaquette counterexample), the correct object is bounded below on the small-field set, the projection class is linear/divergence-free, the ultraviolet statement is renormalized, and the head is H1. It also records the soft-mode and decompactification waves, including the exact periodic constant-mode anchor `1/(2·Cin(2π)) = 0.20511530` and the saturating periodic `k=0` floor `c = 0.22196 ± 0.00011`. **Still not a proof of Clay, and no QED.**

## Honesty conventions (please respect them when reading)

- Every claim is labelled **真定理 / 条件定理 / no-go / 精确缺失步** (true theorem / conditional theorem / no-go / exact missing step).
- **No "QED" on anything unproven.** "Proven", "conditional", and "claimed" are kept distinct.
- Every tool passes a **free-field acid test**: it must *fail* on the massless free field (where `f₁² = 0`), else it is vacuous.
- Citations are checked; corrections are logged in the files (see the `⚠ 更正` banners).

## Reproducing the witnesses

Witness scripts are in `research/unconditional/w7/`. Run one at a time:

```
py -3.12 research/unconditional/w7/<script>.py
```

(Heavy lattice jobs are run *one at a time* — never in parallel — to avoid exhausting memory.)

## Status / license

Research notes shared as-is, for transparency and reproducibility. Not peer-reviewed except where noted (the paper). No warranty. `Clay = 0%`.

---

## 中文摘要

**状态:未解。** 本仓库是 4D 纯 `SU(N)` Yang–Mills **质量隙**(Clay 千禧问题)的**研究台账**,不是证明;`Clay = 0%`。每条结论都带标签(真定理 / 条件定理 / no-go / 精确缺失步),且过自由场酸测。

**真战利品**:Clay 的 G2 归约到单条 `H1`(粗模 Poincaré / 逆磁化率)+ 有限命名输入;G1(存在)已有机器验证的 Lean 构造(仅存在、明确排除隙);强耦隙为 OS'78 / SZZ'23;`P3` 为有限-`L` 强耦残数定理 `[T]` 模 OS/Källén–Lehmann(2026-11-24 已订正);另有 无溶解、`W_geo≥0`、规范不变 PF、`G″=S−D` 等结构定理,以及 A3 行列式恒等式、乘性残数搬运等新工具。

**墙**:隙 = 在交叉点 `α_s~O(1)` 给 `0⁺⁺` 原子一个**非正性、非阿贝尔、无小参数的下界** = 禁闭。`H1` 只是逆磁化率(只上界 `m₁²`),`H*`(边原子 + 隔离)才是隙;原子是 usc、正性工具(lsc)印不出;交叉点无小参数。**缺的那个物**(非线性、非正性、对 IR 支撑敏感、gapped/gapless 取不同值的动力学泛函)尚不存在。

**论文**:`docs/YANG-MILLS-CONDITIONAL-REDUCTION.pdf`(条件归约,Zenodo DOI 10.5281/zenodo.24486316;投稿中;**非 Clay 证明**)。第二版(91 页):`docs/YANG-MILLS-CONDITIONAL-REDUCTION-V2.pdf`,新增 Part III(墙的两个独立坐标:测度/弱-* 与 β/RG/共形)、H1 锐化的机器验证、Savvidy 清理、工具波、格点诊断、验证表与勘误,以及 UV 面专章(强形式 `[NG]`、好集定理、软模、解紧致化、周期常模锚点、H1 七路与载体危机)。聚焦版(17 页):`docs/YANG-MILLS-UV-FACE-V3.pdf`,分类 `(★)`:反例、好集定理、投影/归一化、塔头 = H1,并含周期常模精确锚点 `0.20511530`、`c=0.22196` 饱和、H1 四面等价引理与载体危机。**仍非证明**,Clay = 0%,无 QED。
