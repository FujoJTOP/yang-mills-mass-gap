# Adapting N=4 SUSY line-defect blocks to non-SUSY

> CPB-attack probe (sub-task 10). Verdict: the earlier "4D non-SUSY blocks unavailable" may be
> **too strong**; the 1d blocks are standard; the real gap is the pure-gauge `Z_N` setup.

## SUSY blocks

The 1/2-BPS line in `N=4` SYM gives `OSp(4*|4)` superconformal blocks; `N=2` gives `osp(4*|2)`.
The superconformal algebra binds defect operators into supermultiplets, so the SUSY block is a
multiplet sum.

## Non-SUSY

This agent reports the defect group `SO(2,1)×SO(2)` and the blocks as the **1d conformal
blocks** in explicit hypergeometric form
`Block1D_Δ(χ) = χ^Δ ₂F₁(Δ−Δ₁+Δ₂, Δ+Δ₃−Δ₄; 2Δ; χ)`, labeled by `(Δ, s)` — "standard and
available", not a missing ingredient.

## Deformation?

No: SUSY-breaking destroys the multiplet structure, so the SUSY block is **not** a smooth
limit. Not needed — the non-SUSY block is simpler and independently known (1d hypergeometric).

## Actually missing

Not the blocks but the **setup**: a non-SUSY line-defect bootstrap with a `Z_N` 1-form
symmetry and no flavour symmetry. Non-SUSY line-defect bootstraps exist for `O(2)/O(N)`
magnetic/monodromy defects; the **pure-gauge `Z_N` analogue is the unbuilt piece**.

## Conflict to resolve

`KINEMATICS.md`/`OPE.md` say `SO(2,1)×SO(3)`, 3 cross-ratios, multi-variable blocks **missing**;
this agent says `SO(2,1)×SO(2)`, 1d hypergeometric, **available**. The two are describing
different objects (bulk-2pt-with-defect vs defect-4pt) and/or an outdated literature survey.
**Resolve via arXiv:2501.06900.**

## Honest boundary

The "blocks unavailable" claim is too strong; the real gap is the pure-gauge `Z_N` setup plus
the full 1-form-constraint layer.

Sources: [Universal Constraints for Conformal Line Defects (arXiv:2501.06900)](https://ar5iv.labs.arxiv.org/html/2501.06900);
[Bootstrapping line defects with O(2) symmetry (arXiv:2208.11715)](https://indico.desy.de/event/34520/contributions/129820/attachments/77407/100089/vanVliet.pdf);
[Bootstrapping the half-BPS line defect (arXiv:1806.01862)](https://bib-pubdb1.desy.de/record/405629/files/1806.01862.pdf).
