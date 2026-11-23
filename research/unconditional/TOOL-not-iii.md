# TOOL-not-iii — can the theorem "pure-gauge 4D CFT ⟹ free" be built?

> 2026-11-23 · **Verdict: no existing tool closes ¬(iii).** The four angles are not
> independent escape routes — they are four faces of one unproven statement, and I can
> *prove* the reduction (§0). Among actual *methods*, (d) as the published Wightman-positivity
> no-go of Nikolov–Stanev–Todorov (hep-th/0305200) is closest; completing it needs L1, L3, and
> a no-emergent-SUSY lemma — a new field, not a missing computation. No fabricated theorem.

## 0. The target provably collapses to (b)

**Lemma (¬(iii) ≡ no-interacting-FP).** A pure-gauge IR conformal point is *free* iff it is
Gaussian (`g=0`), and *non-free* iff `g*≠0` = an interacting IR fixed point. In 4D "IR fixed
point" = "IR CFT" (scale ⟹ conformal, Komargodski–Schwimmer 1107.3987). Hence ¬(iii) ⟺
"pure YM has no interacting IR fixed point" = **(b)**. ∎

*Consequence.* The anomaly, `Z_N^{(1)}`, and content clauses add no leverage — they *define*
"pure-gauge CFT." So (a)–(d) are routes to (b), not shortcuts: the real shape of the wall.

## (a) Generalized symmetry / anomaly — dead at θ=0

The one theorem-let: a 4D CFT with a **continuous** 1-form symmetry is free (Hofman–Iqbal
1802.09512; cf. GKSW 1412.5148). Pure YM carries only **discrete** `Z_N^{(1)}` — no conserved
2-form current — so the theorem does not transfer (`13-anomaly-class`). At θ=0 the
purely-electric `Z_N^{(1)}` has **no 't Hooft anomaly**, so anomaly matching has zero force
(`ANOMALY`; GKKS 1703.00501): anomalies constrain phases, they do not make gaps (`B2` §(b)).
Verdict: *(a) ⊬ ¬(iii) at θ=0.*

## (b) β-function — clean statement, no technique

Pure glue: `b₀=11N/3>0` (Gross–Wilczek/Politzer 1973), `b₁=34N²/3>0` (Caswell 1974), so
`g*²=−b₀/b₁<0`: no perturbative IR zero. **No nonperturbative "IR-free" theorem covers 4D
pure glue**; the only 4D monotone is the a-theorem (KS), and `a_IR<a_UV` is *consistent* with
an interacting CFT. Reformulated, (b) is a **one-dimensional** problem — prove the scalar
`β(g)` has no positive zero — the cleanest *statement* of the wall and the farthest from any
*technique*: it needs `β(g)` at all couplings. Verdict: *(b) is the statement; unprovable today.*

## (c)+(d) The analytic-positivity route — closest

The numeric bootstrap is dead ((c)): SDPB yields pointwise **islands**, never global
exclusion, and pure YM has no conformal crossing at all (`III-UNIQUENESS`,
`GLUEBALL-BOOTSTRAP`). Its **analytic** form, though, is a genuine theorem-let — and it forces
the free dimension/coefficient, so it is simultaneously the (d) route:

- **(d1)** conserved-current freeness (HI / Maldacena–Zhiboedov 1105.0473) — blocked by
  discreteness, as in (a).
- **(d2) ESI / Wightman positivity** (Nikolov–Stanev–Todorov hep-th/0305200; hep-th/0504146):
  for a dimension-4 scalar `L=tr F²`, positivity **plus** OPE condition (4.1) gives
  `−3a′≤b≤a′/3`, so `a′=0 ⟹ b=0 ⟹` free Maxwell. **L2 (`a₀=0`) is closed structurally for
  pure YM** — the lowest gauge-invariant scalar is dim 4 (`L2.md`, gated by L3). This is the
  most advanced *published* partial result toward ¬(iii), and it is analytic, not numeric.

**Exactly what completes (d) = ¬(iii):**
- **L1:** extend positivity to 6-/2n-point functions to force `a′=0` — the **non-scalar
  bilocal 6-point** structures are uncomputed (hep-th/0504146 defers them).
- **L3 (decisive):** prove the single-`L` reduction (4.1) for non-abelian SU(N). The abelian
  Furry vanishing fails (`⟨FFF⟩≠0`); the source itself admits "not incompatible with N=4 SYM."
  L3 is *likely false as stated* for non-abelian (`L3.md`).
- **No emergent SUSY:** ESI's surviving branch may *be* N=4-type, so ¬(iii) also needs "pure
  YM cannot flow to an emergent-superconformal IR" (`06-anomaly` extra gap).

## Honest boundary

Close L1 + L3 + no-emergent-SUSY and ¬(iii) follows; with non-triviality (free ⟹ `Λ=0` ⟹
trivial) and the IN-b-type input, the gap follows. But L3 is likely false as stated, L1 is a
hard 4D partial-wave computation, and no-emergent-SUSY is an independent classification
problem. **The completing object is a nonperturbative classification of pure-gauge 4D CFTs —
a field that does not yet exist.** (a),(c) are walls; (b) is the statement; (d)/ESI is the
closest, and it is not a matter of more effort on a known lemma.

*Refs:* GKSW 1412.5148; Komargodski–Schwimmer 1107.3987; Gross–Wilczek/Politzer PRL 30
(1973); Caswell PRL 33 (1974) 244; Hofman–Iqbal 1802.09512; Maldacena–Zhiboedov 1105.0473;
Nikolov–Stanev–Todorov hep-th/0305200 & hep-th/0504146; GKKS 1703.00501. Internal:
III-UNIQUENESS, ANOMALY, B2, IN-B, GLUEBALL-BOOTSTRAP, 06-anomaly, 13-anomaly-class (L1/L2/L3).
