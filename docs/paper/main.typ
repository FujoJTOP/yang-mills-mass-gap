// =====================================================================
//  A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap
//  Typeset source. Compile:  typst compile docs/paper/main.typ out.pdf
// =====================================================================

#set document(
  title: "A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap",
  author: "Yang–Mills campaign dossier",
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, top: 2.2cm, bottom: 2.4cm),
  numbering: "1 / 1",
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 7.8pt, fill: luma(120))
      block(width: 100%, above: 0pt, below: 0pt)[
        #grid(
          columns: (1fr, auto),
          align: (left, right),
          [A conditional reduction of the 4D pure Yang--Mills mass gap],
          [#counter(page).display()],
        )
        #v(-6pt)
        #line(length: 100%, stroke: 0.4pt + luma(185))
      ]
    }
  },
)

#set text(font: ("New Computer Modern", "Libertinus Serif"), size: 9.7pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 0.95em, first-line-indent: (amount: 1.05em, all: true))
#show par: set block(spacing: 0.95em)

#set math.equation(numbering: "(1)", supplement: [Eq.])
#set math.equation(font: "New Computer Modern Math")

#set heading(numbering: "1.1")
#show heading.where(level: 1): it => block(above: 1.4em, below: 0.7em)[
  #set text(size: 12.5pt, weight: "bold")
  #context counter(heading).display("1")#h(0.6em)#it.body
]
#show heading.where(level: 2): it => block(above: 1.05em, below: 0.5em)[
  #set text(size: 10.4pt, weight: "bold")
  #context counter(heading).display("1.1")#h(0.55em)#it.body
]
#show heading.where(level: 3): it => block(above: 0.9em, below: 0.4em)[
  #set text(size: 9.8pt, weight: "bold", style: "italic")
  #context counter(heading).display("1.1.1")#h(0.5em)#it.body
]

#show figure.caption: set text(size: 8.3pt)
#show figure: set block(breakable: false)
#show table: set text(size: 8.4pt)

// ---------- label chips ----------
#let cT = rgb("#1d7a3a")
#let cC = rgb("#1f5fa8")
#let cM = rgb("#b45f06")
#let cN = rgb("#a61b1b")
#let cE = rgb("#6a3d9a")
#let labcol(l) = (T: cT, C: cC, M: cM, N: cN, E: cE).at(l)
#let L(l) = box(
  baseline: -0.18em,
  fill: labcol(l).lighten(90%),
  stroke: 0.4pt + labcol(l).lighten(25%),
  radius: 2pt,
  inset: (x: 3pt, y: 0.8pt),
  text(fill: labcol(l), weight: "bold", size: 6.8pt, l),
)

// ---------- theorem-like environments ----------
#let claim-counter = counter("claim")
#let claim(kind: "Theorem", name: none, lab: "T", body) = {
  claim-counter.step()
  block(
    width: 100%,
    breakable: false,
    inset: (left: 9pt, right: 6pt, y: 5.5pt),
    radius: 2pt,
    stroke: (left: 2.4pt + labcol(lab)),
    fill: labcol(lab).lighten(94%),
    above: 0.9em,
    below: 0.9em,
  )[
    #set text(size: 9.35pt)
    #set par(first-line-indent: 0pt, justify: true, spacing: 0.75em)
    #text(weight: "bold")[#kind #context claim-counter.display()]
    #if name != none [ #text(weight: "bold")[(#name)]]
    #h(0.35em)
    #body
  ]
}
#let thm(name: none, lab: "T", body) = claim(kind: "Theorem", name: name, lab: lab, body)
#let lem(name: none, body) = claim(kind: "Lemma", name: name, lab: "T", body)
#let cor(name: none, body) = claim(kind: "Corollary", name: name, lab: "T", body)

#let proof(body) = block(
  width: 100%,
  inset: (left: 9pt, y: 4pt),
  stroke: (left: 1pt + luma(160)),
  above: 0.25em,
  below: 0.85em,
)[
  #set text(size: 9.35pt)
  #set par(first-line-indent: 0pt, spacing: 0.75em)
  #emph[Proof.] #body #h(1fr) $square$
]

#let src(body) = block(width: 100%, above: 0.15em, below: 0.55em)[
  #set text(size: 8.6pt, fill: luma(70))
  *Source.* #body
]

#let quote-block(body) = block(
  width: 100%,
  inset: (left: 13pt, right: 10pt, y: 6pt),
  stroke: (left: 2pt + luma(150)),
  fill: luma(97%),
  above: 0.7em,
  below: 0.8em,
)[
  #set par(first-line-indent: 0pt, justify: true, spacing: 0.7em)
  #set text(size: 9.2pt)
  #body
]

// ---------- styled table helper ----------
#let tbl(caption, cols, aligns, header-row, ..rows) = figure(
  table(
    columns: cols,
    align: aligns,
    inset: (x: 4.5pt, y: 3.2pt),
    stroke: (x, y) => (
      top: if y == 0 { 0.9pt } else { none },
      bottom: if y == 0 { 0.5pt } else if y == rows.pos().len() { 0.9pt } else { 0.2pt + luma(212) },
      left: none,
      right: none,
      rest: none,
    ),
    fill: (x, y) => if y == 0 { luma(92%) } else if calc.odd(y) { luma(98.5%) } else { white },
    table.header(..header-row),
    ..rows,
  ),
  caption: caption,
  kind: table,
  supplement: [Table],
)

#let note(body) = block(width: 100%, above: 0.55em, below: 0.75em)[
  #set text(size: 8.8pt)
  #set par(first-line-indent: 0pt)
  #body
]

#import "figs.typ": *

// =====================================================================
//  FRONT MATTER
// =====================================================================

#align(center)[
  #v(0.9cm)
  #text(size: 18.5pt, weight: "bold", "A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap")
  #v(0.4em)
  #text(size: 11.5pt, style: "italic", "and a record of what was actually proven")
  #v(1.0em)
  #line(length: 52%, stroke: 0.6pt + luma(120))
  #v(0.8em)
  #text(size: 9.3pt, fill: luma(95))[Yang–Mills campaign · Lean 4 formalization + lattice numerics]
  #v(0.15em)
  #text(size: 9.3pt, fill: luma(95))[23 November 2026]
  #v(0.6em)
]

#block(
  width: 100%,
  fill: luma(96.5%),
  stroke: 0.5pt + luma(180),
  inset: 10pt,
  radius: 3pt,
  above: 1.1em,
  below: 0.9em,
)[
  #set par(first-line-indent: 0pt, justify: true, spacing: 0.7em)
  #set text(size: 8.8pt)
  *Status of this document.* This is a *conditional reduction and an audit*, not a proof of
  the Clay problem. The Clay mass gap for 4D pure SU(N) Yang–Mills remains open. The only
  defensible headline is: _the problem has been rewritten onto one precisely named,
  acid-tested object, and every honest route to it terminates there._

  #v(0.35em)
  *Label legend used throughout.* Each substantive claim carries one of:
  #L("T") #emph[theorem] — proved (any imported standard result is named);
  #L("C") #emph[conditional theorem] — a proved implication from explicitly named hypotheses;
  #L("N") #emph[no-go] — a theorem that a whole class of approaches cannot work;
  #L("M") #emph[exact missing step] — a precisely stated open inequality, with no slack hidden;
  #L("E") #emph[numerical evidence] — measured, not proved.

  #v(0.35em)
  *Provenance rule.* Every theorem below names the file that contains its proof. No claim is
  asserted without a source. Citations in the reference list appear verbatim in the corpus;
  where the corpus itself flags a bibliographic uncertainty, this is reproduced rather than
  silently repaired.
]

#block(width: 100%, inset: (x: 6pt), above: 0.9em, below: 0.6em)[
  #set par(first-line-indent: 0pt, justify: true, spacing: 0.78em)
  #set text(size: 9.15pt)
  #align(center)[#text(size: 10pt, weight: "bold")[Abstract]]
  #v(0.4em)

  We report the outcome of a sustained campaign — Lean 4 formalization plus lattice numerics —
  aimed at the 4D pure SU(N) Yang–Mills mass gap. *We do not prove the mass gap.* We prove,
  instead, a reduction and a ledger.

  The reduction has three parts. *(i) A rewriting of the problem.* The gap is equivalent to the
  statement that the Källén–Lehmann spectral measure $ν$ of the gauge-invariant field-strength
  operator $tr F^2$ has its lightest atom $m_1^2$ at positive mass with positive weight
  $f_1^2 = ν({m_1^2}) = abs(⟨0|tr F^2|G⟩)^2$, and that this atom is isolated. *(ii) An exact
  spectral squeeze.* For any such $ν$ with susceptibility $χ = op("Var")(a) = ∫ s^(-1) dif ν$,
  one has $f_1^2 slash χ ≤ inf op("supp") ν ≤ M_0 slash χ$. It follows that a lower bound on
  the mass requires *both* a residue floor $f_1^2 ≥ r_0$ and a susceptibility ceiling
  $χ ≤ 1\/λ$; the often-quoted condition H1 ($χ ≤ 1\/λ$) is *necessary but not sufficient*,
  because it is blind to a thin atom. The correct sufficient condition is the pair
  $(f_1^2 ≥ r_0) and (χ ≤ 1\/λ)$, which gives $m_1^2 ≥ r_0 λ$. *(iii) A convergence.* Every
  route we could construct — and we constructed on the order of ninety — ends at *one species*:
  a nonperturbative, non-circular, sign-controlling *upper bound on a conditioned variance*,
  appearing once in the ultraviolet (the L3 node of existence theory) and once in the infrared
  (H1). The same inequality, drawn at two scales, is the remaining gap between the campaign and
  Clay.

  The ledger has a positive half and a negative half. On the positive side we record the genuine
  theorems: an unconditional strong-coupling *residue* theorem ($f_1^2 ≥ r_0(β) > 0$ for
  $β < 1\/48$, volume- and $L$-uniform — a statement about the atom weight I1, *not* its
  isolation I2); a kinematic *no-dissolution* theorem for the lightest $0^(+ +)$; the
  *non-negativity of the Haar geometric curvature* $W_("geo") ≥ 0$ by Prékopa; a *gauge-invariant
  Perron–Frobenius* theorem on the physical Hilbert space; the *strong-coupling coarse Poincaré*
  lemma (Lemma SC); the abstract *Feshbach–Bogoliubov two-block lemma*; and the *exact variance
  decomposition* $G'' = S − D$. Six new tools pass the free-field acid test (they must fail on
  the massless free field, where $f_1^2 = 0$ and $G'' = 0$). On the negative side we record the
  no-gos: the Schur-convex master obstruction, Goldstone uniqueness at zero mass only, the
  periodic, nested-cone, rigidity, Feshbach-angle and Krein obstructions.

  Numerically, the residue is seen directly where it is resolvable: $f_1^2 > 0$ at $β = 2.2$ at
  seventeen sigma on SU(2) ($L = 6$); at $β = 1.0$ it dissolves into the cutoff, exactly as the
  analytic strong-coupling picture predicts, not as a failure. A $k → 0$ curvature probe on
  $L = 8$ finds the coarse-mode $G''$ positive, $W_("geo")$-carried, and volume-stable — but it
  is measured on a *gauge-dependent gluon mode*, not the $0^(+ +)$ object, and its
  rigorous-style one-sided interval still includes zero.

  The verdict is sober: *Clay is 0% — open.* The wall did not fall. It became transparent,
  ruled, and instrumented. The tool that would cross it — a nonperturbative machine transporting
  the strong-coupling gap and residue across the Savvidy crossing and taking $a → 0$ — *does not
  exist*, and naming it precisely is the campaign's real result.
]

#pagebreak()
#block(above: 0pt)[
  #set text(size: 10pt)
  #outline(title: [Contents], indent: auto, depth: 2)
]
#pagebreak()

#include "body1.typ"
#include "body2.typ"
#include "body3.typ"
#include "body4.typ"
