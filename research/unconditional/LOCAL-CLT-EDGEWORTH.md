# LOCAL-CLT-EDGEWORTH — species: local CLT / Edgeworth on the coarse marginal `p(a)`. The leading Gaussian term IS the whole corridor; the Edgeworth correction is a 2% decoration and is *provably not the gap* — the **massive free field is exactly Gaussian (`D=0`) yet gapped (`G″=m²>0`)**, so "gap = non-Gaussianity" is refuted by the species' own acid test.

> **判词 (verdict).** Expand the `k→0` marginal as a local CLT with Edgeworth refinement, `p_Z(z)=φ(z)[1+(κ₃/6σ³)H₃(z)+(κ₄/24σ⁴)H₄(z)+…]`, `Z=(a−ā)/σ`, `σ²=Var(a)=1/S`. The exact identity `G″=S−D` (SCBL, INVENT-VARIANCE-INEQ §1) is then literally "Gaussian variance minus non-Gaussianity": the leading term gives `S=1/κ₂=χ`, the corrections `{κ₃/σ³,κ₄/σ⁴,…}` ARE the defect `D`, and `D/S=1−G″Var(a)`. **The acid test passes in the strongest form of the family** — massless free is *exactly* Gaussian (`κ_{n≥3}=0`, `D=0`), the tool's signal is *zero* there (free is the minimizer, not — as in E2/QFB/FOURIER — the maximizer), and `G″=S−D=0−0=0`. **But the same exactness refutes the premise:** the *massive* free field is also exactly Gaussian (`D=0`) yet `σ²=1/m²<∞ ⇒ S=m²>0 ⇒ G″=m²>0`. So `D=0` does **not** give `G″=0`; the gap is `S−D`, `98%` Gaussian variance at β=2.5 (`D/S≈2%`). The correction is the 2%, *not* the gap. **#1,#4 pass (and #4 is the cleanest in the family); #2,#3,#5 fail. Not 5/5.**

## 1 Forge and exact dictionary
`a=⟨v,X⟩`, `|v|₂=1`, `μ∝e^{−V}dX`, `σ²=Var(a)=1/S=χ=G_T(0)` (D1). Cumulants `K(t)=log E e^{t(a−ā)}=½σ²t²+κ₃t³/6+κ₄t⁴/24+…`; the local CLT inverts `e^{K(it)}`. The leading term reads `σ²=1/κ₂` — **the Edgeworth "Gaussian term" *is* `S`**, the very susceptibility `S≥λ` of H1. The corrections are the non-Gaussianity; `G″=S−D` reads `D/S=1−G″Var(a)`, and `D` and `{κ₃/σ³,κ₄/σ⁴}` vanish together iff `p` Gaussian — one defect, two costumes (cf. E2's Fisher/dynamical/MMSE unification).

## 2 Acid test — passed, and it refutes the premise
Massless free: `p` exactly Gaussian ⇒ `κ_{n≥3}=0`, Edgeworth `≡0`, `D=0`; and `σ²=c/k²→∞` ⇒ `S→0`; hence `G″=S−D=0`. ✓ **Genuine pass, strongest of the family:** the free field *minimizes* the signal (correction 0), unlike E2/QFB/FOURIER which peaked on free. Massive free: `p` exactly Gaussian (`D=0`) but `σ²=1/m²` ⇒ `S=m²>0` ⇒ `G″=m²>0` — **gapped with zero non-Gaussianity** (INVENT-VARIANCE-INEQ §4 verbatim). So exact Gaussianity does *not* imply gaplessness; the check `D=0 ⇒ G″=0` holds only because *also* `S=0`, and conflating the two is the error. **The gap is the leading Gaussian variance, not the Edgeworth correction.**

## 3 Sign — the correction is the wrong side, the defect is small and positive
To bound `G″=S−D≥λ` needs `S≥λ` (lower, H1) *and* `D≤S−λ` (upper on defect). Edgeworth controls `|D|≤C` (the upper side, right *form*, but for the 2% term). Measured `D/S≈+2%>0` (Savvidy indefinite fiber), no universal sign (C8); the leading correction is a signed skewness (no symmetry forces `κ₃=0` for nonabelian SU(N)). Thus `G″≈0.98·S`: `λ` must come from `S`, untouched by Edgeworth. FAIL.

## 4 Non-circularity — the uniform remainder is clustering
A V-uniform Edgeworth/Berry–Esseen remainder `sup_a|p_Z−φ|≤C/√V` with bounded `C` requires mixing `ξ<∞` (Dedecker–Doukhan; `C∝` the mixing coefficient). So a uniform bound `|D/S|≤C<1` ⟺ `ξ<∞` ⟺ clustering ⟺ gapped = H1. And `S≥λ` is H1 outright. Two H1s.

## 5 Five-point
1 **PASS** (exact cumulant/SCBL, all β). 2 **FAIL** (Edgeworth bounds defect from above; `S≥λ`=H1; `D>0`). 3 **FAIL** (`|D/S|≤C` uniform ⟺ `ξ<∞` ⟺ clustering; `S≥λ`=H1). 4 **PASS** (free massless exactly Gaussian, tool silent there — cleanest free pass of the family; non-SUSY). 5 **FAIL** (no λ; uniform remainder = mixing = gap).

## Exact failing estimate
```
G″=S−D,  S=1/κ₂=Var(a)^{-1}=χ,   D/S=1−G″Var(a)=1−G″/S
Edgeworth: p_Z=φ[1+(κ₃/6σ³)H₃+(κ₄/24σ⁴)H₄+…] ;  leading term σ²=1/S
massless free: exact Gaussian ⇒ κ_{n≥3}=D=0  AND  σ²=c/k²→∞ ⇒ S=0  ⇒ G″=0     [PASS]
massive   free: exact Gaussian ⇒ D=0  BUT  σ²=1/m² ⇒ S=m²>0  ⇒ G″=m²>0       [REFUTES "gap=D"]
β=2.5:  D/S≈2%  ⇒  G″≈0.98·S   (gap 98% Gaussian variance, 2% Edgeworth correction)
sign :  need S≥λ (lower)=H1  +  uniform |D/S|≤C<1 ⟺ ξ<∞ ⟺ clustering = H1
```

**Sharpest line.** The local CLT is the family's sharpest *diagnostic*: its leading term reads off `S=1/Var(a)=χ` exactly, and — uniquely among the species — its Edgeworth correction is *minimized* on the free massless field (exactly Gaussian), so `free ⇒ G″=0` passes cleanly. But the same exactness is the refutation: the massive free field is exactly Gaussian (`D=0`) and gapped (`G″=m²`), so **"gap = non-Gaussianity" is false**, and the β=2.5 defect `D/S≈2%` is a decoration on a gap that is 98% Gaussian variance `S`. That yields a bankable **one-step no-go**: *no Edgeworth-correction / skewness / kurtosis mechanism can source the gap* — it must come from `S≥λ`, i.e. `Var(a_{k→0})≤1/λ`, i.e. H1. The uniform remainder that could force `|D/S|≤C` is mixing = clustering = the gap itself. Missing step unchanged: an upper bound `Var(a_{k→0})≤1/λ` not derived from `Var(a)`.

*Refs.* INVENT-VARIANCE-INEQ (`G″=S−D`; §4 massive-free `D=0,G″>0`); E2-INFORMATION (`G″Var=1` iff Gaussian; `J`/`D`/MMSE one defect); D1-BLOCK-MARTINGALE (`Var(a)=χ=1/Γ_T(0)`); CONTINUUM-VARIANCE-DEFECT (`D/S≈2%`, `D<S`=H1); FOURIER-UNIFORMITY (char-fn neighbour; Stam sign-reversal); k0-LARGE-DEVIATION (tails, distinct regime); C8-SPECIES-NOGO (no universal sign); ABC-NEW-INEQ (Ward `E[af]=1`); PROVE-Gpp-positive (β=2.5 `G″≈1/Var`, `D/S` measured). External: Edgeworth 1905; Cramér 1928; Gnedenko–Kolmogorov, *Limit Distributions for Sums of Independent Random Variables* (1954); Petrov, *Sums of Independent Random Variables* (1975); Bhattacharya–Ranga Rao, *Normal Approximation and Asymptotic Expansions* (1976/2010); Esseen 1945 / Berry 1941; Rosenblatt, *Markov Processes* (1971); Dedecker et al., *Weak Dependence* (2007).
