# Method 1 — intermediate-coupling cluster expansion (rigorous reach)

> CPB-attack probe. Verdict: **DEAD as a direct route.**

## Rigorous β ranges (4D SU(N))

| Result | Rigorous β range | Constant verified? |
|---|---|---|
| Osterwalder–Seiler 1978 | `β < β_c` (strong coupling), `M(β) ≥ c|log β|` | **No** — abstract only says `β_c>0`, group-dependent |
| SZZ23 (CMP 400, 805–851; arXiv:2204.12737) | SU(N): `β < 1/(16(d−1))`; d=4 → **`β < 1/48 ≈ 0.021`** (refined ≈0.027) | **Yes**, N-independent |
| IJGMP 2026 (Faizal–Shabir) | SU(N), `0<β<β*(N)`, fixed a, volume-uniform | **No** — `β*(N)` not in accessible abstract; "Part (1)", 0 citations |
| CNS25a (Cao–Nissim–Sheffield, arXiv:2509.04688) | 't Hooft / strong-coupling; area law | **No** — numeric β threshold not surfaced |

## Mapping

Project `β = 2N_c h^∨/λ`, SU(3): `β = 18/λ`. Crossover `λ≈7.3` → **`β≈2.47`**.

**Distance** (convention matters — SZZ β may be `1/g²` not `2N/g²`):
- if `1/g²`: crossover `0.41` vs `0.021–0.027` → **~15–20×**;
- if `2N/g²`: `2.47` vs `0.027` → **~90–120×**.

Either way, 1–2 orders of magnitude outside.

## Stronger structural obstruction

All results are for the **bare Wilson action at small bare β**; the crossover needs the
**blocked effective action** at `λ≈7.3` — not Wilson. The theorems do not apply verbatim
even at matching β. So the gap is **structural**, not just numeric. Known refinements
(1/48→0.027) are ~30%, not orders of magnitude.

## Verdict

**DEAD** (direct route). **Needs-X:** a fundamentally new expansion or a non-strong-
coupling mechanism — not in current literature.

## Honest boundary

Reach assessment, not a proof. `β_c` (OS78) and `β*(N)` (IJGMP26) not verified from
primary sources — re-check if larger. IJGMP2026 is 0-citation new work; recorded, not
endorsed. The numeric comparison uses bare coupling; the correspondence to the project's
RT effective coupling is unestablished.
