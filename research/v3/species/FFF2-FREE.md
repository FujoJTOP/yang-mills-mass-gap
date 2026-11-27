# FFF2-FREE — 每个固定 `(L,a)` 上 `Δ_a>0` 是否免费

> 2026-11-27 · 靶：`H_a=−(1/a)log T_a` on `1^⊥` · 只判「每一步是否 `[T]`」· 零新 MC · 不改既有文件

## §0 判词

固定 `(L,a)` 的 `Δ_a>0` = `[T]`；四步中「有限维」= `[NG]`（`D02-DIRECT.md:13,41`：无穷维）；缺的是「紧」：有界核⇒HS⇒紧已 `[T]`。其余三条 `[T]`，PF 给 `λ₁<λ₀`。免费，内容全在一致下界（`[NG]`）；下一手见 §3。

## §1 四条的出处（逐字）

| # | 步骤 | 逐字出处 | 判 |
|---|---|---|---|
| 1 | 有限维 | `research/mass-spectrum/REPORT.md:150`：「时间切片上的 Hilbert 空间 `ℋ_Λ`（有限维）」；同用法 `research/v3/species/J06.md:72`（标 `[T]`）、`research/v3/species/IMP-ARROWS.md:139`（记 ✓4）；本文件判 `[NG]` | `[NG]` |
| 1′ | 反驳 | `research/v3/species/D02-DIRECT.md:13`：「`T` 在 `L=2` 上的维数 **∞**（可分）。载体 `L²(SU(2)²⁴)`」；`:41`：「算子的维数是 ∞……可分无穷维 Hilbert 空间」；一般 `L` 同（`:47–52` 表：`3L³` 个 link） | `[T]` |
| 1″ | 替代 | `D02-DIRECT.md:59`：「有界核 ⇒ Hilbert–Schmidt ⇒ **紧**」；`PROOF-GAUSS-INVARIANT-PF.md:7`：「`T̄` is Hilbert–Schmidt」 | `[T]` |
| 2 | 正性 | `research/gluing/REPORT.md:350`：「`T = e^{−aH}` 存在、`T ≥ 0`（Wilson 作用下为正算子，见 O–S 1978）」；`research/confinement/PROOF.md:171`：「存在正自伴转移矩阵 `T=e^{−aH} ≥ 0`」 | `[T]` |
| 2′ | 收缩 | `research/glue-b-blocking/REPORT.md:267`：「`0 ≤ T_a ≤ 1`，`T_a 1 = 1`」；`research/lit-v3/LIT-ERA-04.md:116`：「`0 ≤ T ≤ 1`」 | `[T]` |
| 3 | 自伴 | `research/confinement/PROOF.md:171`：「正自伴转移矩阵」；`PROOF-GAUSS-INVARIANT-PF.md:7`：「self-adjoint」 | `[T]` |
| 4 | 真空单重 | `PROOF-GAUSS-INVARIANT-PF.md:7`：「Jentzsch/Krein–Rutman: λ₀=1 simple with Ω≡1 strictly positive, and 0≤λ₁<λ₀」；`D02-DIRECT.md:60`：「`λ₁` 单重、`λ₁>0`、本征向量 `φ>0`、`|λ_j|<λ₁`」 | `[T]` |
| 4′ | 严格正核（前提） | `PROOF-GAUSS-INVARIANT-PF.md:5`：「`T(U,V)>0` for **all** U,V … Perron–Frobenius at *every* β>0, finite V」；`D02-DIRECT.md:58`：「核严格正」 | `[T]` |
| 5 | `λ₀=1` | `glue-b-blocking/REPORT.md:267` 的 `T_a 1=1`：真空能平移的**归一化**；判据只用 `λ₀/λ₁`，与归一化无关 | `[T]`（约定） |

检索记录：`grep` OK；式 1 `有限维.{0,40}(转移|Transfer)|(转移|transfer).{0,40}有限维|finite-dimensional.{0,60}(transfer|Hilbert)` 命中 25；式 2 `一致下界|均匀下界|uniform lower bound` 命中 311（显示 250，封顶）。断言 4D 格点 `T` 有限维的条目 ≥10；判其为假的 7 处：`D02-DIRECT.md:13,17,41,126`、`LOOP-RUNBOOK.md:7461,7465,7468`。无 TOOL-FAIL。

## §2 PF 给什么

PF 给 `λ₁<λ₀` 需要：**不可约（本处由严格正核给出）+ 紧（`λ₀` 孤立、次大本征值取到）**；只给「非负自伴收缩」不够（`T=id`，或两块各含一个真空 ⇒ `λ₁=λ₀`）。两条前提本处都 `[T]`（`PROOF-GAUSS:5,7`；`D02:58–60`）。但 PF 只给严格不等号、**不给数**（`D02:62`；`LOOP-RUNBOOK:7472`），常数随 `V` 形如 `e^{−cV}`（`PROOF-GAUSS:7(ii)`）。下一手：把「单重」换成带区间的 `λ₂/λ₁` 证书（截断+严格尾界），不再向 PF 要一致常数。

一致下界的具名出处：`D07-G.md:29`（逐 `(L,a)` 有隙 = 目标本身）、`J06.md:73`（有限体积隙只是上界）、`K0-BLIND.md:36`（困难在 `inf` 一致下界）、`docs/YANG-MILLS-CONDITIONAL-REDUCTION.md:36`（“trivially true by compactness”）。

## §3 判定

`Δ_a>0` = `[T]`；`[NG]` 只在有限维，换紧后全 `[T]`。一致下界 = `[NG]`，内容全在（`D07-G:29`、`J06:73`）。窗：固定 `(L,a)`、Wilson、`β>0`、周期。下一手：① `TransferData` 换紧算子版；② `λ₂/λ₁` 截断+严格尾界证书。

—— 只新建本文件；既有文件未改；未 commit；零新 MC；位置图由 lead 汇总。
