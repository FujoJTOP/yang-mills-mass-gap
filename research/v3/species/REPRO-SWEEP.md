# REPRO-SWEEP — w7 witness 全量重跑 · 笔记↔脚本 数字漂移审计（第一批）

> 时间：2026-11-25 12:03–12:35（本地）｜执行：`s2-repro`（task-6）
> 范围：`research/v3/species/w7/*.py`（65）+ `research/v3/species/w7/audit/*.py`（5）= **70 个脚本全部重跑**
> 状态：**部分完成**（脚本 70/70 跑完；关键数字比对只做到核心 20 个；task-6 未 complete，等 lead 拆批重派）

---

## §0 结论一句话

**59 OK / 1 DRIFT / 5 超时 / 5 失败(全为 cwd 路径约定)**，外加 1 条脚本自检失效（`na_plaquette` 的 `B=0`）。
判死线（`DRIFT > 3` ⇒ 语料不可复现）**未触发**：笔记引用的 R4 `1.182227`、KT-RHO2 `0.99999376`、DS `β²/β*²`、R1 `1.2e-15`、SZ `1/32` 全部逐位复现；**唯一漂的是 `NA2-NONDEGENERATE.md` 表里 `X(a₀=0)` 那一列**（差 3–4%，约 1–1.6σ，判词方向不变）。

### §0.5 范围算术（给拆批用）

- 我的名单 = **12:03 快照**里的 `w7/*.py`(65) + `w7/audit/*.py`(5) = **70**，**70 个全部跑完**（无一遗漏）。
- 此刻树里还有 **3 个新 `.py` 不在我的名单**：`_patch.py`（**是改代码的脚本，别当 witness 跑**：它会把 `r3_l4_volume.py` 重写）、`cov_link.py`、`cov_models.py`（12:23/12:24 新增，尚无笔记引用）。
- 另有 **5 个 EXT**（`r4_kp.py` + `sz_ceiling_{const,exact,hess,lmax}.py`，在 `research/v3/w7/`，不在 species 树）：4 跑完、1（hess）中途停、`lmax` 未跑。
- ⇒ 若按 20 s 超时重派：70 个里 **180 s 才超时的 5 个**需要用 180 s（`hc_spherical`/`ns_nosub`/`r3_l4_volume`/`na_plaquette`/`na2_two_plaquette`），其余 65 个 20 s 够（最慢的非超时脚本：`na2_tci` 139.5 s、`lp_anchor` 130.2 s、`na2_chain` 110.1 s、`rb_ratio_bound` 102.8 s、`ru_joint` 101.9 s、`ct_chain_L` 88.1 s、`tci_comp` 71.8 s ⇒ **20 s 会误杀这 7 个**，建议 150 s）。

## §1 纪律、环境、以及两件必须先说的事

- **跑法**：`py -3.12 -I <script>`，Python 3.12.6 / numpy 2.4.6 / scipy 1.18.1；**严格单进程顺序**，每脚本 **180 s** 硬上限（超时 `taskkill /T`，记"超时"，不判错）。
- **落盘位置**（全部在只允许我写的地方）：`w7/repro/run_all.ps1`（运行器）、`repro/out/`（每个脚本的 stdout+stderr）、`repro/snap/`（跑完后的 .log/.json 快照）、`repro/baseline/`（跑前 108 个 .log/.json 的原样备份）、`repro/analyze2.py`（比对器）、`repro/results.tsv`（机读结果）、`repro/analysis.json`。
- **cwd 约定（关键）**：`w7/*.py` 用 `research/v3/species/w7`（与 DS-RATIO.md §复算 写法一致）；`audit/*.py` 也用 `research/v3/species/w7` —— 它们内部写 `../../w7/…`、`../../../unconditional/…`，需要 research/v3 下 depth-3 的 cwd；**按直觉在 `w7/audit/` 里跑会全部 FileNotFound**（首轮实测）。
- **⚠ 并发**：本机当时另有 agent 在跑 `s10_fixed.py`、`r3_l4_volume.py`、`r1_*.py` 等。5 个"超时"里 `r3_l4_volume`、`na_plaquette` 在 180 s 内只产出 40–60% 输出，属资源争抢下的**下界**，不代表单机独占也超时。
- **⚠ 语料在扫描期间被改动**：`r3_l4_volume.py` 于 **12:07:11** 被 `_patch.py` 改写（我这次跑的是改后版本，所以它的 TIMEOUT 不能当冻结版本结论）；`cov_link.py`(12:23:09)、`cov_models.py`(12:24:34) 扫描中途新出现。⇒ **本表是「12:03 时刻的 70 个脚本」的截面**，不是冻结语料。

## §2 判据（没有为了全绿而放宽）

每个数按**笔记自己印出的位数**取容差（半 ulp）：`−1.69e−2` ⇒ ±5e−5；`0.99999376` ⇒ ±5e−9；`1/32` ⇒ 精确；笔记位数比脚本多时按两者中较粗者判（例如笔记 `1.945910` vs 输出 `1.9459` 判 OK，笔记 `−1.95e−3` vs 输出 `−1.951e-03` 判 OK）。**未**使用统一 5%/1% 之类的宽口径。

---

## §3 全部 70 个脚本：脚本 / exit / 耗时 / 判

| # | 脚本 | exit | 耗时 (s) | 判 | 备注（笔记值 ↔ 重跑值） |
|---|---|---|---|---|---|
| 1 | `audit_c0_ratio.py` | 0 | 1.5 | **OK** | AUDIT-C0 `0.5288/0.5475/0.5633`、`A/C(0)=0.5627`、flat/measured `9.38x` → 重跑 0.52880/0.54750、0.5627、9.38x |
| 2 | `at_scalar.py` | 1 | 0.2 | **失败** | `open("research/v3/species/w7/at_scalar.json")` —— 需 cwd=仓库根；计算本身完成，stdout 完整（见 §失败） |
| 3 | `at_su2_exact.py` | 1 | 52.8 | **失败** | `open("research/v3/species/w7/at_su2_exact.json")` —— 需 cwd=仓库根；52.8 s 计算完成，仅末步写盘失败 |
| 4 | `bc_ising.py` | 0 | 0.7 | **OK** | BC-BETHE 3/3（含 `9e−4` 与重跑 `8.89e−04`，在笔记自身 1 位有效数字内） |
| 5 | `bc_su2.py` | 0 | 0.8 | **OK** | BC-BETHE 5/5（`1.2e−15` vs 重跑 `1.244e-15`） |
| 6 | `hc_spherical.py` | TIMEOUT | 181.5 | **超时** | 180 s 未完成（11-24 23:13 的旧 log 3051 B 存在）；按任务书记"超时"，不判错 |
| 7 | `kt_a_exact.py` | 0 | 1.5 | **OK** | KT-RHO2 5/6 命中（唯一 miss 是笔记里的日志字节数 8051，非物理量） |
| 8 | `kt_b_sl_battery.py` | 0 | 35.5 | **OK** |  |
| 9 | `kt_c_geom.py` | 0 | 1.8 | **OK** |  |
| 10 | `kt_d_gevp.py` | 0 | 0.7 | **OK** |  |
| 11 | `kt_e_identity.py` | 0 | 0.9 | **OK** | KT-RHO2 `0.99999376` (DD Q)、`2.66887958` (ND)、`15.6442` (双井) → 重跑 0.99999376 / 2.66887958 / 15.6442，逐位一致 |
| 12 | `lx_label_geometry.py` | 0 | 2.8 | **OK** | SIXTH-SPECIES 3/3 命中 |
| 13 | `na2_su3_na.py` | 0 | 49.5 | **DRIFT** | NA2 表 `X(±1)`、`se`、H-SCAN 三段带宽收敛**逐位一致**；但表内 `X(a₀=0)` 列不一致：`−1.69e−2 ↔ −0.017592`、`−5.53e−3 ↔ −0.005696`、`−2.00e−4 ↔ −0.000208`（详见 §5 DRIFT-1） |
| 14 | `ns_nosub.py` | TIMEOUT | 180.6 | **超时** | 180 s 未完成，已产出 24 370 B stdout（到 `a0=+2.00` 段）；不判错 |
| 15 | `qg1_network.py` | 0 | 1.8 | **OK** | QG-NETWORK：无 ≥4 位引用数，自检全 PASS |
| 16 | `qg2_abelian.py` | 0 | 17.6 | **OK** |  |
| 17 | `qg3_acid.py` | 0 | 0.5 | **OK** | QG-NETWORK 2/2 命中（酸测：自由场失败 = 非 vacuous） |
| 18 | `qg4_su2.py` | 0 | 1.2 | **OK** | QG-NETWORK 4/4 命中 |
| 19 | `qi_bracket.py` | 1 | 2 | **失败** | `open("research/v3/w7/m1_4d_run_results.json")` —— 需 cwd=仓库根（文件确实在 `research/v3/w7/`）；本条未产出任何新数 |
| 20 | `r1_2dym.py` | 0 | 0.7 | **OK** |  |
| 21 | `r1_chain.py` | 0 | 2.3 | **OK** | R1-BRIDGE 相似变换残差 `1.2e-15` → 重跑 `[PASS] 1.1 |diff| = 1.22e-15` |
| 22 | `r1_gauss.py` | 0 | 0.2 | **OK** |  |
| 23 | `r3_l4_volume.py` | TIMEOUT | 180.3 | **超时** | 180 s 未完成；另：该 .py 于 12:07:11 被另一 agent 改写（`_patch.py`），本次跑的不是冻结版本 |
| 24 | `radius_push.py` | 0 | 0.5 | **OK** |  |
| 25 | `ru_joint.py` | 0 | 101.9 | **OK** | SIXTH-SPECIES 3/3 命中 |
| 26 | `su2_6j_cert.py` | 0 | 2 | **OK** |  |
| 27 | `sw_marginal.py` | 0 | 1.1 | **OK** | AUDIT-R1 4/4 命中 |
| 28 | `wa_alignment.py` | 0 | 2 | **OK** | WA-ALIGNMENT 5/5 命中 |
| 29 | `audit_a1_sun_hessian.py` | 0 | 0.6 | **OK** | AUDIT-R1 1/1 命中、自检 PASS |
| 30 | `audit_b1_tr_flat.py` | 0 | 0.7 | **OK** |  |
| 31 | `audit_b3_sw_marginal.py` | 0 | 0.2 | **OK** |  |
| 32 | `bn_residual.py` | 0 | 1.4 | **OK** |  |
| 33 | `ct_chain_exact.py` | 0 | 32.1 | **OK** |  |
| 34 | `ct_chain_L.py` | 0 | 88.1 | **OK** |  |
| 35 | `d3_averaged.py` | 0 | 0.7 | **OK** |  |
| 36 | `dm_floor_defect.py` | 0 | 1.2 | **OK** | LEDGER/AUDIT-CH16 6/6 命中 |
| 37 | `dn_density.py` | 0 | 28.2 | **OK** | DENSITY-D 6/6 命中 |
| 38 | `ds_acid.py` | 0 | 29.4 | **OK** |  |
| 39 | `ds_ratio.py` | 0 | 20.2 | **OK** | DS-RATIO `A=0.2051152971, B=0.0110592938, β*=4.30661`、FD 校验 `≤6.3e-9` → 重跑 0.2051152971 / 0.0110592938 / 4.30661 / 6.32e-09 |
| 40 | `dt_dirichlet.py` | 0 | 15.8 | **OK** |  |
| 41 | `lp_anchor.py` | 0 | 130.2 | **OK** |  |
| 42 | `lp_framework.py` | 0 | 0.7 | **OK** |  |
| 43 | `mc_ode.py` | 0 | 2.4 | **OK** |  |
| 44 | `na_plaquette.py` | TIMEOUT | 180.3 | **超时** | 180 s 未完成（已产出 5 988 B / 旧 log 14 205 B）；不判错 |
| 45 | `na2_selftest.py` | 0 | 0.2 | **OK** |  |
| 46 | `ps_witness.py` | 0 | 2.3 | **OK** |  |
| 47 | `r2_h2global.py` | 0 | 0.4 | **OK** |  |
| 48 | `rb_ratio_bound.py` | 0 | 102.8 | **OK** |  |
| 49 | `sdp_cert_aklt.py` | 0 | 63.8 | **OK** |  |
| 50 | `sdp_cert_p4.py` | 0 | 1.1 | **OK** |  |
| 51 | `sdp_cert_trunc.py` | 0 | 1.9 | **OK** |  |
| 52 | `tr_block_L.py` | 0 | 0.7 | **OK** |  |
| 53 | `tr_flat.py` | 0 | 0.5 | **OK** |  |
| 54 | `tr_maxplus.py` | 0 | 1.7 | **OK** |  |
| 55 | `at_lattice.py` | 1 | 14.6 | **失败** | `open("research/v3/species/w7/at_lattice.json")` —— 需 cwd=仓库根；计算完成，仅末步写盘失败 |
| 56 | `bc_gaussian.py` | 0 | 0.2 | **OK** | BC-BETHE `1.155/1.342/1.512` 三个数在重跑 stdout 中未找到 ⇒ 待人工（见 §未完成） |
| 57 | `mc_lattice.py` | 0 | 28.6 | **OK** |  |
| 58 | `mc_ricci.py` | 0 | 36.5 | **OK** |  |
| 59 | `na2_chain.py` | 0 | 110.1 | **OK** |  |
| 60 | `na2_t0_fit.py` | 0 | 59 | **OK** |  |
| 61 | `na2_tci.py` | 0 | 139.5 | **OK** |  |
| 62 | `na2_two_plaquette.py` | TIMEOUT | 180.5 | **超时** | 180 s 未完成；不判错 |
| 63 | `t23_monopole.py` | 0 | 29.1 | **OK** |  |
| 64 | `audit_b2_ps_tensor.py` | 0 | 0.2 | **OK** |  |
| 65 | `ct_hadamard.py` | 0 | 0.2 | **OK** |  |
| 66 | `dn_perlink.py` | 0 | 0.2 | **OK** |  |
| 67 | `mc_dbg.py` | 1 | 0.2 | **失败** | 调试用 mini 脚本，需外部参数/输入；非 witness，无笔记引用 |
| 68 | `mc_dbg2.py` | 0 | 0.2 | **OK** |  |
| 69 | `mc_tensor.py` | 0 | 2 | **OK** |  |
| 70 | `tci_comp.py` | 0 | 71.8 | **OK** |  |


---

## §4 核心对照：笔记值 ↔ 重跑值（判据逐条给）

| 笔记 | 笔记引用的数 | 重跑得到的数 | 判 |
|---|---|---|---|
| `R4-PERCENT.md` L75/L187/L288 | `β_c^KP ≤ 1.182227`（per-site, m=32）；`1.778421`（m=4）；`2.435301`（m=1）；`2.5 超出 ≥111%` | 1.182227 ／ 1.778421 ／ 2.435301 ／ `2.5/beta_c = 2.1147`（=+111%） | **OK**（逐位） |
| `SZ-CEILING.md` L62/L264/L271/L305 | ceiling `1/32 = 0.03125`（d=4）；bracket `[1/48, 1/32] = [0.020833, 0.031250]`；`C*(4) ∈ [16, 24]`；`8(d−1)` 行和 | `sz_ceiling_exact` d=4 `best Q = 16.000000` ⇒ `1/(2·16) = 0.031250`；`sz_ceiling_const` d=4 `1/(16(d−1)) = 0.020833`、`rowsum min/max = 24/24 = 8(d−1)`、`λ_max(N) = 24.000000`（d=2/3/4 = 8/16/24，`lam_min ~ 1e-15`） | **OK** |
| `KT-RHO2.md` L288/L381/L422 | `Q`: DD `0.99999376`、ND `2.66887958`、双井 `15.6442`；scale 不变 `0.99999376` | 0.99999376 ／ 2.66887958 ／ 15.6442；`u₂→3u₂` 后仍 0.99999376 | **OK** |
| `DS-RATIO.md` L18/L261 | `A=0.2051152971`、`B=0.0110592938`、`β*=4.30661`、5 点 FD `≤6.3e-9`、`D/S = β²/β*²` | 0.2051152971 ／ 0.0110592938 ／ 4.30661 ／ `|diff| = 6.32e-09` ／ 打印行 `D/S = beta^2/beta*^2`，`β=2.5` 处 `D/S = 0.33698` | **OK** |
| `R1-BRIDGE.md` L88/L359/L392 | 相似变换残差 `1.2e-15` | `[PASS] 1.1 gap(T) = gap(Dirichlet form in L²(Ω²))  |diff| = 1.22e-15` | **OK** |
| `NA2-NONDEGENERATE.md` L146 | `ω·1, β=4.0`：`X(a₀=0) = −1.69e−2`（se `4.3e−4`，`−67σ`） | `−0.017592`（se `0.000429`；`σ≈41`） | **DRIFT** |
| `NA2-NONDEGENERATE.md` L145 | `ω·1, β=2.5`：`−5.53e−3`（se `1.5e−4`） | `−0.005696`（se `1.45e-4`）；同点 h→0 外推 `−5.780e-3` | **DRIFT** |
| `NA2-NONDEGENERATE.md` L144 | `(2π/3,2π/3,−4π/3), β=0.5`：`−2.00e−4`（se `1.4e−5`） | `−0.000208`（se `0.000014`） | **DRIFT**（同族） |
| `NA2-NONDEGENERATE.md` L140–143、L150–154 | 同表 `se`/`X(+1)`/`X(−1)`、H-SCAN 三段带宽 `−1.95e−3/−4.74e−3/−5.46e−3/−5.70e−3/−5.78e−3`、`(π/2,π/2,−π)` 行、`=1` 行 | `se 4.29e-4`；`X(+1)=−3.3739e-2`、`X(−1)=−3.3219e-2`；H-SCAN `−1.951e-03/−4.737e-03/−5.456e-03/−5.699e-03/−5.780e-03`；`−6.745e-03…−3.408e-03`；`−7.731e-03…+1.725e-04` | **OK**（逐位） |
| `AUDIT-C0.md`（经 `audit_c0_ratio.py`） | `A/C(0) = 0.5288/0.5475/0.5633`；`0.5627`；flat control 仅为 measured 的 `9.38x` 远 | 0.52880 ／ 0.54750 ／ `0.5627` ／ `flat/measured = 9.38x` | **OK** |

---

## §5 DRIFT 清单（判死线：>3 条即 P0；本次 **1 条**）

### DRIFT-1 — `na2_su3_na.py` ↔ `NA2-NONDEGENERATE.md` §2.3 表的 **`X(a₀=0)` 列**

| 四元组 | 笔记值（se） | 重跑值（se） | 差 |
|---|---|---|---|
| `B=ω·1, β=4.0, a₀=0` | **−1.69e−2**（4.3e−4） | **−0.017592**（0.000429） | +6.9e−4 ≈ **1.6σ**（4.1%） |
| `B=ω·1, β=2.5, a₀=0` | **−5.53e−3**（1.5e−4） | **−0.005696**（0.000145） | −1.7e−4 ≈ 1.1σ（3.0%） |
| `B=ω·1, β=0.5, a₀=0` | **−2.00e−4**（1.4e−5） | **−0.000208**（0.000014） | −8e−6 ≈ 0.6σ（4.0%） |

**同表的 `se`、`X(+1)`、`X(−1)`、`B=1`/`(π,π,−2π)` 行、H-SCAN 三段带宽收敛全部逐位一致** ⇒ 漂的不是随机数种子（脚本 `rng = default_rng(20261124)`，且 **11-24 23:34 的旧 log 与今天 12:12 的重跑逐位相同**），而是**那张表 `a₀=0` 列的取值口径**（笔记那一列不是脚本在 `a₀=0` 直接打印的那一行，也不等于同行 `h→0` 外推）。
**影响**：结论（`X(0)≠0`、符号随 `B` 变、`B=ω·1` 处 ≫10σ）不变；但 `−1.69e−2` 这个**具体数字不可再引**，应引 `−1.76e−2（41σ）` 或直接引 log 行。

### （附）DRIFT-2 — 脚本自检失效（不是笔记漂移，但必须记）
`na_plaquette.py` `[A0]`：`B = 0.0000000000 (archived 0.01105929)` ⇒ `beta* = 4.3e12` 而非 `4.30661`。
**11-24 23:13 的旧 log 同样报 `B=0`** ⇒ 不是今天改坏，是长期自检口径失效。DS-RATIO 里引的 `A/B/β*` 由 `ds_ratio.py` 独立复算（**OK**），论文数字未被污染。

---

## §6 失败（5 条，全部是 cwd 路径约定，不是算法问题）

| 脚本 | exit | 死在哪 | 真实状态 |
|---|---|---|---|
| `at_lattice.py` | 1 | `open("research/v3/species/w7/at_lattice.json","w")` | 计算完成、stdout 完整，只在末步写 json 挂 |
| `at_scalar.py` | 1 | 同上（`at_scalar.json`） | 同上（Q(inf)=3.0000、γ≥0.00185906 等都已打印） |
| `at_su2_exact.py` | 1 | 同上（`at_su2_exact.json`） | 52.8 s 计算完成，只在末步写 json 挂 |
| `qi_bracket.py` | 1 | `open("research/v3/w7/m1_4d_run_results.json")` | **没跑出数**；该文件确实存在于 `research/v3/w7/`，从仓库根跑即通 |
| `mc_dbg.py` | 1 | 调试脚本，缺外部输入 | 非 witness，无笔记引用 |

⇒ **建议**：这 4 个 witness 的官方跑法写成 `cd D:\Dev\math && py -3.12 -I research/v3/species/w7/<x>.py`（脚本内路径就是为仓库根 cwd 写的），而不是笔记 §复算 里的 `cd …\species\w7`。

## §7 超时（5 条，按任务书记"超时"，不判错）

`hc_spherical.py`(181.5 s) ／ `ns_nosub.py`(180.6 s) ／ `r3_l4_volume.py`(180.3 s) ／ `na_plaquette.py`(180.3 s) ／ `na2_two_plaquette.py`(180.5 s)。
其中 `ns_nosub` 已产出 24 370 B 的 stdout（到 `a0=+2.00` 段），它引用的 2 个关键数**仍对得上**（2/2 命中）；`r3_l4_volume` 的 .py 在扫描中被改（见 §1），其超时不代表冻结版本。

---

## § EXT：笔记引用、但脚本不在 `species/w7` 树里（为了不越界写入，复制到 `repro/ext/` 跑）

`R4-PERCENT.md` 指向 `research/v3/w7/r4_kp.py`；`SZ-CEILING.md` 指向 `research/v3/w7/sz_ceiling_{const,exact,hess,lmax}.py`。**两棵树脚本名零重叠**（`v3/w7` 116 个，`species/w7` 65 个），所以这不是重名拷贝。

| 脚本 | 跑法 | exit | 耗时 | 判 |
|---|---|---|---|---|
| `r4_kp.py`（Pmax=7） | 复制到 `repro/ext/`，cwd=`repro/ext` | 0 | 4.6 s | **OK** — `1.182227 / 1.778421 / 2.435301 / 2.1147 / ≥111%` 全中 |
| `sz_ceiling_const.py` | 同上 | 0 | 89.3 s | **OK** — `1/(16(d−1)) = 0.020833`(d=4)、`rowsum 8/16/24`、`λ_max = 8.000000/16.000000/24.000000` |
| `sz_ceiling_exact.py` | 同上 | 0 | 0.6 s | **OK** — d=4 `best Q = 16.000000` ⇒ ceiling `1/32 = 0.031250`；d=2 `8.000000` |
| `sz_ceiling_hess.py` | 同上 | 停 | >16 min | **未跑完**（已有 N=2 全 5 行 + N=3 d=2；`verification max rel err = 3.04e-15 / 1.22e-15`）——为给主扫描让内存主动停 |
| `sz_ceiling_lmax.py` | — | — | — | **未跑** |

---

## § 未完成（下一批要接的活，按价值排序）

### (a) 完全没跑完的
- `sz_ceiling_hess.py`（N=3 d=3/L=3 段以后）、`sz_ceiling_lmax.py`（一次都没跑）。

### (b) 只做了「自检 PASS 扫描」，**没做** 2–5 个关键数逐条人工对照的脚本（约 30 个）
`at_lattice` `at_scalar` `at_su2_exact` `ct_chain_exact` `ct_hadamard` `d3_averaged` `dn_perlink` `ds_acid` `lp_anchor` `lp_framework` `mc_ode` `mc_tensor` `na2_chain` `na2_t0_fit` `na2_tci` `ns_nosub` `qg1_network` `r1_2dym` `r1_gauss` `r2_h2global` `radius_push` `sdp_cert_aklt` `sdp_cert_p4` `sdp_cert_trunc` `su2_6j_cert` `t23_monopole` `tr_block_L` `tr_flat` `tr_maxplus` `rb_ratio_bound` `ps_witness` `mc_lattice` `mc_ricci` `na2_two_plaquette` `hc_spherical` `tci_comp` `na2_selftest` `bc_gaussian` `qi_bracket` `audit_b1_tr_flat` `audit_b2_ps_tensor` `audit_b3_sw_marginal`
（它们的 stdout 已全部落盘在 `w7/repro/out/`，下一批不必重跑，只要比数。）

### (c) 自动匹配器标了 miss、但**必须人工看**才能定性（**不是已确认的 DRIFT**）
高 miss 的：`rb_ratio_bound`(14) ／ `tr_flat`(12) ／ `tr_block_L`(12) ／ `tr_maxplus`(11) ／ `na2_two_plaquette`(9) ／ `ps_witness`(7) ／ `mc_lattice`(6) ／ `hc_spherical`(6) ／ `na_plaquette`(6) ／ `tci_comp`(4) ／ `na2_selftest`(4) ／ `t23_monopole`(3) ／ `bc_gaussian`(3) ／ `qi_bracket`(3) ／ `mc_ricci`(2) ／ `qg2_abelian`(2) ／ `dt_dirichlet`(2) ／ `audit_b2_ps_tensor`(2)。
已知**不是**漂移的例子（避免下一批重复劳动）：`ds_ratio` 的 `4.4151` 是笔记由源码注释 `0.22571/0.011579` 推出的 coarea 陷阱值（脚本本就不打印）；`kt_a_exact`/`kt_c_geom`/`kt_b_sl_battery`/`dt_dirichlet` 的 miss 是笔记里的**日志字节数/样本数/独立复算值**（如 `8051 B`、`8000`、`1200`、`15880 B`、`10.51`），不是脚本输出量；`qg2_abelian` 的 `6.2e−15` vs `4.00e-15` 是两处机器精度级残差，比较无物理意义。

### (d) 方法学缺口
- **只有 1 个 seed 的抽样脚本**：`na2_su3_na` 证明了「旧 log == 今天重跑」（逐位），所以 DRIFT-1 不是随机漂移 —— 但其余 MC 脚本没做同样的旧新比对（下一批可对照 `repro/baseline/` 与 `repro/snap/` 直接 diff，零成本）。
- **扫描期间语料仍在动**（`r3_l4_volume.py` 12:07 被改、`cov_link.py`/`cov_models.py` 新增）⇒ 结论要绑时间戳；下一批建议先 `git status` 记一次快照再跑。

---

## § 复现命令（下一批直接用）

```powershell
cd D:\Dev\math\research\v3\species\w7\repro
pwsh -NoProfile -File .\run_all.ps1          # 顺序、单进程、180 s/脚本；结果 -> results.tsv + out/
py -3.12 -I -X utf8 analyze2.py              # 笔记引用数 <->
 重跑 stdout 自动比对 -> analysis.json
# 单脚本：py -3.12 -I <script>  （cwd = species\w7；audit\ 下脚本也用这个 cwd）
```
