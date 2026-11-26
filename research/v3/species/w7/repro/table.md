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
| 13 | `na2_su3_na.py` | 0 | 49.5 | **OK** | NA2 表 `X(±1)`、`se`、H-SCAN 三段带宽收敛逐位一致；但 `X(a₀=0)` 列不一致（见 §DRIFT-1） |
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

