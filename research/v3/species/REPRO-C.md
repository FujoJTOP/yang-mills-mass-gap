# REPRO-C · C 批 25 脚本复现（`py -3.12 -I`，单脚本 20 s 超时）

| 脚本 | exit | 秒 | 输出里的第一个数字 | 判 |
|---|---|---|---|---|
| ps_witness.py | 0 | 2.7 | 2 | OK |
| qg1_network.py | 0 | 2 | 1 | OK |
| qg2_abelian.py | - | 20 | (timeout-kill) | SKIP-slow |
| qg3_acid.py | 0 | 1.1 | 3 | OK |
| qg4_su2.py | 0 | 1.3 | 4 | OK |
| qi_bracket.py | 1 | 2 | -1 | FAIL |
| r1_2dym.py | 0 | 0.7 | 1 | OK |
| r1_chain.py | 0 | 2.5 | 1 | OK |
| r1_gauss.py | 0 | 0.2 | 1 | OK |
| r2_h2global.py | 0 | 0.4 | 2 | OK |
| r3_l4_volume.py | - | 20.1 | (timeout-kill) | SKIP-slow |
| radius_push.py | 0 | 0.6 | 2 | OK |
| rb_ratio_bound.py | - | 20 | (timeout-kill) | SKIP-slow |
| ru_joint.py | - | 20 | (timeout-kill) | SKIP-slow |
| sdp_cert_aklt.py | - | 20 | (timeout-kill) | SKIP-slow |
| sdp_cert_p4.py | 0 | 1.5 | 1 | OK |
| sdp_cert_trunc.py | 0 | 2.5 | 3 | OK |
| su2_6j_cert.py | 0 | 3.7 | 3 | OK |
| sw_marginal.py | 0 | 1.4 | 0.0 | OK |
| t23_monopole.py | - | 20 | (timeout-kill) | SKIP-slow |
| tci_comp.py | - | 20 | (timeout-kill) | SKIP-slow |
| tr_block_L.py | 0 | 0.7 | 2 | OK |
| tr_flat.py | 0 | 0.6 | 3 | OK |
| tr_maxplus.py | 0 | 1.9 | 4 | OK |
| wa_alignment.py | 0 | 2.4 | 2.4 | OK |

## 汇总

- 跑完 25 / 25（分 5 批，每批 5 个落盘一次）。
- **OK 17 · FAIL 1 · SKIP-slow 7**。
- SKIP-slow（20 s 内未结束，已 `Stop-Process`/taskkill 终止，不计失败）：`qg2_abelian.py`、`r3_l4_volume.py`、`rb_ratio_bound.py`、`ru_joint.py`、`sdp_cert_aklt.py`、`t23_monopole.py`、`tci_comp.py`。
- FAIL 唯一：`qi_bracket.py`（exit 1，2 s）— 首行错误：`FileNotFoundError: [Errno 2] No such file or directory: 'research/v3/w7/m1_4d_run_results.json'`（第 225 行硬编码相对路径，指向 `research/v3/w7/`，本目录无此文件）。
- `su2_6j_cert.py` 本次 3.7 s 正常跑完（未 SKIP）。
- 未比对任何笔记数字；未改任何 `.py`/`.md`/`.tex`；未 git commit。
