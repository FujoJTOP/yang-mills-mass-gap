# Zenodo 元数据(逐项复制)

## Title

```
A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap — Second Edition: The Wall in Two Coordinates
```

## Version

```
2.0
```

## Publication date

```
2026-11-24
```

## Resource type

```
Publication / Preprint
```

## Creators

```
Jiang, Yuxuan
```

(与 v1 记录保持一致;不添加未确认的机构信息。)

## Language

```
English
```

## Keywords(Zenodo 用分号或逗号分隔,逐项粘贴)

```
Yang–Mills theory; mass gap; Clay Millennium Problem; lattice gauge theory; conformal window; Osterwalder–Schrader axioms; Lean 4; spectral measure; renormalization group; conformal bootstrap
```

## Description

见 `description.html`(整段 HTML 粘贴)。

## Files(本记录包含)

| 文件 | 说明 |
|---|---|
| `YANG-MILLS-CONDITIONAL-REDUCTION-V2.pdf` | 正文,81 页 |
| `ym-v2-lean.zip` | Lean 4 形式化源码补充包(40 个 `.lean`、构建脚本、`lakefile.lean`、toolchain pin、`MANIFEST.sha256`;不含 Mathlib/OSforGFF) |

Zenodo 的 Files 区上传这两个文件。Lean 包的构建入口与依赖说明见包内 `README.md`
和 `scripts/build-core.sh`;论文 Appendix A 的两个模块只依赖 Mathlib。

## Related identifiers

| Relation | Identifier | Scheme |
|---|---|---|
| `isVersionOf` | v1 记录页的 **Concept DOI**(若只显示一个 DOI,则用 `10.5281/zenodo.24486316`) | DOI |
| `isSupplementedBy` | `https://github.com/FujoJTOP/yang-mills-mass-gap` | URL |

## License

**保持与 v1 记录一致。** 不要在未确认的情况下改动许可;新版本默认继承 v1 的许可设置。

## Notes(可选,填入 Notes 字段)

```
Second edition (v2), 81 pages. Supersedes the 31-page v1 (DOI 10.5281/zenodo.24486316).
This is a conditional reduction and an audit, not a proof of the Clay problem: Clay remains open.
New in v2: Part III draws the obstruction in two independent coordinates (measure/weak-* and
beta/RG/conformal); machine-checked H1 sharpening (H1 alone does not imply an atom); Savvidy
"false obstruction" removed; tool-wave classification; lattice diagnostics; verification map;
errata (free-vector a/c corrected to 31/18). Clay = 0%; no QED.
```

## 与 v1 的差异(给审稿/读者的一页说明)

- **新增 Part III(约 20 页)**:弱-* 坐标(原子 = 边界密度的弱-* 不连续;有限数据连续;ONLY-EXIT)
  与 β/RG/共形坐标((★) 为真但化归 (G1);共形窗口真实;8/8 (G1) 路失败)。
- **H1 锐化的 Lean 机器验证**:`H1Sharpening.lean`,0 sorry;证明 `H1` + squeeze 只给
  `m ≥ λf`,而 `H1` 本身不推出 `f > 0`。
- **Savvidy 清理**:假障碍(负模 ⇒ 发散)被定理移除;余下边际近临界 / IR-marginal = H1。
- **工具波**:8 个复分析工具中的 7 个为盲/只读;Hardy 泛函只探测"存在"不探测"权重"。
- **格点诊断**:P3 逼近律、单侧隔离界、Re-slope、有限尺寸律;均标注为"读"而非"证"。
- **附录**:Lean、数值协议、后审计、验证表、测量台账、波次索引、勘误。
- **勘误**:自由矢量 `a/c = 31/18`(此前 `31/36` 为归一化滑档);结论不变。
- **附带 Lean 源码包** `ym-v2-lean.zip`:论文 Appendix A 的 `Reduction.lean` 与
  `H1Sharpening.lean`,以及 `Conditional.lean` 和 glue/handoff 开发;toolchain pin
  `v4.33.0-rc1`,Mathlib commit `79d0395a…`;不含第三方 `OSforGFF` 与 Mathlib 源码。

## 一句话摘要(用于社交媒体/邮件,可选)

```
v2 of the conditional-reduction paper: the Yang–Mills mass gap is drawn in two independent
coordinates (measure/weak-* and beta/RG/conformal); 81 pp; Clay still open (0%).
```
