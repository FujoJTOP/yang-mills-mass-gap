# Zenodo 第二版上载材料

> **目的**:把 81 页的第二版作为**现有记录的新版本**(new version)发布,而不是新建一条记录。
> 这样概念 DOI(concept DOI)保持不变,引用链连续;新版本会得到一个新的 version DOI。

## 0. 待上载的文件

| 项 | 路径 | 说明 |
|---|---|---|
| PDF | `docs/YANG-MILLS-CONDITIONAL-REDUCTION-V2.pdf` | 81 页,1.45 MB |
| 校验 | `docs/zenodo-v2/SHA256SUMS.txt` | 可核对是否传错文件 |
| 描述 | `docs/zenodo-v2/description.html` | 粘贴到 Zenodo 的 Description 字段 |
| 元数据 | `docs/zenodo-v2/zenodo-metadata.md`(人读)/ `zenodo.json`(机读) | 逐项填入 |

**不要上传**:`docs/paper-v2/` 的源码、`main.aux/log/toc`、任何 `out*.txt`。Zenodo 记录只需要 PDF;
源码留在 GitHub 仓库即可。

## 1. 前置检查(30 秒)

- PDF 页数 = **81**(`pdfinfo` 或打开文件看最后一页)。
- PDF 里没有 `??`(未解析引用)。
- 文件哈希与 `SHA256SUMS.txt` 一致:
  `c793ff0d3c87bcd8e6da65a3f6e562aabe156769c91d836bf01c297adb223826`
- 标题页写明 "Second Edition — The Wall in Two Coordinates"、"version 2"、"supersedes the 31-page v1"。

## 2. 在 Zenodo 上发布新版本

1. 登录 Zenodo,打开 v1 记录:**https://doi.org/10.5281/zenodo.24486316**
2. 在记录页面点 **New version**(不要点 New upload)。
3. 在 Files 区:把 v1 的 PDF 替换为 `docs/YANG-MILLS-CONDITIONAL-REDUCTION-V2.pdf`
   (或保留旧文件并新增,二选一;推荐**替换**,记录里只留当前版本)。
4. **Version** 填 `2.0`;**Publication date** 填 `2026-11-24`。
5. **Title** 用:
   `A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap — Second Edition: The Wall in Two Coordinates`
6. **Description**:打开 `description.html`,整段粘贴(HTML 可直接渲染)。
7. **Creators**:`Jiang, Yuxuan`(与 v1 一致)。
8. **Keywords**:见 `zenodo-metadata.md` 的列表,逐项粘贴。
9. **License**:保持与 v1 记录一致(不要在没有确认的情况下改动许可)。
10. **Related identifiers**:
    - `isVersionOf` → v1 记录页显示的 **Concept DOI**;若页面只显示一个 DOI,就用
      `10.5281/zenodo.24486316`。
    - `isSupplementedBy` → `https://github.com/FujoJTOP/yang-mills-mass-gap`(URL)。
11. 检查预览,然后 **Publish**。Zenodo 会给出新的 version DOI;概念 DOI 不变。

## 3. 发布之后

- 把新的 version DOI 记到仓库 `README.md` 的论文小节(以及 `docs/YANG-MILLS-MASTER.md` 若有需要)。
- 若期刊(Rev. Math. Phys.)要求提交更新稿,用同一个 PDF;附一页 response,说明相对 v1 的新增内容
  (Part III 两坐标、H1 机器验证、Savvidy 清理、工具波、诊断、勘误)。
- 不要删除 v1:它是历史版本,概念 DOI 会指向最新版,但 v1 的 version DOI 仍可引用。

## 4. 一句话状态(给审稿人/读者)

**这不是 Clay 的证明。** 第二版仍然是条件归约 + 审计:把墙画成两个独立坐标(测度/弱-* 与
β/RG/共形),把已证定理、no-go、精确缺失步和数值证据逐条标注。`Clay = 0%`,无 QED。
