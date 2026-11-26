# 参考文献核验报告

核验范围：`refs.tex` 中 12 条指定 `\bibitem`（其余条目未动）。
核验日期：2026-11-25。所有书目数据均取自本轮实际抓取的页面；未抓到的字段一律不填。

**判定含义**：`OK` = 条目所述数据全部正确（可能偏简，建议补全，不构成错误）；`PARTIAL` = 部分数据错误 / 有误导性归属；`WRONG` = 关键识别符（DOI / arXiv）指向别的文献。

| 键 | 判定 | 正确书目 | 证据 URL | 检索式 | 工具状态 |
|---|---|---|---|---|---|
| `szz` | **OK**（缺标题与终止页；`1/48` 归属**正确**，但"log-Sobolev 半径"措辞需改） | H. Shen, R. Zhu, X. Zhu, *A stochastic analysis approach to lattice Yang–Mills at strong coupling*, Comm. Math. Phys. **400** (2023), no. 2, 805–851, DOI 10.1007/s00220-022-04609-1, arXiv:2204.12737 [math.PR] | https://arxiv.org/abs/2204.12737 · https://api.crossref.org/works/10.1007/s00220-022-04609-1 · https://api.semanticscholar.org/graph/v1/paper/DOI:10.1007/s00220-022-04609-1 | web_fetch arXiv abs 2204.12737；web_fetch Crossref DOI 10.1007/s00220-022-04609-1；web_fetch SemanticScholar 同 DOI；web_search `Shen Zhu Zhu lattice Yang-Mills "log-Sobolev" "1/48"`；web_search `"1/48" strong coupling lattice Yang-Mills spectral gap` | web_fetch OK(200)×3，命中 3；web_search OK，返回 10 源/10 源。附带失败：web_fetch `https://doi.org/10.1007/s00220-022-04609-1` → **TOOL-FAIL**（跨域重定向 idp.springer.com，未跟随）；web_fetch `https://link.springer.com/article/10.1007/s00220-022-04609-1` → **TOOL-FAIL**（同上）。已用 Crossref + Semantic Scholar 两条独立路补齐，结论不依赖失败路 |
| `dksst` | **OK**（配对**没有**错，怀疑不成立） | A. Dymarsky, Z. Komargodski, A. Schwimmer, S. Theisen, *On Scale and Conformal Invariance in Four Dimensions*, JHEP **10** (2015) 171, DOI 10.1007/JHEP10(2015)171, arXiv:1309.2921 [hep-th] | https://arxiv.org/abs/1309.2921 | web_fetch arXiv abs 1309.2921 | web_fetch OK(200)，命中 1。页面 "Journal reference: JHEP10(2015)171" 与 "Related DOI: 10.1007/JHEP10%282015%29171" 均由 arXiv 元数据直接给出：同一篇 |
| `balaban` | **PARTIAL**（I/II 与 Large field I/II 四条卷页**全对**；`119 (1988) 243--285` 的归属**错**） | T. Bałaban, *Renormalization group approach to lattice gauge field theories. I*, Comm. Math. Phys. **109** (1987) 249–301, DOI 10.1007/BF01215223；*II. Cluster expansions*, Comm. Math. Phys. **116** (1988) 1–22, DOI 10.1007/BF01239022；*Large field renormalization. I. The basic step of the R operation*, Comm. Math. Phys. **122** (1989) 175–202；*II. Localization, exponentiation, and bounds for the R operation*, Comm. Math. Phys. **122** (1989) 355–392。另**独立**一篇：*Convergent renormalization expansions for lattice gauge theories*, Comm. Math. Phys. **119** (1988) 243–285, DOI 10.1007/BF01217741 —— **不属于** Large field 系列 | https://reo.nii.ac.jp/oja/vol_issues?vol_id=ISS0010021632 · https://reo.nii.ac.jp/oja/vol_issues?vol_id=ISS0010021633 · https://api.crossref.org/works?query.bibliographic=Balaban+Convergent+renormalization+expansions+for+lattice+gauge+theories · https://api.crossref.org/works?filter=doi:10.1007/BF01215223,doi:10.1007/BF01211762 | web_fetch NII-REO CMP v122#2（命中 Large field I, p.175-202）；web_fetch NII-REO CMP v122#3（命中 Large field II, p.355-392）；web_fetch Crossref 检索（命中 BF01217741 = CMP 119, 243-285；BF01239022 = CMP 116, 1-22）；web_fetch Crossref filter（命中 BF01215223 = CMP 109, 249-301）；web_search×3 | web_fetch OK(200)×4，命中 5 条书目记录；web_search OK，返回 10/10/10 源。另由 Magnen–Rivasseau–Sénéor 论文自身参考文献表独立复核该系列卷页：`95,17; 96,223; 98,17; 99,75; 99,389; 102,277; 109,249; 116,1; 119,243; 122,175; 122,355` |
| `sw` | **PARTIAL**（"as cited in the Clay statement" 不是可引书目；年份 1964 正确） | R. F. Streater, A. S. Wightman, *PCT, Spin and Statistics, and All That*, W. A. Benjamin, New York, 1964（1999/2000 由 Princeton University Press 收入 *Princeton Landmarks in Physics* 重印） | https://search.worldcat.org/fr/title/930068 · https://nbinet3.ncl.edu.tw/search~S10*chi?/cQC174.45+S77 · https://www.science.org/doi/10.1126/science.145.3631.475.a | web_search `Streater Wightman "PCT, Spin and Statistics, and All That" 1964 W.A. Benjamin publisher`；web_search `Wightman axioms Clay Millennium Yang-Mills Streater Wightman reference` | web_search OK，返回 10 源/10 源，命中 4 条独立书目记录（WorldCat、NBINet 台湾国图书目、加拿大 FSL 书目、*Science* 1964 年书评 145(3631):475）。书评给出 "Benjamin, New York, 1964"；重印版由检索片段给出，标注为"重印"而非替代版 |
| `hofman-li` | **OK** | D. M. Hofman, D. Li, D. Meltzer, D. Poland, F. Rejón-Barrera, *A Proof of the Conformal Collider Bounds*, JHEP **06** (2016) 111, DOI 10.1007/JHEP06(2016)111, arXiv:1603.03771 [hep-th] | https://arxiv.org/abs/1603.03771 | web_fetch arXiv abs 1603.03771 | web_fetch OK(200)，命中 1。作者全名、标题、arXiv 号、DOI（=JHEP06(2016)111）与 `refs.tex` 一致 |
| `fo` | **OK**（缺标题） | J. Feldman, K. Osterwalder, *The Wightman axioms and the mass gap for weakly coupled $(\phi^4)_3$ quantum field theories*, Ann. Physics **97** (1976) 80–135 | https://par.nsf.gov/servlets/purl/10442426 · https://www.e-periodica.ch/cntmng?pid=hpa-001:1980:53::695 · https://www.mathnet.ru/links/a7789f6c1b329c9a0eace5f9d07602b6/intv32.pdf | web_search `Feldman Osterwalder "The Wightman axioms and the mass gap for weakly coupled phi^4" Annals of Physics 97 1976` | web_search OK，返回 10 源，命中 3 条独立引用（NSF-PAR 论文参考文献、Helv. Phys. Acta 1980 论文、Math-Net.Ru），均为 `Ann. Physics 97 (1976) 80-135` |
| `bgdf` | **OK**（缺标题与终止页 487） | D. C. Brydges, C. Giffen, B. Durhuus, J. Fröhlich, *Surface representations of Wilson loop expectations in lattice gauge theory*, Nucl. Phys. B **275** (1986), no. 3, 459–487, DOI 10.1016/0550-3213(86)90609-7 | https://api.crossref.org/works/10.1016/0550-3213(86)90609-7 · https://www.sciencedirect.com/journal/nuclear-physics-b/vol/275/issue/3 | web_fetch Crossref DOI 10.1016/0550-3213(86)90609-7；web_search `Brydges Giffen Durhuus Frohlich "Non-Gaussian fixed points of block-spin transformations" Nuclear Physics B 275 1986 459` | web_fetch OK(200)，命中 1（Crossref 权威记录：v275, i3, 459-487, 1986-11）；web_search OK，10 源，`Pages 459-487` 二次确认。附带：我先前猜的标题 "Non-Gaussian fixed points…" **不是**本篇，已弃用 |
| `gopfert` | **PARTIAL**（卷 82 / 页 545–606 正确；**年份应为 1982**） | M. Göpfert, G. Mack, *Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant*, Comm. Math. Phys. **82** (1982), no. 4, 545–606, DOI 10.1007/BF01961240 | https://reo.nii.ac.jp/oja/vol_issues?vol_id=ISS0010021475 · https://www.scopus.com/pages/publications/0000492747 | web_fetch NII-REO CMP v82#4 目录页；web_search `Göpfert Mack Comm. Math. Phys. 82 1981 545 "Proof of confinement of static quarks"` | web_fetch OK(200)，命中 1（卷头明写 "Volume 82, Issue 4 p.457-606 **January 1982**"，该文 p.545-606）；web_search OK，10 源，Scopus 记 `82, Issue 4, Pages 545-606`。年份差异说明：Springer 卷期为 1982，部分文献按收稿年写 1981 |
| `luscher77` | **OK** | M. Lüscher, *Construction of a selfadjoint, strictly positive transfer matrix for Euclidean lattice gauge theories*, Comm. Math. Phys. **54** (1977) 283–292, DOI 10.1007/BF01614090；*Some analytic results concerning the mass spectrum of Yang–Mills gauge theories on a torus*, Nucl. Phys. B **219** (1983), no. 1, 233–261 | https://ui.adsabs.harvard.edu/abs/1983NuPhB.219..233L · https://mathscinet.ams.org/mathscinet/relay-station?mr=706403 · https://api.crossref.org/works?filter=doi:10.1007/BF01215223,doi:10.1007/BF01211762（本页给出的 BF02097397 参考文献表中含 `Lüscher, CMP 54, 283 (1977)`） | web_search `Luscher "A semiclassical formula for the mass gap" Comm. Math. Phys. 54 1977 283`；web_search `Lüscher "Some analytic results concerning the mass spectrum of Yang-Mills gauge theories on a torus" Nuclear Physics B 219 1983` | web_search OK，10 源/10 源，命中：CMP 54, 283 + DOI BF01614090（Regensburg 学位论文、zenodo 引文表）；NPB 219 (1983) no.1, 233–261（ADS bibcode `1983NuPhB.219..233L` + MathSciNet MR0706403）。终止页 292 由 MRS 论文参考文献表与多条独立引文确认 |
| `kp` | **OK**（缺标题） | R. Kotecký, D. Preiss, *Cluster expansion for abstract polymer models*, Comm. Math. Phys. **103** (1986), no. 3, 491–498, DOI 10.1007/BF01211762 | https://api.crossref.org/works?filter=doi:10.1007/BF01215223,doi:10.1007/BF01211762 · https://www.scopus.com/pages/publications/0000639657 | web_fetch Crossref filter（DOI 精确取记录）；web_search `Kotecky Preiss "Cluster expansion for abstract polymer models" Comm. Math. Phys. 103 1986 491` | web_fetch OK(200)，命中 1（Crossref：v103, 491-498, 1986-09）；web_search OK，10 源，Scopus `103, Issue 3, Pages 491-498` + DOI BF01211762 二次确认 |
| `osterwalder-seiler` | **OK**（缺标题与终止页 471） | K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. **110** (1978), no. 2, 440–471, DOI 10.1016/0003-4916(78)90039-8 | https://www.scopus.com/pages/publications/33744753218 · https://api.crossref.org/works/10.1007/bf02097397（MRS 参考文献表 item 11） | web_search `Osterwalder Seiler "Gauge field theories on a lattice" Annals of Physics 110 1978 440`；web_fetch Crossref DOI 10.1007/bf02097397（取其参考文献表） | web_search OK，10 源，Scopus `110, Issue 2, Pages 440-471`；web_fetch OK(200)，Crossref 参考文献表第 11 条给全 DOI + 卷页 + 年 |
| `mrs` | **WRONG**（卷/页/年**全对**，但 **DOI 是错的**：`BF02079397` 指向一篇流体力学论文） | J. Magnen, V. Rivasseau, R. Sénéor, *Construction of $\mathrm{YM}_4$ with an infrared cutoff*, Comm. Math. Phys. **155** (1993), no. 2, 325–383, DOI **10.1007/BF02097397**（正确）；原文写的 `10.1007/BF02079397` 实为 *Groundwater flow in a medium with periodic inclusions*, Fluid Dynamics **30** (1995) 758–766, A. R. Kasimov & Yu. V. Obnosov | https://api.crossref.org/works/10.1007/BF02079397（错误 DOI 的真身） · https://api.crossref.org/works/10.1007/bf02097397（正确记录） · https://reo.nii.ac.jp/oja/ART1000345241 | web_fetch Crossref 错误 DOI `10.1007/BF02079397`；web_fetch Crossref 书目检索 `Construction of YM4 with an infrared cutoff Magnen Rivasseau Seneor`；web_fetch Crossref 正确 DOI `10.1007/bf02097397`；web_search `Magnen Rivasseau Seneor "Construction of YM4 with an infrared cutoff" Comm. Math. Phys. 155 1993` | web_fetch OK(200)×3 全部有效：误 DOI 命中 1（Fluid Dynamics 30(5) 758-766, 1995，作者 Kasimov/Obnosov）；书目检索命中 1（`10.1007/bf02097397`, CMP 155, 325-383, 1993-07）；正确 DOI 复核命中 1；web_search OK，10 源，NII-REO 记 `Volume155, Issue2, July 1993, p.325-383`。**这是 DOI 数字换位（797 ↔ 973），不是"查不到"** |

## szz 专项结论（`1/48` 溯源）

回答题目三个子问题：

**(a) arXiv:2204.12737 到底对应哪篇论文？**
就是 Shen–Zhu–Zhu 那篇。标题 *A stochastic analysis approach to lattice Yang–Mills at strong coupling*，作者 Hao Shen、Rongchan Zhu、Xiangchan Zhu，math.PR，43 页，2022-04-27 提交，Related DOI `10.1007/s00220-022-04609-1`。**arXiv 号、作者、期刊、卷、起页全部与 `refs.tex` 相符**。

**(b) 是否存在 "Shen, Zhu, Zhu" 的格点规范理论 log-Sobolev 论文？**
存在，且就是同一篇。摘要原文（arXiv abs 页，本轮抓取）：假设
`|β| < (N-2)/(32(d-1)N)`（SO(N)）或 `|β| < 1/(16(d-1))`（SU(N)），则无穷体积极限的 YM 测度满足 **Log-Sobolev 与 Poincaré 不等式**，并有正质量间隙。$d=4$ 时 $1/(16\cdot 3) = 1/48$。

**(c) `1/48` 出处是否另有其人？**
**没有另一个人 —— 就是本篇。** 但论文正文的措辞必须改：

- `1/48` 是**定理的假设窗口**（对 $|\beta|$ 的上界，$d=4$、SU(N)、t'Hooft 标度 $\beta N$），**不是**某个"log-Sobolev 半径"的数值。该文在窗口内证明 LSI 成立；把 `1/48` 说成"log-Sobolev 半径"是把**假设**当成了**结论常数**。
- 建议正文改写为："在 $|\beta|<1/48$（$d=4$, $SU(N)$）的强耦合窗口内，Shen–Zhu–Zhu 证明了无穷体积 YM 测度满足 log-Sobolev 与 Poincaré 不等式"。
- 未能在本轮证据中确证：该文正文是否另给一个形如 $c\cdot e^{-\mathrm{const}\cdot\beta}$ 的**显式 LSI 常数**。arXiv 全文 HTML 未抓取（未做逐字核）。此点标 **UNVERIFIED**，请勿在论文里写具体常数表达式。

## 修正建议

以下 `\bibitem` 可直接替换 `refs.tex` 中同键条目（LaTeX 转义已按原文件风格处理）。

```latex
\bibitem{szz} H.~Shen, R.~Zhu, X.~Zhu, \emph{A stochastic analysis approach to lattice
  Yang--Mills at strong coupling}, Comm.\ Math.\ Phys.\ \textbf{400} (2023) 805--851,
  DOI 10.1007/s00220-022-04609-1, arXiv:2204.12737 [math.PR].
```

```latex
\bibitem{dksst} A.~Dymarsky, Z.~Komargodski, A.~Schwimmer, S.~Theisen, \emph{On scale and
  conformal invariance in four dimensions}, JHEP \textbf{10} (2015) 171,
  DOI 10.1007/JHEP10(2015)171, arXiv:1309.2921 [hep-th].
```
（内容与原文一致，仅补 DOI。原文的 arXiv/JHEP 配对**是对的**，无需改号。）

```latex
\bibitem{balaban} T.~Ba{\l}aban, \emph{Renormalization group approach to lattice gauge field
  theories.\ I}, Comm.\ Math.\ Phys.\ \textbf{109} (1987) 249--301,
  DOI 10.1007/BF01215223; \emph{II.\ Cluster expansions}, Comm.\ Math.\ Phys.\
  \textbf{116} (1988) 1--22, DOI 10.1007/BF01239022; \emph{Large field renormalization.\
  I.\ The basic step of the R operation}, Comm.\ Math.\ Phys.\ \textbf{122} (1989)
  175--202; \emph{II.\ Localization, exponentiation, and bounds for the R operation},
  Comm.\ Math.\ Phys.\ \textbf{122} (1989) 355--392.
\bibitem{balaban-conv} T.~Ba{\l}aban, \emph{Convergent renormalization expansions for
  lattice gauge theories}, Comm.\ Math.\ Phys.\ \textbf{119} (1988) 243--285,
  DOI 10.1007/BF01217741.
```
（原来把 `119 (1988) 243--285` 塞进 "Large field renormalization I--II" 的括号里，属于**归属错误**。它是一篇独立论文；若正文要引，请用新键 `balaban-conv`，并在正文相应位置改写引用键。）

```latex
\bibitem{sw} R.~F.~Streater, A.~S.~Wightman, \emph{PCT, Spin and Statistics, and All That},
  W.~A.~Benjamin, New York, 1964; reprinted in \emph{Princeton Landmarks in Physics},
  Princeton University Press, 2000.
```

```latex
\bibitem{hofman-li} D.~M.~Hofman, D.~Li, D.~Meltzer, D.~Poland, F.~Rej\'on-Barrera, \emph{A
  proof of the conformal collider bounds}, JHEP \textbf{06} (2016) 111,
  DOI 10.1007/JHEP06(2016)111, arXiv:1603.03771 [hep-th].
```
（仅补 DOI；其余原文已正确。）

```latex
\bibitem{fo} J.~Feldman, K.~Osterwalder, \emph{The Wightman axioms and the mass gap for
  weakly coupled $(\phi^4)_3$ quantum field theories}, Ann.\ Physics \textbf{97} (1976)
  80--135.
```

```latex
\bibitem{bgdf} D.~C.~Brydges, C.~Giffen, B.~Durhuus, J.~Fr\"ohlich, \emph{Surface
  representations of Wilson loop expectations in lattice gauge theory}, Nucl.\ Phys.\ B
  \textbf{275} (1986) 459--487, DOI 10.1016/0550-3213(86)90609-7.
```

```latex
\bibitem{gopfert} M.~G\"opfert, G.~Mack, \emph{Proof of confinement of static quarks in
  3-dimensional $U(1)$ lattice gauge theory for all values of the coupling constant},
  Comm.\ Math.\ Phys.\ \textbf{82} (1982) 545--606, DOI 10.1007/BF01961240.
```
（年份 `1981` → `1982`。）

```latex
\bibitem{luscher77} M.~L\"uscher, \emph{Construction of a selfadjoint, strictly positive
  transfer matrix for Euclidean lattice gauge theories}, Comm.\ Math.\ Phys.\ \textbf{54}
  (1977) 283--292, DOI 10.1007/BF01614090; \emph{Some analytic results concerning the mass
  spectrum of Yang--Mills gauge theories on a torus}, Nucl.\ Phys.\ B \textbf{219} (1983)
  233--261.
```

```latex
\bibitem{kp} R.~Koteck\'y, D.~Preiss, \emph{Cluster expansion for abstract polymer models},
  Comm.\ Math.\ Phys.\ \textbf{103} (1986) 491--498, DOI 10.1007/BF01211762.
```

```latex
\bibitem{osterwalder-seiler} K.~Osterwalder, E.~Seiler, \emph{Gauge field theories on a
  lattice}, Ann.\ Phys.\ \textbf{110} (1978) 440--471, DOI 10.1016/0003-4916(78)90039-8.
```

```latex
\bibitem{mrs} J.~Magnen, V.~Rivasseau, R.~S\'en\'eor, \emph{Construction of
  $\mathrm{YM}_4$ with an infrared cutoff}, Comm.\ Math.\ Phys.\ \textbf{155} (1993)
  325--383, DOI 10.1007/BF02097397.
```
（**只改 DOI 的两位数字**：`BF02079397` → `BF02097397`。卷、页、年、标题、作者原样保留。）

### 优先级

1. **`mrs` 的 DOI** —— 现在会把读者送到一篇地下水流的论文，是最硬的一处错。（一个换位数字）
2. **`gopfert` 的年份** —— `1981` → `1982`。
3. **`balaban` 的 `119 (1988) 243--285` 归属** —— 它不是 Large field 系列。
4. **`sw`** —— 需要一条真书目，不能写 "as cited in the Clay statement"。
5. **`szz` 正文措辞** —— `1/48` 是窗口阈值，不是 log-Sobolev 半径；条目本身不用改号。

### 遗留 UNVERIFIED

- `szz` 正文中是否存在**显式** LSI 常数（形如 $c\,e^{-\kappa\beta}$）：未抓全文，**UNVERIFIED**。
- Balaban I 的副标题（"Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions"）本轮只见到检索片段，未取到权威全称，故修正建议里**只写主标题 + `I`**，不补副标题。
- Streater–Wightman 的页码数（x + 181 / viii + 181）两个来源不一致，故上条书目**不写页数**。
