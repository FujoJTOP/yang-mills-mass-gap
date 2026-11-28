// ============================================================
//  从已证结论到 Clay：一条公式链
//  4 维纯 SU(2) Yang–Mills 质量隙 · 第 366 轮 · 总装件
//  build:  typst compile docs/chain-to-clay.typ docs/CHAIN-TO-CLAY.pdf --root .
// ============================================================

#set page(paper: "a4", margin: (x: 2.15cm, y: 2.35cm), numbering: "1")
#set text(font: ("New Computer Modern", "Noto Serif SC"), size: 10.2pt, lang: "zh", fallback: true)
#set par(justify: true, leading: 0.82em, spacing: 1.05em)
#set math.equation(numbering: "(1)")
#set heading(numbering: "1.")

// ---------- 颜色 ----------
#let cT = rgb("#1a7f37")
#let cC = rgb("#8a6100")
#let cU = rgb("#5b6470")
#let cF = rgb("#a4101c")
#let cM = rgb("#16558f")
#let cAcc = rgb("#2b5f9e")
#let cNote = rgb("#6b7280")

// ---------- 判级 ----------
#let T = text(fill: cT, weight: "bold")[[T]]
#let C = text(fill: cC, weight: "bold")[[C]]
#let U = text(fill: cU, weight: "bold")[[U]]
#let F = text(fill: cF, weight: "bold")[[F]]
#let M = text(fill: cM, weight: "bold")[[M]]
#let E = text(fill: rgb("#5b3a99"), weight: "bold")[[E]]

// ---------- 连接词（左）+ 判级/出处（右） ----------
#let c(w, note: none) = block(width: 100%, above: 1.05em, below: -0.55em, breakable: false)[
  #text(size: 8.8pt, style: "italic", weight: "bold", fill: cAcc)[#w]
  #if note != none [
    #h(0.8em)
    #text(size: 7.8pt, fill: cNote)[· #note]
  ]
]
// ---------- 定理框 ----------
#let thmbox(kind, body) = block(
  width: 100%, inset: (x: 0.8em, y: 0.6em), radius: 2pt,
  stroke: (left: 2.4pt + cAcc), fill: rgb("#f4f7fb"), above: 1.25em, below: 1.2em
)[
  #text(weight: "bold", fill: rgb("#1b2a3a"))[#kind] #h(0.3em) #body
]

// ============================================================
#align(center)[
  #text(size: 17pt, weight: "bold")[从已证结论到 Clay：一条公式链]
  #v(0.4em)
  #text(size: 10pt, fill: rgb("#4a525c"))[4 维纯 SU(2) Yang–Mills 质量隙（Clay）· 第 366 轮 · 总装件]
]

#v(1em)

#block(
  inset: (x: 0.95em, y: 0.75em), stroke: (left: 2.5pt + cAcc), fill: rgb("#f7f9fc"), width: 100%
)[
  #text(size: 9.2pt)[
    *摘要。* 本文把本项目中全部已证结论按“因为 — 那么 — 所以 — 得到 — 设”接成一条公式链，终点是
    Clay 命题的等价两环形式
    $ "Clay" <=> ("A3") and ("SCALE") , $
    其中 $("A1")$ 是质量比 $m_"lat" slash a (beta)$ 的一致下界，$("A2")$ 是其一致上界，
    $("A3")$ 是沿渐近自由轨迹的 Osterwalder–Schrader 收敛性。链上除这三条外每一环均为已证。
    文中另给出两段已证的反证：$chi$ 类界只能给出质量的上界（方向定理），以及单链结果到格点的桥为假。
  ]
]

#v(0.45em)

#text(size: 8.8pt, fill: rgb("#4a525c"))[
  判级：#T 已证 · #C 条件定理（条件逐条写明） · #U 未证 · #F 已证否 · #M 数值证据 · #E 外部文献。
  每条箭头只有两种：$arrow.r.double$（已证蕴含）与 $arrow.r.double.not$（已证不蕴含），无第三种。
]

#v(0.3em)

= 目标与记号

#thmbox[定义 1（Clay）][对 $G = "SU"(2)$、$d = 4$：
  $ "Clay" := ("C1") and ("C2") , $ 其中
  $ ("C1") quad exists "非平凡 Osterwalder–Schrader 理论" mu "于" R^4 , $
  $ ("C2") quad exists Delta > 0 : quad "spec" H subset {0} union [Delta, infinity) . $
]

#c("记")
$ Lambda_L = {0, dots, L-1}^4, quad U_e in "SU"(2), quad q_P = 1/2 "tr" U_P, quad beta = 4 slash g^2 $

$ cal(S) = sum_P (1 - q_P) >= 0 , quad S_beta = beta cal(S) , quad T_(beta,L) = e^(-H_(beta,L)) , quad m_"lat" := -ln lambda_1 (T) $

#c("记", note: [#T 谱表示])
$ hat(C)_L (k) = ⟨ D_k, cal(H)^(-1) D_k ⟩ = integral_((0,infinity)) (dif nu_k (lambda)) / lambda , quad nu_k >= 0 , quad lambda_min (k) := min "supp" nu_k $

#c("设", note: [#T 无条件])
$ mu_(beta,L) := Z_(beta,L)^(-1) e^(-beta cal(S)) product_e dif U_e , quad Z_(beta,L) in (0, infinity) . $

= 主链：从 Clay 的定义到等价形式

#c("因为")
$ "Clay" := ("C1") and ("C2") $

#c("那么", note: [#E OS 重构])
$ mu_(beta,L) " 反射正 " quad ( "Wilson 作用量" plus "Haar 测度" ) $

#c("所以", note: [#T])
$ arrow.r.double quad exists ( cal(H)_(beta,L), T_(beta,L) = e^(-H_(beta,L)) ) : quad 0 <= T <= 1 , quad H >= 0 , quad Omega " 唯一" $

#c("所以", note: [#T 谱定理])
$ m_"lat" (beta,L) := - ln lambda_1 (T_(beta,L)) = inf ( "spec" H_(beta,L) without {0} ) $

#c("得到", note: [#T + #E 渐近自由])
$ m_"phys" (beta,L) = m_"lat" (beta,L) slash a (beta) , quad a(beta) = Lambda^(-1) e^(-beta slash (2 b_0)) (1 + o(1)) arrow.r 0 $

#c("设", note: [#U 下界])
$ ("A1") quad exists c_1 > 0, beta_0 : quad forall beta >= beta_0, forall L >= L_0 : quad m_"lat"(beta,L) slash a (beta) >= c_1 $

#c("设", note: [#U 上界])
$ ("A2") quad exists c_2 < infinity, beta_0 : quad forall beta >= beta_0, forall L >= L_0 : quad m_"lat"(beta,L) slash a (beta) <= c_2 $

#c("那么", note: [#T 同一比值的两头夹])
$ ("SCALE") := ("A1") and ("A2") : quad c_1 <= m_"lat"(beta,L) slash a (beta) <= c_2 $

#c("所以", note: [#T])
$ ("A1") and ("A3") quad arrow.r.double quad ("C2") : quad "取" Delta := c_1 , quad "spec" H subset {0} union [c_1, infinity) $

#c("所以", note: [#T 非平凡])
$ ("A2") and ("A3") quad arrow.r.double quad exists O, R_0 : quad lim_(R -> infinity) (⟨ O(0) O(R) ⟩)_c != 0 $

#c("设", note: [#U + #E Balaban 部分])
$ ("A3") quad "沿" beta_n -> infinity , a_n -> 0 , L_n -> infinity : quad S_n -> S_infinity , quad S_infinity " 满足 OS" $

#c("那么", note: [#T])
$ ("A1") and ("A2") and ("A3") quad arrow.r.double quad ("C1") and ("C2") = "Clay" $

#c("得到", note: [#T])
$ "Clay" quad arrow.l.double quad ("A3") and ("SCALE") $

#c("因为", note: [#T 谱的收敛])
$ "Clay" quad arrow.r.double quad "沿子列" m_"lat" slash a -> m_"phys" in (0, infinity) $

#c("所以", note: [#T 等价两环形式])
$ "Clay" quad arrow.l.r.double quad ("A3") and ("SCALE") $

= 输入引理：已证的那一半

#c("因为", note: [#T])
$ T_beta = M_(V^(1 slash 2)) K M_(V^(1 slash 2)) >= 0 quad ( I_(2j+1)(beta) > 0 ) , quad B_kappa^* = B_(-kappa) $

#c("那么", note: [#T])
$ W_mu (n; kappa) = sum_(z_perp) Gamma_L (0, (n, z_perp)) e^(i kappa z_perp) >= 0 , quad w(n) := W_mu (n; 0) >= 0 $

#c("所以", note: [#T 全四维动量])
$ ("Q") : quad hat(C)_L (k) <= hat(C)_L (0) quad forall k in (2 pi slash L) Z_L^4 $

#c("得到", note: [#T])
$ norm(Gamma_L)_(2 -> 2) = max_k hat(C)_L (k) = hat(C)_L (0) = chi_L^"fixed" $

#c("因为", note: [#T 恒等式])
$ dif / (dif beta) ⟨ cal(S) ⟩ = "Var"(cal(S)) >= 0 , quad ⟨ cal(S) ⟩_(beta=0) = 6 L^4 , quad ⟨ cal(S) ⟩_(beta=infinity) = 0 $

#c("那么", note: [#T 无条件])
$ integral_0^infinity chi_L^"clover" (beta) dif beta = 1 slash 6 quad ( "一切" L >= 2 ) $

#c("所以", note: [#T])
$ chi_L^"clover" (beta) <= 1 slash 6 quad forall L >= 2 , forall beta > 0 $

#c("得到", note: [#C 非 Abel GKS 未证])
$ chi_L^"fixed" <= 1 , quad norm(Gamma_L)_(2 -> 2) <= 1 quad arrow.l.double quad v_2 >= 0 $

= 为什么 $chi$ 类进不了主链

#c("设")
$ chi " 类候选：" quad hat(C)_L (k) <= C quad ( C " 与 " L, beta " 无关" ) $

#c("因为", note: [#T 无假设])
$ hat(C)_L (k) = integral (dif nu_k (lambda)) / lambda , quad lambda^(-1) <= 1 slash lambda_min (k) quad "on" "supp" nu_k $

#c("那么", note: [#T])
$ hat(C)_L (k) >= (norm(nu_k)) / (lambda_min (k)) quad arrow.r.double quad lambda_min (k) <= norm(nu_k) slash hat(C)_L (k) $

#c("所以", note: [#T 而 (A1) 要下界])
$ lambda_min (k) <= (norm(nu_k)) / C <= (norm(D_k)^2) / C <= 4 slash C quad arrow.r.double quad "上界" $

#c("得到", note: [#T 方向相反])
$ chi " 类 " quad arrow.r.double.not quad ("A1") $

#c("设", note: [#T 可容许谱测度])
$ dif nu (lambda) = lambda dot 1_((0,1]) (lambda) dif lambda $

#c("那么", note: [#T 锐性])
$ hat(C) = integral_0^1 lambda^(-1) dot lambda dif lambda = 1 , quad norm(nu) = integral_0^1 lambda dif lambda = 1 slash 2 , quad lambda_min = 0 $

#c("所以", note: [#T])
$ chi (beta) -> 0 quad ( beta -> infinity , "由" integral chi = 1 slash 6 ) quad arrow.r.double quad "任何 " chi >= norm(nu) slash m_0 " 型下界必失败" $

= 为什么转移矩阵线进不了主链

#c("因为", note: [#T 两条独立])
$ I_1 (z) I_3 (z) < I_2 (z)^2 quad forall z > 0 quad ( "级数 + Gosper 证书" slash "Riccati" ) $

#c("那么", note: [#T])
$ ("R") : quad t_2 (lambda) slash t_1 (lambda) <= I_2 (beta) slash I_1 (beta) quad forall lambda >= 0 , quad beta " 固定" $

#c("因为", note: [#T V 的结构])
$ "宇称超选" , quad B = 1/2 bb(1)_4 , quad "高度阶梯约化" , quad c_1 = lambda_max (A(G)) slash 2 $

#c("那么", note: [#T])
$ V " 的效应瞬态 " quad arrow.l.r.double quad L_s = 1 , quad r(0) = r(infinity) = u := I_2 slash I_1 $

#c("所以", note: [#T 全程幂次])
$ L_s = 1 : quad "隙" = 1 - r(beta_s) >= 1 - u approx 3 slash (2 beta) $

#c("得到", note: [#M 4D 环面读数])
$ m_"lat" = O(1) quad ( beta_s -> infinity ) quad arrow.r.double quad m_"lat" slash a (beta) -> infinity $

#c("所以", note: [#T])
$ "违 (A2)" quad arrow.r.double quad "该扇区的连续极限超局部" quad ( "C1 平凡" ) $

#c("设")
$ "桥" : quad V <= 1 quad ( "单链结果延拓到多链" ) $

#c("那么", note: [#F 差 0.185])
$ r_"nonproduct" = 0.692 > 0.507 = I_2 slash I_1 quad arrow.r.double quad "桥为假" $

#c("所以", note: [#F + #T])
$ "转移矩阵线" : quad "只到 " L_s = 1 quad arrow.r.double quad "不进 " ("A1") "，不进 " ("A2") $

= Clay 的等价形式

#thmbox[定理 1][
  $ "Clay" quad arrow.l.r.double quad ("A3") and ("SCALE") , $
  其中 $("A3")$ 与 $("SCALE")$ 如上。“$arrow.l.double$”由 $("A1") and ("A2") and ("A3") arrow.r.double "Clay"$ 给出；
  “$arrow.r.double$”由 Clay $arrow.r.double$ 沿子列 $m_"lat" slash a -> m_"phys" in (0, infinity)$ 给出。
]

#v(0.15em)
#text(size: 9.4pt)[
  *已证。* 测度与 OS 重构、隙 $=$ 谱隙、四条蕴含 $("A1") and ("A3") arrow.r.double ("C2")$、
  $("A2") and ("A3") arrow.r.double$ 非平凡、以及等价性；另有两段反证（方向定理与锐性例；桥的反例）。

  *未证。* $("A1")$、$("A2")$、$("A3")$ 三条。#text(weight: "bold")[Clay 未证]。
]

= 五环与五洞

#figure(
  table(
    columns: (auto, 1fr, auto, auto), stroke: 0.4pt + rgb("#c9d0d8"),
    align: (left, left, center, left), inset: 5pt,
    table.header([*环*], [*公式*], [*判*], [*洞*]),
    [R1], [$mu_(beta,L) -> (cal(H), T = e^(-H)), quad H >= 0$], T, [—],
    [R2], [$m_"lat" = -ln lambda_1 (T) = inf "spec" H without {0}$], T, [—],
    [R3], [$exists c_1 forall beta forall L : quad m_"lat" slash a >= c_1$], U, [D1：方向定理封 $chi$ 类],
    [R4], [$exists c_2 forall beta forall L : quad m_"lat" slash a <= c_2$], U, [D2：相关长度上界],
    [R5], [$S_n -> S_infinity$ 且满足 OS], U, [D3：Balaban 部分],
  ),
  caption: [五环。$0^(++)$ 是最低道，故其隙只是上界（#T）；D5（模型 $arrow.r$ 格点的桥）已 #F。],
)

= 死路台账

#figure(
  table(
    columns: (2.2fr, auto, 2.3fr), stroke: 0.4pt + rgb("#c9d0d8"),
    align: (left, center, left), inset: 5pt,
    table.header([*路*], [*判*], [*数*]),
    [HS / 谱 / 求和规则（$chi$ 类）], T, [$lambda_min <= norm(nu) slash hat(C)$；$integral chi = 1 slash 6 arrow.r.double chi -> 0$],
    [单链延拓多链 $V <= 1$], F, [$r = 0.692 > 0.507$],
    [$V$ 的秩截断], F, [谱填满 $[0.81, 1.50]$],
    [特征标低阶], F, [$c_J tilde beta^(-3 slash 2)$],
    [$V$ 的微扰（物理区）], F, [半径 $approx 2.5 < beta_s$],
    [朴素 $op("Tr")_2$ / 加权 block-spin], F, [残差 $>= 0.89$],
    [$L >= 2$ 精确对角化], F, [$2.8 times 10^(11)$ 维],
    [两链截断 Peter–Weyl], F, [$0.5692 -> 0.8802 -> 0.6179$],
    [非 Abel GKS $v_2 >= 0$], U, [文献只到 Abel / Villain；数值 $+6.24 times 10^(-3)$],
  ),
)

= 限定

- $G = "SU"(2)$、$d = 4$；一般单李群未逐条重跑。
- 求和规则与 $chi^"clover" <= 1 slash 6$：一切 $beta > 0$、$L >= 2$；$("R")$ 的数值核验 $beta in [0.2, 400]$（11 个值）。
- $("Q")$ 与 $norm(Gamma) = chi$：固定取向 plaquette；clover 版 $("Q")$ 未证 #U。
- 极限次序：先 $L -> infinity$ 再 $a -> 0$；两个极限不可交换 #T。
- 外部输入：OS 重构（标准）；渐近自由（标准）；Balaban RG（部分）；Ginibre 对非 Abel *不可用*（已核文献）。
- 大 $beta_s$ 数值必须用 `ive` 归一化（第 223 轮更正）。

= 下一手

#c("设", note: [攻 D1])
$ ("P1") quad L_s = 2 : quad "gap"(beta_s) >= c dot e^(-(3 pi^2 slash 11) beta_s) quad <=> quad m_"lat" slash a >= c $

#c("设", note: [攻 D2])
$ ("P2") quad lambda_1 (T_(beta,L)) >= e^(-c_2 a (beta)) quad <=> quad xi >= 1 slash c_2 $

#c("设", note: [攻 D5；$xi$ 为物理相关长度])
$ ("P3") quad V <= 1 " 不修；改找比较对象：按 " h = sum_x j_x " 分层，或 " L_s = 2 " 直接接格点键" $

#v(0.6em)
#block(inset: 8pt, stroke: (left: 2pt + rgb("#c79a2e")), fill: rgb("#fffaf0"), width: 100%)[
  #text(size: 9pt)[
    *报告纪律。* 不以百分比度量 Clay（项目宪章 §2）。可报告的三件事：① R1–R2、四条蕴含与等价性已 #T；
    ② 缺 R3、R4、R5；③ 缺的东西带量词：$exists c_1 forall beta >= beta_0 forall L$ 的下界、
    $exists c_2 forall beta >= beta_0 forall L$ 的上界、紧性。
  ]
]

#v(1.1em)
#text(size: 8.4pt, fill: rgb("#4a525c"))[
  出处：`docs/THEOREM-machine-closed.md`、`docs/CORRECTION-ginibre-not-available.md`、
  `docs/NOGO-spectral-route-cannot-give-gap.md`、`docs/CLAY-CHAIN-MAP.md`、`docs/FINAL-STATE.md`、
  `docs/ROUND259-4d-first.md`、`docs/VERDICT-C-undecidable.md`、`docs/position-current.html`。
]




