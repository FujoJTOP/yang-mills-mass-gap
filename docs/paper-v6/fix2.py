import io, glob, sys

files = ["main.tex"] + sorted(glob.glob("sec0*.tex"))

pairst = [
    # ---- window numbers ----
    ("in the strong-coupling window ($\\beta<0.0736$)", "in the strong-coupling window ($\\beta<0.018394$)"),
    ("\\beta<0.073576=3.53\\times(1/48)", "\\beta<0.018394=0.8829\\times(1/48)"),
    ("\\beta<0.0736=3.53\\times(1/48)", "\\beta<0.018394=0.8829\\times(1/48)"),
    ("\\beta<0.0736", "\\beta<0.018394"),
    ("0.073576", "0.018394"),
    ("0.0736", "0.018394"),
    ("3.53\\times(1/48)", "0.8829\\times(1/48)"),
    ("3.532\\times(1/48)", "0.8829\\times(1/48)"),
    # ---- claim tail: five-name identification -> under adjudication ----
    ("The residual gap is the single statement $\\Delta_0>0$, which is the same proposition as $H1$, $I2$, $(E\\text{-}I2)$ and $H^*$ under different names.",
     "The residual gap is the single statement $\\Delta_0>0$. \\textbf{Its status is under adjudication}: "
     "(i) $\\Delta_0$ sits on $G''$ while $\\Lambda_\\star$ sits on the spectral measure $\\nu$, so the identification with $I2$\\,/\\,$(E\\text{-}I2)$\\,/\\,$H^*$ is \\emph{not} established; "
     "(ii) on 2D YM $G''=0$ while $m=3g^2/8>0$, so if $\\Delta_{\\mathrm{dyn}}=a^2\\Delta_0^2$ and $W_{\\mathrm{full}}\\ge N/12>0$ then $\\Delta_0^2<0$, i.e. $\\Delta_0>0$ fails on a model that does have a gap; "
     "(iii) $\\Delta_0$ may not be defined beyond its single occurrence in \\code{app.tex:190--191}. "
     "\\textbf{We therefore do not assert the five-name identification.}"),
    # ---- Chinese five-name assertions -> provisional ----
    ("残余缺口只有\\textbf{一条}陈述 $\\Delta_0>0$（$=H1=I2=(E\\text{-}I2)=H^*$，\\codep{g03}/\\codep{m05}）。",
     "残余缺口的候选是\\textbf{一条}陈述 $\\Delta_0>0$；它与 $H1$、$I2$、$(E\\text{-}I2)$、$H^*$ 的等价\\textbf{正在裁决}（$\\Delta_0$ 在 $G''$ 上、$\\Lambda_\\star$ 在 $\\nu$ 上），本文\\textbf{不宣称}（\\S\\ref{sec:assembly}）。"),
    ("而是\\textbf{一条}陈述 $\\Delta_0>0$（$=H1=I2=(E\\text{-}I2)=H^*$）。",
     "而是候选陈述 $\\Delta_0>0$（它与 $H1$、$I2$、$(E\\text{-}I2)$、$H^*$ 的等价\\textbf{正在裁决}，本文不宣称）。"),
    ("\\textbf{残余缺口只有一条陈述} $\\Delta_0>0$，它与 $H1$、$I2$、$(E\\text{-}I2)$、$H^*$ 是同一命题的五个名字（\\code{g03}；$\\delta:=a^2\\Delta_0^2+W_{\\mathrm{full}}$，\\code{m05}）。",
     "\\textbf{残余缺口的候选陈述}是 $\\Delta_0>0$；它与 $H1$、$I2$、$(E\\text{-}I2)$、$H^*$ 的等价\\textbf{正在裁决}（三条怀疑：$G''$ vs $\\nu$、2D YM 上 $G''=0$ 而 $m>0$、\\code{app.tex:190--191} 的悬空记号），本文\\textbf{不宣称}。"),
    ("\\textbf{所以主力陈述收成一条：}$\\Delta_0>0$。它与 $H1$、$I2$、$(E\\text{-}I2)$、$H^*$ 是\\textbf{同一命题的五个名字}（\\code{g03}），",
     "\\textbf{所以主力陈述收成一条候选：}$\\Delta_0>0$。它与 $H1$、$I2$、$(E\\text{-}I2)$、$H^*$ 的等价\\textbf{正在 12 路裁决}（三条硬怀疑见 \\S\\ref{sec:kl}），本文\\textbf{不宣称}；"),
    # ---- D.3 wording keeps the new number (string already swapped above) ----
]

for fn in files:
    s = io.open(fn, encoding="utf-8").read()
    n0 = s
    for a, b in pairst:
        s = s.replace(a, b)
    if s != n0:
        io.open(fn, "w", encoding="utf-8", newline="").write(s)
        print("patched", fn)
print("done")
