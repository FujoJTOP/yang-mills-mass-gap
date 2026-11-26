import io, re, glob

NEW = (r"\textbf{The residual gap is the uniform statement $\forall\beta\in[1/48,2.5]\ \forall L,a:\ G''_{k=0}/a^2\ge\delta^2>0$.} "
       r"An earlier draft named this gap $\Delta_0>0$ and identified it with $H1$, $I2$, $(E\text{-}I2)$ and $H^*$. "
       r"\textbf{That identification is withdrawn}, on three independent grounds: "
       r"(i) \code{n12}: $\Delta_0>0$ fails the four-model test on three of four models (on 2D YM $G''=0$ with $W_{\mathrm{full}}\ge N/12>0$ forces $\Delta_{\mathrm{dyn}}=-W_{\mathrm{full}}<0$; "
       r"on M6, $\Phi=0\Rightarrow W_{\mathrm{full}}=0\Rightarrow\Delta_{\mathrm{dyn}}=2\omega>0$ on a gapless model; the single plaquette gives $\Delta_{\mathrm{dyn}}=-0.0691$); "
       r"(ii) \code{n06}: $\Delta_0$ and $\Lambda_\star$ live on different carriers ($G''$ vs $\operatorname{supp}\nu$), are of different type (position vs distance), "
       r"and have different quantifier domains, so they are \textbf{incomparable}; "
       r"(iii) \code{n09}, Theorem N09-1: on M6 the only degenerate scale is the spectral bottom, hence "
       r"\textbf{any quantity passing the four-model test must be a function of $m$ itself} --- the target is unique and no proxy exists. "
       r"\textbf{The correct equivalence class is three, not five: $\{\Lambda_\star\ \text{lower bound}\equiv I2\equiv(E\text{-}I2)\}$; $H^*$ and $\Delta_0$ do not belong to it.}")

TAIL = (r"\emph{Consequently this paper does \textbf{not} propose any new proxy for the gap. Five proof classes are closed "
        r"(averaging, extremal spectrum, transport, similarity ratios, triviality); the only one not refuted is \textbf{constructive RG}, "
        r"whose status is \code{[O]} (not read through), \textbf{not} \code{[NG]} (refuted).}")

files = ["main.tex"] + sorted(glob.glob("sec0*.tex"))
pat = re.compile(r"The residual gap is the single statement \$\\Delta_0>0\$.*?five-name identification\.", re.S)

for fn in files:
    s = io.open(fn, encoding="utf-8").read()
    n0 = s
    s, k = pat.subn(lambda m: NEW, s)
    # drop the now contradictory "provisional" sentence
    s = s.replace(r" \textbf{The withdrawn-transport section (\S\ref{sec:transport-withdrawn}) stands; only the identification of the residual gap with $\Delta_0>0$ is provisional.}", "")
    s = s.replace(r"\textbf{The withdrawn-transport section (\S\ref{sec:transport-withdrawn}) stands; only the identification of the residual gap with $\Delta_0>0$ is provisional.} ", "")
    if k or s != n0:
        io.open(fn, "w", encoding="utf-8", newline="").write(s)
        print("patched", fn, "blocks:", k)
print("done")

