$ErrorActionPreference = 'Stop'
function S4([double]$a) { $e = [math]::Exp(-$a); return [math]::Pow(1.0 / (1.0 - $e), 4) - 1.0 }
function Fsum([double]$w, [int]$N) {
  $s = 0.0
  for ($i = 2; $i -le $N; $i++) {
    $s = $s + [math]::Pow([double]$i, [double]($i - 2)) * [math]::Pow($w, [double]($i - 1))
  }
  return $s
}
"--- S4(alpha) = (1/(1-e^-a))^4 - 1 ---"
foreach ($a in @(4.0, 3.75, 3.5, 3.25, 3.0, 2.5, 2.328)) {
  "alpha=" + $a.ToString("F3") + "  S4=" + (S4 $a).ToString("F10")
}
$w = S4 4.0
"w(A=1,alpha=4) = " + $w.ToString("F10")
"1/e = " + (1.0 / [math]::Exp(1.0)).ToString("F10") + "   A_conv = " + ((1.0 / [math]::Exp(1.0)) / $w).ToString("F6")
"--- partial sums at w = S4(4), A=1 ---"
$cum = 0.0
for ($i = 2; $i -le 14; $i++) {
  $tr = [math]::Pow([double]$i, [double]($i - 2))
  $t = $tr * [math]::Pow($w, [double]($i - 1))
  $cum = $cum + $t
  "n=" + $i.ToString() + "  trees=" + $tr.ToString("F0") + "  term=" + $t.ToString("E4") + "  cum=" + $cum.ToString("F6")
}
"F(w,60) = " + (Fsum $w 60).ToString("F8")
"F(w,200) = " + (Fsum $w 200).ToString("F8")
$G = 0.90537044375382710143
"gate G = " + $G.ToString("F10")
"margin G/F = " + ($G / (Fsum $w 200)).ToString("F3")
# A threshold: F(A*w) = G
$lo = 0.000001; $hi = 6.0
for ($k = 0; $k -lt 120; $k++) {
  $m = 0.5 * ($lo + $hi)
  if ((Fsum ($m * $w) 120) -lt $G) { $lo = $m } else { $hi = $m }
}
"A_thr (F(A*w)=G) = " + $lo.ToString("F6")
# alpha threshold at A=1
$lo2 = 3.0; $hi2 = 8.0
for ($k = 0; $k -lt 120; $k++) {
  $m = 0.5 * ($lo2 + $hi2)
  if ((Fsum (S4 $m) 120) -gt $G) { $lo2 = $m } else { $hi2 = $m }
}
"alpha_min = " + $hi2.ToString("F6") + "   S4(alpha_min)=" + (S4 $hi2).ToString("F6") + "   1/alpha_min=" + (1.0 / $hi2).ToString("F6")
"n=2 only: A*S4(4) = " + $w.ToString("F8") + "   G/w = " + ($G / $w).ToString("F4")
# per-n shares at A=1
"--- share of n=2,3,4 in total ---"
$t2 = $w; $t3 = 3.0 * [math]::Pow($w, 2); $t4 = 16.0 * [math]::Pow($w, 3)
"t2=" + $t2.ToString("F6") + " t3=" + $t3.ToString("F6") + " t4=" + $t4.ToString("F6") + " tot=" + (Fsum $w 200).ToString("F6")
