$ErrorActionPreference = 'Stop'
# t = e^{-alpha}
function S4full([double]$a) { $t = [math]::Exp(-$a); return [math]::Pow((1.0 + $t) / (1.0 - $t), 4) - 1.0 }   # sum over Z^4 \ {0}
function S4half([double]$a) { $t = [math]::Exp(-$a); return [math]::Pow(1.0 / (1.0 - $t), 4) - 1.0 }        # sum over N^4 \ {0}
function Fsum([double]$w, [int]$N) { $s = 0.0; for ($i = 2; $i -le $N; $i++) { $s = $s + [math]::Pow([double]$i,[double]($i-2))*[math]::Pow($w,[double]($i-1)) }; return $s }
function Ffact([double]$w, [int]$N) {
  $s = 0.0; $f = 1.0
  for ($i = 2; $i -le $N; $i++) { $f = $f * $i; $s = $s + [math]::Pow([double]$i,[double]($i-2))*[math]::Pow($w,[double]($i-1))/$f }
  return $s
}
$G = 0.90537044375382710143
"G (Gaussian-side gate) = " + $G.ToString("F10")
foreach ($a in @(4.0, 3.5, 3.0, 2.52, 2.5, 2.328)) {
  "alpha=" + $a.ToString("F3") + "  S4_full=" + (S4full $a).ToString("F8") + "  S4_half=" + (S4half $a).ToString("F8")
}
$wf = S4full 4.0; $wh = S4half 4.0
"--- frozen alpha=4, A=1: n=2 term ---"
"full Z^4 : w=" + $wf.ToString("F8") + "  gate/w=" + ($G/$wf).ToString("F4")
"half N^4 : w=" + $wh.ToString("F8") + "  gate/w=" + ($G/$wh).ToString("F4")
"--- Cayley WITHOUT 1/n! : partial sums (half-line w) ---"
$cum = 0.0
for ($i = 2; $i -le 16; $i++) { $t = [math]::Pow([double]$i,[double]($i-2))*[math]::Pow($wh,[double]($i-1)); $cum = $cum + $t; "n=" + $i + " term=" + $t.ToString("E4") + " cum=" + $cum.ToString("F4") }
"Cayley cum at n=16 = " + $cum.ToString("E4") + "   (F(60)=" + (Fsum $wh 60).ToString("E3") + ")"
"--- WITH 1/n! (KP/set counting) ---"
$cum2 = 0.0
for ($i = 2; $i -le 12; $i++) { $f=[double]1.0; for($k=1;$k -le $i;$k++){$f=$f*$k}; $t = [math]::Pow([double]$i,[double]($i-2))*[math]::Pow($wh,[double]($i-1))/$f; $cum2 = $cum2 + $t; "n=" + $i + " term=" + $t.ToString("E4") + " cum=" + $cum2.ToString("F6") }
"Ffact(w_half,40) = " + (Ffact $wh 40).ToString("F8")
"Ffact(w_full,40) = " + (Ffact $wf 40).ToString("F8")
"--- alpha thresholds for the n=2 term = G, A=1 ---"
$lo=0.5; $hi=9.0
for ($k=0;$k -lt 150;$k++){ $m=0.5*($lo+$hi); if ((S4full $m) -gt $G) { $lo=$m } else { $hi=$m } }
"alpha_min(full) = " + $hi.ToString("F6") + "   1/alpha_min = " + (1.0/$hi).ToString("F6")
$lo=0.5; $hi=9.0
for ($k=0;$k -lt 150;$k++){ $m=0.5*($lo+$hi); if ((S4half $m) -gt $G) { $lo=$m } else { $hi=$m } }
"alpha_min(half) = " + $hi.ToString("F6") + "   1/alpha_min = " + (1.0/$hi).ToString("F6")
"--- gate rewritten: Bc* <= 1/alpha_min ---"
"full: Bc* <= " + (1.0/$hi).ToString("F6") + " ; B=4,c*=1/16 -> Bc*=" + (4.0/16.0).ToString("F6")
