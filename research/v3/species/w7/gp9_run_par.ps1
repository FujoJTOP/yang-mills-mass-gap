# 32-thread parallel, gap=1 (cost/cfg ~1 sweep), 8 procs x 4 chains x 4 threads, 2 waves x 100k = 1.6M cfg
$WAVE = 2; $PROC = 8; $N = 100000; $GAP = 1
$T00 = Get-Date
for ($w = 1; $w -le $WAVE; $w++) {
  $jobs = @()
  for ($p = 1; $p -le $PROC; $p++) {
    $seed = 20261126 + 100*$w + $p
    $tag  = "x$($w)_$($p)"
    $jobs += Start-Process -PassThru -NoNewWindow -FilePath "py" `
      -ArgumentList "-3.12","gh3_probe.py","--Ls","6","--Lt","20","--n","$N","--gap","$GAP","--therm","400","--chains","4","--threads","4","--seed","$seed","--tag","$tag" `
      -RedirectStandardOutput "gp9w_$tag.log" -RedirectStandardError "gp9e_$tag.log"
  }
  "wave $w : launched $($jobs.Count) procs  " + (Get-Date -Format HH:mm:ss)
  $jobs | Wait-Process
  "wave $w : all done in $([math]::Round(((Get-Date)-$T00).TotalSeconds,0))s"
  & py -3.12 gp9_acc_par.py $w
  Remove-Item gh3_x${w}_*.npz -Force -ErrorAction SilentlyContinue
}
"ALL WAVES DONE total $([math]::Round(((Get-Date)-$T00).TotalSeconds,0))s"
