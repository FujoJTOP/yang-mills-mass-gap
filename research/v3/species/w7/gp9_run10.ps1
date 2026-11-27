for ($k=1; $k -le 10; $k++) {
  $seed = 20261126 + $k
  $t0 = Get-Date
  "=== chunk $k/10 seed=$seed  " + (Get-Date -Format HH:mm:ss) + " ==="
  & py -3.12 gh3_probe.py --Ls 6 --Lt 20 --n 160000 --gap 4 --therm 400 --chains 8 --threads 8 --seed $seed *>> gp9_10x.log
  if ($LASTEXITCODE -ne 0) { "chunk $k FAILED exit=$LASTEXITCODE"; continue }
  $o = & py -3.12 gp9_acc.py
  $o | Add-Content gp9_10x.log
  Rename-Item gp9_block_acc.npy ("gp9_acc_{0}.npy" -f $k) -Force
  Copy-Item gp9_vac.npy ("gp9_vac_{0}.npy" -f $k) -Force
  Remove-Item gh3_probe.npz -Force -ErrorAction SilentlyContinue
  "chunk $k done in $([math]::Round(((Get-Date)-$t0).TotalSeconds,0))s"
}
"ALL 10 CHUNKS DONE"
