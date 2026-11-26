# run_all.ps1 -- sequential repro sweep for research/v3/species/w7
# one process at a time, 180 s timeout per script, stdout/stderr captured under repro/out
# NOTE: audit scripts are run with cwd = species/w7 (their ../../w7/... paths need depth 3 under research/v3)
$ErrorActionPreference = 'Continue'
$root = 'D:\Dev\math\research\v3\species\w7'
$out  = Join-Path $root 'repro\out'
$snap = Join-Path $root 'repro\snap'
$tsv  = Join-Path $root 'repro\results.tsv'
New-Item -ItemType Directory -Force -Path $out, $snap | Out-Null

$rows = @()
$w7 = Get-ChildItem $root -Filter *.py | Sort-Object Name
$au = Get-ChildItem (Join-Path $root 'audit') -Filter *.py | Sort-Object Name

$list = @()
foreach ($f in @($w7) + @($au)) {
    $txt = Get-Content $f.FullName -Raw
    $hasLog  = [bool]($txt -match '\.log')
    $hasJson = [bool]($txt -match '\.json')
    $rank = if ($hasLog -and $hasJson) { 0 } elseif ($hasLog) { 1 } elseif ($hasJson) { 2 } else { 3 }
    $isAudit = ($f.DirectoryName -like '*\audit')
    $list += [pscustomobject]@{
        File  = $f.Name
        Dir   = if ($isAudit) { $root } else { $root }   # audit: depth-3 cwd (see NOTE)
        Src   = $f.FullName
        Group = if ($isAudit) { 'audit' } else { 'w7' }
        Rank  = $rank
    }
}
$list = $list | Sort-Object Rank, Group, File

foreach ($s in $list) {
    $stem = $s.Group + '__' + [IO.Path]::GetFileNameWithoutExtension($s.File)
    $o = Join-Path $out ($stem + '.out.txt')
    $e = Join-Path $out ($stem + '.err.txt')
    $sw = [Diagnostics.Stopwatch]::StartNew()
    $p = Start-Process -FilePath 'py' -ArgumentList '-3.12', '-I', $s.Src `
         -WorkingDirectory $s.Dir -RedirectStandardOutput $o -RedirectStandardError $e `
         -NoNewWindow -PassThru
    $finished = $p.WaitForExit(180000)
    if (-not $finished) {
        & taskkill /F /T /PID $p.Id 2>&1 | Out-Null
        try { $p.WaitForExit(15000) | Out-Null } catch {}
    }
    $sw.Stop()
    $code = if ($finished) { $p.ExitCode } else { 'TIMEOUT' }
    # snapshot whatever the script wrote into its own directory
    $sd = Join-Path $snap $stem
    New-Item -ItemType Directory -Force -Path $sd | Out-Null
    $wrote = @()
    foreach ($f2 in (Get-ChildItem -Path (Join-Path $s.Dir '*') -Include *.log, *.json -File)) {
        try { Copy-Item $f2.FullName (Join-Path $sd $f2.Name) -Force; $wrote += $f2.Name } catch {}
    }
    $rows += [pscustomobject]@{
        script   = $s.File
        group    = $s.Group
        exit     = $code
        secs     = [math]::Round($sw.Elapsed.TotalSeconds, 1)
        wrote    = ($wrote -join ';')
        outbytes = if (Test-Path $o) { (Get-Item $o).Length } else { 0 }
        errbytes = if (Test-Path $e) { (Get-Item $e).Length } else { 0 }
    }
    Write-Host ("[{0,3}] {1,-26} exit={2,-7} {3,6}s  out={4}B err={5}B" -f $rows.Count, $s.File, $code, [math]::Round($sw.Elapsed.TotalSeconds,1), $rows[-1].outbytes, $rows[-1].errbytes)
    $rows | Export-Csv -Path $tsv -NoTypeInformation -Delimiter "`t" -Encoding UTF8
}
Write-Host ("DONE " + $rows.Count + " scripts -> " + $tsv)
