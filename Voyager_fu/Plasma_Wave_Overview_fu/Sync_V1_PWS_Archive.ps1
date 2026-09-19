param([switch]$ResumeFromInventory,[string]$StartDate='19900101',[string]$StopDate='20250630',[int]$Workers=8)
# Official source transfer only; scientific computation is MATLAB.
$ErrorActionPreference='Stop'
$root='Z:\SPART-WORK\Data\Voyager\voyager1\pws'
$metadata=Join-Path $root 'source_verification\overview_1990_20250630'
$waveURL='https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager1/wave_spectra_pws/spectrum_analyzer_cdf/'
$iowaURL='https://space.physics.uiowa.edu/plasma-wave/voyager/data/voyager-1-pws-sa/data/'
$densityURL='https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/'
$densityRoot=Join-Path $root 'derived\electron_density\native\PDS_release_20260910'
if (-not (Test-Path -LiteralPath $root)) { throw 'Authorized PWS database unavailable.' }
New-Item -ItemType Directory -Force -Path $metadata | Out-Null
$utf8=[Text.UTF8Encoding]::new($false)
function Write-Json($Path,$Value) { [IO.File]::WriteAllText($Path,($Value | ConvertTo-Json -Depth 10),$utf8) }
function Read-Index($Html,$Prefix) {
    $selected=@{}
    foreach ($m in [regex]::Matches($Html,('href="('+[regex]::Escape($Prefix)+'(\d{8})_v([\d.]+)\.cdf)"'))) {
        $day=$m.Groups[2].Value
        if ($day -lt $StartDate -or $day -gt $StopDate) { continue }
        $version=[version]$m.Groups[3].Value
        if (-not $selected.ContainsKey($day) -or $version -gt $selected[$day].Version) {
            $selected[$day]=[pscustomobject]@{Name=$m.Groups[1].Value;Version=$version}
        }
    }
    return $selected
}
$densityNames=@('bundle_voyager-pws-vlism-density.lblx','readme.md','data/collection_data.csv',
    'data/collection_data.lblx','data/vg1-vlism-density-2012-2025.csv','data/vg1-vlism-density-2012-2025.lblx')
if (-not $ResumeFromInventory) {
    foreach ($name in $densityNames) {
        $target=Join-Path $densityRoot $name
        New-Item -ItemType Directory -Force -Path (Split-Path $target) | Out-Null
        $temporary=$target+'.part-pwsh-'+[guid]::NewGuid().ToString('N')
        Invoke-WebRequest -Uri ($densityURL+$name) -OutFile $temporary -TimeoutSec 60
        if (Test-Path -LiteralPath $target) {
            $oldHash=(Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
            if ($oldHash -ne (Get-FileHash -LiteralPath $temporary -Algorithm SHA256).Hash) {
                $previous=Join-Path (Split-Path $target) ('previous_versions\'+$oldHash)
                New-Item -ItemType Directory -Force -Path $previous | Out-Null
                Copy-Item -LiteralPath $target -Destination (Join-Path $previous (Split-Path $target -Leaf)) -Force
            }
        }
        Move-Item -LiteralPath $temporary -Destination $target -Force
    }
    $inventory=[Collections.Generic.List[object]]::new();$comparisons=[Collections.Generic.List[object]]::new()
    foreach ($year in ([int]$StartDate.Substring(0,4))..([int]$StopDate.Substring(0,4))) {
        $spdf=(Invoke-WebRequest -Uri ($waveURL+$year+'/') -TimeoutSec 60).Content
        $iowa=(Invoke-WebRequest -Uri ($iowaURL+$year+'/') -TimeoutSec 60).Content
        [IO.File]::WriteAllText((Join-Path $metadata ('SPDF_'+$year+'_index.html')),$spdf,$utf8)
        [IO.File]::WriteAllText((Join-Path $metadata ('Iowa_'+$year+'_index.html')),$iowa,$utf8)
        $s=Read-Index $spdf 'vg1_pws_lr_';$i=Read-Index $iowa 'vg1pws_lr_'
        $extraIowa=@($i.Keys | Where-Object { -not $s.ContainsKey($_) })
        $extraSPDF=@($s.Keys | Where-Object { -not $i.ContainsKey($_) })
        $versions=@($i.Keys | Where-Object { $s.ContainsKey($_) -and $i[$_].Version -ne $s[$_].Version })
        $comparisons.Add([pscustomobject]@{Year=$year;IowaDays=$i.Count;SPDFDays=$s.Count;
            AdditionalIowaDates=$extraIowa;AdditionalSPDFDates=$extraSPDF;VersionDifferences=$versions})
        foreach ($date in ($s.Keys | Sort-Object)) {
            $name=$s[$date].Name
            $inventory.Add([pscustomobject]@{URL=$waveURL+$year+'/'+$name;File=(Join-Path $root "calibrated\spectrum_analyzer\native\$year\$name")})
        }
        Write-Host "Official inventory $year : $($s.Count) CDFs"
    }
    $match=@($comparisons | Where-Object { $_.AdditionalIowaDates.Count -or $_.AdditionalSPDFDates.Count -or $_.VersionDifferences.Count }).Count -eq 0
    Write-Json (Join-Path $metadata 'SPDF_Iowa_comparison.json') ([ordered]@{
        CheckedUTC=[DateTime]::UtcNow.ToString('o');Start=$StartDate;StopInclusive=$StopDate;CompleteMatch=$match;Years=$comparisons.ToArray()})
    Write-Json (Join-Path $metadata 'selected_official_files.json') $inventory.ToArray()
    if (-not $match) { throw 'Official mirrors differ; inspect comparison.' }
}
$inventory=@(Get-Content -LiteralPath (Join-Path $metadata 'selected_official_files.json') -Raw | ConvertFrom-Json)
$comparison=Get-Content -LiteralPath (Join-Path $metadata 'SPDF_Iowa_comparison.json') -Raw | ConvertFrom-Json
if (-not $comparison.CompleteMatch) { throw 'Official mirror comparison did not pass.' }
foreach ($folder in ($inventory.File | ForEach-Object { Split-Path $_ } | Sort-Object -Unique)) {
    New-Item -ItemType Directory -Force -Path $folder | Out-Null
}
$bag=[Collections.Concurrent.ConcurrentBag[object]]::new();$total=$inventory.Count
Write-Host "Checking $total official CDFs; downloading only missing files."
$inventory | ForEach-Object -Parallel {
    $ErrorActionPreference='Stop';$item=$_;$records=$using:bag
    $allowedRoot=[IO.Path]::GetFullPath($using:root).TrimEnd('\')+'\'
    $file=[IO.Path]::GetFullPath([string]$item.File);$temporary=$null
    try {
        if (-not $file.StartsWith($allowedRoot,[StringComparison]::OrdinalIgnoreCase)) { throw 'Target outside authorized archive.' }
        if (([uri]$item.URL).Host -ne 'spdf.gsfc.nasa.gov') { throw 'Unexpected host.' }
        $status='existing'
        if (-not [IO.File]::Exists($file)) {
            $temporary=$file+'.part-pwsh-'+[guid]::NewGuid().ToString('N')
            for ($attempt=1;$attempt -le 4;$attempt++) {
                try {
                    $response=Invoke-WebRequest -Uri $item.URL -OutFile $temporary -PassThru -TimeoutSec 30
                    $bytes=[IO.File]::ReadAllBytes($temporary);$length=$response.Headers['Content-Length']
                    if ($length -and $bytes.Length -ne [long]@($length)[0]) { throw 'HTTP content length mismatch.' }
                    if ($bytes.Length -lt 8 -or [BitConverter]::ToString($bytes,0,4) -notin @('CD-F3-00-01','CD-F2-60-02')) { throw 'Response is not CDF.' }
                    [IO.File]::Move($temporary,$file);$temporary=$null;$status='downloaded';break
                } catch {
                    if ($attempt -eq 4) { throw }
                    Start-Sleep -Seconds ([Math]::Pow(2,$attempt-1))
                }
            }
        } else { $bytes=[IO.File]::ReadAllBytes($file) }
        if ($bytes.Length -lt 8 -or [BitConverter]::ToString($bytes,0,4) -notin @('CD-F3-00-01','CD-F2-60-02')) { throw 'Existing file is not CDF.' }
        $hash=[Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($bytes)).ToLowerInvariant()
        $records.Add([pscustomobject]@{URL=$item.URL;File=$file;Bytes=$bytes.Length;SHA256=$hash;Status=$status;CheckedUTC=[DateTime]::UtcNow.ToString('o')})
    } catch {
        $records.Add([pscustomobject]@{URL=$item.URL;File=$file;Status='failed';Error=$_.Exception.Message})
        if ($temporary -and [IO.File]::Exists($temporary)) { Remove-Item -LiteralPath $temporary }
    }
    $done=$records.Count
    if ($done % 250 -eq 0 -or $done -eq $using:total) { Write-Host "CDF checked $done/$using:total" }
} -ThrottleLimit $Workers
$records=@($bag.ToArray() | Sort-Object File);$failed=@($records | Where-Object Status -eq 'failed')
foreach ($name in $densityNames) {
    $file=Join-Path $densityRoot $name
    $records+=[pscustomobject]@{URL=$densityURL+$name;File=$file;Bytes=(Get-Item -LiteralPath $file).Length;
        SHA256=(Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant();Status='existing';CheckedUTC=[DateTime]::UtcNow.ToString('o')}
}
$manifest=Join-Path $metadata 'download_manifest.jsonl'
[IO.File]::WriteAllLines($manifest,@($records | ForEach-Object { $_ | ConvertTo-Json -Compress -Depth 4 }),$utf8)
$summary=[ordered]@{Start=$StartDate;StopInclusive=$StopDate;CreatedUTC=[DateTime]::UtcNow.ToString('o');
    SelectedCDFs=$total;Records=$records.Count;Failures=$failed;Downloaded=@($records | Where-Object Status -eq 'downloaded').Count;
    Bytes=($records | Measure-Object Bytes -Sum).Sum;Manifest=$manifest;Transport='PowerShell HTTPS, original bytes'}
Write-Json (Join-Path $metadata 'download_summary.json') $summary
$summary | ConvertTo-Json -Depth 6
if ($failed.Count) { throw 'Some CDFs could not be downloaded; inspect source report.' }
