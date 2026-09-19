param([int]$StartYear=1990,[int]$StopYear=2024)
# Official inventory, original CDF transfer and source documentation only.
$ErrorActionPreference='Stop'
$root='Z:\SPART-WORK\Data\Voyager\voyager2'
$base='https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/'
$metadata=Join-Path $root 'source_verification\plasma_overview_1990_2024'
New-Item -ItemType Directory -Force -Path $metadata | Out-Null
$utf8=[Text.UTF8Encoding]::new($false)
$listings=[Collections.Generic.List[object]]::new()
$files=[Collections.Generic.List[object]]::new()
function Get-SourceText($url,$name) {
    $r=Invoke-WebRequest -UseBasicParsing -Uri $url -TimeoutSec 60
    $file=Join-Path $metadata $name
    [IO.File]::WriteAllText($file,$r.Content,$utf8)
    $listings.Add([pscustomobject]@{URL=$url;File=$file;SHA256=(Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash;CheckedUTC=[DateTime]::UtcNow.ToString('o')})
    return $r.Content
}
function Select-CDF($html,$prefix) {
    $selected=@{}
    foreach($m in [regex]::Matches($html,('href="('+[regex]::Escape($prefix)+'(\d{8})_v(\d+)\.cdf)"'))) {
        $key=$m.Groups[2].Value;$v=[int]$m.Groups[3].Value
        $year=[int]$key.Substring(0,4)
        if($year -lt $StartYear -or $year -gt $StopYear){ continue }
        if(-not $selected.ContainsKey($key) -or $v -gt $selected[$key].Version) {
            $selected[$key]=[pscustomobject]@{Name=$m.Groups[1].Value;Date=$key;Year=$year;Version=$v}
        }
    }
    return @($selected.Values | Sort-Object Date)
}
function Ensure-CDF($url,$file,$product) {
    $existed=Test-Path -LiteralPath $file
    if(-not $existed) {
        New-Item -ItemType Directory -Force -Path (Split-Path $file) | Out-Null
        $temp=$file+'.part-'+[guid]::NewGuid().ToString('N')
        Invoke-WebRequest -UseBasicParsing -Uri $url -OutFile $temp -TimeoutSec 120
        $stream=[IO.File]::OpenRead($temp);$header=New-Object byte[] 4
        try { [void]$stream.Read($header,0,4) } finally { $stream.Dispose() }
        if([BitConverter]::ToString($header) -notin @('CD-F3-00-01','CD-F2-60-02')) { throw "Invalid downloaded CDF: $url" }
        Move-Item -LiteralPath $temp -Destination $file
        Write-Host "Downloaded original CDF: $file"
    }
    $stream=[IO.File]::OpenRead($file);$header=New-Object byte[] 4
    try { [void]$stream.Read($header,0,4) } finally { $stream.Dispose() }
    if([BitConverter]::ToString($header) -notin @('CD-F3-00-01','CD-F2-60-02')) { throw "Invalid archived CDF: $file" }
    $files.Add([pscustomobject]@{Product=$product;URL=$url;File=$file;ExistedBefore=$existed;
        Bytes=(Get-Item -LiteralPath $file).Length;SHA256=(Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash})
}
$cohoRoot=Get-SourceText ($base+'coho1hr_magplasma/') 'COHO_root_index.html'
$cohoYears=@([regex]::Matches($cohoRoot,'href="(\d{4})/"') | ForEach-Object {[int]$_.Groups[1].Value})
foreach($year in @($cohoYears | Where-Object {$_ -ge $StartYear -and $_ -le $StopYear})) {
    $url=$base+'coho1hr_magplasma/'+$year+'/'
    $html=Get-SourceText $url ('COHO_'+$year+'_index.html')
    foreach($entry in (Select-CDF $html 'voyager2_coho1hr_merged_mag_plasma_')) {
        $folder=Join-Path $root ('coho\1hr\l2\merged_mag_plasma\'+$year+'\'+$entry.Date.Substring(4,2))
        Ensure-CDF ($url+$entry.Name) (Join-Path $folder $entry.Name) 'COHO'
    }
    Write-Host "COHO source inventory: $year"
}
$plsRoot=Get-SourceText ($base+'plasma_cdaweb/hires_plasma/') 'PLS_root_index.html'
$plsYears=@([regex]::Matches($plsRoot,'href="(\d{4})/"') | ForEach-Object {[int]$_.Groups[1].Value})
foreach($year in @($plsYears | Where-Object {$_ -ge $StartYear -and $_ -le $StopYear})) {
    $url=$base+'plasma_cdaweb/hires_plasma/'+$year+'/'
    $html=Get-SourceText $url ('PLS_'+$year+'_index.html')
    foreach($entry in (Select-CDF $html 'voyager2_pls_hires_plasma_data_')) {
        $folder=Join-Path $root ('pls\hires\l3\solar_wind\'+$year)
        Ensure-CDF ($url+$entry.Name) (Join-Path $folder $entry.Name) 'PLS_solar_wind'
    }
}
$url=$base+'plasma_cdaweb/hires_plasma/heliosheath/'
$html=Get-SourceText $url 'PLS_heliosheath_index.html'
foreach($entry in (Select-CDF $html 'voyager2_pls_hires_plasma_data_hsh_')) {
    Ensure-CDF ($url+$entry.Name) (Join-Path $root ('pls\hires\l3\heliosheath\'+$entry.Year+'\'+$entry.Name)) 'PLS_heliosheath'
}
$url=$base+'magnetic_fields_cdaweb/vim_48secmag/'
$html=Get-SourceText $url 'MAG48s_index.html'
foreach($entry in (Select-CDF $html 'voyager2_48s_mag-vim_')) {
    Ensure-CDF ($url+$entry.Name) (Join-Path $root ('mag\48s\reviewed_vim\'+$entry.Year+'\'+$entry.Name)) 'MAG48s'
}
$extra=@(
    'plasma/hires/','plasma/hires/heliosheath/','plasma/hour/','plasma/daily/',
    'plasma_cdaweb/ions/l/','plasma_cdaweb/ions/m/',
    'plasma/hires/vy2pla_hires_fmt.txt','plasma/hires/heliosheath/00readme.txt',
    'plasma/hour/vy2pla_1h_fmt.txt')
foreach($part in $extra) {
    $name=($part -replace '/','_').TrimEnd('_')
    if(-not $name.EndsWith('.txt')) { $name+='_index.html' }
    [void](Get-SourceText ($base+$part) $name)
}
$result=[ordered]@{CheckedUTC=[DateTime]::UtcNow.ToString('o');StartYear=$StartYear;StopYear=$StopYear;
    LatestOfficialCOHOYear=($cohoYears | Measure-Object -Maximum).Maximum;
    DownloadedCDFs=@($files | Where-Object {-not $_.ExistedBefore}).Count;
    SelectedCDFs=$files.Count;Files=$files.ToArray();Listings=$listings.ToArray();
    Policy='Original CDFs retained. COHO and reviewed MAG used for overview; independent PLS files checked for archive completeness. No additional fitting or synthetic extension.'}
[IO.File]::WriteAllText((Join-Path $metadata 'official_source_inventory.json'),($result | ConvertTo-Json -Depth 8),$utf8)
Write-Host "Source inventory complete: $($files.Count) CDFs; $($result.DownloadedCDFs) new."
