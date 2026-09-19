# Original official files only; no scientific conversion.
$ErrorActionPreference='Stop'
$base='https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/'
$root='Z:\SPART-WORK\Data\Voyager\voyager2\pws\derived\electron_density\native\PDS_release_20260910'
$verify='Z:\SPART-WORK\Data\Voyager\voyager2\source_verification\plasma_overview_density_extension'
New-Item -ItemType Directory -Force -Path $root,$verify | Out-Null
$names=@('readme.md','bundle_voyager-pws-vlism-density.lblx','data/collection_data.csv','data/collection_data.lblx',
'data/vg2-vlism-density-2019-2025.csv','data/vg2-vlism-density-2019-2025.lblx')
$manifest=@()
foreach($name in $names){
    $file=Join-Path $root $name
    New-Item -ItemType Directory -Force -Path (Split-Path $file) | Out-Null
    $temp=$file+'.download'
    Invoke-WebRequest -Uri ($base+$name) -OutFile $temp -UseBasicParsing -TimeoutSec 45
    if(Test-Path -LiteralPath $file){
        $prior=(Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash
        if($prior -ne (Get-FileHash -LiteralPath $temp -Algorithm SHA256).Hash){
            $history=Join-Path $root ('previous_versions\'+$prior)
            New-Item -ItemType Directory -Force -Path $history | Out-Null
            Copy-Item -LiteralPath $file -Destination (Join-Path $history (Split-Path $file -Leaf))
        }
    }
    Copy-Item -LiteralPath $temp -Destination $file -Force
    Remove-Item -LiteralPath $temp
    $manifest+=[pscustomobject]@{URL=$base+$name;File=$file;SHA256=(Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash;Bytes=(Get-Item -LiteralPath $file).Length;CheckedUTC=[DateTime]::UtcNow.ToString('o')}
}
$pages=@(
'https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/',
'https://web.mit.edu/space/www/voyager/voyager_data/voyager_data.html',
'https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma_cdaweb/hires_plasma/heliosheath/',
'https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma/hour/',
'https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma/daily/',
'https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/wave_spectra_pws/')
for($k=0;$k -lt $pages.Count;$k++){
    $file=Join-Path $verify ('official_index_'+$k+'.html')
    Invoke-WebRequest -Uri $pages[$k] -OutFile $file -UseBasicParsing -TimeoutSec 45
    $manifest+=[pscustomobject]@{URL=$pages[$k];File=$file;SHA256=(Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash;Bytes=(Get-Item -LiteralPath $file).Length;CheckedUTC=[DateTime]::UtcNow.ToString('o')}
}
[IO.File]::WriteAllText((Join-Path $verify 'official_sources.json'),($manifest | ConvertTo-Json -Depth 4),[Text.UTF8Encoding]::new($false))
$csv=Join-Path $root 'data\vg2-vlism-density-2019-2025.csv'
$label=Join-Path $root 'data\vg2-vlism-density-2019-2025.lblx'
$rows=Import-Csv -LiteralPath $csv
$xml=[IO.File]::ReadAllText($label)
$checksum=[regex]::Match($xml,'<md5_checksum>([^<]+)').Groups[1].Value
if($checksum -ne (Get-FileHash -LiteralPath $csv -Algorithm MD5).Hash){throw 'Official MD5 does not match downloaded CSV'}
[pscustomobject]@{Rows=@($rows).Count;First=$rows[0];Last=$rows[-1];MD5=$checksum;SourceFiles=$manifest.Count} | ConvertTo-Json -Depth 4
