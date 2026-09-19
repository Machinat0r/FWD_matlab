# 只归档官方源文件与来源元数据；不执行密度反演或版本合并。
$ErrorActionPreference='Stop'
$base='Z:\SPART-WORK\Data\Voyager\voyager1\pws\derived\electron_density\native'
$check=Join-Path $base 'source_check_20260915'
$official='https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/'
$annex='https://pds-ppi.igpp.ucla.edu/annex/voyager-pws-vlism-density/'
New-Item -ItemType Directory -Force -Path $check | Out-Null
$items=@(
    @{URL=$official;Name='bundle_index.html'},
    @{URL=$official+'data/';Name='data_index.html'},
    @{URL=$official+'document/';Name='document_index.html'},
    @{URL=$official+'browse/';Name='browse_index.html'},
    @{URL=$official+'data/vg1-vlism-density-2012-2025.csv';Name='vg1-vlism-density-2012-2025.csv'},
    @{URL=$official+'data/vg1-vlism-density-2012-2025.lblx';Name='data/vg1-vlism-density-2012-2025.lblx'},
    @{URL=$official+'bundle_voyager-pws-vlism-density.lblx';Name='bundle_voyager-pws-vlism-density.lblx'},
    @{URL=$official+'document/errata.md';Name='document/errata.md'},
    @{URL=$official+'document/errata.lblx';Name='document/errata.lblx'},
    @{URL='https://space.physics.uiowa.edu/voyager/data/';Name='iowa_data_index.html'},
    @{URL='https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager1/wave_spectra_pws/';Name='spdf_pws_index.html'},
    @{URL='https://cdaweb.gsfc.nasa.gov/misc/NotesV.html';Name='cdaweb_notes_v.html'}
)
$manifest=[Collections.Generic.List[object]]::new()
foreach($item in $items) {
    $path=Join-Path $check $item.Name
    New-Item -ItemType Directory -Force -Path (Split-Path $path) | Out-Null
    # Existing snapshots are immutable; a refresh uses a timestamped sibling.
    $temp=$path+'.check-'+[guid]::NewGuid().ToString('N')
    $response=Invoke-WebRequest -UseBasicParsing -Uri $item.URL -OutFile $temp -PassThru -TimeoutSec 60
    $newHash=(Get-FileHash -LiteralPath $temp -Algorithm SHA256).Hash
    if(Test-Path -LiteralPath $path) {
        $oldHash=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
        if($oldHash -eq $newHash) {
            Remove-Item -LiteralPath $temp
        } else {
            $path=$path+'.snapshot-'+[DateTime]::UtcNow.ToString('yyyyMMddTHHmmssZ')
            Move-Item -LiteralPath $temp -Destination $path
        }
    } else { Move-Item -LiteralPath $temp -Destination $path }
    $manifest.Add([pscustomobject]@{URL=$item.URL;File=$path;SHA256=$newHash;
        Bytes=(Get-Item -LiteralPath $path).Length;CheckedUTC=[DateTime]::UtcNow.ToString('o');
        HTTPStatus=[int]$response.StatusCode})
}
foreach($name in @('data/vg1-vlism-density-2012-2025.csv','data/vg1-vlism-density-2012-2025.lblx','bundle_voyager-pws-vlism-density.lblx')) {
    $path=Join-Path (Join-Path $base 'PDS_annex_release_1') $name
    if(-not (Test-Path -LiteralPath $path)) {
        New-Item -ItemType Directory -Force -Path (Split-Path $path) | Out-Null
        Invoke-WebRequest -UseBasicParsing -Uri ($annex+$name) -OutFile $path -TimeoutSec 60
    }
    $manifest.Add([pscustomobject]@{URL=$annex+$name;File=$path;SHA256=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash;
        Bytes=(Get-Item -LiteralPath $path).Length;CheckedUTC=[DateTime]::UtcNow.ToString('o');HTTPStatus=$null})
}
$current=Join-Path $base 'PDS_release_20260910/data/vg1-vlism-density-2012-2025.csv'
$retrieved=@($manifest | Where-Object URL -eq ($official+'data/vg1-vlism-density-2012-2025.csv'))[0]
$matches=(Get-FileHash -LiteralPath $current -Algorithm SHA256).Hash -eq $retrieved.SHA256
$result=[ordered]@{CheckedUTC=[DateTime]::UtcNow.ToString('o');CurrentArchiveMatchesWebsite=$matches;
    CurrentCSV=$current;CurrentMD5=(Get-FileHash -LiteralPath $current -Algorithm MD5).Hash;
    Files=$manifest.ToArray()}
$out=Join-Path $check 'source_verification_manifest.json'
[IO.File]::WriteAllText($out,($result | ConvertTo-Json -Depth 8),[Text.UTF8Encoding]::new($false))
if(-not $matches) { throw 'Official current CSV changed. New source retained; review before plotting.' }
Write-Host "Official density verification passed. Source manifest: $out"
