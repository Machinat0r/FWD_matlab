# 官网图片归档：按仪器/卫星/模式分类；不创建日期文件夹。
$ErrorActionPreference='Stop'
$cache=Join-Path $env:TEMP 'MMS_plots_20260919'
$root=[IO.Path]::GetFullPath('C:\Users\Administrator\Documents\KH\SMILE_MMS')
$worker='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\SMILE_MMS_download_plots_20260919.ps1'
$stateFile=Join-Path $cache 'controller_status.json'
$manifest=[Collections.Generic.List[object]]::new()
foreach($kind in @('burst','quicklook')){
 $name=if($kind -eq 'burst'){'mms_burst_plots.txt.selected'}else{'mms_ql_plots.txt.selected'}
 foreach($line in [IO.File]::ReadLines((Join-Path $cache $name))){
  if($line -notmatch '^\./([a-z0-9_]+)/(\d{4})/(\d{2})/(\d{2})/([a-z0-9_]+\.png)$'){throw 'Unexpected index path'}
  $plot=$Matches[1];$date=$Matches[2]+'-'+$Matches[3]+'-'+$Matches[4];$file=$Matches[5]
  if($date -lt '2026-07-20' -or $date -gt '2026-09-19'){throw 'Date outside requested snapshot'}
  $inst=if($plot.StartsWith('all_')){'综合图'}else{($plot -split '_')[0].ToUpperInvariant()}
  $sat=if($plot -match 'mms([1-4])'){'MMS'+$Matches[1]}else{'多星综合'}
  $path=[IO.Path]::GetFullPath((Join-Path $root ($inst+'\'+$sat+'\'+$kind+'\'+$file)))
  if(-not $path.StartsWith($root+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Unsafe output path'}
  $manifest.Add([pscustomobject]@{kind=$kind;instrument=$inst;spacecraft=$sat;plot=$plot;date=$date;name=$file;path=$path})
 }
}
$manifest=@($manifest | Sort-Object path -Unique)
[IO.File]::WriteAllText((Join-Path $cache 'expected_manifest.json'),($manifest | ConvertTo-Json -Depth 4 -Compress),[Text.UTF8Encoding]::new($false))
function Write-State($phase,$pass,$valid,$missing,$bytes,$errorText=''){
 $s=[ordered]@{phase=$phase;complete=($phase -eq 'complete');pid=$PID;pass=$pass;expected=$manifest.Count;valid=$valid;missing=$missing;bytes=$bytes;updatedUTC=[datetime]::UtcNow.ToString('o');destination=$root;error=$errorText}
 [IO.File]::WriteAllText($stateFile,($s | ConvertTo-Json -Compress),[Text.UTF8Encoding]::new($false))
}
try {
 for($pass=1;$pass -le 6;$pass++){
  Write-State 'downloading' $pass 0 $manifest.Count 0
  & $worker -Mode all -Workers 4
  Write-State 'verifying' $pass 0 $manifest.Count 0
  $missing=[Collections.Generic.List[object]]::new();$valid=0;$bytes=0L;$good=[Collections.Generic.List[object]]::new()
  foreach($entry in $manifest){
   $ok=$false;$f=$null
   try{
    if([IO.File]::Exists($entry.path)){
     $f=[IO.File]::OpenRead($entry.path)
     if($f.Length -gt 32){
      $head=[byte[]]::new(8);$tail=[byte[]]::new(12)
      $f.ReadExactly($head,0,8)
      $f.Seek(-12,[IO.SeekOrigin]::End) | Out-Null
      $f.ReadExactly($tail,0,12)
      $ok=([Convert]::ToHexString($head) -eq '89504E470D0A1A0A' -and [Convert]::ToHexString($tail) -eq '0000000049454E44AE426082')
      if($ok){$valid++;$bytes+=$f.Length;$good.Add($entry)}
     }
    }
   }catch{$ok=$false}finally{if($f){$f.Dispose()}}
   if(-not $ok){$missing.Add($entry)}
  }
  [IO.File]::WriteAllText((Join-Path $cache 'missing_images.json'),(ConvertTo-Json -InputObject $missing.ToArray() -Depth 4 -Compress),[Text.UTF8Encoding]::new($false))
  if($missing.Count -eq 0){
   $byInstrument=@($good | Group-Object instrument,spacecraft,kind | ForEach-Object {[pscustomobject]@{group=$_.Name;count=$_.Count}})
   $summary=[ordered]@{completedUTC=[datetime]::UtcNow.ToString('o');expected=$manifest.Count;downloaded=$valid;missing=0;bytes=$bytes;byInstrument=$byInstrument;sourceIndexes=@('https://lasp.colorado.edu/mms/sdc/public/data/sdc/mms_burst_plots.txt','https://lasp.colorado.edu/mms/sdc/public/data/sdc/mms_ql_plots.txt');scope='2026-07-20 inclusive through website inventory retrieved 2026-09-19';classification='instrument / spacecraft / burst or quicklook; original filenames'}
   [IO.File]::WriteAllText((Join-Path $cache 'download_complete.json'),($summary | ConvertTo-Json -Depth 6),[Text.UTF8Encoding]::new($false))
   Write-State 'complete' $pass $valid 0 $bytes
   exit 0
  }
  Write-State 'retrying_missing' $pass $valid $missing.Count $bytes
  for($wait=0;$wait -lt 10;$wait++){Start-Sleep -Seconds 30}
 }
 Write-State 'needs_attention' 6 $valid $missing.Count $bytes 'Some indexed images remain unavailable after six passes; see missing_images.json and all_download.jsonl.'
 exit 2
}catch{
 Write-State 'failed' 0 0 $manifest.Count 0 $_.Exception.Message
 throw
}
