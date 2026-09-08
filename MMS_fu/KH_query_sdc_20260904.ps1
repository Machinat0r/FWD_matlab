$ErrorActionPreference='Stop'
$cacheDir='Z:\SPART-WORK\Data\MMS\derived\KH\catalog_audit_20260903\sdc_metadata'
$bundle=Get-Content -LiteralPath (Join-Path $env:TEMP 'KH_catalog_20260903\input_bundle.json') -Raw | ConvertFrom-Json
$dates=@($bundle.catalog | ForEach-Object {$_.StartUTC.Substring(0,10);$_.EndUTC.Substring(0,10)})+@('2018-11-06','2015-10-02')
$dates=@($dates | Sort-Object -Unique)
$specs=@(@('B','brst','fgm',''),@('Vi','brst','fpi','dis-moms'),@('Ve','brst','fpi','des-moms'),@('E','brst','edp','dce'),@('B','srvy','fgm',''),@('Vi','fast','fpi','dis-moms'),@('Ve','fast','fpi','des-moms'),@('E','fast','edp','dce'))
$jobs=foreach($day in $dates){foreach($spec in $specs){[pscustomobject]@{day=$day;product=$spec[0];mode=$spec[1];instrument=$spec[2];descriptor=$spec[3]}}}
$results=$jobs | ForEach-Object -Parallel {
 $j=$_;$folder=$using:cacheDir
 $key=@($j.day,$j.product,$j.mode,$j.instrument,$j.descriptor)-join '_'
 $fp=Join-Path $folder ($key+'.json')
 if(Test-Path -LiteralPath $fp){
  $existing=Get-Content -LiteralPath $fp -Raw | ConvertFrom-Json
  if(-not $existing.error){return [pscustomobject]@{key=$key;error='';cached=$true}}
 }
 $end=([datetime]::ParseExact($j.day,'yyyy-MM-dd',$null)).AddDays(1).ToString('yyyy-MM-dd')
 $u='https://lasp.colorado.edu/mms/sdc/public/files/api/v1/file_info/science?sc_id=mms1,mms2,mms3,mms4&instrument_id='+$j.instrument+'&data_rate_mode='+$j.mode+'&data_level=l2&start_date='+$j.day+'&end_date='+$end
 if($j.descriptor){$u+='&descriptor='+$j.descriptor}
 $data=$null
 for($attempt=0;$attempt -lt 2;$attempt++){
  try{
   $data=Invoke-RestMethod -Uri $u -TimeoutSec 25 -ErrorAction Stop
   $data | Add-Member -NotePropertyName url -NotePropertyValue $u -Force
   $data | Add-Member -NotePropertyName checkedUTC -NotePropertyValue ([datetime]::UtcNow.ToString('s')+'Z') -Force
   $data | Add-Member -NotePropertyName error -NotePropertyValue '' -Force
   break
  }catch{$data=[pscustomobject]@{url=$u;checkedUTC=[datetime]::UtcNow.ToString('s')+'Z';error=$_.Exception.Message;files=@()}}
 }
 [IO.File]::WriteAllText($fp,($data | ConvertTo-Json -Depth 5 -Compress),[Text.UTF8Encoding]::new($false))
 [pscustomobject]@{key=$key;error=$data.error;cached=$false}
} -ThrottleLimit 4
$results | Where-Object error | ConvertTo-Json -Compress
'SDC cached requests: '+$results.Count
