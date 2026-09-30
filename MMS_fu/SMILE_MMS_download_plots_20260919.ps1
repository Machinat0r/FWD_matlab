param([ValidateSet('burst','quicklook','all')][string]$Mode='all',[int]$Workers=4)
$ErrorActionPreference='Stop'
$cache=Join-Path $env:TEMP 'MMS_plots_20260919'
$root=[IO.Path]::GetFullPath('C:\Users\Administrator\Documents\KH\SMILE_MMS')
$base='https://lasp.colorado.edu/mms/sdc/public/data/sdc/'
$items=[Collections.Generic.List[object]]::new()
$kinds=if($Mode -eq 'all'){@('burst','quicklook')}else{@($Mode)}
foreach($kind in $kinds){
 $name=if($kind -eq 'burst'){'mms_burst_plots.txt.selected'}else{'mms_ql_plots.txt.selected'}
 foreach($line in [IO.File]::ReadLines((Join-Path $cache $name))){
  if($line -notmatch '^\./([a-z0-9_]+)/(\d{4})/(\d{2})/(\d{2})/([a-z0-9_]+\.png)$'){throw ('Unexpected index path '+$line)}
  $plot=$Matches[1];$year=$Matches[2];$month=$Matches[3];$day=$Matches[4];$file=$Matches[5]
  $remote=if($kind -eq 'burst'){'burst'}else{'ql'}
  $relative=$line.Substring(2)
  $instrument=if($plot.StartsWith('all_')){'综合图'}else{($plot -split '_')[0].ToUpperInvariant()}
  $sat=if($plot -match 'mms([1-4])'){'MMS'+$Matches[1]}else{'多星综合'}
  $dest=[IO.Path]::GetFullPath((Join-Path $root ($instrument+'\'+$sat+'\'+$kind+'\'+$file)))
  if(-not $dest.StartsWith($root+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Unsafe destination'}
  $items.Add([pscustomobject]@{Kind=$kind;Plot=$plot;Date=($year+'-'+$month+'-'+$day);Url=($base+$remote+'/'+$relative);Path=$dest})
 }
}
# 用户要求：先下载全部仪器、全部卫星的120分钟quicklook，其余时长随后继续。
$items=@($items | Sort-Object Url -Unique | Sort-Object @{Expression={if($_.Kind -eq 'quicklook' -and $_.Url.EndsWith('_0120.png')){0}else{1}}},@{Expression={if($_.Plot.StartsWith('all_')){0}else{1}}},Date,Plot,Kind,Url)
$log=Join-Path $cache ($Mode+'_download.jsonl')
$progressFile=Join-Path $cache ($Mode+'_progress.json')
$writer=[IO.StreamWriter]::new($log,$true,[Text.UTF8Encoding]::new($false));$writer.AutoFlush=$true
$control=[Collections.Concurrent.ConcurrentDictionary[string,long]]::new()
# 全部线程共用请求间隔；旧版只降速不恢复，新版在连续成功后逐渐恢复。
$control['retryAfter']=0
$control['nextRequest']=0
$control['spacingMs']=2000
$control['successStreak']=0
$control['lastAdjustment']=[datetime]::UtcNow.Ticks
$control['networkDownloaded']=0
$control['rateLimitResponses']=0
$throttleFile=Join-Path $cache 'throttle.json'
if(Test-Path -LiteralPath $throttleFile){
 $saved=Get-Content -LiteralPath $throttleFile -Raw | ConvertFrom-Json
 $control['retryAfter']=[long]$saved.retryAfter
 if($saved.version -eq 2){$control['spacingMs']=[Math]::Max(2000,[long]$saved.spacingMs)}
}
[IO.File]::WriteAllText($throttleFile,(@{version=2;spacingMs=$control['spacingMs'];retryAfter=$control['retryAfter']} | ConvertTo-Json -Compress))
$done=0;$success=0;$failed=0;$bytes=0L;$started=[datetime]::UtcNow;$lastReport=[datetime]::UtcNow
[pscustomobject]@{mode=$Mode;total=$items.Count;workers=$Workers;startUTC=$started.ToString('o')} | ConvertTo-Json -Compress
try {
 $items | ForEach-Object -Parallel {
  $item=$_;$ctrl=$using:control;$throttlePath=$using:throttleFile
  if(-not $script:mmsImageClient){
   $handler=[Net.Http.SocketsHttpHandler]::new();$handler.MaxConnectionsPerServer=2
   $script:mmsImageClient=[Net.Http.HttpClient]::new($handler)
   $script:mmsImageClient.Timeout=[TimeSpan]::FromSeconds(45)
   $script:mmsImageClient.DefaultRequestHeaders.UserAgent.ParseAdd('MMS-Research-Image-Archive/1.0')
  }
  $status='failed';$err='';$size=0;$hash='';$attempt=0;$http=0
  for($attempt=1;$attempt -le 4;$attempt++){
   $response=$null
   try {
    if([IO.File]::Exists($item.Path)){
     $buf=[IO.File]::ReadAllBytes($item.Path)
     if($buf.Length -gt 32 -and [Convert]::ToHexString($buf[0..7]) -eq '89504E470D0A1A0A' -and [Convert]::ToHexString($buf[($buf.Length-12)..($buf.Length-1)]) -eq '0000000049454E44AE426082'){
      $status='existing';$size=$buf.LongLength;$hash=[Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($buf));break
     }
    }
    # 每次发请求前重新检查全局冷却，避免其他线程在收到429后继续发出排队请求。
    while($true){
     $granted=$false
     [Threading.Monitor]::Enter($ctrl)
     try {
      $now=[datetime]::UtcNow.Ticks
      $eligible=[Math]::Max($ctrl['retryAfter'],$ctrl['nextRequest'])
      if($now -ge $eligible){
       $ctrl['nextRequest']=$now+[TimeSpan]::FromMilliseconds($ctrl['spacingMs']).Ticks
       $granted=$true
      }else{$wait=[int][Math]::Min(1000,[Math]::Max(1,[Math]::Ceiling(($eligible-$now)/10000.0)))}
     }finally{[Threading.Monitor]::Exit($ctrl)}
     if($granted){break}
     Start-Sleep -Milliseconds $wait
    }
    $requestStarted=[datetime]::UtcNow
    $http=0
    $response=$script:mmsImageClient.GetAsync($item.Url,[Net.Http.HttpCompletionOption]::ResponseContentRead).GetAwaiter().GetResult()
    $http=[int]$response.StatusCode
    if($http -eq 429 -or $http -eq 503){
     $delay=60
     if($response.Headers.RetryAfter.Delta){$delay=[Math]::Max(1,[Math]::Ceiling($response.Headers.RetryAfter.Delta.TotalSeconds))}
     elseif($response.Headers.RetryAfter.Date){$delay=[Math]::Max(1,[Math]::Ceiling(($response.Headers.RetryAfter.Date.UtcDateTime-[datetime]::UtcNow).TotalSeconds))}
     [Threading.Monitor]::Enter($ctrl)
     try {
      $now=[datetime]::UtcNow.Ticks
      if($now -ge $ctrl['retryAfter']){$ctrl['spacingMs']=[Math]::Min(12000,[Math]::Ceiling($ctrl['spacingMs']*1.25))}
      $ctrl['retryAfter']=[Math]::Max($ctrl['retryAfter'],[datetime]::UtcNow.AddSeconds($delay).Ticks)
      $ctrl['lastAdjustment']=$now
      $ctrl['successStreak']=0
      $ctrl['rateLimitResponses']++
      [IO.File]::WriteAllText($throttlePath,(@{version=2;spacingMs=$ctrl['spacingMs'];retryAfter=$ctrl['retryAfter']} | ConvertTo-Json -Compress))
     }finally{[Threading.Monitor]::Exit($ctrl)}
     throw ('HTTP '+$http+'; server retry after '+$delay+' seconds')
    }
    if(-not $response.IsSuccessStatusCode){throw ('HTTP '+$http)}
    $buf=$response.Content.ReadAsByteArrayAsync().GetAwaiter().GetResult()
    if($response.Content.Headers.ContentLength -and $buf.LongLength -ne $response.Content.Headers.ContentLength){throw 'Content-Length mismatch'}
    if($buf.Length -le 32 -or [Convert]::ToHexString($buf[0..7]) -ne '89504E470D0A1A0A' -or [Convert]::ToHexString($buf[($buf.Length-12)..($buf.Length-1)]) -ne '0000000049454E44AE426082'){throw 'Invalid or incomplete PNG'}
    [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($item.Path)) | Out-Null
    $partial=$item.Path+'.part'
    [IO.File]::WriteAllBytes($partial,$buf)
    [IO.File]::Move($partial,$item.Path,$true)
    $size=$buf.LongLength;$hash=[Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($buf))
    [Threading.Monitor]::Enter($ctrl)
    try {
     $ctrl['networkDownloaded']++
     if($requestStarted.Ticks -ge $ctrl['retryAfter']){$ctrl['successStreak']++}
     if($ctrl['successStreak'] -ge 40 -and ([datetime]::UtcNow.Ticks-$ctrl['lastAdjustment']) -ge [TimeSpan]::FromMinutes(2).Ticks){
      $ctrl['spacingMs']=[Math]::Max(2000,[Math]::Floor($ctrl['spacingMs']*0.85))
      $ctrl['successStreak']=0
      $ctrl['lastAdjustment']=[datetime]::UtcNow.Ticks
      [IO.File]::WriteAllText($throttlePath,(@{version=2;spacingMs=$ctrl['spacingMs'];retryAfter=$ctrl['retryAfter']} | ConvertTo-Json -Compress))
     }
    }finally{[Threading.Monitor]::Exit($ctrl)}
    $status='downloaded';$err='';break
   }catch{
    $err=$_.Exception.Message
    Write-Warning ($item.Url+' attempt '+$attempt+': '+$err)
    if($http -eq 404){break}
    if($attempt -lt 4){Start-Sleep -Seconds ([Math]::Min(20,3*$attempt))}
   }finally{if($response){$response.Dispose()}}
  }
  [pscustomobject]@{kind=$item.Kind;plot=$item.Plot;date=$item.Date;url=$item.Url;path=$item.Path;status=$status;size=$size;sha256=$hash;http=$http;attempts=$attempt;error=$err;checkedUTC=[datetime]::UtcNow.ToString('o')}
 } -ThrottleLimit $Workers | ForEach-Object {
  $rec=$_;$writer.WriteLine(($rec | ConvertTo-Json -Compress))
  $done++;if($rec.status -eq 'failed'){$failed++}else{$success++;$bytes+=$rec.size}
  if(($done % 100 -eq 0) -or (([datetime]::UtcNow-$lastReport).TotalSeconds -ge 20) -or $done -eq $items.Count){
   $p=[pscustomobject]@{mode=$Mode;total=$items.Count;done=$done;success=$success;failed=$failed;bytes=$bytes;networkDownloaded=$control['networkDownloaded'];rateLimitResponses=$control['rateLimitResponses'];spacingMs=$control['spacingMs'];retryAfterUTC=([datetime]::new($control['retryAfter'],[DateTimeKind]::Utc)).ToString('o');elapsedSec=[int]([datetime]::UtcNow-$started).TotalSeconds;updatedUTC=[datetime]::UtcNow.ToString('o');last=$rec.path}
   $json=$p | ConvertTo-Json -Compress
   [IO.File]::WriteAllText($progressFile,$json,[Text.UTF8Encoding]::new($false))
   $json;$lastReport=[datetime]::UtcNow
  }
 }
} finally {$writer.Dispose()}
