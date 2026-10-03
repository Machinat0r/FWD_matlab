param([ValidateSet('before','after')][string]$Phase='before')
$ErrorActionPreference='Stop'
$work=Join-Path $env:TEMP 'MMS1_marked1_20260930'
$manifest=Get-Content -Raw -LiteralPath (Join-Path $work 'build\manifest.json') | ConvertFrom-Json
$renderDir=Join-Path $work $Phase
[void](New-Item -ItemType Directory -Path $renderDir -Force)
if ($Phase -eq 'before') {$source=$manifest.snapshot} else {$source=Join-Path $work 'output\MMS1_marked1.pptx'}
$app=New-Object -ComObject PowerPoint.Application
$beforeCount=$app.Presentations.Count
$deck=$null
try {
    $deck=$app.Presentations.Open($source,-1,0,0)
    for($i=0;$i -lt $manifest.selected_pages.Count;$i++) {
        if($Phase -eq 'before') {$page=$manifest.selected_pages[$i]} else {$page=$i+1}
        $dest=Join-Path $renderDir ('slide_'+($i+1)+'.png')
        $deck.Slides.Item([int]$page).Export($dest,'PNG',1700,919)
    }
    Write-Output ($Phase+' rendered: '+$manifest.selected_pages.Count)
} finally {
    if ($null -ne $deck) {$deck.Close(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck)}
    if ($beforeCount -eq 0 -and $app.Presentations.Count -eq 0) {$app.Quit()}
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
