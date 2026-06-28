$ErrorActionPreference = 'Stop'

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$prefDir = Join-Path $env:TEMP ("matlab_pref_ace20260119_{0}" -f ([guid]::NewGuid().ToString('N')))
New-Item -ItemType Directory -Force -Path $prefDir | Out-Null
$env:MATLAB_PREFDIR = $prefDir

try {
    matlab -batch "cd('$scriptDir'); out=plot_ace_20260119_event('Z:\SPART-WORK\Data\ACE'); disp(out.velocity_note); disp(out.png);"
}
finally {
    Remove-Item -LiteralPath $prefDir -Recurse -Force -ErrorAction SilentlyContinue
}

$pythonExe = 'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
if (!(Test-Path -LiteralPath $pythonExe)) {
    $pythonExe = 'python'
}

$cacheDir = Join-Path $scriptDir 'hapi_cache'
$mfiMat = Join-Path $cacheDir 'AC_H3_MFI_20260119_Magnitude_BGSEc.mat'
$sweMat = Join-Path $cacheDir 'AC_K0_SWE_20260119_Vp_Tpr.mat'
$epmMat = Join-Path $cacheDir 'AC_H3_EPM_20260119T1200_20260120T0000_P1_P2_P3_P4_P5_P6_P7_P8_DE1_DE2_DE3_DE4.mat'
$plot0 = '740001.5'
$plot1 = '740002'
$tBmax = '740001.7914199884'

& $pythonExe (Join-Path $scriptDir 'plot_ace_overview_pillow.py') `
    (Join-Path $scriptDir 'ace_20260119_overview.png') `
    $plot0 $plot1 $tBmax $mfiMat $sweMat $epmMat

& $pythonExe (Join-Path $scriptDir 'plot_ace_overview_vector_pdf.py') `
    (Join-Path $scriptDir 'ace_20260119_overview_editable.pdf') `
    $plot0 $plot1 $tBmax $mfiMat $sweMat $epmMat
