param(
    [string]$DatabaseRoot = 'Z:\SPART-WORK\Data\Juno',
    [ValidateRange(329,362)]
    [int]$StartDay = 329,
    [ValidateRange(329,362)]
    [int]$EndDay = 362,
    [string]$StagingRoot = ''
)

$ErrorActionPreference = 'Stop'
$officialRoot = 'https://pds-ppi.igpp.ucla.edu/data/JNO-J_SW-JAD-5-CALIBRATED-V1.0/DATA/2024'
if ([string]::IsNullOrWhiteSpace($StagingRoot)) {
    $StagingRoot = Join-Path ([System.IO.Path]::GetTempPath()) `
        'Juno_JADE_PJ67_PJ68_staging'
}
if ($EndDay -lt $StartDay) {
    throw 'EndDay must be greater than or equal to StartDay.'
}
$days = $StartDay..$EndDay

function Get-LabelInteger {
    param(
        [string]$Text,
        [string]$Name
    )
    $match = [regex]::Match($Text, "(?m)^\s*$Name\s*=\s*(\d+)")
    if (-not $match.Success) {
        throw "Missing $Name in PDS label."
    }
    return [int64]$match.Groups[1].Value
}

foreach ($day in $days) {
    $dateCode = '2024{0:D3}' -f $day
    $targetDir = Join-Path $DatabaseRoot "JADE_data\2024\$dateCode\ELECTRONS"
    if (-not (Test-Path -LiteralPath $targetDir)) {
        New-Item -ItemType Directory -Path $targetDir -Force | Out-Null
    }

    $stem = "JAD_L50_LRS_ELC_ANY_DEF_${dateCode}_V02"
    $labelName = "$stem.LBL"
    $dataName = "$stem.DAT"
    $remoteDir = "$officialRoot/$dateCode/ELECTRONS"
    $labelPath = Join-Path $targetDir $labelName
    $dataPath = Join-Path $targetDir $dataName
    $databasePart = "$dataPath.part"
    $databaseCopy = "$dataPath.copying"
    $stagingDir = Join-Path $StagingRoot $dateCode
    if (-not (Test-Path -LiteralPath $stagingDir)) {
        New-Item -ItemType Directory -Path $stagingDir -Force | Out-Null
    }
    $dataPart = Join-Path $stagingDir "$dataName.part"

    if (-not (Test-Path -LiteralPath $labelPath)) {
        $labelPart = "$labelPath.part"
        Write-Host "[$dateCode] downloading label"
        & curl.exe --fail --location --retry 5 --retry-delay 3 `
            --output $labelPart "$remoteDir/$labelName"
        if ($LASTEXITCODE -ne 0) {
            throw "curl failed for $labelName (exit $LASTEXITCODE)."
        }
        Move-Item -LiteralPath $labelPart -Destination $labelPath
    }

    $labelText = [System.IO.File]::ReadAllText($labelPath)
    $recordBytes = Get-LabelInteger -Text $labelText -Name 'RECORD_BYTES'
    $fileRecords = Get-LabelInteger -Text $labelText -Name 'FILE_RECORDS'
    $expectedBytes = $recordBytes * $fileRecords

    if (Test-Path -LiteralPath $dataPath) {
        $actualBytes = (Get-Item -LiteralPath $dataPath).Length
        if ($actualBytes -eq $expectedBytes) {
            if (Test-Path -LiteralPath $databasePart) {
                Remove-Item -LiteralPath $databasePart -Force
            }
            if (Test-Path -LiteralPath $databaseCopy) {
                Remove-Item -LiteralPath $databaseCopy -Force
            }
            if (Test-Path -LiteralPath $dataPart) {
                Remove-Item -LiteralPath $dataPart -Force
            }
            Write-Host "[$dateCode] verified existing DAT ($actualBytes bytes)"
            continue
        }
        $stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
        $backupPath = "$dataPath.invalid_$stamp"
        Write-Warning "[$dateCode] existing DAT has $actualBytes bytes; expected $expectedBytes. Moving it to $backupPath"
        Move-Item -LiteralPath $dataPath -Destination $backupPath
    }

    if (-not (Test-Path -LiteralPath $dataPart) -and `
            (Test-Path -LiteralPath $databasePart)) {
        Write-Host "[$dateCode] staging existing database partial file"
        Copy-Item -LiteralPath $databasePart -Destination $dataPart
    }
    Write-Host "[$dateCode] downloading DAT ($expectedBytes bytes)"
    & curl.exe --fail --location --retry 5 --retry-delay 3 `
        --continue-at - --output $dataPart "$remoteDir/$dataName"
    if ($LASTEXITCODE -ne 0) {
        throw "curl failed for $dataName (exit $LASTEXITCODE)."
    }

    $downloadedBytes = (Get-Item -LiteralPath $dataPart).Length
    if ($downloadedBytes -ne $expectedBytes) {
        throw "Downloaded $dataName has $downloadedBytes bytes; expected $expectedBytes."
    }
    $copied = $false
    for ($attempt = 1; $attempt -le 5 -and -not $copied; $attempt++) {
        try {
            if (Test-Path -LiteralPath $databaseCopy) {
                Remove-Item -LiteralPath $databaseCopy -Force
            }
            Write-Host "[$dateCode] copying verified DAT to database (attempt $attempt)"
            Copy-Item -LiteralPath $dataPart -Destination $databaseCopy
            $copiedBytes = (Get-Item -LiteralPath $databaseCopy).Length
            if ($copiedBytes -ne $expectedBytes) {
                throw "Database copy has $copiedBytes bytes; expected $expectedBytes."
            }
            $copied = $true
        }
        catch {
            if ($attempt -eq 5) {
                throw
            }
            Write-Warning "[$dateCode] database copy attempt $attempt failed: $($_.Exception.Message)"
            Start-Sleep -Seconds 5
        }
    }
    if (Test-Path -LiteralPath $dataPath) {
        $concurrentBytes = (Get-Item -LiteralPath $dataPath).Length
        if ($concurrentBytes -eq $expectedBytes) {
            Remove-Item -LiteralPath $databaseCopy -Force
        }
        else {
            throw "A concurrent database DAT has $concurrentBytes bytes; expected $expectedBytes."
        }
    }
    else {
        Move-Item -LiteralPath $databaseCopy -Destination $dataPath
    }
    Remove-Item -LiteralPath $dataPart -Force
    if (Test-Path -LiteralPath $databasePart) {
        Remove-Item -LiteralPath $databasePart -Force
    }
    Write-Host "[$dateCode] saved and verified"
}

Write-Host 'All PJ67-to-PJ68 JADE Level-5 files are present and record-length verified.'
