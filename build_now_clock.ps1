# build_now_clock.ps1
# Usage: run this from project root in PowerShell (not MSYS2).
# It invokes MSYS2's MinGW bash to compile cmd/main/now_clock.c and produce now_clock.dll
# Then copies the DLL to project root and runs moon clean/build/test.

$proj = (Resolve-Path ".").Path
Write-Host "Project path: $proj"

# Locate MSYS2 bash
$defaultBash = 'C:\msys64\usr\bin\bash.exe'
$bashPath = $null
if (Test-Path $defaultBash) { $bashPath = $defaultBash }
else {
  $where = & where.exe bash 2>$null
  if ($where) { $bashPath = $where[0] }
}

if (-not $bashPath) {
  Write-Host "MSYS2 bash not found. Please install MSYS2 (https://www.msys2.org/) and ensure 'C:\msys64\usr\bin\bash.exe' exists." -ForegroundColor Red
  exit 1
}

# Convert Windows path to MSYS path like /d/path/to/project
if ($proj -match '^([A-Za-z]):\\(.*)') {
  $drive = $matches[1].ToLower()
  $rest = $matches[2] -replace '\\','/'
  $msysProj = "/$drive/$rest"
} else {
  $msysProj = $proj -replace '\\','/'
}

Write-Host "Using MSYS2 bash: $bashPath"
Write-Host "MSYS project path: $msysProj"

# Build commands to run under MSYS2 MinGW 64-bit
$gccCompileCmd = "gcc -O2 -c '$msysProj/cmd/main/now_clock.c' -o '$msysProj/cmd/main/now_clock.o'"
$gccDllCmd = "gcc -shared -O2 -o '$msysProj/cmd/main/now_clock.dll' '$msysProj/cmd/main/now_clock.c' -Wl,--export-all-symbols -Wl,--enable-auto-import"
$fullCmd = "$gccCompileCmd && $gccDllCmd"

Write-Host "Running in MSYS2: $fullCmd"

# Invoke bash -lc
& $bashPath -lc $fullCmd
if ($LASTEXITCODE -ne 0) {
  Write-Host "MSYS2 gcc build failed with exit code $LASTEXITCODE" -ForegroundColor Red
  exit $LASTEXITCODE
}

# Copy DLL to project root for runtime discovery
$srcDll = Join-Path $proj 'cmd\main\now_clock.dll'
if (Test-Path $srcDll) {
  Copy-Item $srcDll -Destination $proj -Force
  Write-Host "Copied now_clock.dll to project root"
} else {
  Write-Host "DLL not found at $srcDll" -ForegroundColor Yellow
}

# Run moon commands
Write-Host "Running: moon clean; moon build; moon test"
& moon clean
& moon build
& moon test

Write-Host "Done."